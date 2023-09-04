'''
Module: **pydft2kp/rotatebasis.py**
'''

from numpy import array, allclose, trace, eye, kron, \
                  vstack, zeros, exp, diag, \
                  hstack, sum, abs, \
                  append, ones, pi, argwhere
from numpy.linalg import norm
from scipy.linalg import null_space, lstsq
from scipy.optimize import minimize
from sympy import lambdify
from numpy.random import default_rng
rng = default_rng()
from .constants import QSkeys, DFTkeys
from .util import convert_units_coeffs

class basis_transform():
    '''
    Reads the qsymm and irrep data to identify
    the matching symmetries and find the transformation
    matrix U that rotates the bases from DFT into the 
    representation informed into qsymm.

    Parameters
    ---------
    qsymm : object of our qsymm class
        Model built with our qsymm class
    irrep : object of our irrep class
        DFT data read by the irrep package
    nullspace_thresh : float
        Threshold for the null space calculation
    diagonal : bool
        If True, the U matrix enforces that H is diagonal at k=0.

    Attributes
    ----------
    qs2irrep : array
        Table relating indices of the qsymm symmetries to
        the corresponding ones in the irrep data.
    num_phases : int
        Number of free phases in the definition of the matrix U.
        These are relative phases between the different irreps that
        compose U.
    U : array
        The unitary transformation matrix.
    Q : array
        The Q matrix, built as in https://mathoverflow.net/q/391653/147408
    US, USconj : list of arrays
        The nullspace of Q in the form of partial matrices that compose U
    report : output of minimize
        Full output of the minimize call, stored for testing purposes
    coeffs : array
            Values for each of the cn QSYMM model coefficients.
    keys : array
        Label for the powers of k related to each coefficient.
    Heff : callable
        A function Heff(kx,ky,kz) for the effective Hamiltonian
        as a function of k.
    '''
    def __init__(self, qsymm, irrep, nullspace_thresh=1e-4, diagonal=True):
        # match qsymm and irrep symmetries
        self.qs2irrep = self.match_qsymm_irrep(qsymm, irrep)
        
        # check traces
        TRQS = []
        TRDFT = []
        for i,j in self.qs2irrep:
            TRQS += [trace(qsymm.symms[i].U)]
            TRDFT += [trace(irrep.GammaDFT[j])]
        TRQS = array(TRQS)
        TRDFT = array(TRDFT)
        if allclose(TRQS, TRDFT) == False:
            # check also the anti-unitary
            raise Exception("QSYMM and DFT traces do not match!")

        # find U
        results = self.get_basis_transform(qsymm, irrep, nullspace_thresh)
        # extract results
        self.num_phases = results[0]
        self.US = results[1]
        self.USconj = results[2]
        self.Q = results[3]

        # check if size of null space match deggree of freedom from irreps
        if self.num_phases != irrep.deg_of_freedom:
            print(f'WARNING: nullity ({self.num_phases}) does not match irreps degree of freedom ({irrep.deg_of_freedom}). \
                  There might be symmetry constraints missing. Please check if you are using all symmetry group generators.')

        # minimize residues to find linear combination
        # of the null space and get the optimized U
        residue = 1
        count = 0
        while (residue > 1e-6) and (count < 50):
            self.report, self.U = self.findU(qsymm, irrep, diagonal)
            residue = self.report.fun
            count += 1
        if residue > 1e-6:
            raise Exception(f'WARNING: findU did not converge after {count} trials.')

        # apply rotation U
        self.coeffs, self.keys, self.Heff = self.getHeff(qsymm, irrep)

    def match_qsymm_irrep(self, qsymm, irrep):
        '''
        Finds the matching symmetries in qsymm and irrep.

        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package

        Returns
        -------
        table : array
            Matrix where each line refers to a matching symmetry
            and the columns indicate their indices in qsymm and irrep.
        '''
        table = [] # conversion table
        for i in range(len(qsymm.symm_pg)):
            for j in range(len(irrep.symm_pg)):
                # if labels (1,I,S,R,M) don't match, keep searching..
                if qsymm.symm_pg[i][0] != irrep.symm_pg[j][0]: continue
                # if angles don't match, keep searching...
                if qsymm.symm_pg[i][1] != irrep.symm_pg[j][1]: continue
                # if labels and angle match, check axis
                qs_axis = array(qsymm.symm_pg[i][2], dtype=float)
                qe_axis = array(irrep.symm_pg[j][2], dtype=float)
                # normalize axis
                qs_axis = qs_axis / norm(qs_axis)
                qe_axis = qe_axis / norm(qe_axis)
                if (qsymm.symm_pg[i][0] == '1') or (qsymm.symm_pg[i][0] == 'I'):
                    # (1,I): axis is irrelevant, always [0,0,1]
                    table += [[i,j]]
                elif qsymm.symm_pg[i][1] == 180: # both M and C2 have angle = 180
                    # M: sign of axis (normal) does not matter, check +/-
                    # C2 = R(180): sign of axis does not matter, check +/-
                    if allclose(qs_axis, +qe_axis) or allclose(qs_axis, -qe_axis):
                        table += [[i,j]]
                else: # R or S with angle != 180
                    # (R,S): sign is important: R(-theta, axis) = R(theta, -axis)
                    if allclose(qs_axis, +qe_axis):
                        table += [[i,j]]
        
        # number of qsymm irreps without TRS and chiraliy
        # must match number of found symms
        nqsymm = 0
        for S in qsymm.symms:
            if (S.conjugate == False) and (S.antisymmetry == False):
                nqsymm += 1
        if len(table) != nqsymm:
            raise Exception('Qsymm symmetries not found in DFT.')

        return table


    def get_basis_transform(self, qsymm, irrep, nullspace_thresh=1e-4):
        '''
        Finds the set of transformation matrices US that compose the 
        transformation U as a linear combination.
        
        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package
        nullspace_thresh : float
            Threshold to identify the nullspace, or zeros of the SVD.
        
        Returns
        -------
        nullity : int
            Number of matrices US
        US : list of arrays
            The matrices US
        USconj : list of arrays
            The conjugate of US. Returned for testing purposes only.
        Q : array
            The Q matrix, built as in https://mathoverflow.net/q/391653/147408
        '''        

        if irrep.antiU != []:
            # if there are anti-unitary symmetries, call the generalized version
            US, USconj, Q = self.get_basis_transform_with_trs(qsymm, irrep, nullspace_thresh)
            nullity = US.shape[0]
        
        else: # there there are no anti-unitary, calculation is simpler

            A = [] # list from DFT   = original basis
            B = [] # list from QSYMM = desired  basis
            for g1,g2 in self.qs2irrep:
                A.append(irrep.GammaDFT[g2])
                B.append(qsymm.symms[g1].U)

            M = len(A) # number of blocks
            N = len(A[0]) # size of each A, B, U
            idN = eye(N) # identity on hilbert space
            
            Q = kron(idN, A[0].T) - kron(B[0], idN)
            for i in range(1, M):
                Qi = kron(idN, A[i].T) - kron(B[i], idN)
                Q = vstack([Q, Qi])
            
            US = null_space(Q, nullspace_thresh)
            if US.size == 0:
                raise Exception("Empty nullspace")
            
            nullity = US.shape[1]
            aux = zeros((nullity, N, N), dtype=complex)
            for i in range(nullity):
                aux[i] = US[:,i].reshape(N,N)
            US = aux
            USconj = US.conj() # for compatibility with TRS case

        # END IF
        # return common from with and without TRS cases
        return nullity, US, USconj, Q

    def get_basis_transform_with_trs(self, qsymm, irrep, nullspace_thresh=1e-4):
        '''
        Finds the set of transformation matrices US that compose the 
        transformation U as a linear combination.
        
        This is the generalized
        routine for cases where there are anti-unitary symmetries.

        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package
        nullspace_thresh : float
            Threshold to identify the nullspace, or zeros of the SVD.
        
        Returns
        -------
        U1 : list of arrays
            The matrices US
        U2 : list of arrays
            Should be U2 = U1*
        Q : array
            The Q matrix, built as in https://mathoverflow.net/q/391653/147408
        '''
        A = []
        B = []
        for g1,g2 in self.qs2irrep:
            A.append(irrep.GammaDFT[g2])
            B.append(qsymm.symms[g1].U)

        NB = len(A) # number of blocks
        N = len(A[0]) # order of the irrep. Size of A,B,U = NxN
        idN = eye(N) # identity on hilbert space
        
        Q11 = kron(idN, A[0].T) - kron(B[0], idN) # U
        Q12 = zeros([N**2,N**2]) # U*
        Q21 = zeros([N**2,N**2]) # U
        Q22 = kron(idN, A[0].T.conj()) - kron(B[0].conj(), idN) # U*
        Q1 = hstack([Q11, Q12])
        Q2 = hstack([Q21, Q22])
        for i in range(1, NB):
            Q11 = kron(idN, A[i].T) - kron(B[i], idN) # U
            Q12 = zeros([N**2,N**2]) # U*
            Q21 = zeros([N**2,N**2]) # U
            Q22 = kron(idN, A[i].T.conj()) - kron(B[i].conj(), idN)  # U*
            Q1 = vstack([Q1, hstack([Q11, Q12])])
            Q2 = vstack([Q2, hstack([Q21, Q22])])

        for A,B in irrep.antiU:
            Q11 = -kron(B, idN) # U
            Q12 =  kron(idN, A.T) # U*
            Q21 =  kron(idN, A.T.conj()) # U
            Q22 = -kron(B.conj(), idN)  # U*
            Q1 = vstack([Q1, hstack([Q11, Q12])])
            Q2 = vstack([Q2, hstack([Q21, Q22])])

        Q = vstack([Q1, Q2])

        Uall = null_space(Q, nullspace_thresh)
        if Uall.size == 0:
            raise Exception("Empty nullspace")

        nullity = Uall.shape[1]
        # Uall shape is (2*N², nullity)
        # flip it to put nullity first
        U1 = zeros((nullity, N, N), dtype=complex)
        U2 = zeros((nullity, N, N), dtype=complex)
        for i in range(nullity):
            U1[i] = Uall[:(N**2),i].reshape(N,N) # U
            U2[i] = Uall[(N**2):,i].reshape(N,N) # U*

        return U1, U2, Q
        
    def getU(self, pars):
        '''
        Builds the matrix U from the linear combination
        set by the pars coefficients.

        Parameters
        ----------
        pars : array
            Intensities and phases of the linear combination of US matrices

        Returns
        -------
        U : array
            The U matrix.
        '''
        nu, N = self.US.shape[:2]
        U = zeros([N,N], dtype=complex)
        for i in range(nu):
            U += pars[i]*exp(1j*pars[i+nu]) * self.US[i]
        return U

    def errorU(self, pars, qsymm, irrep, diagonal=True):
        '''
        Calculates the error for a given matrix U(pars).

        Parameters
        ----------
        pars : array
            Set of parameters that define U
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package
        diagonal : bool
            If True, the U matrix enforces that H is diagonal at k=0.

        Returns
        -------
        residue : float
            The error estimative
        '''
        U = self.getU(pars)
        N = U.shape[0]
        
        # U must be unitary
        residue = sum(abs(U.T.conj() @ U - eye(N)))**2
        
        # B = U.A.U† (for unitary symmetries A and B)
        for g1,g2 in self.qs2irrep:
            A = irrep.GammaDFT[g2]
            B = qsymm.symms[g1].U
            residue += sum(abs(B - U @ A @ U.T.conj()))**2
        
        # B = U*.A.U† (for anti-unitary symmetries A and B)
        if irrep.antiU != []:
            for A,B in irrep.antiU:
                # B = U†.A.U, A,B antiunitary
                residue += sum(abs(B - U.conj() @ A @ U.T.conj()))**2
                
        # H0 must remain diagonal (optional)
        if diagonal:
            emin = irrep.energies[irrep.setA][0]
            emax = irrep.energies[irrep.setA][-1] # assuming ordered
            # scale energies to 0..1
            energies = (irrep.energies[irrep.setA] - emin) / (emax-emin)
            H0 = U @ diag(energies) @ U.T.conj()
            residue += sum(abs(H0 - diag(diag(H0))))**2
        
        # total error
        return residue

    def findU(self, qsymm, irrep, diagonal=True):
        '''
        Minimizes the error to find the optimal matrix U.

        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package

        Returns
        -------
        sol : solution
            The full solution from the minimize call
        U : array
            The optimal matrix U
        '''
        # nullity and Hilbert space size
        nu, N = self.US.shape[:2]
        # initial guess with intensity 1 and random phases
        x0 = append(ones(nu), rng.uniform(0, 2*pi, nu))
        # minimizes the error
        sol = minimize(self.errorU, x0, tol=1e-12, args=(qsymm, irrep, diagonal))
        # returns the solution and U
        return sol, self.getU(sol.x)

    def getHeff(self, qsymm, irrep):
        '''
        Compares DFT and QSYMM to extract the numerical
        values of the coefficients and build the effective
        Hamiltonian.

        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package

        Returns
        -------
        coeffs : array
            Values for each of the cn QSYMM model coefficients.
        keys : array
            Label for the powers of k related to each coefficient.
        Heff : callable
            A function Heff(kx,ky,kz) for the effective Hamiltonian
            as a function of k.
        '''
        coeffs, keys = self.get_coeffs(qsymm, irrep)
        tosub = []
        for n in range(len(coeffs)):
            tosub += [('c'+str(n), coeffs[n])]
        return coeffs, keys, lambdify(['k_x', 'k_y', 'k_z'], qsymm.Hmodel.subs(tosub))
    
    def get_coeffs(self, qsymm, irrep):
        '''
        Compares the DFT data to the QSYMM model
        and extracts the numerical value for the coefficients cn.

        Parameters
        ----------
        qsymm : qsymm object
            Model built with our qsymm class
        irrep : irrep object
            DFT data read by the irrep package
        
        Returns
        -------
        coeffs : array
            Values for each of the cn QSYMM model coefficients.
        keys : array
            Label for the powers of k related to each coefficient.
        '''
        N = len(self.U)
        # compare DFT and QSYMM
        # builds and solve a system of equations
        equations = []
        values = []
        for i in range(N):
            for j in range(N):
                for qsk, qek in zip(QSkeys, DFTkeys):
                    line = [q[qsk][i,j] for q in qsymm.model]
                    Haux = self.U @ irrep.Hdict[qek] @ self.U.T.conj()
                    if not allclose(line, 0):
                        equations += [line]
                        values += [Haux[i,j]]
        equations = array(equations)
        sol = lstsq(equations, values)
        # coeffs are real by construction
        coeffs = sol[0].real

        # identify the labels of k-powers on each coefficient
        ncoeffs = len(qsymm.model)
        keys = {}
        for n in range(ncoeffs):
            keys[n] = [] # init empty
        for qsk, qek in zip(QSkeys, DFTkeys):
            if qek not in irrep.Hdict.keys():
                continue
            for n in range(ncoeffs):
                ij = argwhere(abs(qsymm.model[n][qsk]) > 1e-4)
                if len(ij) > 0:
                    keys[n] += [qek]
        # return value and label of the coefficients
        return coeffs, keys

    def print_report(self, sigdigits=5):
        """
        Convert coefficient values from a.u. to eV and nm.

        Parameters
        ----------
        sigdigits : int
            Number of significant digits printed in the report, by default 5

        Returns
        -------
        array
            Coefficients with units converted to eV and nm.
        """
        return convert_units_coeffs(self.coeffs, self.keys, True, sigdigits)
