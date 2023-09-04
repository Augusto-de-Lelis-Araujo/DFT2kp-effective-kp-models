'''
Doc for the **pydft2kp/irrepwrapper.py** module.

The class **irrep** is a wrapper to call
the routines from the `IrRep package <https://github.com/stepan-tsirkin/irrep/>`_.

The routines **p_matrix**, **symm_matrix**, and **sigma_matrix** 
are modified versions of the `symm_matrix` routine from the irrep package file *gvectors.py*.
'''

import numpy as np
from irrep.gvectors import transformed_g 
import irrep.bandstructure as bs
from .util import cwd, R_to_spin, R_to_bvec
from .constants import Ry, a0, sx, sy, sz, alpha
from .qe_aux import read_espresso, read_kp_dat
from .lowdin import getHpowers, H_of_k

class irrep():
    """
    Uses the irrep package to read the DFT data and calculate
    the matrix elements of the symmetry operations and the
    matrix elements of the momentum operator.
    
    Parameters from the **irrep package**, class bandstructure, file bandstructure.py
    
    - fWAV 
    - fWFK 
    - fPOS 
    - Ecut 
    - IBstart 
    - IBend 
    - spinor 
    - code 
    - EF 
    - onlysym 
    - spin_channel 
    - refUC 
    - shiftUC 
    - identify_irreps 
    
    Parameters
    ----------
    dftdir : str
        Path to the directory where the DFT data is stored.
    outdir : str
        Name of the directory informed as outdir to QE.
    prefix : str
        Prefix of the DFT data on QE
    kpt : int
        Single number indicating the k point for the expansion.
        Equivalent to the kplist from the irrep package, but here
        it must be a single point.
    kname : str
        Label of the kpoint as defined in the irrep package.
    degen_thresh : float
        Threshold for the energy degeneracy.
        
    Attributes
    ----------
    dftdir : string 
        Path to the directory where the DFT data is stored.
    prefix : string
        QE prefix parameter
    outdir : string
        QE outdir parameter
    kindex : int
        index for the k point along the k path    
    alat : float
        The lattice constant in Bohr units.
    fermi : float
        The Fermi energy in Rydberg units.
    bandstr : irrep object
        The main irrep package object that stores the information read from the DFT data.
    energies : array
        Band energies in Rydberg units.
    symm_pg : list
        Poing group (S) components of the symmetry operation {S,T} (Seitz notation).
    symm_translation : list
        Translation (T) components of the symmetry operation {S,T} (Seitz notation).
    irreps : list
        List of degenerate bands with their irrep identification. Each entry contains [band indices, irreps list, degeneracy].
    antiU : list
        Matrix representations of the symmetry operators.
        Each line correspond to an operator. The first column
        refers to the matrix representation calculated from DFT,
        and the second column the one read from the QSYMM object. 
    Hdict : dict
        Dictionary with the matrices that multiply the powers of momentum.
    setA : list or slice
        Stores the informed set A.
    num_irreps : int
        Number of irreps that compose set A.
    deg_of_freedom : int
        Number of free parameters that define the
        unitary transformation matrix U
    GammaDFT : list
        List of matrices for each symmetry operation
    p[x|y|z] : array
        Matrix of the momentum operator along the [x|y|z]-axis in Bohr units.
    sigma_[x|y|z] : array
        Matrix of the spin operator along the [x|y|z]-axis.
    psoc_[x|y|z] : array
        Matrix of the SOC momentum operator along the [x|y|z]-axis in Bohr units.
    """
    def __init__(self,
            dftdir='.',
            outdir='.',
            fWAV=None,
            fWFK=None,
            prefix=None,
            fPOS=None,
            Ecut=None,
            IBstart=None,
            IBend=None,
            kpt=None, # single kpt instead of kplist
            spinor=None,
            code='espresso',
            EF='auto',
            onlysym=False,
            spin_channel=None,
            refUC = None,
            shiftUC = None,
            identify_irreps = False,
            kname='',
            degen_thresh=1e-4
        ):

        kplist = np.ascontiguousarray(kpt)
        assert len(kplist) == 1, 'kpt should refer to a single k point.'

        # store paths and kpoint
        # ----------------------
        self.dftdir = dftdir 
        self.prefix = prefix
        self.outdir = outdir
        self.kindex = kpt

        if code == 'espresso':
            self.alat, self.fermi = read_espresso(dftdir, prefix, outdir)
        else:
            raise Exception("Code not ready for " + code)

        # runs the irrep code on the outdir directory
        with cwd(dftdir + '/' + outdir):
            self.bandstr = bs.BandStructure(fWAV, fWFK, 
                                    prefix, fPOS, Ecut,
                                    IBstart, IBend, kplist,
                                    spinor, code, EF, onlysym, 
                                    spin_channel, refUC,
                                    shiftUC, identify_irreps)
        
        # extracts band energies and stores in a.u.
        self.energies = self.bandstr.kpoints[0].Energy / Ry - self.fermi
        # identify the symmetry operations
        self.symm_pg, self.symm_translation = self.identify_symmetries()
        # identify irreps
        self.irreps = self.get_irreps(kname, degen_thresh)
        # init list of anti-unitary symmetries
        self.antiU = []

    def fold_down_H(self, NB=None, maxorder=2):
        """
        Uses Löwdin partitioning to calculate a dictionary
        with the matrices for each power of k.

        Example: key 'xx' refers to Hxx * k_x^2.
            
        The calculation is done in the crude, original QE basis.

        Parameters
        ----------
        NB : int or None
            Number of bands to consider the set B above set A.
        maxorder : int
            Maximum power of momentum.
        """
        self.Hdict = getHpowers(self, NB, maxorder)
    
    def build_H_of_k(self, all_bands=False):
        """
        Builds a callable function H(kx, ky, kz, [maxorder=2]).

        Parameters
        ----------
        all_bands : bool
            If True, builds a first order model with all bands.
            If False, builds the folded model for set A only up to maxorder.

        Returns
        -------
        callable
            H(kx, ky, kz, [maxorder=2]).
        """

        if all_bands:
            return H_of_k(self)
        else:
            return H_of_k(self, self.Hdict)


    def identify_symmetries(self):
        """
        Reads the symmetry operations from the k point
        and stores it as human readable labels.

        Seitz notation {S,T}, where S is the point group part, 
        and T the translation part of the symmetry operation.

        Returns
        -------
        symm_pg : list
            Poing group (S) component of the
            symmetry operation {S,T} (Seitz notation). See Notes.
        symm_translation : list of arrays(float)
            Translation component of each symmetry operation
            in catesian coordinates.
            
        Notes
        -----
        
        Each element in **symm_pg** describes the point group part of the
        symmetry operation.Each entry has three components [S, theta, axis]:
        
        [S=1/I/S/C/M] 
            The first component S is a string that identifies the type of operation:
            1 for identity, I for inversion, S for proto-rotation,
            C for rotation, M for mirror. 
        
        [theta] 
            The second entry (int) is the angle theta (degrees). 
            For I and 1 we use angle = 0, and for M the angle is 180. 
        
        [axis] 
            The third entry is the axis (array, shape (3,), float).
            For I and 1 we use [0,0,1]. For C and S it is the
            rotation axis, and for M it is the normal to the mirror plane.
            In all cases the axis are in cartesian coordinates.
        """
        symm_pg = [] # label point group elements
        symm_translation = [] # label translations
        for op in self.bandstr.kpoints[0].symmetries:
            theop = [] # [pg, deg, axis], pg=(1,I,S,R,M)
            angle = round(np.rad2deg(op.angle))
            if angle == 0:
                if op.inversion:
                    theop = ['I', 0, [0,0,1]]
                else:
                    theop = ['1', 0, [0,0,1]]
            elif np.abs(angle) != 180:
                if op.inversion: # proto-rotation
                    theop = ['S', angle, list(np.round(op.axis, 2))]
                else: # rotation
                    theop = ['R', angle, list(np.round(op.axis, 2))]
            else: # angle = 180
                if op.inversion: # mirror
                    theop = ['M', 180, list(np.round(op.axis, 2))]
                else: # C2 rotation
                    theop = ['R', angle, list(np.round(op.axis, 2))]
            # store
            symm_pg += [theop]
            symm_translation += [np.round(op.translation, 2)]
        # return lists
        return symm_pg, symm_translation

    def get_irreps(self, kname='', degen_thresh=1e-4):
        """
        Identifies the irreps of each set of degenerate bands.

        Parameters
        ----------
        kname : str
            Label for the k point.
        degen_thresh : float
            Energy threshold (in eV) to identify degeneracies.

        Returns
        -------
        list
            List of degenerate bands with their irrep identification.
            Each entry contains [band indices, irreps list, degeneracy].
        """
        # get irrep tables for the k point
        irrep_tables = self.bandstr.spacegroup.get_irreps_from_table(kname, self.bandstr.kpoints[0].K)
        # run write_characters to get the jdata dict with the irreps
        # this code also creates the file irreps.dat, but we don't use it.
        _, _, _, jdata = self.bandstr.kpoints[0].write_characters(
            degen_thresh = degen_thresh,
            symmetries = [op.ind for op in self.bandstr.kpoints[0].symmetries],
            irreptable = irrep_tables)

        # reads jdata and stores only the necessary information
        # each entry has: band indices, irreps list, degeneracy
        irreps = []
        idb = 0
        for bnd,dim in zip(jdata['irreps'], jdata['dimensions']):
            txt = ''
            for key in bnd.keys():
                if bnd[key][0] > 0.3:
                    txt += '('+key+')'
            irreps.append([list(range(idb, idb+dim)), txt, dim])
            idb += dim

        return irreps
    
    def define_set_A(self, setA, verbose=True, NB=None, maxorder=2):
        """
        Verifies if the chosen set A is composed by full sets of irreps.
        If not, raises an error. If successful, defines set A and applies fold down.

        Parameters
        ----------
        setA : list or slice
            List of band indices to be used in the calculation.
        verbose : bool, optional
            If True, prints information about the irreps from set A.
        NB : int or None, optional
            Number of bands to consider the set B above set A.
        maxorder : int, optional
            Maximum power of momentum.

        Attributes
        ----------
        setA : list or slice
            Stores the informed set A.
        num_irreps : int
            Number of irreps that compose set A.
        """
        # print space group number and name
        if verbose:
            print('Space group ',
                    str(self.bandstr.spacegroup.number),
                    ':', self.bandstr.spacegroup.name)
            print('Group of the k-vector: <code not ready>')
            print('Verifying set A:', setA)

        # print report and store band indices
        # that match setA
        bandindices = []
        self.num_irreps = 0
        irreps_list = ''
        for ir in self.irreps: # band indices, irreps list, degeneracy
            if set(ir[0]) & set(setA):
                bandindices += ir[0]
                irreps_list += ir[1] # concatenate as string
                # irreps list has the structure: (-K7)(-K9)
                # the split of )( help us count the number of irreps in the list
                self.num_irreps += len(ir[1].split(')('))
                if verbose:
                    print('Band indices:', ir[0], 'Irreps:', ir[1], 'Degeneracy:', ir[2])

        # count occurrences of each irrep and 
        # calculate number of degrees of freedom (dog)
        # this will be the nullity in the basis transformation
        irreps_list = irreps_list[1:-1].split(')(') # remove (...) and split into list
        irreps_count = {} # init empty dict, entries are irreps
        for ir in irreps_list:
            if ir in irreps_count.keys():
                irreps_count[ir] += 1 # increment number of that irrep
            else:
                irreps_count[ir] = 1 # init if it is the first found
        # now calculate the dog = number of paramters in U(n) = n²
        # where n is the number of repeated irreps
        self.deg_of_freedom = 0
        for key in irreps_count.keys():
            self.deg_of_freedom += irreps_count[key]**2    
        
        # check if setA is complete
        if set(bandindices) != set(setA):
            raise Exception('Set A does not match a full set of irreps.')
        
        # store setA
        self.setA = setA

        # apply folding down
        self.fold_down_H(NB, maxorder)

    def get_symm_matrices(self, setA=None, store=True):
        """
        Calculates the symmetry matrices of all operations of the k point,
        for the bands in set A.

        Parameters
        ----------
        setA : list, slice, optional
            List of band indices to be used in the calculation.
        store : bool, optional
            If True, stores the matrices as an attribute of the object.

        Returns
        -------
        GammaDFT : list
            List of matrices for each symmetry operation
        """
        if setA is None:
            setA = self.setA

        GammaDFT = []
        for op in self.bandstr.kpoints[0].symmetries:
            # calls our modified symm_matrix instead of the original
            # on of the irrep package
            GammaDFT += [symm_matrix(
                            self.bandstr.kpoints[0].K,
                            self.bandstr.kpoints[0].RecLattice, 
                            self.bandstr.kpoints[0].WF, 
                            self.bandstr.kpoints[0].ig, 
                            op.rotation, 
                            op.spinor_rotation, 
                            op.translation, 
                            op.spinor,
                            setA)]
        if store:
            self.GammaDFT = GammaDFT
        return GammaDFT

    def get_p_matrices(self, setA=slice(None), SOC=False, qekp=''):
        """
        Calculates the matrices of the momentum operator in Bohr units
        for the bands in set A.

        Parameters
        ----------
        setA : list, slice, optional
            List of band indices to be used in the calculation.
        SOC : bool, optional
            If True, estimates p_soc and calculates sigma matrices
        qekp : str, optional
            Name of the QE file with the matrix elements of p
            
        Notes
        -----
        If qekp is informed, reads the full matrix elements of
        velocity operator from the QE data generated by our patch.
        
        Otherwise, calculates matrix elements of p without SOC or
        PAW corrections, and estimates values for the p_soc and 
        sigma (spin) matrix elements (testing purposes only)
        """

        if qekp != '': # read p from kp.dat
            kpdat = self.dftdir + '/' + qekp
            aux = read_kp_dat(kpdat)
            self.px = np.copy(aux.p1)
            self.py = np.copy(aux.p2)
            self.pz = np.copy(aux.p3)
            del aux
        else: # do no use kp.dat, calculate p from plane-waves
            # reciprocal lattice vector in Bohr units
            # RecLattice is in 1/Angstrom, 
            # so we multiply by 10*a0 = 0.529177249 [Angstrom]
            bvec = self.bandstr.RecLattice * (10*a0)
            self.px = p_matrix(self.bandstr.kpoints[0].K, bvec, self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, 0, setA)
            self.py = p_matrix(self.bandstr.kpoints[0].K, bvec, self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, 1, setA)
            self.pz = p_matrix(self.bandstr.kpoints[0].K, bvec, self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, 2, setA)

        if SOC: # estimate SOC and sigma matrices
            # matrix elements of sigma <m|sigma|n>
            self.sigma_x = sigma_matrix(self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, sx, setA)
            self.sigma_y = sigma_matrix(self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, sy, setA)
            self.sigma_z = sigma_matrix(self.bandstr.kpoints[0].WF, self.bandstr.kpoints[0].ig, self.bandstr.spinor, sz, setA)

            # aux function: <m|s_mu . V_nu|n>
            sV = lambda s, p: 1j*(np.einsum('n,mj,jn->mn',self.energies, s, p) - np.einsum('j,mj,jn->mn',self.energies, s, p))
            # <m|psoc_x|n> = (α²/8) (σy.Vz - σz.Vy)
            # <m|psoc_y|n> = (α²/8) (σz.Vx - σx.Vz)
            # <m|psoc_z|n> = (α²/8) (σx.Vy - σy.Vx)
            self.psoc_x = ((alpha**2)/8) * (sV(self.sigma_y, self.pz) - sV(self.sigma_z, self.py))
            self.psoc_y = ((alpha**2)/8) * (sV(self.sigma_z, self.px) - sV(self.sigma_x, self.pz))
            self.psoc_z = ((alpha**2)/8) * (sV(self.sigma_x, self.py) - sV(self.sigma_y, self.px))

    def add_antiunitary_symm(self, QS, T, bands=None):
        '''
        Calculates the matrix representation of the anti-unitary symmetry operator
        S_mn = <Psi_m|{A|T}|Psi_n>, where A is the anti-unitary operator, assuming 
        that the symmetry is anti-unitary and applies the complex conjugate on the 
        <Psi_m|.

        Parameters
        ----------
        QS : QSYMM object
            The corresponding symmetry as a QSYMM object.
        T : array, shape=(3,)
            Translational part of the symmetry operation, in terms of the basis 
            vectors of the unit cell.
        bands : slice, list, array, shape=(N,), optional
            Selects which bands (m,n) will be used in the calculation. If None, all 
            bands will be used. Defaults to None.
        '''
        if bands == None:
            bands = self.setA

        # convert R to reciprocal space and to spin space
        A = R_to_bvec(QS.R, self.bandstr.kpoints[0].RecLattice)
        S = R_to_spin(QS.R)
        # multiplies by TRS
        S = 1j*sy @ S
        # calculate U
        U = symm_matrix(self.bandstr.kpoints[0].K,
                        self.bandstr.kpoints[0].RecLattice, 
                        self.bandstr.kpoints[0].WF, 
                        self.bandstr.kpoints[0].ig, 
                        A, S, T, self.bandstr.spinor, bands, True)
        # add to list of antiU
        self.antiU += [[U, QS.U]]


##################################################
# MODIFIED VERSIONS OF THE symm_matrix ROUTINE   #
# FROM THE gvectors.py FILE OF THE IRREP PACKAGE #
##################################################

def p_matrix(K, RecLattice, WF, igall, spinor, xyz, bands=slice(None)):
    '''
    Calculates the matrix elements of the momentum operator <m|p_nu|n> along
    direction nu=xyz.

    This routine is based on the symm_matrix routine from the irrep package (gvectors.py).
    
    Parameters
    ----------
    K : array, shape=(3,)
        Direct coordinates of the k-point.
    RecLattice : array, shape=(3,3)
        Each row contains the cartesian coordinates of a basis vector forming 
        the unit-cell in reciprocal space.
    WF : array
        `WF[i,j]` contains the coefficient corresponding to :math:`j^{th}`
        plane-wave in the expansion of the wave-function in :math:`i^{th}`
        band. It contains only plane-waves if energy smaller than `Ecut`.
    igall : array
        Returned by `__sortIG`.
        Every column corresponds to a plane-wave of energy smaller than 
        `Ecut`. The number of rows is 6: the first 3 contain direct 
        coordinates of the plane-wave, the fourth row stores indices needed
        to short plane-waves based on energy (ascending order). Fitfth 
        (sixth) row contains the index of the first (last) plane-wave with 
        the same energy as the plane-wave of the current column.
    spinor : bool
        `True` if wave functions are spinors, `False` if they are scalars.
    xyz : int
        Values 0, 1, 2 refer to the cartesian coordinates x, y, z.
    bands : slice, list, array, shape=(N,)
        Selects which bands (m,n) will be used in the calculation.

    Returns
    -------
        Matrix of the momentum operator in the basis of eigenstates of the 
        Bloch Hamiltonian :math:`H(k)`.
    '''
    npw = igall.shape[1]
    WF = WF[bands] # local copy if sliced, original (copy by reference) if not sliced
    if spinor:
        WF = np.stack([WF[:, :npw], WF[:, npw:]], axis=2)
        return np.einsum("mgs,ngs,ig,i->mn", WF.conj(), WF, igall[:3,:] + K[:, None], RecLattice[:,xyz])
    else:
        return np.einsum("mg ,ng ,ig,i->mn", WF.conj(), WF, igall[:3,:] + K[:, None], RecLattice[:,xyz])

def symm_matrix(K, RecLattice, WF, igall, A, S, T, spinor, bands=slice(None), TRS=False):
    """
    Modifies the original symm_matrix routine from the irrep package (gvectors.py).
    The original routine computes the matrix S_mn = <Psi_m|{A|T}|Psi_n>.
    Here we add two parameters: bands and TRS.

    The bands parameter controls which bands will be used in the calculation.
    The default is to use all bands.

    If TRS == True the code assumes that the symmetry operation is anti-unitary,
    adding the complex conjugation to the symmetry operation, acting to the left.

    Parameters
    ----------
    K : array, shape=(3,)
        Direct coordinates of the k-point.
    RecLattice : array, shape=(3,3)
        Each row contains the cartesian coordinates of a basis vector forming 
        the unit-cell in reciprocal space.
    WF : array
        `WF[i,j]` contains the coefficient corresponding to :math:`j^{th}`
        plane-wave in the expansion of the wave-function in :math:`i^{th}`
        band. It contains only plane-waves if energy smaller than `Ecut`.
    igall : array
        Returned by `__sortIG`.
        Every column corresponds to a plane-wave of energy smaller than 
        `Ecut`. The number of rows is 6: the first 3 contain direct 
        coordinates of the plane-wave, the fourth row stores indices needed
        to short plane-waves based on energy (ascending order). Fitfth 
        (sixth) row contains the index of the first (last) plane-wave with 
        the same energy as the plane-wave of the current column.
    A : array, shape=(3,3)
        Matrix describing the tranformation of basis vectors of the 
        reciprocal space vectors b_i
    S : array, shape=(2,2)
        Matrix describing how spinors transform under the symmetry.
    T : array, shape=(3,)
        Translational part of the symmetry operation, in terms of the basis 
        vectors of the unit cell.
    spinor : bool
        `True` if wave functions are spinors, `False` if they are scalars.
    bands : slice, list, array, shape=(N,), optional
        Selects which bands (m,n) will be used in the calculation.
    TRS : bool, default = False
        If True, adds the complex conjugation to the symmetry operation.
    
    Returns
    -------
    array
        Matrix of the symmetry operation in the basis of eigenstates of the 
        Bloch Hamiltonian :math:`H(k)`.
    """
    npw = igall.shape[1]
    TRS_sign = (-1)**TRS # 1 if False, -1 if True
    multZ = np.exp(-1.0j * (2 * np.pi * A.dot(T).dot(igall[:3, :] + K[:, None])))
    igrot = transformed_g(K, igall, RecLattice, TRS_sign*A)
    WF = WF[bands] # local copy if sliced, original (copy by reference) if not sliced
    if spinor:
        WF1 = np.stack([WF[:, igrot], WF[:, igrot + npw]], axis=2).conj()
        WF2 = np.stack([WF[:, :npw], WF[:, npw:]], axis=2)
        if TRS:
            WF1 = WF1.conj()
        return np.einsum("mgs,ngt,g,st->mn", WF1, WF2, multZ, S)
    else:
        if TRS:
            return np.einsum("mg,ng,g->mn", WF[:, igrot], WF, multZ)
        else:
            return np.einsum("mg,ng,g->mn", WF[:, igrot].conj(), WF, multZ)


def sigma_matrix(WF, igall, spinor, sigma, bands=slice(None)):
    '''
    Calculates the matrix elements of the momentum operator <m|sigma_nu|n> along
    direction nu=xyz.

    This routine is based on the symm_matrix routine from the irrep package (gvectors.py).

    Parameters
    ----------
    WF : ndarray
        An array containing the coefficients for the plane-wave expansion of the wave-function
        in the Bloch basis. `WF[i,j]` contains the coefficient corresponding to the j-th plane-wave
        in the expansion of the wave-function in the i-th band. It contains only plane-waves if energy
        is smaller than `Ecut`.
    igall : ndarray
        Array returned by `__sortIG`. Every column corresponds to a plane-wave of energy smaller than 
        `Ecut`. The number of rows is 6: the first 3 contain direct coordinates of the plane-wave, the
        fourth row stores indices needed to sort plane-waves based on energy (ascending order). The
        fifth (sixth) row contains the index of the first (last) plane-wave with the same energy as the
        plane-wave of the current column.
    spinor : bool
        A flag indicating if the wave functions are spinors (True) or scalars (False).
    sigma : ndarray
        A Pauli matrix sigma_x/y/z.
    bands : slice, list, ndarray, optional
        A slice, list or array indicating which bands (m,n) will be used in the calculation.

    Returns
    -------
    sigma_mn : ndarray
        The matrix of the momentum operator in the basis of eigenstates of the Bloch Hamiltonian H(k).
        The matrix elements of the momentum operator are given by <m|sigma_nu|n> along the direction nu=xyz.
    '''
    npw = igall.shape[1]
    WF = WF[bands] # local copy if sliced, original (copy by reference) if not sliced
    if spinor:
        WF = np.stack([WF[:, :npw], WF[:, npw:]], axis=2)
        return np.einsum("mgs,st,ngt->mn", WF.conj(), sigma, WF)
    else:
        raise Exception('Sigma matrix elements not defined for spinles case.')