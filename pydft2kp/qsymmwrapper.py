'''
Doc for the **pydft2kp/qsymmwrapper.py** module.
'''

from numpy import fromstring, append
from qsymm.groups import pretty_print_pge
from qsymm import continuum_hamiltonian, display_family
from qsymm import inversion, rotation, mirror, time_reversal, PointGroupElement
from qsymm.model import _commutative_momenta
from qsymm.hamiltonian_generator import hamiltonian_from_family
from .util import convertpi

class qsymm():
    '''
    Uses the qsymm package to build a kp model.
    
    Most of the input parameters of follow from the qsymm routines:
    
    From qsymm.continuum_hamiltonian(...)::
    
    - symmetries: list of PointGroupElement
    - total_power: allowed total powers of the momentum variables
    - all_powers: allowed total powers of the momentum variables
    - dim: number of spatial dimensions
    - momenta: names of momentum variables
    - sparse_linalg: use sparse matrix
    - prettify: sympy rounding
    - num_digits: significant digits to prettify

    From qsymm.display_family(...)::
    
    - summed: return summed H or individual terms
    - coeffs: string/sympy list with coeffs names
    - nsimplify: whether to use sympy.nsimplify on the output or not

    Parameters
    ----------
    print_model : bool
        If the user wants to call qsymm.display_family(...) when generating the model.


    Attributes
    ----------
    symms : list
        The symmetries built with qsymm
    model : list
        List of qsymm/Model objects with the matrices
        that define the kp model
    symm_pg : list
        Poing group (S) component of the
        symmetry operation {S,T} (Seitz notation). See Notes.
    Hmodel : matrix
        The Hamiltonian.

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
    '''
    def __init__(self,
                 # parameters from continuum_hamiltonian(...)
                 symmetries, 
                 total_power, 
                 all_powers=None, 
                 dim=3,
                 momenta=_commutative_momenta,
                 sparse_linalg=False, 
                 prettify=True, 
                 num_digits=10,
                 # extra parameter: print model or not?
                 print_model=True,
                 # parameters from display_family(...)
                 summed=True,
                 coeffs=None,
                 nsimplify=True
                 ):


        self.symms = symmetries
        self.model = continuum_hamiltonian(symmetries,
                                           dim, 
                                           total_power, 
                                           all_powers, 
                                           momenta,
                                           sparse_linalg, 
                                           prettify, 
                                           num_digits)
        
        self.Hmodel = hamiltonian_from_family(self.model)
        
        self.symm_pg = self.identify_symmetries()
        
        if print_model:
            display_family(self.model, summed, coeffs, nsimplify)



    def identify_symmetries(self):
        '''
        Reads the symmetry operations and identifies their
        point group parts. There's no information about
        non-symmorphic translations, since these are implied
        in the U matrices.

        Returns
        -------
        symm_pg : list
            Point group elements [S, theta, axis]. See Notes above.
        '''
        # [(1,I,S,R,M), deg, axis]
        # if conjugate, adds T label
        symm_pg = []
        for n in range(len(self.symms)):
            opstr = pretty_print_pge(self.symms[n])
            TRS = ''
            if self.symms[n].conjugate:
                TRS = 'T'
            if (opstr[0] == '1') or (opstr[0] == 'I'):
                theop = [TRS+opstr[0], 0, [0,0,1.]]
            elif (opstr[0] == 'R') or (opstr[0] == 'S'):
                theop = [TRS+opstr[0]] # label R or S
                if ',' in opstr: # 3D case as R(angle, axis)
                    opstr = opstr[2:-1].split(',') # angle,axis
                    angle = convertpi(opstr[0])
                    theax = opstr[1].strip()
                    theax = fromstring(theax[1:-1], sep=' ')                        
                    theop += [angle, list(theax)]
                else: # 2D case as R(angle), while axis is assumed [0,0,1]
                    angle = convertpi(opstr[2:-1])
                    theop += [angle, [0,0,1]]
            else: # M(axis)
                theax = opstr[2:-1].strip()
                theax = fromstring(theax[1:-1], sep=' ')
                # add zero if 2D
                if len(theax) == 2: theax = append(theax, 0)
                theop = [TRS+opstr[0], 180, list(theax)]
            # store
            symm_pg += [theop]
        # return list
        return symm_pg
    
