'''
Auxiliary methods used along the code.

Module: **pydft2kp/util.py**
'''

###############################################
# CONVERT REPRESENTATION MATRICES
###############################################
def R_to_spin(R):
    """
    Calculates the spin inv/rotation matrix from the O(3) inv/rotation matrix in r-space.

    Parameters
    ----------
    R : array (3,3)
        The rotation matrix in r-space

    Returns
    -------
    S : array (3,3)
        The rotation matrix in spin 1/2 space
    """

    from scipy.linalg import null_space, expm
    from numpy import array, eye, trace, arccos, cos, sin, outer, allclose
    from numpy.linalg import det
    # Pauli matrices
    sx = array([[0, 1],[1, 0]])
    sy = array([[0, -1j],[1j, 0]])
    sz = array([[1, 0],[0, -1]])

    # detR = +/- 1
    hasinv = det(R)
    R = hasinv * R # make sure R is a rotation, no inversion
    # calculate u and theta
    i3 = eye(3)
    u = null_space(R - i3)[:,0]
    trR = trace(R)
    theta = arccos((trR-1)/2)
    # cross product matrix
    uc = array([[0, -u[2], u[1]], 
                   [u[2], 0, -u[0]], 
                   [-u[1], u[0], 0]])
    # reconstruct R and check sign of theta
    checkR = cos(theta)*i3 + sin(theta)*uc + (1-cos(theta))*outer(u, u)
    if allclose(checkR, R) == False:
        theta = -theta
    # build spin matrix
    su = u[0]*sx + u[1]*sy + u[2]*sz
    S = hasinv * expm(-1j*(theta/2)*su)
    return S

def R_to_bvec(R, bs):
    """
    Calculates the rotation matrix in the basis of the reciprocal vectors.

    Parameters
    ----------
    R : array (3,3)
        Rotation matrix in r-space
    bs: array (3,3)
        Reciprocal vectors, one per line

    Returns
    -------
    RotB : array (3,3)
        Rotation matrix in reciprocal space
    """
    # import locally for compatibility
    from numpy import array
    from numpy.linalg import inv
    invB = inv(bs.T)
    RotB = array([invB @ R @ bs[i] for i in range(3)])
    return RotB


###############################################
# CONVERT UNITS FROM AU TO meV AND nm
###############################################
def convert_units_coeffs(coeffs, keys, verbose=True, sigdigits=5):
    """
    Convert coefficient values from a.u. to eV and nm.

    Parameters
    ----------
    coeffs : array
        List of coefficients cn from qsymm
    keys : array
        Labels of the k-powers of each coefficient
    verbose : bool, optional
        If true, prints a report, by default True
    sigdigits : int, optional
        Number of significant digits printed in the report, by default 5

    Returns
    -------
    array
        Coefficients with units converted to eV and nm.
    """
    # import locally for compatibility
    from .constants import Ry, a0
    from numpy import round as npround # to avoid confusion with round -> int
    from tabulate import tabulate

    # all must have energy units
    coeffs_with_units = coeffs * Ry
    # if term of k^n, multiply by a0^n
    for n in range(len(coeffs)):
        kpow = 0
        if keys[n] != [0]:
            kpow = len(keys[n][0])
            coeffs_with_units[n] *= a0**kpow
    
    if verbose:
        cs = []
        units = []
        strkeys = []
        for n in range(len(coeffs)):
            kpow = 0
            if keys[n] != [0]:
                kpow = len(keys[n][0])
            # coeffs
            cs += ['c' + str(n)]
            # units
            unit = 'eV'
            if kpow == 1:
                unit += '.nm'
            elif kpow == 2:
                unit += '.nm²'
            elif kpow > 2:
                unit += '.nm^'+str(kpow)
            units += [unit]
            # k powers
            strkey = ''
            for key in keys[n]:
                strkey += str(key)+','
            strkeys += [strkey[:-1]]
        # organize data to print
        toprint = [cs, 
                   [f'{coeffs[n]:.{sigdigits}}' for n in range(len(coeffs))],
                   [f'{coeffs_with_units[n]:.{sigdigits}}' for n in range(len(coeffs))],
                   units, strkeys]
        # tranpose
        toprint = list(map(list, zip(*toprint)))
        head = ['cn', 'a.u. (Ry, a0)', 'with (eV, nm)', 'units', 'k powers']
        print(tabulate(toprint, headers=head))

    if not verbose:
        return coeffs_with_units   


###############################################
# CONVERT STRING TO ANGLE DEGREES
###############################################
def convertpi(thestr):
    """
    Uses sympy to convert a string (e.g. 2pi/3) into a numerical angle in degrees (e.g. 120).

    Parameters
    ----------
    thestr : str
        A string representing an angle in terms of pi (e.g. 2pi/3)

    Returns
    -------
    float
        The angle in degrees
    """
    # import locally for compatibility
    from sympy import symbols
    from sympy.parsing.sympy_parser import standard_transformations, parse_expr, implicit_multiplication
    pi = symbols("π")
    transformation = standard_transformations+(implicit_multiplication,)
    return round(parse_expr(thestr,transformations=transformation).subs(pi,180))


###############################################
# DEFINE A CONTEXT TO LOCALY CHANGE DIR
###############################################
import os
from contextlib import contextmanager

#based on: https://stackoverflow.com/questions/299446/
@contextmanager
def cwd(path):
    """
    Defines a context manager to locally change the current working directory.

    Parameters
    ----------
    path : str
        The path to the new working directory
    """
    
    try:
        oldpwd = os.getcwd()
        os.chdir(path)
        yield
    finally:
        os.chdir(oldpwd)
