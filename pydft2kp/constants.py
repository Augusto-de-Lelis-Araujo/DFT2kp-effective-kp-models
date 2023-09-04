'''
This module lists constants, matrices and dict keys of commom use along the code.

Module: **pydft2kp/constants.py**
'''

from numpy import array

# global constants
Ry = (13.6056980659) #: Rydberg energy in [eV]
a0 = 0.0529177249 #: Bohr radius in [nm]
hbar = 0.658211951e-3 #: Reduced Planck's constant in [eV.ps]
m0 = (hbar**2)/(2*Ry*a0**2) #: Bare electron mass in [eV.ps²/nm²]
alpha = 1/137 #: Fine-structure constant
h2m = Ry*a0**2 #: Auxiliary :math:`\hbar^2/2m` in [eV.nm²]

# Pauli matrices
s0 = array([[1,0],[0,1]])    #: Pauli matrix, identity
sx = array([[0,1],[1,0]])    #: Pauli matrix, :math:`\sigma_x`
sy = array([[0,-1j],[1j,0]]) #: Pauli matrix, :math:`\sigma_y`
sz = array([[1,0],[0,-1]])   #: Pauli matrix, :math:`\sigma_z`

# DICTIONARY KEYS TO READ THE DFT AND QSYMM DATA
# DFT   data: 0, 'x', 'y', 'z', 'xx', 'xy', ...
# QSYMM data: 1, 'k_x', 'k_y', 'k_z', 'k_x**2', 'kx*k_z', ...
#                 kx**2 = kx*kx (sympy simplifies it)

strKs = ['k_x', 'k_y', 'k_z']
#: List of strings labeling the dict keys to read the QSymm families up to order 3.
QSkeys = []
QSkeys += [1] # order 0
QSkeys += strKs # order 1
QSkeys += [strKs[i]+'*'+strKs[j] for i in range(3) for j in range(i,3)] # order 2
QSkeys += [strKs[i]+'*'+strKs[j]+'*'+strKs[l] for i in range(3) for j in range(i,3) for l in range(j,3)] # order 3

strKs = ['x', 'y', 'z']
#: List of strings labeling the dict keys to read the QE data up to order 3.
DFTkeys = []
DFTkeys += [0]
DFTkeys += strKs
DFTkeys += [strKs[i]+strKs[j] for i in range(3) for j in range(i,3)]
DFTkeys += [strKs[i]+strKs[j]+strKs[l] for i in range(3) for j in range(i,3) for l in range(j,3)]