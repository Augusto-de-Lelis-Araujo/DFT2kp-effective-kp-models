'''
Doc for the **pydft2kp/lowdin.py** module.
'''

import numpy as np

def getHpowers(irrep, NB=None, maxorder=2):
    """
    Returns a dictionary representing H terms that multiply powers of the momentum.

    Parameters
    ----------
    irrep : irrep object
        DFT data read by the irrep package
    NB : int or None, optional
        Number of bands to consider the set B above set A. Defaults to None.
    maxorder : int, optional
        Maximum power of momentum. Defaults to 2.

    Returns
    -------
    dict
        Dictionary with the matrices that multiply the powers of momentum.

    Examples
    --------
    >>> # get the H matrices that multiply k_x^2, k_y^2 and k_z^2, up to second order
    >>> getHpowers(kp, maxorder=2)
    {'xx': array([...]), 'yy': array([...]), 'zz': array([...])}

    Notes
    -----
    
    The key 'xx' refers to the matrix :math:`H_{xx}` that contributues to the
    total Hamiltonian with :math:`H_{xx} k_x^2`. Similarly, 'xy' refers to the
    :math:`H_{xy}` term in the contribution :math:`H_{xy} k_x k_y`, and so on.
    """
    H0 = np.diag(irrep.energies)
    # factor 2 below due to H = H0 + 2k.p + k²
    return lowdin(irrep.setA, H0, 2*irrep.px, 2*irrep.py, 2*irrep.pz, NB, maxorder)



# Build Hamiltonians
def H_of_k(irrep, Hpow=None):
    '''
    Builds and return a function that calculates H(k).

    Using Rydberg units (:math:`\hbar` = 1, m = 1/2). 
    The kp Hamiltonian is :math:`H = H_0 + 2 k \cdot p + k^2`,
    where :math:`H_0` contains the QE eigenstates.

    Parameters
    ----------
    irrep: irrep object
        DFT data read by the irrep package
    Hpow : dict, optional
        Dictionary with the matrices that multiply the powers of momentum.
        Generated by getHpowers(...). Defaults to None.

    Returns
    -------
    H(kx, ky, kz, [maxorder=2]) : callable 
        A function that calculates the Hamiltonian H(k) up to a maximum order
        in the Löwdin expansion if Hpow is provided. Otherwise, it returns
        the full Hamiltonian H.

    Notes
    -----
    If Hpow is None, the function calculates the full Hamiltonian H, which is
    a sum of the kinetic energy, H0, the dot product between momentum and the
    position operator, 2k.p, and the momentum squared, k².

    If Hpow is not None, the function calculates the Hamiltonian H(k) using a
    Löwdin expansion truncated up to the maximum order specified by `maxorder`.
    The matrices that multiply the powers of momentum are provided in the
    `Hpow` dictionary.
    '''

    if Hpow is None: # returns full H
        H0 = np.diag(irrep.energies)
        I0 = np.eye(len(irrep.energies))
        def H(kx=0, ky=0, kz=0):
            k2 = kx**2 + ky**2 + kz**2
            # building H = H0 + 2k.p + k²
            return (H0 + k2*I0 + 2*(irrep.px*kx + irrep.py*ky + irrep.pz*kz))
        return H
    else: # returns H reduced to set A
        def H(kx=0, ky=0, kz=0, maxorder=2):
            '''
            Returns H(kx,ky,kz) up to maxorder in the Löwdin expansion.
            '''
            # order 0
            h  = Hpow[0] + 0j # trick: add 0j to make sure h is complex
            # order 1
            if maxorder >= 1:
                h += Hpow['x']*kx + Hpow['y']*ky + Hpow['z']*kz
            # order 2
            if maxorder >= 2:
                h += Hpow['xx']*kx**2 + Hpow['xy']*kx*ky + Hpow['xz']*kx*kz
                h += Hpow['yy']*ky**2 + Hpow['yz']*ky*kz + Hpow['zz']*kz*kz
            # order 3
            if maxorder >= 3:
                h += Hpow['xxx']*kx**3 + Hpow['yyy']*ky**3 + Hpow['zzz']*kz**3
                h += Hpow['xxy']*kx*kx*ky + Hpow['xxz']*kx*kx*kz + Hpow['xyy']*kx*ky*ky
                h += Hpow['xyz']*kx*ky*kz + Hpow['xzz']*kx*kz*kz + Hpow['yyz']*ky*ky*kz
                h += Hpow['yzz']*ky*kz*kz
            if maxorder >= 4:
                h += Hpow['xxxx']*kx**4
                h += Hpow['yyyy']*ky**4
                h += Hpow['zzzz']*kz**4
                h += Hpow['xxxy']*kx**3*ky
                h += Hpow['xxxz']*kx**3*kz
                h += Hpow['xyyy']*ky**3*kx
                h += Hpow['yyyz']*ky**3*kz
                h += Hpow['xzzz']*kz**3*kx
                h += Hpow['yzzz']*kz**3*ky
                h += Hpow['xxyy']*kx**2*ky**2
                h += Hpow['xxzz']*kx**2*kz**2
                h += Hpow['yyzz']*ky**2*kz**2
                h += Hpow['xxyz']*kx**2*ky*kz
                h += Hpow['xyyz']*ky**2*kx*kz
                h += Hpow['xyzz']*kz**2*kx*ky
            return h
        return H

# the order 2 term in Löwdin
def order2(a, a1, B, e0, H1, H2):
    """
    Calculates the second-order term in the Löwdin expansion.

    Parameters
    ----------
    a : int
        Index of the bra.
    a1 : int
        Index of the ket.
    B : list, array
        List of indices for the virtual states.
    e0 : array
        Eigenvalues.
    H1 : array
        Matrix for first element
    H2 : array
        Matrix for second element

    Returns
    -------
    o2 : float
        The calculated second-order Löwdin term.
    """

    return np.sum([0.5*H1[a,b]*H2[b,a1]*(1/(e0[a]-e0[b]) + 1/(e0[a1]-e0[b])) for b in B])

# the order 3 term in Löwdin
def order3(a, a1, A, B, e0, H1, H2, H3):
    """
    Calculates the third order term in the Löwdin expansion.

    Parameters
    ----------
    a : int
        Index of the bra.
    a1 : int
        Index of the ket.
    A : list
        List of the band indices in set A.
    B : list
        List of the band indices in set B.
    e0 : array
        Eigenenergies
    H1 : array
        Matrix for first element
    H2 : array
        Matrix for second element
    H3 : array
        Matrix for third element

    Returns
    -------
    o3 : float
        Third order term in the Löwdin expansion.
    """

    o3  = np.sum([-0.5*H1[a,b] *H2[b,a2]*H3[a2,a1]/((e0[a1]-e0[b])*(e0[a2]-e0[b] )) for b in B for a2 in A])
    o3 += np.sum([-0.5*H1[a,a2]*H2[a2,b]*H3[b,a1] /((e0[a] -e0[b])*(e0[a2]-e0[b] )) for b in B for a2 in A])
    o3 += np.sum([+0.5*H1[a,b] *H2[b,b1]*H3[b1,a1]/((e0[a] -e0[b])*(e0[a] -e0[b1])) for b in B for b1 in B])
    o3 += np.sum([+0.5*H1[a,b] *H2[b,b1]*H3[b1,a1]/((e0[a1]-e0[b])*(e0[a1]-e0[b1])) for b in B for b1 in B])
    return o3

# the order 4 term in Löwdin
# VERY SLOW, REQUIRES MULTIPROCESSING: DO NOT USE!!!
def order4(a, a1, A, B, e0, H1, H2, H3, H4):
    """
    Calculates the order 4 term in the Löwdin expansion.
    
    Parameters
    ----------
    a : int
        Index of the bra
    a1 : int
        Index of the ket
    A : list
        List with the indices of the bands that belong to set A.
    B : list
        List with the indices of the bands that belong to set B.
    e0 : array
        Eigenvalues
    H1 : array
        Matrix for first element
    H2 : array
        Matrix for second element
    H3 : array
        Matrix for third element
    H4 : array
        Matrix for fourth element
    
    Returns
    -------
    o4 : float
        The value of the order 4 term in the Löwdin expansion.
        
    Notes
    -----
    
    This function is very slow and requires multiprocessing.
    """
    d = lambda i,j: 1/(e0[i]-e0[j])
    o4  = np.sum([ +0.5 *H1[a,a2] *H2[a2,a3] *H3[a3,b]  *H4[b,a1]  *d(a2,b) *d(a3,b)  *d(a,b)   for b in B for a2 in A for a3 in A ])
    o4 += np.sum([ +0.5 *H1[a,b]  *H2[b,a2]  *H3[a2,a3] *H4[a3,a1] *d(a2,b) *d(a3,b)  *d(a1,b)  for b in B for a2 in A for a3 in A ])
    #------------------
    o4 += np.sum([ -0.5 *H1[a,b]  *H2[b,b1]  *H3[b1,a2] *H4[a2,a1] *d(a1,b) *d(a2,b1) *d(a2,b)  for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -0.5 *H1[a,b]  *H2[b,b1]  *H3[b1,a2] *H4[a2,a1] *d(a1,b) *d(a2,b1) *d(a1,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -0.5 *H1[a,a2] *H2[a2,b]  *H3[b,b1]  *H4[b1,a1] *d(a,b1) *d(a2,b)  *d(a2,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -0.5 *H1[a,a2] *H2[a2,b]  *H3[b,b1]  *H4[b1,a1] *d(a,b1) *d(a2,b)  *d(a,b)   for b in B for b1 in B for a2 in A ])
    #------------------
    o4 += np.sum([ -(8/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a,b)  *d(a,b1)  *d(a2,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -(8/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a1,b) *d(a1,b1) *d(a2,b)  for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -(4/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a,b1) *d(a2,b)  *d(a,b)   for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -(4/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a,b1) *d(a2,b)  *d(a2,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -(4/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a1,b) *d(a2,b1) *d(a1,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ -(4/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a1,b) *d(a2,b1) *d(a2,b)  for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ +(1/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a2,b) *d(a2,b1) *d(a,b)   for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ +(1/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a2,b) *d(a2,b1) *d(a1,b1) for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ +(3/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a,b)  *d(a1,b1) *d(a2,b)  for b in B for b1 in B for a2 in A ])
    o4 += np.sum([ +(3/24) *H1[a,b]  *H2[b,a2]  *H3[a2,b1] *H4[b1,a1] *d(a,b)  *d(a1,b1) *d(a2,b1) for b in B for b1 in B for a2 in A ])
    #------------------
    o4 += np.sum([ +0.5 *H1[a,b] *H2[b,b1] *H3[b1,b2] *H4[b2,a1] *d(a,b)  *d(a,b1)  *d(a,b2)  for b in B for b1 in B for b2 in B ])
    o4 += np.sum([ +0.5 *H1[a,b] *H2[b,b1] *H3[b1,b2] *H4[b2,a1] *d(a1,b) *d(a1,b1) *d(a1,b2) for b in B for b1 in B for b2 in B ])
    #--------
    return o4

def lowdin(A, H0, Hx, Hy, Hz, NB=None, maxorder=2):
    '''
    Folds down H into the selected setA
    
    Parameters
    ----------
        A: int list/array
            set of states considered as Löwdin's set A
        H0: NxN array
            Diagonal entries.
        Hx, Hy, Hz: NxN array
            Perturbation terms proportional to kx, ky, kz
        NB: int
            Number of bands in set B, above A
        maxorder: int
            Calculate the expansion up to this order
    
    Returns
    -------
        h : dict
            folded down h = h0 + hx.kx + hy.ky + hz.kz + hxx.kx² + hxy.kx.ky + ...
    '''
    # get size of Hs (next version: use asserts)
    N = len(H0)
    NA = len(A)
    
    # define set B
    fullset = np.arange(N)
    A = np.array(A)
    if NB is None:
        B = np.array([F for F in fullset if F not in A])
    else:
        firstA = np.where(fullset == A[0])[0][0]
        lastA = np.where(fullset == A[-1])[0][0]
        B = fullset[0:firstA] # all until firtA
        B = np.append(B, fullset[(lastA+1):(lastA+1+NB)]) # NB above lastA
    
    # eigenstates, assuming E0 from H0 diagonal
    e0 = np.diag(H0)
    
    # order 0
    h0 = H0[A,:][:,A]
    
    # order 1
    hx = np.zeros([NA,NA], dtype=complex)
    hy = np.zeros([NA,NA], dtype=complex)
    hz = np.zeros([NA,NA], dtype=complex)
    if maxorder >= 1:
        hx = Hx[A,:][:,A]
        hy = Hy[A,:][:,A]
        hz = Hz[A,:][:,A]
    
    # order 2
    hxx = np.zeros([NA,NA], dtype=complex)
    hxy = np.zeros([NA,NA], dtype=complex)
    hxz = np.zeros([NA,NA], dtype=complex)
    hyy = np.zeros([NA,NA], dtype=complex)
    hyz = np.zeros([NA,NA], dtype=complex)
    hzz = np.zeros([NA,NA], dtype=complex)
    if maxorder >= 2:
        for ia in range(NA):
            for ja in range(NA):
                a = A[ia]
                a1 = A[ja]
                
                # (X+Y+Z)² = X² + Y² + Z² + (X.Y + Y.X) + (X.Z + Z.X) + (Y.Z + Z.Y)
                
                hxx[ia,ja]  = order2(a, a1, B, e0, Hx, Hx)
                hyy[ia,ja]  = order2(a, a1, B, e0, Hy, Hy)
                hzz[ia,ja]  = order2(a, a1, B, e0, Hz, Hz)
                
                hxy[ia,ja]  = order2(a, a1, B, e0, Hx, Hy)
                hxy[ia,ja] += order2(a, a1, B, e0, Hy, Hx)
                
                hxz[ia,ja]  = order2(a, a1, B, e0, Hx, Hz)
                hxz[ia,ja] += order2(a, a1, B, e0, Hz, Hx)
                
                hyz[ia,ja]  = order2(a, a1, B, e0, Hy, Hz)
                hyz[ia,ja] += order2(a, a1, B, e0, Hz, Hy)
   
    # order 3
    hxxx = np.zeros([NA,NA], dtype=complex)
    hyyy = np.zeros([NA,NA], dtype=complex)
    hzzz = np.zeros([NA,NA], dtype=complex)
    hxxy = np.zeros([NA,NA], dtype=complex)
    hxxz = np.zeros([NA,NA], dtype=complex)
    hxyy = np.zeros([NA,NA], dtype=complex)
    hxyz = np.zeros([NA,NA], dtype=complex)
    hxzz = np.zeros([NA,NA], dtype=complex)
    hyyz = np.zeros([NA,NA], dtype=complex)
    hyzz = np.zeros([NA,NA], dtype=complex)
    if maxorder >= 3:
        for ia in range(NA):
            for ja in range(NA):
                a = A[ia]
                a1 = A[ja]
                
                # (X+Y+Z)³ = X³ + Y³ + Z³
                #          + XXY + XYX + YXX 
                #          + XXZ + XZX + ZXX 
                #          + XYY + YXY + YYX 
                #          + XYZ + XZY + YXZ + YZX + ZXY + ZYX 
                #          + XZZ + ZXZ + ZZX
                #          + YYZ + YZY + ZYY
                #          + YZZ + ZYZ + ZZY
                
                hxxx[ia,ja]  = order3(a, a1, A, B, e0, Hx, Hx, Hx)
                hyyy[ia,ja]  = order3(a, a1, A, B, e0, Hy, Hy, Hy)
                hzzz[ia,ja]  = order3(a, a1, A, B, e0, Hz, Hz, Hz)
                
                hxxy[ia,ja]  = order3(a, a1, A, B, e0, Hx, Hx, Hy)
                hxxy[ia,ja] += order3(a, a1, A, B, e0, Hx, Hy, Hx)
                hxxy[ia,ja] += order3(a, a1, A, B, e0, Hy, Hx, Hx)
                
                hxxz[ia,ja]  = order3(a, a1, A, B, e0, Hx, Hx, Hz)
                hxxz[ia,ja] += order3(a, a1, A, B, e0, Hx, Hz, Hx)
                hxxz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hx, Hx)
                
                hxyy[ia,ja]  = order3(a, a1, A, B, e0, Hy, Hy, Hx)
                hxyy[ia,ja] += order3(a, a1, A, B, e0, Hy, Hx, Hy)
                hxyy[ia,ja] += order3(a, a1, A, B, e0, Hx, Hy, Hy)
                
                hxyz[ia,ja]  = order3(a, a1, A, B, e0, Hx, Hy, Hz)
                hxyz[ia,ja] += order3(a, a1, A, B, e0, Hx, Hz, Hy)
                hxyz[ia,ja] += order3(a, a1, A, B, e0, Hy, Hx, Hz)
                hxyz[ia,ja] += order3(a, a1, A, B, e0, Hy, Hz, Hx)
                hxyz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hx, Hy)
                hxyz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hy, Hx)
                
                hxzz[ia,ja]  = order3(a, a1, A, B, e0, Hz, Hz, Hx)
                hxzz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hx, Hz)
                hxzz[ia,ja] += order3(a, a1, A, B, e0, Hx, Hz, Hz)
                
                hyyz[ia,ja]  = order3(a, a1, A, B, e0, Hy, Hy, Hz)
                hyyz[ia,ja] += order3(a, a1, A, B, e0, Hy, Hz, Hy)
                hyyz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hy, Hy)
                
                hyzz[ia,ja]  = order3(a, a1, A, B, e0, Hz, Hz, Hy)
                hyzz[ia,ja] += order3(a, a1, A, B, e0, Hz, Hy, Hz)
                hyzz[ia,ja] += order3(a, a1, A, B, e0, Hy, Hz, Hz)

    # order 4   
    hxxxx = np.zeros([NA,NA], dtype=complex)
    hyyyy = np.zeros([NA,NA], dtype=complex)
    hzzzz = np.zeros([NA,NA], dtype=complex)
    hxxxy = np.zeros([NA,NA], dtype=complex)
    hxxxz = np.zeros([NA,NA], dtype=complex)
    hxyyy = np.zeros([NA,NA], dtype=complex)
    hyyyz = np.zeros([NA,NA], dtype=complex)
    hxzzz = np.zeros([NA,NA], dtype=complex)
    hyzzz = np.zeros([NA,NA], dtype=complex)
    hxxyy = np.zeros([NA,NA], dtype=complex)
    hxxzz = np.zeros([NA,NA], dtype=complex)
    hyyzz = np.zeros([NA,NA], dtype=complex)
    hxxyz = np.zeros([NA,NA], dtype=complex)
    hxyyz = np.zeros([NA,NA], dtype=complex)
    hxyzz = np.zeros([NA,NA], dtype=complex)
    if maxorder >= 4:
        for ia in range(NA):
            for ja in range(NA):
                a = A[ia]
                a1 = A[ja]        

                # (X+Y+Z)^4 = x^4 + y^4 + z^4 
                #           + 4x^3y + 4x^3z + 4y^3x + 4y^3z + 4z^3x + 4z^3y 
                #           + 6x^2y^2 + 6x^2z^2 + 6y^2z^2 
                #           + 12x^2yz + 12y^2xz + 12z^2xy

                # XXXX | YYYY | ZZZZ
                hxxxx[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hx, Hx)
                hyyyy[ia,ja]  = order4(a, a1, A, B, e0, Hy, Hy, Hy, Hy)
                hzzzz[ia,ja]  = order4(a, a1, A, B, e0, Hz, Hz, Hz, Hz)

                # XXXY + XXYX + XYXX + YXXX
                hxxxy[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hx, Hy)
                hxxxy[ia,ja] += order4(a, a1, A, B, e0, Hx, Hx, Hy, Hx)
                hxxxy[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hx, Hx)
                hxxxy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hx, Hx)
                
                # XXXZ + XXZX + XZXX + ZXXX 
                hxxxz[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hx, Hz)
                hxxxz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hx, Hz, Hx)
                hxxxz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hx, Hx)
                hxxxz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hx, Hx)
                
                # YYYX + YYXY + YXYY + XYYY
                hxyyy[ia,ja]  = order4(a, a1, A, B, e0, Hy, Hy, Hy, Hx)
                hxyyy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hy, Hx, Hy)
                hxyyy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hy, Hy)
                hxyyy[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hy, Hy)
                
                # YYYZ + YYZY + YZYY + ZYYY
                hyyyz[ia,ja]  = order4(a, a1, A, B, e0, Hy, Hy, Hy, Hz)
                hyyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hy, Hz, Hy)
                hyyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hy, Hy)
                hyyyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hy, Hy)
                
                # ZZZX + ZZXZ + ZXZZ + XZZZ
                hxzzz[ia,ja]  = order4(a, a1, A, B, e0, Hz, Hz, Hz, Hx)
                hxzzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hz, Hx, Hz)
                hxzzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hz, Hz)
                hxzzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hz, Hz)
                
                # ZZZY + ZZYZ + ZYZZ + YZZZ
                hyzzz[ia,ja]  = order4(a, a1, A, B, e0, Hz, Hz, Hz, Hy)
                hyzzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hz, Hy, Hz)
                hyzzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hz, Hz)
                hyzzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hz, Hz)
                
                # XXYY + XYYX + YYXX + YXXY + XYXY + YXYX
                hxxyy[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hy, Hy)
                hxxyy[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hy, Hx)
                hxxyy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hy, Hx, Hx)
                hxxyy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hx, Hy)
                hxxyy[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hx, Hy)
                hxxyy[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hy, Hx)
                
                # XXZZ + XZZX + ZZXX + ZXXZ + XZXZ + ZXZX
                hxxzz[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hz, Hz)
                hxxzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hz, Hx)
                hxxzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hz, Hx, Hx)
                hxxzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hx, Hz)
                hxxzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hx, Hz)
                hxxzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hz, Hx)
                
                # YYZZ + YZZY + ZZYY + ZYYZ + YZYZ + ZYZY
                hyyzz[ia,ja]  = order4(a, a1, A, B, e0, Hy, Hy, Hz, Hz)
                hyyzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hz, Hy)
                hyyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hz, Hy, Hy)
                hyyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hy, Hz)
                hyyzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hy, Hz)
                hyyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hz, Hy)
                
                # (XXYZ + XYZX + YZXX + YXXZ + XYXZ + YXZX) + (XXZY + XZYX + ZYXX + ZXXY + XZXY + ZXYX)
                hxxyz[ia,ja]  = order4(a, a1, A, B, e0, Hx, Hx, Hy, Hz)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hz, Hx)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hx, Hx)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hx, Hz)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hx, Hz)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hz, Hx)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hx, Hz, Hy)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hy, Hx)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hx, Hx)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hx, Hy)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hx, Hy)
                hxxyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hy, Hx)
                
                # (YYXZ + YXZY + XZYY + XYYZ + YXYZ + XYZY) + (YYZX + YZXY + ZXYY + ZYYX + YZYX + ZYXY)
                hxyyz[ia,ja]  = order4(a, a1, A, B, e0, Hy, Hy, Hx, Hz)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hz, Hy)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hy, Hy)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hy, Hz)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hy, Hz)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hz, Hy)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hy, Hz, Hx)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hx, Hy)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hy, Hy)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hy, Hx)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hy, Hx)
                hxyyz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hx, Hy)
                
                # (ZZXY + ZXYZ + XYZZ + XZZY + XZYZ + ZXZY) + (ZZYX + ZYXZ + YXZZ + YZZX + YZXZ + ZYZX)
                hxyzz[ia,ja]  = order4(a, a1, A, B, e0, Hz, Hz, Hx, Hy)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hy, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hy, Hz, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hz, Hy)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hx, Hz, Hy, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hx, Hz, Hy)                        
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hz, Hy, Hx)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hx, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hx, Hz, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hz, Hx)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hy, Hz, Hx, Hz)
                hxyzz[ia,ja] += order4(a, a1, A, B, e0, Hz, Hy, Hz, Hx)
            
    # return as a dictionary
    h = {}
    # order 0
    h[0] = h0
    # order 1
    h['x'] = hx
    h['y'] = hy
    h['z'] = hz
    # order 2
    h['xx'] = hxx + np.eye(hxx.shape[0])
    h['xy'] = hxy
    h['xz'] = hxz
    h['yy'] = hyy + np.eye(hyy.shape[0])
    h['yz'] = hyz
    h['zz'] = hzz + np.eye(hzz.shape[0])
    # order 3
    h['xxx'] = hxxx
    h['yyy'] = hyyy
    h['zzz'] = hzzz
    h['xxy'] = hxxy
    h['xxz'] = hxxz
    h['xyy'] = hxyy
    h['xyz'] = hxyz
    h['xzz'] = hxzz
    h['yyz'] = hyyz
    h['yzz'] = hyzz
    # order 4
    h['xxxx'] = hxxxx
    h['yyyy'] = hyyyy
    h['zzzz'] = hzzzz
    h['xxxy'] = hxxxy
    h['xxxz'] = hxxxz
    h['xyyy'] = hxyyy
    h['yyyz'] = hyyyz
    h['xzzz'] = hxzzz
    h['yzzz'] = hyzzz
    h['xxyy'] = hxxyy
    h['xxzz'] = hxxzz
    h['yyzz'] = hyyzz
    h['xxyz'] = hxxyz
    h['xyyz'] = hxyyz
    h['xyzz'] = hxyzz

    return h
