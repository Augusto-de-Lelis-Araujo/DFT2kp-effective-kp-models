'''
Doc for the **pydft2kp/qe_aux.py** module.
'''

from numpy import loadtxt, pi, array, where, fromstring, unique, argwhere, copy
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from .constants import Ry, a0

class qe_plotter:
    '''
    A class that reads and plots Quantum Espresso band structures.

    Parameters
    ----------
    kp : kp object
        A `kp` object that contains the information about the DFT calculations.
    bandsgnu : str
        The name of the file that contains the bands data in gnuplot format.
    kpath : str
        The name of the file that contains the k-path data.

    Attributes
    ----------
    en_at_k : ndarray
        An array with the energy levels at each k point.
    kpaths : list of str
        A list of strings with the names of the k paths.
    klabels : ndarray
        An array with the distances along the k-path for each k point.
    kdist : ndarray
        An array with the distances along the k-path.
    bands : ndarray
        A 2D array with the energy levels for each band and each k-point.
    k3D : ndarray
        A 2D array with the Cartesian coordinates of the k-points.

    Notes
    -----
    This class is designed to work with Quantum Espresso DFT calculations that
    provide the band structure data in gnuplot format and k-path data in a
    separate file. It assumes that the band structure data is stored in a file
    named `bandsgnu` in the working directory, and the k-path data is stored in
    a file named `kpath` in the same directory. The `kp` object is used to
    extract other relevant information, such as the lattice parameter, the
    Fermi level, and the number of k-points.
    
    Examples
    --------
    
    >>> from qeplot import qe_plotter
    >>> kpdata = read_kpdata('kp.xml')
    >>> plotter = qe_plotter(kpdata, 'bands.gnu', 'kpath')
    >>> fig, ax = plt.subplots()
    >>> ax.plot(plotter.kdist, plotter.bands[:,0], label='Band 1')
    >>> ax.plot(plotter.kdist, plotter.bands[:,1], label='Band 2')
    >>> plotter.set_labels_and_limits(ax)
    >>> ax.legend()
    >>> plt.show()
    
    '''

    def __init__(self, kp, bandsgnu, kpath):
        # extract info from kp object
        dftdir = kp.dftdir
        outdir = kp.outdir
        prefix = kp.prefix
        alat = kp.alat
        fermi = kp.fermi
        kpt = kp.kindex
        
        # energy list at the k point
        self.en_at_k = kp.energies

        # kpaths and labels for k axis
        # and bands for plotting
        # kpt -> kpt-1 because python starts at 0
        self.kpaths, self.klabels, data = qe_bands_fix_units(dftdir, bandsgnu, kpath, kpt-1, alat, fermi)
        self.kdist = unique(data[:,0]) # the k axis as distances = 1D path
        nk = len(self.kdist)
        ne = data.shape[0]//nk
        self.bands = data[:,1].reshape([ne,nk]).T # the energy data

        # read the 3D version of k axis
        xmlpath = dftdir + '/' + outdir + '/' + prefix + '.xml'
        xmlroot = ET.parse(xmlpath).getroot()
        kpts = xmlroot.find('input').find('k_points_IBZ').findall('k_point')
        self.k3D = []
        for each in kpts:
            self.k3D += [fromstring(each.text, sep=' ')]
        self.k3D = array(self.k3D) * 2*pi/alat # fix units
        self.k3D -= self.k3D[kpt-1] # shift central k to origin
    
    def set_labels_and_limits(self, ax, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set axis labels and limits for a given matplotlib axis.

        Parameters
        ----------
        ax : matplotlib axis
            The axis for which to set the labels and limits.
        xmin : float, optional
            The minimum x-axis limit. If not provided, the minimum k-distance is used.
        xmax : float, optional
            The maximum x-axis limit. If not provided, the maximum k-distance is used.
        ymin : float, optional
            The minimum y-axis limit. If not provided, -1 Ry is used.
        ymax : float, optional
            The maximum y-axis limit. If not provided, +1 Ry is used.

        Notes
        -----
        This method sets the y-axis label to "Energy [Ry]". If xmin, xmax, ymin, or ymax
        are not provided, the method uses default values. It then sets the x-axis limits
        and labels based on the k-paths and labels of the `qe_plotter` object.

        If xmin or xmax fall outside the range of the k-paths, this method will shift
        the labels of the nearest k-point to indicate the limit. If the limit is to the
        left of the k-path, the label will have an arrow pointing left. If the limit is
        to the right, the label will have an arrow pointing right.

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> plotter.set_labels_and_limits(ax, xmin=0, xmax=1, ymin=-0.5, ymax=0.5)
        """
        # set y label
        ax.set_ylabel(R"Energy [Ry]")
        # set limits
        if xmin is None: xmin = self.kdist[0]
        if xmax is None: xmax = self.kdist[-1]
        if ymin is None: ymin = -1
        if ymax is None: ymax = +1
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        klabels = list(copy(self.klabels))
        kpaths = list(copy(self.kpaths))

        # adjust xmin labels if needed
        idx = argwhere(array(kpaths) < xmin)
        if len(idx) > 0:
            idx = idx[-1][0]
            klabels[idx] = klabels[idx] + r'$\leftarrow$'
            kpaths[idx] = xmin  
        # the same for xmax
        idx = argwhere(array(kpaths) > xmax)
        if len(idx) > 0:
            idx = idx[0][0]
            klabels[idx] = r'$\rightarrow$' + klabels[idx]
            kpaths[idx] = xmax
        
        # set x labels
        ax.set_xticks(kpaths, klabels)
        ticklabels = ax.get_xticklabels()
        ticklabels[ 0].set_ha("left")
        ticklabels[-1].set_ha("right")


def qe_bands_fix_units(dftdir, bandsgnu, kpath, kpt, alat, fermi):
    '''
    Reads bands in gnuplot format. The bands data is adjusted to set the Fermi 
    energy at zero, and express the k points in Bohr units (2pi/alat).

    Parameters
    ----------
    dftdir : str
        Directory where the Quantum ESPRESSO (QE) data is stored.
    bandsgnu : str
        Name of the file with the bands in gnuplot format.
    kpath : array, shape (2, 3)
        An array containing the number of points along each path and label of 
        the k points used in QE 'bands' calculation. The first line contains 
        the number of points in each section of the k path, and the second 
        line contains the labels of these k points.
    kpt : int
        The index of the central k point selected for the kp expansion.
    alat : float
        The lattice parameter in Bohr units.
    fermi : float
        The Fermi energy in Ry units.
    
    Returns
    -------
    kpts : array, shape (n,)
        The indexes of the k points along the path.
    klabels : array, shape (n,)
        The label of the k points along the path.
    data : array, shape (n*m, 2)
        The k points and energies from the gnuplot file.

    Examples
    --------
    To understand the paths informed here, consider the following three::
    
        .
        ├── your_python_code.ipynb
        ├── ...
        └──dftdata
            ├─── bands.gnu
            ├── ...
            └─── outdir
                    ├── ...
                    └──── graphene.save
                            ├── ...
                            └──── data-file-schema.xml

    This implies these parameters:
    
    - dftdir = 'dftdata'
    - prefix = 'graphene'
    - outdir = 'outdir'
    - bandsgnu = 'bands.gnu'

    To understand the kpath parameter, consider that the QE bands calculation
    is set with the following K_POINTS namelist::

        K_POINTS crystal_b
        3
        +0.0000000000 0.000000000 0.000000000 30 ! G
        +0.3333333333 0.333333333 0.000000000 40 ! K
        +0.0000000000 0.500000000 0.000000000  1 ! M

    The path from G to K contains 30 points, and from K to M 40 points.
    The last number (1) in the M line is irrelevant. The first line of kpath 
    must inform these numbers: [30, 40, 1]. The second line of kpath must 
    have the labels for these k points. It can be simple strings, as ['G', 
    'K', 'M'], or Latex-formatted strings as [R'$\Gamma$', R'$K$', R'$M$'].
    So, here we get::

        kpath = [[30, 40, 1], [R'$\Gamma$', R'$K$', R'$M$']]
    '''

    bandspath = dftdir + '/' + bandsgnu
    
    # example: kpath = [[30, 30, 30], [R'$\Gamma$', R'$K$', R'$M$']]
    klabels = kpath[1] # example: [R'$\Gamma$', R'$K$', R'$M$']
    kids = [0]
    for klen in kpath[0][:-1]:
        kids.append(kids[-1]+klen) # example: [0, 30, 60]

    # READ DFT BANDS
    data = loadtxt(bandspath)
    data[:,0] -= data[kpt,0]
    data[:,0] *= 2*pi/alat # fix units
    data[:,1] /= Ry
    data[:,1] -= fermi # sets Fermi energy at 0
    # set k path and labels for the plots
    kpts = [data[k,0] for k in kids]

    return kpts, klabels, data

def read_espresso(dftdir, prefix, outdir):
    '''
    Reads QE data: alat, Fermi energy.

    Parameters
    ----------
    dftdir : str
        Directory where the QE data is stored
    prefix : str
        Prefix used in the QE calculation
    outdir : str
        Outdir used in the QE calculation
    
    Returns
    -------
    alat : float
        The lattice parameter in Bohr units.
    fermi : float
        The Fermi energy in Ry units.
    
    Raises
    ------
    FileNotFoundError
        If the specified XML file cannot be found.

    Examples
    --------

    To understand the paths informed here, consider the following three::
    
        .
        ├── your_python_code.ipynb
        ├── ...
        └──dftdata
            ├─── bands.gnu
            ├── ...
            └─── outdir
                    ├── ...
                    └──── graphene.save
                            ├── ...
                            └──── data-file-schema.xml

    This implies these parameters:
    
    - dftdir = 'dftdata'
    - prefix = 'graphene'
    - outdir = 'outdir'
    - bandsgnu = 'bands.gnu'

    Notes
    --------
    This function assumes that the Quantum Espresso data has been generated
    with the `save` option in the input file. The lattice parameter and Fermi
    energy are extracted from the data-file-schema.xml file located in the
    specified output directory.

    The Fermi energy is multiplied by a factor of 2 because Quantum Espresso
    reports it in Hartree units, but the function returns it in Rydberg units.

    If the Fermi energy is not found in the XML file, a warning message is
    printed and a value of 0 is returned instead.
    '''

    xmlpath = dftdir + '/' + outdir + '/' + prefix + '.save/data-file-schema.xml'
    # verify if file exists

    # read alat and Fermi
    mytree = ET.parse(xmlpath)
    myroot = mytree.getroot()
    alat = float(myroot.find('input').find('atomic_structure').attrib['alat'])
    # factor 2 due to Hartree to Rydberg conversion
    try:
        fermi = 2 * float(myroot.find('output').find('band_structure').find('fermi_energy').text)
    except:
        fermi = 0
        print('''
###################################################
# WARNING: Fermi level not found in the XML file. #
#                                                 #
#       Using fermi = 0 instead                   #
#                                                 #
# Are you using an old version of QE?             #
###################################################
''')
    
    return alat, fermi


class read_kp_dat():
    """
    Reads the kp.dat file and builds an object with its properties.

    Parameters
    ----------
    kp_path : str
        Path to the kp.dat file.

    Attributes
    ----------
    nbnd : int
        Number of bands.
    p1 : ndarray
        Matrix of the momentum operator in the x-direction.
    p2 : ndarray
        Matrix of the momentum operator in the y-direction.
    p3 : ndarray
        Matrix of the momentum operator in the z-direction.
    """
    def __init__(self, kp_path):
        lines = open(kp_path, 'r').read().split('\n')
        # read nbnd and nocc from line 0
        self.nbnd = int(lines[0].split(',')[0].split('nbnd=')[1])
        # identify lines that mark labels
        line1 = where(array(lines) == '  1')[0][0]
        line2 = where(array(lines) == '  2')[0][0]
        line3 = where(array(lines) == '  3')[0][0]
        # read p1
        p1 = []
        for line in range(line1+1, line2):
            p1 += list(array(lines[line].split(), dtype='float'))
        p1 = array(p1)
        # read p2
        p2 = []
        for line in range(line2+1, line3):
            p2 += list(array(lines[line].split(), dtype='float'))
        p2 = array(p2)
        # read p3
        p3 = []
        for line in range(line3+1, len(lines)):
            p3 += list(array(lines[line].split(), dtype='float'))
        p3 = array(p3)

        # converts matrix elements into matrices
        # mix real and imag parts
        p1 = p1[0::2] + 1j*p1[1::2]
        p2 = p2[0::2] + 1j*p2[1::2]
        p3 = p3[0::2] + 1j*p3[1::2]
        # reshape into matrix and fix units
        self.p1 = p1.reshape((self.nbnd, self.nbnd)) / (20*a0)
        self.p2 = p2.reshape((self.nbnd, self.nbnd)) / (20*a0)
        self.p3 = p3.reshape((self.nbnd, self.nbnd)) / (20*a0)
