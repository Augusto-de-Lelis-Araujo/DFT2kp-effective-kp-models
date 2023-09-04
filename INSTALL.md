# Installation

The installation requires two steps. First, one needs to install the python package using *pip*. Second, dowload the QE patch and recompile QE with our modifications.

## Compatibility

The table below lists the compatibility of each version of our DFT2kp code with the [Quantum Espresso](https://gitlab.com/QEF/q-e), [IrRep](https://github.com/stepan-tsirkin/irrep), and [Qsymm](https://github.com/quantum-tinkerer/qsymm) codes.

| DFT2kp version | QE               | IrRep | Qsymm | 
| :------------: | :--------------: | :---: | :---: |
| 0.0.1          |  7.0 / 7.1 / 7.2 | 1.7.1 | 1.3.0 |

## Installing the python package

To install the latest version, run:

```bash
pip install dft2kp
```

To specify a specific realease (e.g. version 0.0.1), run:

```bash
pip install dft2kp==0.0.1
```

### Development version

To install directly from the git repos, run:

```bash
pip install git+https://gitlab.com/dft2kp/dft2kp.git
```

## Patch the QE code

Download the appropriate patch from our [gitlab repos' patch folder](https://gitlab.com/dft2kp/dft2kp), and the desired version of the [QE code](https://gitlab.com/QEF/q-e). 

To apply the patch, go the root folder of the QE code and run:

```bash
git apply path/to/qe2kp-<version>.patch
```

Then proceed with the compilation of QE. A minimal compilation with MPI and/or OPENMP looks like:

```bash
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DQE_ENABLE_OPENMP=[OFF/ON] ..
make -jN
```