# Citing DFT2kp

The main paper describing the code is

> J. V. V. Cassiano, A. L. Araújo, P. E. Faria Junior, G. J. Ferreira. "**DFT2kp: effective kp models from ab-initio data**" [arXiv]

## Other relevant references to cite

If you use our code you will certanly cite QE and the other relevant DFT references. But additionally, our code relies on the python packages:

- [IrRep](https://github.com/stepan-tsirkin/irrep): to identify the space group, and calculate the symmetry representation matrices. Additionally, we use modified versions of routines from this package to calculate the matrix elements of the momentum operator.

- [Qsymm](https://github.com/quantum-tinkerer/qsymm): to build effective models using the theory of invariants.

Therefore, if we use our code, please cite these references:

> Dániel Varjas, Tómas Ö Rosdahl, and Anton R Akhmero. "**Qsymm: algorithmic symmetry finding and symmetric Hamiltonian generation**", [New J. Phys. 20 093026 (2018)](https://doi.org/10.1088/1367-2630/aadf67)
---
> M. Iraola, J. L. Mañes, B. Bradlyn, M. K. Horton, T. Neupert, M. G. Vergniory and S. S. Tsirkin. "**IrRep: Symmetry eigenvalues and irreducible representations of ab initio band structures**", [Computer Physics Communications 272, 108226 (2022)](https://doi.org/10.1016/j.cpc.2021.108226)
---
> L. Elcoro, B. Bradlyn, Z. Wang, M. G. Vergniory, J. Cano, C. Felser, B. A. Bernevig, D. Orobengoa, G. de la Flor and M. I. Aroyo. "**Double crystallographic groups and their representations on the Bilbao Crystallographic Server**", [J. of Appl. Cryst. (2017). 50, 1457-1477](https://doi.org/10.1107/S1600576717011712)

