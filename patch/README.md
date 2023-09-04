# Patches

This folder contains patches to the QE code. These modify the files `.../PP/src/bands.f90` and `/PP/src/write_p_avg.f90` to add a new flag `lpall` to the `bands.x` code. 

Instructions on **how to apply the patch** and compile QE are presented in [INSTALL.md](../INSTALL.md) file and in the code documentation.

## The new flag

The patches add a new flag `lpall` to the namelist set of paramters for `bands.x`. This flag is only relevant if the flag `lp = True`. 

If `lpall=False`, the code `bands.x` will follow its default behaviour. It calculates the matrix elements of the velocity operator `v = i*m*[H, x]` between valence and conduction bands only, and stores `|v|²` in the output file `filp`.

If `lpall=True`, it computes the matrix elements of `v = i*m*[H, x]` for **all bands** and stores `v` into `filp`. For our kp model, we need `lpall=True`.

See the files `bands.in` in the Examples.

## Naming convention and compatibility

The patches are generated on top of a stable version of QE. For instance, `qe2kp-7.1.patch` was based on the QE code version 7.1, which is available in their [QE gitlab repository](https://gitlab.com/QEF/q-e) under the tag `qe-7.1`. 

The patch may work on other versions of QE, but it is not guaranteed.

To check if a patch will work with your QE version, check the `PP/src/bands.f90` and `PP/src/write_p_avg.f90` on your version match the ones in the supported QE version. For instance, to check if the patch built for `qe-7.1` will work with the `qe-7.2` tag, run

```bash
git clone https://gitlab.com/QEF/q-e.git

cd q-e/PP/src/

git diff qe-7.1 qe-7.2 -- bands.f90 write_p_avg.f90
```

If the output of `git diff` is blank, it means that these files match in both branches/tags and you can use the patch meant for one version in the other. If the output is not blank, you'll have to wait or ask us for a new patch.




