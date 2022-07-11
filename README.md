# A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations

> A parallel implementation of the preconditioner suitably interfaced with the PETSc library.

<p align="center"><img src="https://github.com/prj-/aldaas2021robust/raw/main/header.png" height="300"></p>

The code available in this repository can reproduce the results from the following [paper](https://epubs.siam.org/doi/abs/10.1137/21M1434891).
```
@article{aldaas2021robust,
    Author = {Al Daas, Hussam and Jolivet, Pierre and Scott, Jennifer A.},
    Title = {A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations},
    Year = {2022},
    Journal = {SIAM Journal on Scientific Computing},
    Pages = {A1047--A1068},
    Volume = {44},
    Issue = {3},
    Url = {https://github.com/prj-/aldaas2021robust}
}
```

## Getting started
### Dependencies
Make sure you have access to a recent [PETSc](https://petsc.org/) installation (version 3.18.0 or above), configured with the options `--download-slepc --download-hpddm`, using 32-bit indices, double-precision real-valued scalars (`--with-64-bit-indices=0 --with-precision=double --with-scalar-type=real`).  
Then, after setting the appropriate environment variable `${PETSC_DIR}` and `${PETSC_ARCH}` compile the single C source file from the repository with the following command.
```
$ make sparse_ls
```
When using Fortran, the above command should read `make sparse_ls_f` instead.  
When using Python, make sure that petsc4py is visible in your environment variable `${PYTHONPATH}`, e.g., `PYTHONPATH=${PYTHONPATH}:${PETSC_DIR}/${PETSC_ARCH}/lib`.

### Usage example
One should be able to launch the following commands, which solve respectively the linear least-squares problems associated to mesh_deform and lp_stocfor3 from the SuiteSparse Matrix Collection.
```
$ mpirun -np 4 ./sparse_ls -mat_name datafiles/mesh_deform.dat -options_file default.rc -pc_type hpddm
$ mpirun -np 4 ./sparse_ls -mat_name datafiles/lp_stocfor3.dat -options_file default.rc -pc_type hpddm
```
The command line option `-pc_type` may also be set to `asm`, `hypre` (if PETSc has been configured with `--download-hypre`), `gamg`, `qr` (if configured with `--download-suitesparse`), or whatever the user feels like trying out.  
All the other matrices are available at the following URL: `http://joliv.et/aldaas2021robust/mat_name.ext`, where `mat_name` is any of the identifier from the paper, and `ext` is either `dat` (PETSc binary format) or `mat` (MATLAB binary format).  
Here are two examples: http://joliv.et/aldaas2021robust/Hardesty2.dat and http://joliv.et/aldaas2021robust/cont11_l.mat. Matrices supplied through the command-line option `-mat_name` must always be in PETSc binary format. `.mat` files are merely provided for an easier inspection in MATLAB, e.g., using the [spy](https://mathworks.com/help/matlab/ref/spy.html) command.

## Acknowledgements
* HPC resources of [TGCC@CEA](http://www-hpc.cea.fr/index-en.htm) under the allocation A0090607519 made by [GENCI](http://www.genci.fr/en)
* L. Dalcin, V. Hapla, and T. Isaac for their recent contributions to PETSc that made the implementation more flexible
