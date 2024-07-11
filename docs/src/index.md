# CHEASE.jl

Julia wrapper for the [CHEASE](https://gitlab.epfl.ch/spc/chease.git) fixed boundary equilibrium solver

See `examples/run_chease_example.jl` on how to use this package.

## Compile CHEASE executable

```
mamba install -c conda-forge gfortran
git clone https://gitlab.epfl.ch/spc/chease.git
cd chease/src-f90
make chease
```
