# CHEASE.jl

Julia wrapper for the CHEASE fixed boundary equilibrium solver (https://gitlab.epfl.ch/spc/chease.git)

## Compile CHEASE executable

```
mamba install -c conda-forge gfortran
git clone https://gitlab.epfl.ch/spc/chease.git
cd chease/src-f90
make chease
```

## Online documentation
For more details, see the [online documentation](https://projecttorreypines.github.io/CHEASE.jl/dev).

![Docs](https://github.com/ProjectTorreyPines/CHEASE.jl/actions/workflows/make_docs.yml/badge.svg)