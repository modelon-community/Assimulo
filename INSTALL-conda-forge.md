# Building Assimulo with all Fortran solvers in a conda-forge environment

Verified 2026-09-17 on Ubuntu 24.04 (WSL2), micromamba 2.9.0, Python 3.11.16.
Script: `./build_conda_forge.sh [--test]`.

## Why build from source at all

The conda-forge `assimulo` 3.8.0 package (`np2py311` build) is missing every
f2py-wrapped Fortran solver: `assimulo.lib.rodas`, `dopri5`, `odepack`, `odassl`,
`radar5`, `glimda` do not exist, so

```
from assimulo.solvers import RodasODE, Dopri5, LSODAR   # ImportError
```

`CVode`, `IDA`, `Radau5ODE` (C) still work. The cause is `setup.py`: the Fortran
extensions are only built when `numpy.distutils` can be imported (`have_nd`), and
when it cannot, they are dropped silently. `numpy.distutils` no longer exists on
Python >= 3.12 and is broken by setuptools >= 70 on older Pythons, so most
current build environments lose the Fortran solvers without any error. Fixing
this properly (a Meson build) is the reason this fork exists; until then the
build below is how to get a complete Assimulo.

## Requirements

Python 3.11 (not 3.12+), plus from conda-forge:

```
micromamba create -n pyfmi -c conda-forge python=3.11 numpy scipy cython \
    "setuptools=69.1.0" gfortran compilers cmake sundials superlu \
    libblas liblapack matplotlib pytest
```

`setuptools=69.1.0` is the same pin PyFMI's Dockerfile uses; newer versions
break `numpy.distutils`. numpy 2.x is fine on 3.11 (`numpy.distutils` is still
shipped there, with a deprecation warning).

conda-forge `sundials` is 7.8.0 here. `setup.py` reads the version from
`sundials_config.h` and adapts; 7.x works. The package is built with
`SUNDIALS_SUPERLUMT_ENABLED` and ships `lib/libsuperlu_mt_OPENMP.a` and
`include/superlu_mt/*.h`, which is what Assimulo's Radau5 sparse solver and the
SUNDIALS SuperLU_MT linear solver need. (The plain `superlu` package is a
different library and is not what `setup.py` looks for.)

## The SuperLU_MT shim

`setup.py --superlu-home=<home>` expects

```
<home>/include/supermatrix.h, slu_mt_ddefs.h, ...   (or <home>/SRC/)
<home>/lib/libsuperlu_mt_OPENMP.a
```

and links `-lsuperlu_mt_OPENMP -lblas_OPENMP` (the `blas_OPENMP` name comes
from the upstream Docker build, which copies `libblas.so` to
`libblas_OPENMP.so`). conda-forge has the files but not that layout, so the
script creates `$CONDA_PREFIX/superlu_mt_home` with symlinks:

```
P=$CONDA_PREFIX; S=$P/superlu_mt_home
mkdir -p $S/include $S/lib
ln -sfn $P/include/superlu_mt/* $S/include/
ln -sfn $P/lib/libsuperlu_mt_OPENMP.a $S/lib/
ln -sfn $P/lib/libblas.so $S/lib/libblas_OPENMP.so
```

Without `--superlu-home` the build still succeeds, but Radau5 is compiled without
`__RADAU5_WITH_SUPERLU` and 7 tests fail with
"Radau5 solver has not been compiled with superLU enabled".

## Build

```
micromamba activate pyfmi
pip uninstall -y assimulo          # remove the Fortran-less conda-forge package
python setup.py install --sundials-home=$CONDA_PREFIX --blas-home=$CONDA_PREFIX/lib \
    --lapack-home=$CONDA_PREFIX/lib --superlu-home=$CONDA_PREFIX/superlu_mt_home
python -m pytest tests/ -q          # 330 passed, 5 skipped
```

About 10 minutes, nearly all of it gfortran on `thirdparty/`. Things in the log
that look like errors but are not:

- `gcc: error: unrecognized command-line option '-mavx512er'` /
  `-mavx512pf`, and `cpu_avx512_knl.c ... implicit declaration` — these are
  `numpy.distutils` CPU-feature probes for Knights Landing; they fail on a
  modern gcc and are simply skipped.
- `Maybe empty "odepack-f2pywrappers.f"` — f2py noise.
- `sortvarnames: failed to compute dependencies because of cyclic dependencies`
  — f2py signature-file noise, also present in the upstream Docker build.

The installed package reports `assimulo.__version__ == '3.7.0.dev0'` unless
`--version` is passed to `setup.py`; that is cosmetic.

## Checking the result

```
python -c "from assimulo.solvers import CVode, IDA, Radau5ODE, RodasODE, Dopri5, LSODAR, ODASSL, GLIMDA, Radar5ODE, DASP3ODE"
ls $CONDA_PREFIX/lib/python3.11/site-packages/assimulo/lib/*.so
# dasp3dp dopri5 glimda odassl odepack radar5 radau5 radau5ode rodas
```
