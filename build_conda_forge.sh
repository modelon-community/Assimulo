#!/usr/bin/env bash
# Build and install Assimulo from source, with ALL Fortran solvers, inside an
# active conda-forge environment (micromamba/conda). See INSTALL-conda-forge.md
# for the background.
#
# Usage (env must be active, i.e. $CONDA_PREFIX set):
#   ./build_conda_forge.sh            # build + install
#   ./build_conda_forge.sh --test     # ... then run pytest
#
# Prerequisites in the env (all from conda-forge):
#   python=3.11 numpy scipy cython "setuptools=69.1.0" gfortran compilers
#   sundials superlu libblas liblapack matplotlib pytest
# Python must be < 3.12: setup.py relies on numpy.distutils for the f2py-wrapped
# Fortran solvers, and numpy.distutils does not exist on 3.12+.
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "error: activate the conda env first (CONDA_PREFIX is not set)" >&2
    exit 1
fi
cd "$(dirname "$0")"

python - <<'PY'
import sys, numpy
assert sys.version_info < (3, 12), "Python >= 3.12 has no numpy.distutils; the Fortran solvers would be dropped silently"
import numpy.distutils  # noqa: F401  (warns about deprecation, that is fine)
import setuptools
major = int(setuptools.__version__.split(".")[0])
assert major < 70, f"setuptools {setuptools.__version__} breaks numpy.distutils; pin setuptools=69.1.0"
PY

# setup.py looks for a SuperLU_MT prefix laid out as <home>/include/slu_mt_*.h and
# <home>/lib/libsuperlu_mt_OPENMP.a, and links -lsuperlu_mt_OPENMP -lblas_OPENMP.
# conda-forge's sundials package ships exactly that library (it is built with
# SUPERLUMT_ENABLE=ON), but puts the headers under include/superlu_mt/ and has no
# libblas_OPENMP. Build a symlink shim with the layout setup.py wants.
SLU_HOME="$CONDA_PREFIX/superlu_mt_home"
mkdir -p "$SLU_HOME/include" "$SLU_HOME/lib"
ln -sfn "$CONDA_PREFIX"/include/superlu_mt/* "$SLU_HOME/include/"
ln -sfn "$CONDA_PREFIX/lib/libsuperlu_mt_OPENMP.a" "$SLU_HOME/lib/"
ln -sfn "$CONDA_PREFIX/lib/libblas.so" "$SLU_HOME/lib/libblas_OPENMP.so"

# A previous install (e.g. the conda-forge assimulo package, which has no Fortran
# solvers) would shadow this one; remove it and any stale build tree.
pip uninstall -y assimulo >/dev/null 2>&1 || true
rm -rf build

python setup.py install \
    --sundials-home="$CONDA_PREFIX" \
    --blas-home="$CONDA_PREFIX/lib" \
    --lapack-home="$CONDA_PREFIX/lib" \
    --superlu-home="$SLU_HOME"

# Verify that the f2py-wrapped Fortran solvers actually got built.
( cd / && python - <<'PY'
from assimulo.solvers import CVode, IDA, Radau5ODE, RodasODE, Dopri5, LSODAR, ODASSL, GLIMDA, Radar5ODE, DASP3ODE
print("assimulo: all solvers import (Fortran: RodasODE Dopri5 LSODAR ODASSL GLIMDA Radar5ODE DASP3ODE)")
PY
)

if [[ "${1:-}" == "--test" ]]; then
    python -m pytest tests/ -q
fi
