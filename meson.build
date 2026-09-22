project('assimulo', 'c',
  version: '3.8.0',
  meson_version: '>= 1.1.0',
  default_options: ['c_std=c99', 'warning_level=1'])

# ---------------------------------------------------------------- tools and python
py = import('python').find_installation(pure: false)
py_dep = py.dependency()
cc = meson.get_compiler('c')
fs = import('fs')

have_fortran = add_languages('fortran', native: false, required: get_option('fortran'))

np_incdir = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'], check: true).stdout().strip()
np_inc = include_directories(np_incdir)
f2py_srcdir = run_command(py, ['-c', 'import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))'],
                          check: true).stdout().strip()
f2py_inc = include_directories(f2py_srcdir)

# Same flags as the numpy.distutils build (setup.py), so that both builds give identical numerics
c_flags = cc.get_supported_arguments(['-O2', '-fno-strict-aliasing', '-fwrapv', '-fno-merge-constants', '-msse', '-msse2', '-msse3'])
cython_c_flags = c_flags + ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']   # Cython modules only; f2py's C uses the older API
fortran_flags = []
if have_fortran
  fc = meson.get_compiler('fortran')
  fortran_flags = fc.get_supported_arguments(['-O3', '-funroll-loops', '-fno-second-underscore', '-fPIC'] + get_option('extra-fortran-flags'))
endif

# ---------------------------------------------------------------- SUNDIALS
sundials_home = get_option('sundials-home')
if sundials_home == ''
  sundials_home = '/usr/local'
endif
sundials_incdir = sundials_home / 'include'
sundials_libdir = sundials_home / 'lib'
have_sundials = fs.exists(sundials_incdir / 'cvodes' / 'cvodes.h')
sundials_libs = []
sundials_cython_env = []
if have_sundials
  sconf = run_command(py, [meson.current_source_dir() / 'tools' / 'sundials_config.py', sundials_incdir], check: true).stdout().strip().split(' ')
  sundials_version = sconf[0].split('.')
  sundials_index_bits = sconf[1]
  sundials_with_superlu = sconf[2] == '1'
  sundials_rtol_vec = sconf[3] == '1'
  message('SUNDIALS @0@ found in @1@ (index @2@ bit, SuperLU_MT @3@)'.format(sconf[0], sundials_home, sundials_index_bits, sundials_with_superlu))
  sundials_libnames = ['sundials_cvodes', 'sundials_nvecserial', 'sundials_idas', 'sundials_sunlinsoldense',
                       'sundials_sunlinsolspgmr', 'sundials_sunmatrixdense', 'sundials_sunmatrixsparse']
  if sundials_version[0].to_int() >= 7
    sundials_libnames += ['sundials_core']
  endif
  kinsol_libnames = ['sundials_kinsol', 'sundials_nvecserial']
  if sundials_version[0].to_int() >= 7
    kinsol_libnames += ['sundials_core']
  endif
else
  warning('SUNDIALS not found under ' + sundials_home + ': CVode, IDA and Kinsol will not be built (-Dsundials-home=...)')
endif

# ---------------------------------------------------------------- SuperLU_MT
superlu_home = get_option('superlu-home')
have_superlu = false
superlu_libs = []
superlu_inc = []
if superlu_home != ''
  # layouts: a SuperLU_MT source tree (SRC/), an install prefix (include/), or a conda-forge
  # prefix whose sundials package ships the library (include/superlu_mt/, lib/)
  # (slu_mt_ddefs.h is specific to SuperLU_MT; supermatrix.h also exists in plain SuperLU)
  superlu_incdir = superlu_home / 'SRC'
  foreach cand : [superlu_home / 'include', superlu_home / 'include' / 'superlu_mt']
    if not fs.exists(superlu_incdir / 'slu_mt_ddefs.h')
      superlu_incdir = cand
    endif
  endforeach
  if fs.exists(superlu_incdir / 'slu_mt_ddefs.h')
    have_superlu = true
    superlu_inc = include_directories(superlu_incdir)
    superlu_libdir = superlu_home / 'lib'
    superlu_libs = [cc.find_library('superlu_mt_OPENMP', dirs: [superlu_libdir], static: true),
                    cc.find_library('blas_OPENMP', dirs: [superlu_libdir], required: false)]
    message('SuperLU_MT found in ' + superlu_home)
  else
    warning('SuperLU_MT not found under ' + superlu_home + ' (no slu_mt_ddefs.h): sparse solvers disabled')
  endif
endif

# ---------------------------------------------------------------- BLAS / LAPACK (GLIMDA, SuperLU)
blas_dirs = get_option('blas-home') != '' ? [get_option('blas-home')] : []
lapack_dirs = get_option('lapack-home') != '' ? [get_option('lapack-home')] : []
blas_dep = cc.find_library('blas', dirs: blas_dirs, required: false)
lapack_dep = cc.find_library('lapack', dirs: lapack_dirs, required: false)
m_dep = cc.find_library('m', required: false)
gomp_dep = cc.find_library('gomp', required: false)

# ---------------------------------------------------------------- Cython
# The .pyx files cimport "assimulo.ode" etc.: build-support/assimulo/ holds copies of the .pxd
# files under that package name for the Cython include path (the sources live in src/).
# Cython is driven by tools/cythonize_one.py (a compile-time env with a version tuple cannot
# be passed on cython's command line).
subdir('build-support')
cythonize = [py, meson.current_source_dir() / 'tools' / 'cythonize_one.py']
stage_dir = meson.current_build_dir() / 'build-support'
cython_inc = ['-I', stage_dir, '-I', stage_dir / 'assimulo', '-I', stage_dir / 'assimulo' / 'lib',
              '-I', stage_dir / 'assimulo' / 'thirdparty' / 'radau5']    # same as setup.py's include_path
sundials_cython_env = []
if have_sundials
  sundials_cython_env = ['-E', 'SUNDIALS_VERSION=(@0@,@1@,@2@)'.format(sundials_version[0], sundials_version[1], sundials_version[2]),
                         '-E', 'SUNDIALS_WITH_SUPERLU=' + ((sundials_with_superlu and have_superlu) ? 'True' : 'False'),
                         '-E', 'SUNDIALS_VECTOR_SIZE="@0@"'.format(sundials_index_bits),
                         '-E', 'SUNDIALS_CVODE_RTOL_VEC=' + (sundials_rtol_vec ? 'True' : 'False')]
endif

subdir('src')
subdir('thirdparty')
subdir('examples')
