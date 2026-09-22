# ---------------------------------------------------------------- Radau5 (C) -> assimulo.lib.radau5ode
radau5_dir = 'radau5'
radau5_src = files(radau5_dir / 'radau5.c', radau5_dir / 'radau5_io.c')
radau5_c_args = cython_c_flags
radau5_deps = [py_dep, m_dep]
radau5_incs = [np_inc, src_inc, include_directories(radau5_dir)]
if have_superlu
  radau5_src += files(radau5_dir / 'superlu_double.c', radau5_dir / 'superlu_complex.c', radau5_dir / 'superlu_util.c')
  radau5_c_args += ['-D__RADAU5_WITH_SUPERLU']
  radau5_deps += superlu_libs + [blas_dep, gomp_dep]
  radau5_incs += [superlu_inc]
endif
radau5ode_c = custom_target('radau5ode_cython', input: stage['thirdparty/radau5/radau5ode.pyx'], output: 'radau5ode.c',
  command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.thirdparty.radau5.radau5ode'] + cython_inc,   # as setup.py; installed as assimulo.lib.radau5ode
  depend_files: core_pxd + files(radau5_dir / 'radau5ode.pxd'))
py.extension_module('radau5ode', [radau5ode_c] + radau5_src, include_directories: radau5_incs,
  c_args: radau5_c_args, dependencies: radau5_deps, install: true, subdir: pkg / 'lib')
install_data(radau5_dir / 'LICENSE_HAIRER', install_dir: py.get_install_dir() / pkg / 'lib')

# ---------------------------------------------------------------- Fortran solvers (f2py) -> assimulo.lib.*
if have_fortran
  f2py_gen = [py, meson.current_source_dir() / '..' / 'tools' / 'f2py_generate.py']
  fortranobject = files(f2py_srcdir / 'fortranobject.c')

  # name : [directory, signature file, fortran sources, extra fortran args, extra deps]
  fortran_modules = {
    'dopri5':  ['hairer', 'dopri5.pyf', ['dopri5.f'], [], []],
    'rodas':   ['hairer', 'rodas_decsol.pyf', ['rodas_decsol.f'], [], []],
    'radau5':  ['hairer', 'radau_decsol.pyf', ['radau_decsol.f'], [], []],
    'radar5':  ['hairer', 'radar5.pyf', ['contr5.f90', 'radar5_int.f90', 'radar5.f90', 'dontr5.f90', 'decsol.f90', 'dc_decdel.f90'], [], []],
    'odepack': ['odepack', 'odepack.pyf', ['opkdmain.f', 'opkda1.f', 'opkda2.f', 'odepack_aux.f90'], [], []],
    'odassl':  ['odassl', 'odassl.pyf', ['odassl.f', 'odastp.f', 'odacor.f', 'odajac.f', 'd1mach.f', 'daxpy.f', 'ddanrm.f', 'ddatrp.f',
                                         'ddot.f', 'ddwats.f', 'dgefa.f', 'dgesl.f', 'dscal.f', 'idamax.f', 'xerrwv.f'], [], []],
    'dasp3dp': ['dasp3', 'dasp3dp.pyf', ['DASP3.f', 'ANORM.f', 'CTRACT.f', 'DECOMP.f', 'HMAX.f', 'INIVAL.f', 'JACEST.f', 'PDERIV.f',
                                         'PREPOL.f', 'SOLVE.f', 'SPAPAT.f'], ['-fdefault-double-8', '-fdefault-real-8'], []],
  }
  if blas_dep.found() and lapack_dep.found()
    fortran_modules += {'glimda': ['glimda', 'glimda_complete.pyf', ['glimda_complete.f'], [], [lapack_dep, blas_dep]]}
  else
    warning('BLAS/LAPACK not found: GLIMDA will not be built (-Dblas-home=... -Dlapack-home=...)')
  endif

  foreach name, spec : fortran_modules
    dir = spec[0]
    f_src = []
    foreach f : spec[2]
      f_src += files(dir / f)
    endforeach
    gen = custom_target(name + '_f2py',
      input: dir / spec[1],
      output: [name + 'module.c', name + '-f2pywrappers.f', name + '-f2pywrappers2.f90'],
      command: f2py_gen + ['@INPUT@', '@OUTDIR@', name])
    py.extension_module(name, [gen, fortranobject] + f_src,
      include_directories: [np_inc, f2py_inc],
      c_args: c_flags,
      fortran_args: fortran_flags + spec[3],
      dependencies: [py_dep] + spec[4],
      install: true, subdir: pkg / 'lib')
    lic = dir / 'LICENSE_' + (dir == 'hairer' ? 'HAIRER' : dir.to_upper())
    if fs.exists(lic)
      install_data(lic, install_dir: py.get_install_dir() / pkg / 'lib')
    endif
  endforeach
endif
