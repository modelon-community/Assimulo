py.install_sources(['__init__.py', 'dasp3.py', 'glimda.py', 'odassl.py', 'odepack.py', 'radar5.py', 'radau5.py',
                    'rosenbrock.py', 'runge_kutta.py', 'trbdf2.py'], subdir: pkg / 'solvers')

# Euler
euler_c = custom_target('euler_cython', input: stage['solvers/euler.pyx'], output: 'euler.c',
  command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.solvers.euler'] + cython_inc, depend_files: core_pxd)
py.extension_module('euler', euler_c, include_directories: [np_inc, src_inc], c_args: cython_c_flags,
  dependencies: [py_dep], install: true, subdir: pkg / 'solvers')

# SUNDIALS: CVode / IDA and Kinsol
if have_sundials
  sundials_incs = [np_inc, src_inc, include_directories(sundials_incdir)]
  sundials_deps = [py_dep]
  foreach l : sundials_libnames
    sundials_deps += cc.find_library(l, dirs: [sundials_libdir])
  endforeach
  kinsol_deps = [py_dep]
  foreach l : kinsol_libnames
    kinsol_deps += cc.find_library(l, dirs: [sundials_libdir])
  endforeach
  if sundials_with_superlu and have_superlu
    sundials_deps += [cc.find_library('sundials_sunlinsolsuperlumt', dirs: [sundials_libdir])] + superlu_libs
    kinsol_deps += superlu_libs
    sundials_incs += [superlu_inc]
  endif
  sundials_c = custom_target('sundials_cython', input: stage['solvers/sundials.pyx'], output: 'sundials.c',
    command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.solvers.sundials'] + cython_inc + sundials_cython_env,
    depend_files: core_pxd + sundials_pxi)
  py.extension_module('sundials', sundials_c, include_directories: sundials_incs, c_args: cython_c_flags,
    dependencies: sundials_deps, install: true, subdir: pkg / 'solvers')
  kinsol_c = custom_target('kinsol_cython', input: stage['solvers/kinsol.pyx'], output: 'kinsol.c',
    command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.solvers.kinsol'] + cython_inc + sundials_cython_env,
    depend_files: core_pxd + sundials_pxi)
  py.extension_module('kinsol', kinsol_c, include_directories: sundials_incs, c_args: cython_c_flags,
    dependencies: kinsol_deps, install: true, subdir: pkg / 'solvers')
endif
