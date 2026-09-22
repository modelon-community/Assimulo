py.install_sources(['__init__.py', 'radau_core.py'], subdir: pkg / 'lib')
sundials_pxi = files('sundials_includes.pxd', 'sundials_constants.pxi', 'sundials_callbacks.pxi',
                     'sundials_callbacks_ida_cvode.pxi', 'sundials_callbacks_kinsol.pxi')

# TR-BDF2 (C stepper + Cython wrapper) -> assimulo.lib.trbdf2ode
trbdf2ode_c = custom_target('trbdf2ode_cython', input: stage['lib/trbdf2ode.pyx'], output: 'trbdf2ode.c',
  command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.lib.trbdf2ode'] + cython_inc,
  depend_files: core_pxd + files('trbdf2ode.pxd', 'trbdf2.h'))
py.extension_module('trbdf2ode', [trbdf2ode_c, files('trbdf2.c')], include_directories: [np_inc, include_directories('.')],
  c_args: cython_c_flags, dependencies: [py_dep, m_dep], install: true, subdir: pkg / 'lib')
