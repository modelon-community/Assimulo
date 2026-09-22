# ---------------------------------------------------------------- package files
pkg = 'assimulo'
version_txt = configure_file(output: 'version.txt', command: [py, '-c', 'import sys; print("@0@")'.format(meson.project_version())],
                             capture: true)
install_data(version_txt, install_dir: py.get_install_dir() / pkg)
py.install_sources(['__init__.py', 'exception.py', 'problem_algebraic.py',
                    'ode.pxd', 'explicit_ode.pxd', 'implicit_ode.pxd', 'algebraic.pxd', 'problem.pxd', 'support.pxd'],
                   subdir: pkg)
install_data(['../LICENSE', '../CHANGELOG', '../README.md'], install_dir: py.get_install_dir() / pkg)

# ---------------------------------------------------------------- core Cython modules
core_pxd = files('ode.pxd', 'explicit_ode.pxd', 'implicit_ode.pxd', 'algebraic.pxd', 'problem.pxd', 'support.pxd', 'constants.pxi')

src_inc = include_directories('.')       # ode_event_locator.h, declared in explicit_ode.pxd

core_modules = {
  'explicit_ode': ['ode_event_locator.c'],
  'algebraic': [], 'implicit_ode': [], 'ode': [], 'problem': [], 'special_systems': [], 'support': [],
}
foreach name, extra_c : core_modules
  c_src = custom_target(name + '_cython',
    input: stage[name + '.pyx'], output: name + '.c',
    command: cythonize + ['@INPUT@', '@OUTPUT@', '-m', pkg + '.' + name] + cython_inc,
    depend_files: core_pxd)
  py.extension_module(name, [c_src] + extra_c,
    include_directories: [np_inc, src_inc],
    c_args: cython_c_flags,
    dependencies: [py_dep],
    install: true, subdir: pkg)
endforeach

subdir('lib')
subdir('solvers')
