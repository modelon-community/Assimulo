# A mirror of the Cython sources under the package name (what setup.py copies to build/assimulo/),
# so that "from assimulo.ode cimport ODE" resolves and each .pyx is compiled as assimulo.<name>.
stage = {}
foreach f : ['__init__.py', 'constants.pxi',
             'ode.pyx', 'ode.pxd', 'explicit_ode.pyx', 'explicit_ode.pxd', 'implicit_ode.pyx', 'implicit_ode.pxd',
             'algebraic.pyx', 'algebraic.pxd', 'problem.pyx', 'problem.pxd', 'special_systems.pyx',
             'support.pyx', 'support.pxd']
  stage += {f: configure_file(input: '../../src' / f, output: f, copy: true)}
endforeach
subdir('lib')
subdir('solvers')
subdir('thirdparty')
