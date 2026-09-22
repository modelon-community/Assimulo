foreach f : ['__init__.py', 'sundials_includes.pxd', 'sundials_constants.pxi', 'sundials_callbacks.pxi',
             'sundials_callbacks_ida_cvode.pxi', 'sundials_callbacks_kinsol.pxi', 'arkode_includes.pxd', 'trbdf2ode.pyx', 'trbdf2ode.pxd', 'trbdf2.h']
  stage += {'lib/' + f: configure_file(input: '../../../src/lib' / f, output: f, copy: true)}
endforeach
