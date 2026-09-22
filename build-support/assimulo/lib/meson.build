foreach f : ['__init__.py', 'sundials_includes.pxd', 'sundials_constants.pxi', 'sundials_callbacks.pxi',
             'sundials_callbacks_ida_cvode.pxi', 'sundials_callbacks_kinsol.pxi']
  stage += {'lib/' + f: configure_file(input: '../../../src/lib' / f, output: f, copy: true)}
endforeach
