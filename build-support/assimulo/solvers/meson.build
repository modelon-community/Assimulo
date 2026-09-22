foreach f : ['__init__.py', 'euler.pyx', 'sundials.pyx', 'kinsol.pyx']
  stage += {'solvers/' + f: configure_file(input: '../../../src/solvers' / f, output: f, copy: true)}
endforeach
