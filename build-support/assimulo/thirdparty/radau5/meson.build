foreach f : ['radau5ode.pyx', 'radau5ode.pxd']
  stage += {'thirdparty/radau5/' + f: configure_file(input: '../../../../thirdparty/radau5' / f, output: f, copy: true)}
endforeach
