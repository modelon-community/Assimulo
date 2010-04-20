from distutils.core import setup, Extension
import numpy as N
import sys as S

lib_file = ''

if S.platform == 'win32':
    lib_file = 'sundials_core.pyd'
else:
    lib_file = 'sundials_core.so'

setup(name='Assimulo',
      version='b1.02',
      description='A package for solving ordinary differential equations',
      author='Claus Führer and Christian Andersson',
      author_email='claus@maths.lth.se chria@kth.se',
      url='http://wwww.jmodelica.org/assimulo',
      package_dir = {'Assimulo':'src'},
      packages=['Assimulo', 'Assimulo.lib'],
      package_data={'Assimulo.lib':[lib_file]}
     )


