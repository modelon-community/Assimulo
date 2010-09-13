#!/usr/bin/env python 
# -*- coding: utf-8 -*-
from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy as N
import sys as S
import os as O

incdirs = ''
libdirs = ''

if S.platform == 'win32':
    incdirs = ''
    libdirs = ''
else:
    incdirs = '/usr/local/include'
    libdirs = '/usr/local/lib'

copy_args=S.argv[1:]

for x in S.argv[1:]:
    if not x.find('--sundials-home'):
        incdirs = O.path.join(x[16:],'include')
        libdirs = O.path.join(x[16:],'lib')
        copy_args.remove(x)
    if not x.find('--prefix'):
        copy_args[copy_args.index(x)] = x.replace('/',O.sep)

if O.path.exists(O.path.join(O.path.join(incdirs,'cvodes'), 'cvodes.h')):
    
    cordir = O.path.join(O.path.join('src','lib'),'sundials_core.pyx')
    cordir_KINSOL = O.path.join(O.path.join('src','lib'),'sundials_kinsol_core.pyx')

    setup(name='Assimulo',
      version='Assimulo-1.2.x',
      description='A package for solving ordinary differential equations',
      author='Claus Führer and Christian Andersson',
      author_email='claus@maths.lth.se chria@kth.se',
      url='http://wwww.jmodelica.org/assimulo',
      package_dir = {'assimulo':'src'},
      packages=['assimulo', 'assimulo.lib'],
      cmdclass = {'build_ext': build_ext},
      ext_package='assimulo',
      ext_modules = [
        Extension('lib.sundials_core',
            [cordir],
            include_dirs=[incdirs, N.get_include()],
            library_dirs=[libdirs],
            libraries=['sundials_cvodes','sundials_idas','sundials_nvecserial']),
        Extension('lib.sundials_kinsol_core',
            [cordir_KINSOL],
            include_dirs=[incdirs, N.get_include()],
            library_dirs=[libdirs],
            libraries=['sundials_kinsol','sundials_nvecserial'])
            
    ],
    script_args=copy_args
     )

else:
    raise Exception('Could not find Sundials. Recheck Sundials path.')
