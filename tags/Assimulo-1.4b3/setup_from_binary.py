#!/usr/bin/env python 
# -*- coding: utf-8 -*-
from distutils.core import setup, Extension
import numpy as N
import sys as S
import os as O

lib_file = ''

if S.platform == 'win32':
    lib_file = 'sundials_core.pyd'
    lib_file = 'sundials_kinsol_core.pyd'
    suffix = '.pyd'
else:
    lib_file = 'sundials_core.so'
    lib_file = 'sundials_kinsol_core.so'
    suffix = '.so'


setup(name='Assimulo',
      version='1.4b3',
      description='A package for solving ordinary differential equations',
      author='Claus Führer and Christian Andersson',
      author_email='claus@maths.lth.se chria@kth.se',
      url='http://wwww.jmodelica.org/assimulo',
      package_dir = {'assimulo':'src'},
      packages=['assimulo', 'assimulo.lib'],
      package_data={'assimulo':['lib'+O.sep+'*'+suffix]}
     )


