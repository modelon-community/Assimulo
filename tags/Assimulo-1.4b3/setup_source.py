#!/usr/bin/env python 
# -*- coding: utf-8 -*-
from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy as N
import sys as S
import os as O

incdirs = ''
libdirs = ''
SLUdir = ""
BLASdir = ""

if S.platform == 'win32':
    incdirs = ''
    libdirs = ''
elif S.platform == 'win64':
    incdirs = ''
    libdirs = ''
else:
    incdirs = '/usr/local/include'
    libdirs = '/usr/local/lib'
    
BLASname = 'blas'
    
BLASname_t = ""

copy_args=S.argv[1:]

for x in S.argv[1:]:
    if not x.find('--sundials-home'):
        print x[16:]
        incdirs = O.path.join(x[16:],'include')
        libdirs = O.path.join(x[16:],'lib')
        copy_args.remove(x)
    if not x.find('--prefix'):
        copy_args[copy_args.index(x)] = x.replace('/',O.sep)
    if not x.find('--superlu-home'):
        SLUdir = x[15:]
        print SLUdir
        copy_args.remove(x)
    if not x.find('--blas-home'):
        BLASdir = x[12:]
        print BLASdir
        copy_args.remove(x)
    if not x.find('--blas-name'):
        print x[12:]
        BLASname_t = x[12:]
        copy_args.remove(x)

if O.path.exists(O.path.join(O.path.join(incdirs,'cvodes'), 'cvodes.h')):
    
    cordir = O.path.join(O.path.join('src','lib'),'sundials_core.pyx')
    
    cordir_KINSOL_wSLU = O.path.join(O.path.join('src','lib'),'sundials_kinsol_core_wSLU.pyx')
    cordir_KINSOL = O.path.join(O.path.join('src','lib'),'sundials_kinsol_core.pyx')
    
    cordir_KINSOL_jmod_wSLU = O.path.join(O.path.join('src','lib'),'kinsol_jmod_wSLU.c')
    cordir_KINSOL_jmod = O.path.join(O.path.join('src','lib'),'kinsol_jmod.c')
    
    cordir_kinpinv = O.path.join(O.path.join('src','lib'),'kinpinv.c')
    cordir_kinslug = O.path.join(O.path.join('src','lib'),'kinslug.c')
    cordir_reg_routines = O.path.join(O.path.join('src','lib'),'reg_routines.c')

    wSLU = True
    if SLUdir != "":    
        SLUincdir = O.path.join(SLUdir,'SRC')
        SLUlibdir = O.path.join(SLUdir,'lib')
        if not O.path.exists(O.path.join(SLUincdir,'supermatrix.h')):
            wSLU = False
            print "--------------- WARNING ---------------"
            print "Could not find SuperLU at the given path."
            print "usage: --superlu-home=path"
            print "KINSOL will not be compiled with support for SUperLU."
            
        print "SLUinc: ",SLUincdir
        print "SLUlib: ", SLUlibdir

    else:
        print "--------------- WARNING ---------------"
        print "No path to SuperLU supplied, KINSOL will not be compiled with support for SUperLU."
        print "usage: --superlu-home=path"
        print "Note: the path required is to the folder where the folders 'SRC' and 'lib' are found."
        wSLU = False
        
    if BLASname_t != "":
        if BLASname_t.startswith("lib"):
            BLASname = BLASname_t[3:]
        else:
            BLASname = BLASname_t
            BLASname_t = "lib"+BLASname_t
    else:
        BLASname_t = "lib" + BLASname
           
    if BLASdir == "":
        print "--------------- WARNING ---------------"
        print "No path to BLAS supplied, KINSOL will not be compiledwith support for SUperLU."
        print "usage: --blas-home=path"
        print "Note: the path required is to where the static library","lib"+BLASname,"is found"
        wSLU = False
    else:
        if not O.path.exists(O.path.join(BLASdir,BLASname_t+'.a')):
            print "--------------- WARNING ---------------"
            print "Could not find BLAS at the given path."
            print "usage: --blas-home=path"
            print "KINSOL will not be compiled with support for SUperLU."
            wSLU = False
            
        print "BLAS: ", BLASdir + "/" + BLASname_t
        
    if wSLU:
        setup(name='Assimulo',
              version='1.4b3',
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
                Extension('lib.sundials_kinsol_core_wSLU',
                          [cordir_KINSOL_wSLU,cordir_KINSOL_jmod_wSLU,cordir_kinpinv,cordir_kinslug,cordir_reg_routines],
                          include_dirs=[incdirs, N.get_include(),SLUincdir],
                          library_dirs=[libdirs,SLUlibdir,BLASdir],
                          libraries=['sundials_kinsol','sundials_nvecserial','superlu_4.1',BLASname])
            
                          ],
            script_args=copy_args
            )
    else:
        setup(name='Assimulo',
              version='1.4b3',
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
                          [cordir_KINSOL,cordir_KINSOL_jmod,cordir_kinpinv],
                          include_dirs=[incdirs, N.get_include()],
                          library_dirs=[libdirs],
                          libraries=['sundials_kinsol','sundials_nvecserial'])
            
                          ],
            script_args=copy_args
            )

else:
    raise Exception('Could not find Sundials. Recheck Sundials path.')
