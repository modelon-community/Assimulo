#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#from distutils.core import setup, Extension
import numpy as N
import logging as L
import sys as S
import os as O
import shutil as SH
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    raise Exception("Please upgrade to a newer Cython version, >= 0.15.")

#L.basicConfig(format='%(levelname)s:%(message)s')

incdirs = ''
libdirs = ''
SLUdir = ""
BLASdir = ""
BLASname = 'blas'
BLASname_t = ""
debug = False

if S.platform == 'win32':
    incdirs = ''
    libdirs = ''
elif S.platform == 'win64':
    incdirs = ''
    libdirs = ''
else:
    incdirs = '/usr/local/include'
    libdirs = '/usr/local/lib'
    
static_link_gcc = ["-static-libgcc"]
static_link_gfortran = ["-static-libgfortran"]

copy_args=S.argv[1:]

for x in S.argv[1:]:
    if not x.find('--sundials-home'):
        incdirs = O.path.join(x[16:],'include')
        libdirs = O.path.join(x[16:],'lib')
        copy_args.remove(x)
    if not x.find('--prefix'):
        copy_args[copy_args.index(x)] = x.replace('/',O.sep)
    if not x.find('--superlu-home'):
        SLUdir = x[15:]
        copy_args.remove(x)
    if not x.find('--blas-home'):
        BLASdir = x[12:]
        copy_args.remove(x)
    if not x.find('--blas-name'):
        BLASname_t = x[12:]
        copy_args.remove(x)
    if not x.find('--debug'):
        debug = x[8:]
        if x[8:].upper() == "TRUE":
            debug = True
        copy_args.remove(x)
    if not x.find('--static'):
        static = x[9:]
        if x[9:].upper() == "TRUE":
            static = True
        else:
            static = False
        copy_args.remove(x)
    else:
        static = False
    if not x.find('--log'):
        level = x[6:]
        try:
            num_level = getattr(L, level.upper())
        except AttributeError:
            L.warning("No log-level defined for: "+level)
            num_level = 30
        L.basicConfig(level=num_level)
        copy_args.remove(x)
        
def pre_processing():
    join = O.path.join
    
    def create_dir(d):
        try:
            O.mkdir(d) #Create the build directory
        except:
            pass #Directory already exists
    create_dir(O.path.join("build"))
    create_dir(O.path.join("build","assimulo"))
    create_dir(O.path.join(O.path.join("build","assimulo"),"lib"))
    create_dir(O.path.join(O.path.join("build","assimulo"),"solvers"))
    create_dir(O.path.join(O.path.join("build","assimulo"),"examples"))
    create_dir(O.path.join(O.path.join("build","assimulo"),"tests"))
    create_dir(O.path.join(O.path.join("build","assimulo"),"tests","solvers"))
    create_dir(join("build","assimulo","thirdparty"))
    create_dir(join("build","assimulo","thirdparty","hairer"))
    
    fileSrc     = O.listdir("src")
    fileLib     = O.listdir(O.path.join("src","lib"))
    fileSolvers = O.listdir(O.path.join("src","solvers"))
    fileExamples= O.listdir("examples")
    fileMain    = ["setup.py","README","INSTALL","CHANGELOG","MANIFEST.in"]
    fileTests   = O.listdir("tests")
    fileTestsSolvers = O.listdir(O.path.join("tests","solvers"))
    fileThirdPartyHairer = O.listdir(join("thirdparty","hairer"))
    
    curdir = O.path.dirname(O.path.abspath(__file__))
    
    desSrc = O.path.join(curdir,O.path.join("build","assimulo"))
    desLib = O.path.join(curdir,O.path.join(O.path.join("build","assimulo"),"lib"))
    desSolvers = O.path.join(curdir,O.path.join("build","assimulo"),"solvers")
    desExamples = O.path.join(curdir,O.path.join("build","assimulo"),"examples")
    desMain = O.path.join(curdir,"build")
    desTests = O.path.join(curdir,O.path.join("build","assimulo"),"tests")
    desTestsSolvers = O.path.join(curdir,O.path.join("build","assimulo"),"tests","solvers")
    desThirdPartyHairer = join(curdir,"build","assimulo","thirdparty","hairer")

    for f in fileSrc:
        if not O.path.isdir(O.path.join("src",f)):
            SH.copy2(O.path.join("src",f), desSrc)
    for f in fileLib:
        if not O.path.isdir(O.path.join(O.path.join("src","lib"),f)):
            SH.copy2(O.path.join(O.path.join("src","lib"),f), desLib)
    for f in fileSolvers:
        if not O.path.isdir(O.path.join(O.path.join("src","solvers"),f)):
            SH.copy2(O.path.join(O.path.join("src","solvers"),f), desSolvers)
    for f in fileExamples:
        if not O.path.isdir(O.path.join("examples",f)):
            SH.copy2(O.path.join("examples",f), desExamples)
    for f in fileMain:
        if not O.path.isdir(f):
            SH.copy2(f,desMain)
    for f in fileTests:
        if not O.path.isdir(O.path.join("tests",f)):
            SH.copy2(O.path.join("tests",f), desTests)
    for f in fileTestsSolvers:
        if not O.path.isdir(O.path.join("tests","solvers",f)):
            SH.copy2(O.path.join("tests","solvers",f),desTestsSolvers)
    for f in fileThirdPartyHairer:
        if not O.path.isdir(join("thirdparty","hairer",f)):
            SH.copy2(join("thirdparty","hairer",f),desThirdPartyHairer)
        if f == "LICENSE":
            SH.copy2(join("thirdparty","hairer",f),join(curdir,"build","assimulo","lib"))
            
    #Delete OLD renamed files
    delFiles = [("lib","sundials_kinsol_core_wSLU.pxd")]
    for item in delFiles:
        dirDel = desSrc
        for f in item[:-1]:
            dirDel = O.path.join(dirDel, f)
        dirDel = O.path.join(dirDel, item[-1])
        if O.path.exists(dirDel):
            try:
                O.remove(dirDel)
            except:
                L.warning("Could not remove: "+str(dirDel))

def check_extensions():
    
    if static:
        extra_link_flags = static_link_gcc
    else:
        extra_link_flags = [""]
    
    #Cythonize main modules
    ext_list = cythonize(["assimulo"+O.path.sep+"*.pyx"], include_path=[".","assimulo"],include_dirs=[N.get_include()],pyrex_gdb=debug)
    
    #Cythonize Euler
    ext_list = ext_list + cythonize(["assimulo"+O.path.sep+"solvers"+O.path.sep+"euler.pyx"], include_path=[".","assimulo"],include_dirs=[N.get_include()],pyrex_gdb=debug)
    
    for i in ext_list:
        i.include_dirs = [N.get_include()]
        
        #Debug
        if debug:
            i.extra_compile_args = ["-g","-fno-strict-aliasing"]
            i.extra_link_args = ["-g"]
        else:
            i.extra_compile_args = ["-O2", "-fno-strict-aliasing"]
            
    #If Sundials
    if O.path.exists(O.path.join(O.path.join(incdirs,'cvodes'), 'cvodes.h')):
        ext_list = ext_list + cythonize(["assimulo"+O.path.sep+"solvers"+O.path.sep+"sundials.pyx"], include_path=[".","assimulo","assimulo"+O.sep+"lib"],include_dirs=[N.get_include()],pyrex_gdb=debug)
        ext_list[-1].include_dirs = [N.get_include(), "assimulo","assimulo"+O.sep+"lib", incdirs]
        ext_list[-1].library_dirs = [libdirs]
        ext_list[-1].libraries = ["sundials_cvodes", "sundials_nvecserial", "sundials_idas"]
        if debug:
            ext_list[-1].extra_compile_args = ["-g", "-fno-strict-aliasing"]
        else:
            ext_list[-1].extra_compile_args = ["-O2", "-fno-strict-aliasing"]
    
    #Sundials found
    if O.path.exists(O.path.join(O.path.join(incdirs,'cvodes'), 'cvodes.h')):
        cordir = O.path.join(O.path.join('assimulo','lib'),'sundials_core.pyx')
        cordir_KINSOL_wSLU = O.path.join(O.path.join('assimulo','lib'),'sundials_kinsol_core_wSLU.pyx')
        cordir_KINSOL = O.path.join(O.path.join('assimulo','lib'),'sundials_kinsol_core.pyx')
    
        cordir_KINSOL_jmod_wSLU = O.path.join(O.path.join('assimulo','lib'),'kinsol_jmod_wSLU.c')
        cordir_KINSOL_jmod = O.path.join(O.path.join('assimulo','lib'),'kinsol_jmod.c')
    
        cordir_kinpinv = O.path.join(O.path.join('assimulo','lib'),'kinpinv.c')
        cordir_kinslug = O.path.join(O.path.join('assimulo','lib'),'kinslug.c')
        cordir_reg_routines = O.path.join(O.path.join('assimulo','lib'),'reg_routines.c')

        
        wSLU = check_wSLU()
        if wSLU:
            SLUincdir = O.path.join(SLUdir,'SRC')
            SLUlibdir = O.path.join(SLUdir,'lib')
            #ext_list = ext_list + [Extension('assimulo.lib.sundials_kinsol_core_wSLU',
            #              [cordir_KINSOL_wSLU,cordir_KINSOL_jmod_wSLU,cordir_kinpinv,cordir_kinslug,cordir_reg_routines],
            #              include_dirs=[incdirs, N.get_include(),SLUincdir],
            #              library_dirs=[libdirs,SLUlibdir,BLASdir],
            #              libraries=['sundials_kinsol','sundials_nvecserial','superlu_4.1',BLASname])]
            ext_list = ext_list + cythonize([cordir_KINSOL_wSLU], include_path=[".","assimulo","assimulo"+O.sep+"lib"])
            ext_list[-1].sources += [cordir_KINSOL_jmod_wSLU,cordir_kinpinv,cordir_kinslug,cordir_reg_routines]
            ext_list[-1].include_dirs = [N.get_include(), SLUincdir, incdirs]
            ext_list[-1].library_dirs = [libdirs,SLUlibdir,BLASdir]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial", "superlu_4.1",BLASname]
            if debug:
                ext_list[-1].extra_compile_args = ["-g", "-fno-strict-aliasing"]
            else:
                ext_list[-1].extra_compile_args = ["-O2", "-fno-strict-aliasing"]
                
        else:
            #ext_list = ext_list + [Extension('assimulo.lib.sundials_kinsol_core',
            #              [cordir_KINSOL,cordir_KINSOL_jmod,cordir_kinpinv],
            #              include_dirs=[incdirs, N.get_include()],
            #              library_dirs=[libdirs],
            #              libraries=['sundials_kinsol','sundials_nvecserial'])]

            ext_list = ext_list + cythonize([cordir_KINSOL])#, include_path=[".","assimulo","assimulo"+O.sep+"lib"])
            ext_list[-1].sources += [cordir_KINSOL_jmod,cordir_kinpinv]
            ext_list[-1].include_dirs = [N.get_include(), incdirs]
            ext_list[-1].library_dirs = [libdirs]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial"]
            if debug:
                ext_list[-1].extra_compile_args = ["-g", "-fno-strict-aliasing"]
            else:
                ext_list[-1].extra_compile_args = ["-O2", "-fno-strict-aliasing"]
    
    for i in ext_list:
        i.extra_link_args += extra_link_flags
    
    return ext_list

def check_wSLU():
    wSLU = True
    
    global BLASname, BLASname_t
    
    if SLUdir != "":    
        SLUincdir = O.path.join(SLUdir,'SRC')
        SLUlibdir = O.path.join(SLUdir,'lib')
        if not O.path.exists(O.path.join(SLUincdir,'supermatrix.h')):
            wSLU = False
            L.warning("Could not find SuperLU, disabling support. View more information using --log=DEBUG")
            L.debug("Could not find SuperLU at the given path.")
            L.debug("usage: --superlu-home=path")
            L.debug("KINSOL will not be compiled with support for SUperLU.")
            
        L.debug("SLUinc: "+SLUincdir)
        L.debug("SLUlib: "+SLUlibdir)

    else:
        L.warning("No path to SuperLU supplied, disabling support. View more information using --log=DEBUG")
        L.debug("No path to SuperLU supplied, KINSOL will not be compiled with support for SUperLU.")
        L.debug("usage: --superlu-home=path")
        L.debug("Note: the path required is to the folder where the folders 'SRC' and 'lib' are found.")
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
        L.warning("No path to BLAS supplied, disabling support. View more information using --log=DEBUG")
        L.debug("No path to BLAS supplied, KINSOL will not be compiled with support for SUperLU.")
        L.debug("usage: --blas-home=path")
        L.debug("Note: the path required is to where the static library lib"+BLASname+" is found")
        wSLU = False
    else:
        if not O.path.exists(O.path.join(BLASdir,BLASname_t+'.a')):
            L.warning("Could not find BLAS, disabling support. View more information using --log=DEBUG")
            L.debug("Could not find BLAS at the given path.")
            L.debug("usage: --blas-home=path")
            L.debug("KINSOL will not be compiled with support for SUperLU.")
            wSLU = False
            
        L.debug("BLAS: "+BLASdir+"/"+BLASname_t)
    
    return wSLU


def check_fortran_extensions():
    """
    Adds the Fortran extensions using Numpy's distutils extension.
    """
    if static:
        extra_link_flags = static_link_gfortran+static_link_gcc
    else:
        extra_link_flags = [""]
    
    config = Configuration()

    config.add_extension('assimulo.lib.dopri5',
                         sources=['assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'dopri5.f','assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'dopri5.pyf']
                         ,extra_link_args=extra_link_flags)#include_dirs=[N.get_include()])
    
    config.add_extension('assimulo.lib.rodas',
                         sources=['assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'rodas_decsol.f','assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'rodas_decsol.pyf'],
                         include_dirs=[N.get_include()],extra_link_args=extra_link_flags)
    
    config.add_extension('assimulo.lib.radau5',
                         sources=['assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'radau_decsol.f','assimulo'+O.sep+'thirdparty'+O.sep+'hairer'+O.sep+'radau_decsol.pyf'],
                         include_dirs=[N.get_include()],extra_link_args=extra_link_flags)
    
    return config.todict()["ext_modules"]

"""
Pre-processing is necessary due to the setup of the repository. 
"""
if not O.path.isdir("assimulo"):
    pre_processing()
    O.chdir("build") #Change dir
    change_dir = True
else:
    change_dir = False
      
ext_list = check_extensions()
ext_list += check_fortran_extensions()


NAME = "Assimulo"
AUTHOR = "C. Andersson, C. Führer, J. Åkesson, M. Gäfvert"
AUTHOR_EMAIL = "chria@maths.lth.se"
VERSION = "2.1.1"
LICENSE = "GPL"
URL = "http://www.jmodelica.org/assimulo"
DOWNLOAD_URL = "http://www.jmodelica.org/assimulo"
DESCRIPTION = "A package for solving ordinary differential equations and differential algebraic equations."
PLATFORMS = ["Linux", "Windows", "MacOS X"]
CLASSIFIERS = [ 'Programming Language :: Python',
                'Programming Language :: Cython',
                'Programming Language :: C',
                'Programming Language :: Fortran',
                'Operating System :: MacOS :: MacOS X',
                'Operating System :: Microsoft :: Windows',
                'Operating System :: Unix']

LONG_DESCRIPTION = """
Assimulo is a Cython / Python based simulation package that allows for 
simulation of both ordinary differential equations (ODEs), f(t,y), and 
differential algebraic equations (DAEs), f(t,y,yd). It combines a 
variety of different solvers written in C, FORTRAN and Python via a 
common high-level interface.

Assimulo currently supports Explicit Euler, adaptive Runge-Kutta of 
order 4 and Runge-Kutta of order 4. It also wraps the popular SUNDIALS 
(https://computation.llnl.gov/casc/sundials/main.html) solvers CVode 
(for ODEs) and IDA (for DAEs). Ernst Hairer's 
(http://www.unige.ch/~hairer/software.html) codes Radau5, Rodas and 
Dopri5 are also available.

Documentation and installation instructions can be found at: 
http://www.jmodelica.org/assimulo . 

For questions and comments, visit: 
http://www.jmodelica.org/forums/jmodelicaorg-users/assimulo

The package requires Numpy, Scipy and Matplotlib and additionally for 
compiling from source, Cython 0.15 and Sundials 2.4 together with a 
C-compiler and a FORTRAN-compiler.
"""



setup(name=NAME,
      version=VERSION,
      license=LICENSE,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      platforms=PLATFORMS,
      classifiers=CLASSIFIERS,
      package_dir = {'assimulo':'assimulo'},
      packages=['assimulo', 'assimulo.lib','assimulo.solvers','assimulo.examples','assimulo.tests','assimulo.tests.solvers'],
      #cmdclass = {'build_ext': build_ext},
      ext_modules = ext_list,
      package_data={'assimulo': ['thirdparty'+O.sep+'hairer'+O.sep+'LICENSE','lib'+O.sep+'LICENSE']},
      script_args=copy_args)

if change_dir:
    O.chdir("..") #Change back to dir
