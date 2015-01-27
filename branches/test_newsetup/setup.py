#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB  
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#from distutils.core import setup, Extension
import numpy as np
import logging as L
import sys 
import os
import shutil as SH
import ctypes.util
import argparse
import np.distutils


static_link_gcc = ["-static-libgcc"]
static_link_gfortran = ["-static-libgfortran"]
flag_32bit = ["-m32"]


parser = argparse.ArgumentParser(description='Assimulo setup script.')
package_arguments=['plugins','sundials','blas','superlu','lapack']
package_arguments.sort()
for pg in package_arguments:
    parser.add_argument("--{}-home".format(pg), 
           help="Location of the {} directory".format(pg.upper()),type=str,default='')
parser.add_argument("--blas-name", help="name of the blas package",default='blas')   
parser.add_argument("--extra-c-flags", help='Extra C-flags (a list enclosed in " ")',default='')                  
parser.add_argument("--is_static", action="store_true", help="set to true if present",default=False)
parser.add_argument("--debug", action="store_true", help="set to true if present",default=False)
parser.add_argument("--force-32bit", action="store_true", help="set to true if present",default=False)
parser.add_argument("--no-msvcr", action="store_true", help="set to true if present",default=False)
parser.add_argument("--log",choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),default='NOTSET')
parser.add_argument("--log_file",default=None,type=str,help='Path of a logfile')
                                       
args = parser.parse_known_args()


L.basicConfig(level=getattr(L,args[0].log),format='%(levelname)s:%(message)s',filename=args[0].log_file)
L.debug('setup.py called with the following optional args\n %s\n argument parsing completed.',vars(args[0]))
try:
    from subprocess import Popen, PIPE
    _p = Popen(["svnversion", "."], stdout=PIPE)
    revision = _p.communicate()[0].decode('ascii')
except:
    revision = "unknown"
L.debug('Source from svn revision {}.format(revision))

try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    msg="Please upgrade to a newer Cython version, >= 0.15."
    L.error(msg)
    raise Exception(msg)


L.debug('Python version used: {}'.format(sys.version.split()[0]))


thirdparty_methods = ["hairer","voigtmann", "odepack","odassl","dasp3"]


#  Has to be checked with chria
copy_args=sys.argv[1:]

for x in sys.argv[1:]:
    if not x.find('--prefix'):
        copy_args[copy_args.index(x)] = x.replace('/',os.sep)

class Assimulo_setup(object):
# helper functions
    def create_dir(self,d):
        try:
            os.makedirs(d) #Create the build directory
        except OSError:
            pass #Directory already exists
    def copy_file(self,fi, to_dir):
        # copies only files not directories
        if not os.path.isdir(fi):
            SH.copy2(fi, to_dir)
    def copy_all_files(self,file_list, from_dir, to_dir):
        for f in file_list:
            if from_dir:
                copy_file(os.path.join(from_dir,f),to_dir)
            else:
                copy_file(f,to_dir)
    def __init__(self,args,thirdpartymethods):
        # args[0] are optinal arguments given above
        # args[1] are argumenets passed to disutils 
        self.distutil_args=args[1]
        self.SLUdir = args[0].superlu_home
        self.BLASdir = args[0].blas_home 
        self.BLASname_t = args[0].blas_name if args[0].blas_name.startswith('lib') else 'lib'+args[0].blas_name
        self.BLASname = self.BLASname_t[3:]    # the name without "lib"
        self.debug_flag = args[0].debug 
        self.LAPACKdir = args[0].lapack_home
        self.PLUGINSdir = args[0].plugins_home
        self.static = args[0].is_static 
        self.debug_flag = args[0].debug 
        self.force_32bit = args[0].force_32bit 
        self.no_mvscr = args[0].no_msvcr 
        self.extra_c_flags = args[0].extra_c_flags
        self.thirdpartymethods = thirdpartymethods
        
        if args[0].no_msvcr:
        # prevent the MSVCR* being added to the DLLs passed to the linker
            def msvc_runtime_library_mod(): 
                return None
            np.distutils.misc_util.msvc_runtime_library = msvc_runtime_library_mod
            L.debug('numpy.distutils.misc_util.msvc_runtime_library overwritten.')
        
        self.platform = 'linux'
        self.platform = 'win' if 'win' in sys.platform
        self.platform = 'mac' if 'darwin' in sys.platform
        
        self.is_python3 = True if sys.version_info.major >= 3 else False
        L.debug('Platform {}'.format(self.platform))
    
        if args[0].sundials_home:
            self.incdirs = os.path.join(sundials_home,'include')
            self.libdirs = os.path.join(sundials_home,'lib')
        elif 'win' in self.platform:
            self.incdirs = ''
            self.libdirs = ''
        else:
            self.incdirs = '/usr/local/include'
            self.libdirs = '/usr/local/lib'            
        # directory paths
        self.curdir = os.path.dirname(os.path.abspath(__file__))
        # build directories
        self.build_assimulo = os.path.join("build","assimulo")
        self.build_assimulo_thirdparty = os.path.join(build_assimulo,'thirdparty')
        # destination directories
        self.desSrc = os.path.join(self.curdir,self.build_assimulo)
        self.desLib = os.path.join(self.desSrc,"lib")
        self.desSolvers = os.path.join(self.desSrc,"solvers")
        self.desExamples = os.path.join(self.desSrc,"examples")
        self.desMain = os.path.join(self.curdir,"build")
        self.desTests = os.path.join(self.desSrc,"tests")
        self.desTestsSolvers = os.path.join(self.desTests,"solvers")
        self.desThirdParty=dict([(thp,os.path.join(self.curdir,self.build_assimulo_thirdparty,thp)) 
                                          for thp in self.thirdparty_methods])
        
        # filelists
        
        self.fileSrc     = os.listdir("src")
        self.fileLib     = os.listdir(os.path.join("src","lib"))
        self.fileSolvers = os.listdir(os.path.join("src","solvers"))
        self.fileExamples= os.listdir("examples")
        self.fileMain    = ["setup.py","README","INSTALL","CHANGELOG","MANIFEST.in"]
        self.fileTests   = os.listdir("tests")
        self.filelist_thirdparty=dict([(thp,os.listdir(join("thirdparty",thp))) 
                                         for thp in self.thirdparty_methods])
        self.fileTestsSolvers = os.listdir(os.path.join("tests","solvers"))
        
    def create_assimulo_dirs_and_populate(self):
        for subdir in ["lib," "solvers", "examples"]:
            self.create_dir(os.path.join(self.build_assimulo,subdir))
        self.create_dir(os.path.join(build_assimulo, "tests", "solvers"))
        for pck in self.thirdparty_methods:
            self.create_dir(os.path.join(self.build_assimulo_thirdparty, pck))
        
        self.copy_all_files(self.fileSrc, "src", self.desSrc)
        self.copy_all_files(self.fileLib, "lib", self.desLib)
        self.copy_all_files(self.fileSolvers, os.path.join("src","solvers"), self.desSolvers)
        self.copy_all_files(self.fileExamples, "examples", self.desExamples)
        self.copy_all_files(self.fileMain, None, self.desMain)
        self.copy_all_files(self.fileTests, os.path.join("tests", f), self.desTests)
        self.copy_all_files(self.fileTestSolvers, os.path.join("tests","solvers"), self.desTestsSolvers)

        for solver in thirdparty_methods:
            for f in filelist_thirdparty.items():
                self.copy_all_files(f[1],os.path.join("thirdparty", f[0]), self.desThirdParty[f[0]])
                if fi == "LICENSE_{}".format(f[0].upper()):   
                    SH.copy2(join("thirdparty",f[0],fi),self.desLib)

        #Delete OLD renamed files
        delFiles = [("lib","sundials_kinsol_core_wSLU.pxd")]
        for item in delFiles:
            dirDel = desSrc
            for f in item[:-1]:
                dirDel = os.path.join(dirDel, f)
            dirDel = os.path.join(dirDel, item[-1])
            if os.path.exists(dirDel):
                try:
                    os.remove(dirDel)
                except:
                    L.debug("Could not remove: "+str(dirDel))
    
    def check_BLAS(self):
        self.with_BLAS = True
        msg=", disabling support. View more information using --log=DEBUG"
        if self.BLASdir == "":
            L.warning("No path to BLAS supplied" + msg)
            L.debug("usage: --blas-home=path")
            L.debug("Note: the path required is to where the static library lib"+BLASname+" is found")
            self.with_BLAS = False
        else:
            if not os.path.exists(os.path.join(self.BLASdir,self.BLASname_t+'.a')):
                L.warning("Could not find BLAS"+msg)
                L.debug("Could not find BLAS at the given path {}.".format(self.BLASdir)
                L.debug("usage: --blas-home=path")
                self.with_BLAS = False
            else:
                L.debug("BLAS found at "+BLASdir+"/"+BLASname_t)
        return self.with_BLAS
        
    def check_SuperLU(self):
        """
        Check if SuperLU package installed
        """
        self.with_SLU = True and self.check_BLAS()
        kinsol_msg='KINSOL will not be compiled with support for SuperLU'
        
        if !self.with_BLAS:
           L.warning(kinsol_msg+' as BLAS is missing.')
           return self.with_SLU
    
        if self.SLUdir != "":    
            self.SLUincdir = os.path.join(SLUdir,'SRC')
            self.SLUlibdir = os.path.join(SLUdir,'lib')
            if not os.path.exists(os.path.join(self.SLUincdir,'supermatrix.h')):
                self.with_SLU = False
                L.warning("Could not find SuperLU, disabling support. View more information using --log=DEBUG")
                L.debug("Could not find SuperLU at the given path {}.".format(self.SLUdir))
                L.debug("usage: --superlu-home path")
                L.debug(kinsol_msg+'.')
            
            L.debug("SuperLU found in {} and {}: ".format(self.SLUincdir, self.SLUlibdir)
        else:
            L.warning("No path to SuperLU supplied, disabling support. View more information using --log=DEBUG")
            L.debug("No path to SuperLU supplied, KINSOL will not be compiled with support for SUperLU.")
            L.debug("usage: --superlu-home=path")
            L.debug("Note: the path required is to the folder where the folders 'SRC' and 'lib' are found.")
            self.with_SLU = False
            L.debug(kinsol_msg+'.')
        return wSLU
    
    def check_SUNDIALS(self):
        if os.path.exists(O.path.join(O.path.join(incdirs,'cvodes'), 'cvodes.h')):
            self.with_SUNDIALS=True
            L.debug('SUNDIALS found.')
        else:    
            L.warning("Could not find Sundials, check the provided path (--sundials-home={}) "+ 
                    "to see that it actually points to Sundials.".format(self.sundials_home))
            L.debug("Could not find cvodes.h in " + O.path.join(incdirs,'cvodes'))
            self.with_SUNDIALS=False
        return self.with_SUNDIALS    
            

    def cython_extensionlists(self):
        extra_link_flags = []
        if self.static:
            extra_link_flags += self.static_link_gcc
        if self.force_32bit:
            extra_link_flags += self.flag_32bit
    
        #Cythonize main modules
        ext_list = cythonize(["assimulo"+os.path.sep+"*.pyx"], 
                             include_path=[".","assimulo"],
                             include_dirs=[np.get_include()],
                             pyrex_gdb=self.debug_flag)
         #Cythonize Solvers
         # Euler
        ext_list += cythonize(["assimulo"+os.path.sep+"solvers"+os.path.sep+"euler.pyx"], 
                             include_path=[".","assimulo"],
                             include_dirs=[np.get_include()],
                             pyrex_gdb=self.debug_flag)
        for el in ext_list:
            el.include_dirs = [np.get_include()]
            
        # SUNDIALS
        self.check_SUNDIALS()
        if self.with_SUNDIALS:
            #CVode and IDA
            ext_list += cythonize(["assimulo" + os.path.sep + "solvers" + os.path.sep + "sundials.pyx"], 
                                 include_path=[".","assimulo","assimulo" + os.sep + "lib"],
                                 include_dirs=[np.get_include()],
                                 pyrex_gdb=self.debug_flag)
            ext_list[-1].include_dirs = [np.get_include(), "assimulo","assimulo"+os.sep+"lib", incdirs]
            ext_list[-1].library_dirs = [libdirs]
            ext_list[-1].libraries = ["sundials_cvodes", "sundials_nvecserial", "sundials_idas"]
        
            #Kinsol
            ext_list += cythonize(["assimulo"+os.path.sep+"solvers"+os.path.sep+"kinsol.pyx"], 
                        include_path=[".","assimulo","assimulo"+os.sep+"lib"],
                        include_dirs=[np.get_include()],
                        pyrex_gdb=self.debug_flag)
            ext_list[-1].include_dirs = [np.get_include(), "assimulo","assimulo"+os.sep+"lib", incdirs]
            ext_list[-1].library_dirs = [libdirs]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial"]
    
        
            for el in ext_list:
                #Debug
                if self.debug_flag:
                    el.extra_compile_args = ["-g","-fno-strict-aliasing"]
                    el.extra_link_args = ["-g"]
                else:
                    el.extra_compile_args = ["-O2", "-fno-strict-aliasing"]
                if self.platform == "mac":
                    el.extra_compile_args += ["-Wno-error=return-type"]
                if self.force_32bit:
                    el.extra_compile_args += self.flag_32bit
                if self.extra_c_flags:
                    for f in self.extra_c_flags.split(' '):
                        el.extra_compile_args.append(f)
    
    #Sundials found
    if os.path.exists(os.path.join(os.path.join(incdirs,'cvodes'), 'cvodes.h')):
        cordir = os.path.join(os.path.join('assimulo','lib'),'sundials_core.pyx')
        cordir_KINSOL_wSLU = os.path.join(os.path.join('assimulo','lib'),'sundials_kinsol_core_wSLU.pyx')
        cordir_KINSOL = os.path.join(os.path.join('assimulo','lib'),'sundials_kinsol_core.pyx')
    
        cordir_KINSOL_jmod_wSLU = os.path.join(os.path.join('assimulo','lib'),'kinsol_jmod_wSLU.c')
        cordir_KINSOL_jmod = os.path.join(os.path.join('assimulo','lib'),'kinsol_jmod.c')
    
        cordir_kinpinv = os.path.join(os.path.join('assimulo','lib'),'kinpinv.c')
        cordir_kinslug = os.path.join(os.path.join('assimulo','lib'),'kinslug.c')
        cordir_reg_routines = os.path.join(os.path.join('assimulo','lib'),'reg_routines.c')

        
        self.check_SuperLU()
        if self.with_SLU:
            SLUincdir = os.path.join(args[0].superlu_home,'SRC')
            SLUlibdir = os.path.join(args[0].superlu_home,'lib')
            ext_list = ext_list + cythonize([cordir_KINSOL_wSLU], include_path=[".","assimulo","assimulo"+os.sep+"lib"])
            ext_list[-1].sources += [cordir_KINSOL_jmod_wSLU,cordir_kinpinv,cordir_kinslug,cordir_reg_routines]
            ext_list[-1].include_dirs = [np.get_include(), SLUincdir, incdirs]
            ext_list[-1].library_dirs = [libdirs,SLUlibdir,args[0].blas]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial", "superlu_4.1",args[0].blas_name,'gfortran']
            if debug_flag:
                ext_list[-1].extra_compile_args = ["-g", "-fno-strict-aliasing"]
            else:
                ext_list[-1].extra_compile_args = ["-O2", "-fno-strict-aliasing"]
            if self.platform == "mac":
                ext_list[-1].extra_compile_args += ["-Wno-error=return-type"]
            if self.force_32bit:
                ext_list[-1].extra_compile_args += flag_32bit
            for f in extra_c_flags.split(' '):
                    ext_list[-1].extra_compile_args.append(f)
        else:
            ext_list = ext_list + cythonize([cordir_KINSOL])#, include_path=[".","assimulo","assimulo"+os.sep+"lib"])
            ext_list[-1].sources += [cordir_KINSOL_jmod,cordir_kinpinv]
            ext_list[-1].include_dirs = [np.get_include(), incdirs]
            ext_list[-1].library_dirs = [libdirs]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial"]
            if debug_flag:
                ext_list[-1].extra_compile_args = ["-g", "-fno-strict-aliasing"]
            else:
                ext_list[-1].extra_compile_args = ["-O2", "-fno-strict-aliasing"]
            if self.platform == "mac":
                ext_list[-1].extra_compile_args += ["-Wno-error=return-type"]
            if force_32bit:
                ext_list[-1].extra_compile_args += flag_32bit
            for f in extra_c_flags.split(' '):
                    ext_list[-1].extra_compile_args.append(f)
        
            for el in ext_list:
                if is_python3:
                    el.cython_directives = {"language_level": 3} 
                el.extra_link_args += extra_link_flags
        return ext_list




def check_fortran_extensions():
    """
    Adds the Fortran extensions using Numpy's distutils extension.
    """
    extra_link_flags = []
    extra_compile_flags = []
    if static:
        extra_link_flags += static_link_gfortran + static_link_gcc
    if force_32bit:
        extra_link_flags += flag_32bit
        extra_compile_flags += flag_32bit
    if extra_c_flags:
        flags = extra_c_flags.split(' ')
        for f in flags:
            extra_compile_flags.append(f)
    
    from numpy.distutils.misc_util import Configuration
    config = Configuration()

    if force_32bit:
        config.add_extension('assimulo.lib.dopri5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dopri5.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dopri5.pyf']
                             ,extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:])#include_dirs=[np.get_include()])
        
        config.add_extension('assimulo.lib.rodas',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'rodas_decsol.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'rodas_decsol.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:])
        
        config.add_extension('assimulo.lib.radau5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radau_decsol.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radau_decsol.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:])
    
        config.add_extension('assimulo.lib.radar5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'contr5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5_int.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dontr5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'decsol.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dc_decdel.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])#, extra_f90_compile_args=["-O2"])#, extra_f77_compile_args=['-O2']) # extra_compile_args=['--noopt'])
        
        #ODEPACK
        config.add_extension('assimulo.lib.odepack',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkdmain.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkda1.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkda2.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'odepack_aux.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'odepack.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])
        
        #ODASSL
        odassl_dir='assimulo'+os.sep+'thirdparty'+os.sep+'odassl'+os.sep
        odassl_files=['odassl.pyf','odassl.f','odastp.f','odacor.f','odajac.f','d1mach.f','daxpy.f','ddanrm.f','ddatrp.f','ddot.f',
                      'ddwats.f','dgefa.f','dgesl.f','dscal.f','idamax.f','xerrwv.f']
        config.add_extension('assimulo.lib.odassl',
                             sources=[odassl_dir+file for file in odassl_files],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])
    else:
        config.add_extension('assimulo.lib.dopri5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dopri5.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dopri5.pyf']
                             ,extra_link_args=extra_link_flags[:])#include_dirs=[np.get_include()])
        
        config.add_extension('assimulo.lib.rodas',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'rodas_decsol.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'rodas_decsol.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])
        
        config.add_extension('assimulo.lib.radau5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radau_decsol.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radau_decsol.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])
    
        config.add_extension('assimulo.lib.radar5',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'contr5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5_int.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dontr5.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'decsol.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'dc_decdel.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'radar5.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])#, extra_f90_compile_args=["-O2"])#, extra_f77_compile_args=['-O2']) # extra_compile_args=['--noopt'])
        
        #ODEPACK
        config.add_extension('assimulo.lib.odepack',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkdmain.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkda1.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'opkda2.f',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'odepack_aux.f90',
                                      'assimulo'+os.sep+'thirdparty'+os.sep+'opepack'+os.sep+'odepack.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])
        
        #ODASSL
        odassl_dir='assimulo'+os.sep+'thirdparty'+os.sep+'odassl'+os.sep
        odassl_files=['odassl.pyf','odassl.f','odastp.f','odacor.f','odajac.f','d1mach.f','daxpy.f','ddanrm.f','ddatrp.f','ddot.f',
                      'ddwats.f','dgefa.f','dgesl.f','dscal.f','idamax.f','xerrwv.f']
        config.add_extension('assimulo.lib.odassl',
                             sources=[odassl_dir+file for file in odassl_files],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])
        
    #DASP3
    dasp3_f77_compile_flags = ["-fdefault-double-8","-fdefault-real-8"]
    if force_32bit:
        dasp3_f77_compile_flags += flag_32bit
    
    if np.version.version > "1.6.1": #NOTE, THERE IS A PROBLEM WITH PASSING F77 COMPILER ARGS FOR NUMPY LESS THAN 1.6.1, DISABLE FOR NOW
        dasp3_dir='assimulo'+os.sep+'thirdparty'+os.sep+'dasp3'+os.sep
        dasp3_files = ['dasp3dp.pyf', 'DASP3.f', 'ANORM.f','CTRACT.f','DECOMP.f',
                       'HMAX.f','INIVAL.f','JACEST.f','PDERIV.f','PREPOL.f','SOLVE.f','SPAPAT.f']
        config.add_extension('assimulo.lib.dasp3dp',
                              sources=[dasp3_dir+file for file in dasp3_files],
                              include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_f77_compile_args=dasp3_f77_compile_flags[:],extra_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])
    else:
        L.warning("DASP3 requires a numpy > 1.6.1. Disabling...")

    
    #GLIMDA
    #ADD liblapack and libblas
    lapack = False
    blas = False
    if LAPACKdir != "":
        lapack = True
        extra_link_flags += ["-L"+LAPACKdir, "-llapack"]
    else: #Try to see if Lapack exists in PATH
        name = ctypes.util.find_library("lapack")
        if name != None:
            extra_link_flags += ["-l"+name.split(os.path.sep)[-1].replace("lib","").split(".")[0]]
            lapack = True
    if BLASdir != "":
        blas = True
        extra_link_flags += ["-L"+BLASdir, "-lblas"]
    else: #Try to see if Blas exists in PATH
        name = ctypes.util.find_library("blas")
        if name != None:
            extra_link_flags += ["-l"+name.split(os.path.sep)[-1].replace("lib","").split(".")[0]]
            blas = True
    
    if lapack and blas:
        if force_32bit:
            config.add_extension('assimulo.lib.glimda',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'voigtmann'+os.sep+'glimda_complete.f','assimulo'+os.sep+'thirdparty'+os.sep+'voigtmann'+os.sep+'glimda_complete.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:],extra_compile_args=extra_compile_flags[:], extra_f77_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])
        else:
            config.add_extension('assimulo.lib.glimda',
                             sources=['assimulo'+os.sep+'thirdparty'+os.sep+'voigtmann'+os.sep+'glimda_complete.f','assimulo'+os.sep+'thirdparty'+os.sep+'voigtmann'+os.sep+'glimda_complete.pyf'],
                             include_dirs=[np.get_include()],extra_link_args=extra_link_flags[:])

    else:
        L.warning("Could not find Blas or Lapack, disabling support for the solver GLIMDA.")
    

    return config.todict()["ext_modules"]

"""
Pre-processing is necessary due to the setup of the repository. 
"""
if not os.path.isdir("assimulo"):
    pre_processing()
    os.chdir("build") #Change dir
    change_dir = True
else:
    change_dir = False
      
ext_list = check_extensions()

#MAJOR HACK DUE TO NUMPY CHANGE IN VERSION 1.6.2 THAT DOES NOT SEEM TO
#HANDLE EXTENSIONS OF BOTH TYPE (DISTUTILS AND NUMPY DISTUTILS) AT THE
#SAME TIME.
for e in ext_list:
    e.extra_f77_compile_args = []
    e.extra_f90_compile_args = []

ext_list += check_fortran_extensions()


NAME = "Assimulo"
AUTHOR = "C. Andersson, C. Führer, J. Åkesson, M. Gäfvert"
AUTHOR_EMAIL = "chria@maths.lth.se"
VERSION = "trunk"
LICENSE = "LGPL"
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
http://www.jmodelica.org/forums/jmodelicaorg-platform/assimulo

The package requires Numpy, Scipy and Matplotlib and additionally for 
compiling from source, Cython 0.15, Sundials 2.4/2.5, BLAS and LAPACK 
together with a C-compiler and a FORTRAN-compiler.
"""


version_txt = 'assimulo'+os.path.sep+'version.txt'
#If a revision is found, always write it!
if revision != "unknown" and revision!="":
    with open(version_txt, 'w') as f:
        f.write(VERSION+'\n')
        f.write("r"+revision)
else:# If it does not, check if the file exists and if not, create the file!
    if not os.path.isfile(version_txt):
        with open(version_txt, 'w') as f:
            f.write(VERSION+'\n')
            f.write("unknown")

from numpy.distutils.core import setup

license_info=[place+os.sep+pck+os.sep+'LICENSE_{}'.format(pck.upper()) 
               for pck in  thirdparty_methods for place in ['thirdparty','lib']]
L.debug(license_info)
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
      package_data={'assimulo': ['version.txt']+license_info+['examples'+os.sep+'kinsol_ors_matrix.mtx',
                                'examples'+os.sep+'kinsol_ors_matrix.mtx']},
      script_args=distutil_args)

if change_dir:
    os.chdir("..") #Change back to dir
