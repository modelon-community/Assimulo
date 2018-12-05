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
from os import listdir
from os.path import isfile, join

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def remove_prefix(name, prefix):
    if name.startswith(prefix):
        return name[len(prefix):]
    return name

parser = argparse.ArgumentParser(description='Assimulo setup script.')
parser.register('type','bool',str2bool)
package_arguments=['plugins','sundials','blas','superlu','lapack']
package_arguments.sort()
for pg in package_arguments:
    parser.add_argument("--{}-home".format(pg), 
           help="Location of the {} directory".format(pg.upper()),type=str,default='')
parser.add_argument("--blas-name", help="name of the blas package",default='blas')   
parser.add_argument("--extra-c-flags", help='Extra C-flags (a list enclosed in " ")',default='')
parser.add_argument("--with_openmp", type='bool', help="set to true if present",default=False)
parser.add_argument("--is_static", type='bool', help="set to true if present",default=False)
parser.add_argument("--sundials-with-superlu", type='bool', help="(DEPRECATED) set to true if Sundials has been compiled with SuperLU",default=None)
parser.add_argument("--debug", type='bool', help="set to true if present",default=False)
parser.add_argument("--force-32bit", type='bool', help="set to true if present",default=False)
parser.add_argument("--no-msvcr", type='bool', help="set to true if present",default=False)
parser.add_argument("--log",choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),default='NOTSET')
parser.add_argument("--log_file",default=None,type=str,help='Path of a logfile')
parser.add_argument("--prefix",default=None,type=str,help='Path to destination directory')
parser.add_argument("--extra-fortran-link-flags", help='Extra Fortran link flags (a list enclosed in " ")', default='')
parser.add_argument("--extra-fortran-link-files", help='Extra Fortran link files (a list enclosed in " ")', default='')
parser.add_argument("--version", help='Package version number', default='Default')
                                       
args = parser.parse_known_args()
version_number_arg = args[0].version

L.basicConfig(level=getattr(L,args[0].log),format='%(levelname)s:%(message)s',filename=args[0].log_file)
L.debug('setup.py called with the following optional args\n %s\n argument parsing completed.',vars(args[0]))
try:
    from subprocess import Popen, PIPE
    _p = Popen(["svnversion", "."], stdout=PIPE)
    revision = _p.communicate()[0].decode('ascii')
except:
    revision = "unknown"
L.debug('Source from svn revision {}'.format(revision[:-1])) # exclude newline 

#If prefix is set, we want to allow installation in a directory that is not on PYTHONPATH
#and this is only possible with distutils, not setuptools
if not args[0].prefix:
    import setuptools
import numpy.distutils as nd

try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    msg="Please upgrade to a newer Cython version, >= 0.18."
    L.error(msg)
    raise Exception(msg)

#Verify Cython version
import Cython
cython_version = Cython.__version__.split(".")
if not (cython_version[0] > '0' or (cython_version[0] == '0' and cython_version[1] >= '18')):
    msg="Please upgrade to a newer Cython version, >= 0.18."
    L.error(msg)
    raise Exception(msg)

L.debug('Python version used: {}'.format(sys.version.split()[0]))

thirdparty_methods= ["hairer","glimda", "odepack","odassl","dasp3"] 

class Assimulo_prepare(object):
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
        L.debug('fromdir {}  todir {}'.format(from_dir,to_dir))
        for f in file_list:
            if from_dir:
                self.copy_file(os.path.join(from_dir,f),to_dir)
            else:
                self.copy_file(f,to_dir)
    def __init__(self,args, thirdparty_methods):
        # args[0] are optinal arguments given above
        # args[1] are argumenets passed to disutils 
        self.distutil_args=args[1]
        if args[0].prefix:
            self.prefix = args[0].prefix.replace('/',os.sep)   # required in this way for cygwin etc.
            self.distutil_args.append('--prefix={}'.format(self.prefix))
        self.SLUdir = args[0].superlu_home
        self.BLASdir = args[0].blas_home 
        self.sundialsdir = args[0].sundials_home
        self.sundials_with_superlu = args[0].sundials_with_superlu
        self.BLASname_t = args[0].blas_name if args[0].blas_name.startswith('lib') else 'lib'+args[0].blas_name
        self.BLASname = self.BLASname_t[3:]    # the name without "lib"
        self.debug_flag = args[0].debug 
        self.LAPACKdir = args[0].lapack_home
        self.LAPACKname = ""
        self.PLUGINSdir = args[0].plugins_home
        self.static = args[0].is_static 
        self.static_link_gcc = ["-static-libgcc"] if self.static else []
        self.static_link_gfortran = ["-static-libgfortran"] if self.static else []
        self.force_32bit = args[0].force_32bit
        self.flag_32bit = ["-m32"] if self.force_32bit else [] 
        self.no_mvscr = args[0].no_msvcr 
        self.extra_c_flags = args[0].extra_c_flags.split()
        self.extra_fortran_link_flags = args[0].extra_fortran_link_flags.split()
        self.extra_fortran_link_files = args[0].extra_fortran_link_files.split()
        self.thirdparty_methods  = thirdparty_methods
        self.with_openmp = args[0].with_openmp

        if self.sundials_with_superlu is not None:
            L.warning("The option 'sundials_with_superlu' has been deprecated and has no effect. Support for SuperLU using Sundials is automatically checked.")
        
        if self.no_mvscr:
        # prevent the MSVCR* being added to the DLLs passed to the linker
            def msvc_runtime_library_mod(): 
                return None
            nd.misc_util.msvc_runtime_library = msvc_runtime_library_mod
            L.debug('numpy.distutils.misc_util.msvc_runtime_library overwritten.')
        
        self.platform = 'linux'
        if 'win' in sys.platform: self.platform = 'win'
        if 'darwin' in sys.platform: self.platform = 'mac' 
        
        self.is_python3 = True if sys.version_info.major >= 3 else False
        L.debug('Platform {}'.format(self.platform))
        
        if args[0].sundials_home:
            self.incdirs = os.path.join(self.sundialsdir,'include')
            self.libdirs = os.path.join(self.sundialsdir,'lib')
        elif 'win' in self.platform:
            self.incdirs = ''
            self.libdirs = ''
        else:
            self.incdirs = '/usr/local/include'
            self.libdirs = '/usr/local/lib'            
        
        self.assimulo_lib = os.path.join('assimulo','lib')
        
        # check packages
        self.check_BLAS()
        self.check_SuperLU()
        self.check_SUNDIALS()
        self.check_LAPACK()
        
    def _set_directories(self):
        # directory paths
        self.curdir = os.path.dirname(os.path.abspath(__file__))
        # build directories
        self.build_assimulo = os.path.join("build","assimulo")
        self.build_assimulo_thirdparty = os.path.join(self.build_assimulo,'thirdparty')
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
        self.fileMainIncludes = ["README","CHANGELOG", "LICENSE"]
        self.fileTests   = os.listdir("tests")
        self.filelist_thirdparty=dict([(thp,os.listdir(os.path.join("thirdparty",thp))) 
                                         for thp in self.thirdparty_methods])
        self.fileTestsSolvers = os.listdir(os.path.join("tests","solvers"))
        
    def create_assimulo_dirs_and_populate(self):
        self._set_directories()
        
        for subdir in ["lib", "solvers", "examples"]:
            self.create_dir(os.path.join(self.build_assimulo,subdir))
        self.create_dir(os.path.join(self.build_assimulo, "tests", "solvers"))
        for pck in self.thirdparty_methods:
            self.create_dir(os.path.join(self.build_assimulo_thirdparty, pck))
        
        self.copy_all_files(self.fileSrc, "src", self.desSrc)
        self.copy_all_files(self.fileLib, "src/lib", self.desLib)
        self.copy_all_files(self.fileSolvers, os.path.join("src","solvers"), self.desSolvers)
        self.copy_all_files(self.fileExamples, "examples", self.desExamples)
        self.copy_all_files(self.fileMain, None, self.desMain)
        self.copy_all_files(self.fileMainIncludes, None, self.desSrc)
        self.copy_all_files(self.fileTests, "tests", self.desTests)
        self.copy_all_files(self.fileTestsSolvers, os.path.join("tests","solvers"), self.desTestsSolvers)

        for f in self.filelist_thirdparty.items():
            L.debug('Thirdparty method {} file {} copied'.format(f[0],f[1]))
            self.copy_all_files(f[1],os.path.join("thirdparty", f[0]), self.desThirdParty[f[0]])
            try:   
                SH.copy2(os.path.join("thirdparty",f[0],"LICENSE_{}".format(f[0].upper())),self.desLib)
            except IOError:
                L.warning('No license file {} found.'.format("LICENSE_{}".format(f[0].upper())))

        #Delete OLD renamed files
        delFiles = [("lib","sundials_kinsol_core_wSLU.pxd")]
        for item in delFiles:
            dirDel = self.desSrc
            for f in item[:-1]:
                dirDel = os.path.join(dirDel, f)
            dirDel = os.path.join(dirDel, item[-1])
            if os.path.exists(dirDel):
                try:
                    os.remove(dirDel)
                except:
                    L.debug("Could not remove: "+str(dirDel))
        
        if self.extra_fortran_link_files:
            for extra_fortran_lib in self.extra_fortran_link_files:
                path_extra_fortran_lib = ctypes.util.find_library(extra_fortran_lib)
                if path_extra_fortran_lib != None:
                    SH.copy2(path_extra_fortran_lib,self.desSrc)
                else:
                    L.debug("Could not find Fortran link file: "+str(extra_fortran_lib))
    
    def check_BLAS(self):
        """
        Check if BLAS can be found
        """
        self.with_BLAS = True
        msg=", disabling support. View more information using --log=DEBUG"
        if self.BLASdir == "":
            L.warning("No path to BLAS supplied" + msg)
            L.debug("usage: --blas-home=path")
            L.debug("Note: the path required is to where the static library lib is found")
            self.with_BLAS = False
        else:
            if not os.path.exists(os.path.join(self.BLASdir,self.BLASname_t+'.a')):
                L.warning("Could not find BLAS"+msg)
                L.debug("Could not find BLAS at the given path {}.".format(self.BLASdir))
                L.debug("usage: --blas-home=path")
                self.with_BLAS = False
            else:
                L.debug("BLAS found at "+self.BLASdir)
                self.with_BLAS = True
        
    def check_SuperLU(self):
        """
        Check if SuperLU package installed
        """
        self.with_SLU = True
        sundials_msg='SUNDIALS will not be compiled with support for SuperLU'
        
        if self.SLUdir != "":    
            self.SLUincdir = os.path.join(self.SLUdir,'SRC')
            self.SLUlibdir = os.path.join(self.SLUdir,'lib')
            if not os.path.exists(os.path.join(self.SLUincdir,'supermatrix.h')):
                self.with_SLU = False
                L.warning("Could not find SuperLU, disabling support. View more information using --log=DEBUG")
                L.debug("Could not find SuperLU at the given path {}.".format(self.SLUdir))
                L.debug("usage: --superlu-home path")
                L.debug(sundials_msg+'.')
            else:
                L.debug("SuperLU found in {} and {}: ".format(self.SLUincdir, self.SLUlibdir))
            
            potential_files = [remove_prefix(f.rsplit(".",1)[0],"lib") for f in listdir(self.SLUlibdir) if isfile(join(self.SLUlibdir, f)) and f.endswith(".a")]
            potential_files.sort(reverse=True)
            L.debug("Potential SuperLU files: "+str(potential_files))
            
            self.superLUFiles = []
            for f in potential_files:
                if "superlu" in f:
                    self.superLUFiles.append(f)
                #if self.with_BLAS == False and "blas" in f:
                #    self.superLUFiles.append(f)
                if "blas" in f:
                    self.superLUFiles.append(f)
                    
            #if self.with_BLAS:
            #    self.superLUFiles.append(self.BLASname)
            
            L.debug("SuperLU files: "+str(self.superLUFiles))
            
        else:
            L.warning("No path to SuperLU supplied, disabling support. View more information using --log=DEBUG")
            L.debug("No path to SuperLU supplied, SUNDIALS will not be compiled with support for SuperLU.")
            L.debug("usage: --superlu-home=path")
            L.debug("Note: the path required is to the folder where the folders 'SRC' and 'lib' are found.")
            self.with_SLU = False
    
    def check_SUNDIALS(self):
        """
        Check if Sundials installed
        """
        if os.path.exists(os.path.join(os.path.join(self.incdirs,'cvodes'), 'cvodes.h')):
            self.with_SUNDIALS=True
            L.debug('SUNDIALS found.')
            sundials_version = None
            sundials_vector_type_size = None
            sundials_with_superlu = False
            
            try:
                if os.path.exists(os.path.join(os.path.join(self.incdirs,'sundials'), 'sundials_config.h')):
                    with open(os.path.join(os.path.join(self.incdirs,'sundials'), 'sundials_config.h')) as f:
                        for line in f:
                            if "SUNDIALS_PACKAGE_VERSION" in line or "SUNDIALS_VERSION" in line:
                                sundials_version = tuple([int(f) for f in line.split()[-1][1:-1].split(".")])
                                L.debug('SUNDIALS %d.%d found.'%(sundials_version[0], sundials_version[1]))
                                break
                    with open(os.path.join(os.path.join(self.incdirs,'sundials'), 'sundials_config.h')) as f:
                        for line in f:
                            if "SUNDIALS_INT32_T" in line and line.startswith("#define"):
                                sundials_vector_type_size = "32"
                                L.debug('SUNDIALS vector type size %s bit found.'%(sundials_vector_type_size))
                                break
                            if "SUNDIALS_INT64_T" in line and line.startswith("#define"):
                                sundials_vector_type_size = "64"
                                L.debug('SUNDIALS vector type size %s bit found.'%(sundials_vector_type_size))
                                if self.with_SLU:
                                    L.warning("It is recommended to set the SUNDIALS_INDEX_TYPE to an 32bit integer when using SUNDIALS together with SuperLU (or make sure that SuperLU is configured to use the same int size).")
                                    L.warning("SuperLU may not function properly.")
                                break
                    with open(os.path.join(os.path.join(self.incdirs,'sundials'), 'sundials_config.h')) as f:
                        for line in f:
                            if "SUNDIALS_SUPERLUMT" in line and line.startswith("#define"): #Sundials compiled with support for SuperLU
                                sundials_with_superlu = True
                                L.debug('SUNDIALS found to be compiled with support for SuperLU.')
                                break
            except Exception as e:
                if os.path.exists(os.path.join(os.path.join(self.incdirs,'arkode'), 'arkode.h')): #This was added in 2.6
                    sundials_version = (2,6,0)
                    L.debug('SUNDIALS 2.6 found.')
                else:
                    sundials_version = (2,5,0)
                    L.debug('SUNDIALS 2.5 found.')
                
            self.SUNDIALS_version = sundials_version
            self.SUNDIALS_vector_size = sundials_vector_type_size
            self.sundials_with_superlu = sundials_with_superlu
            if not self.sundials_with_superlu:
                L.debug("Could not detect SuperLU support with Sundials, disabling support for SuperLU.")
        else:    
            L.warning(("Could not find Sundials, check the provided path (--sundials-home={}) "+ 
                    "to see that it actually points to Sundials.").format(self.sundialsdir))
            L.debug("Could not find cvodes.h in " + os.path.join(self.incdirs,'cvodes'))
            self.with_SUNDIALS=False
            
    def check_LAPACK(self):
        """
        Check if LAPACK installed
        """
        msg=", disabling support. View more information using --log=DEBUG"
        self.with_LAPACK=False
        if self.LAPACKdir != "":
            if not os.path.exists(self.LAPACKdir):
                L.warning('LAPACK directory {} not found'.format(self.LAPACKdir))
            else:
                L.debug("LAPACK found at "+self.LAPACKdir)
                self.with_LAPACK = True
        else:
            """
            name = ctypes.util.find_library("lapack")
            if name != "":
                L.debug('LAPACK found in standard library path as {}'.format(name))
                self.with_LAPACK=True
                self.LAPACKname = name
            else:
            """
            L.warning("No path to LAPACK supplied" + msg)
            L.debug("usage: --lapack-home=path")
            L.debug("Note: the path required is to where the static library lib is found")
            self.with_LAPACK = False
            
    def cython_extensionlists(self):
        extra_link_flags = self.static_link_gcc + self.flag_32bit
    
        #Cythonize main modules
        ext_list = cythonize(["assimulo"+os.path.sep+"*.pyx"], 
                             include_path=[".","assimulo"])
         #Cythonize Solvers
         # Euler
        ext_list += cythonize(["assimulo"+os.path.sep+"solvers"+os.path.sep+"euler.pyx"], 
                             include_path=[".","assimulo"])
        for el in ext_list:
            el.include_dirs = [np.get_include()]
            
        # SUNDIALS
        if self.with_SUNDIALS:
            compile_time_env = {'SUNDIALS_VERSION': self.SUNDIALS_version,
                                'SUNDIALS_WITH_SUPERLU': self.sundials_with_superlu and self.with_SLU,
                                'SUNDIALS_VECTOR_SIZE': self.SUNDIALS_vector_size}
            #CVode and IDA
            ext_list += cythonize(["assimulo" + os.path.sep + "solvers" + os.path.sep + "sundials.pyx"], 
                                 include_path=[".","assimulo","assimulo" + os.sep + "lib"],
                                 compile_time_env=compile_time_env, force=True)
            ext_list[-1].include_dirs = [np.get_include(), "assimulo","assimulo"+os.sep+"lib", self.incdirs]
            ext_list[-1].library_dirs = [self.libdirs]
            
            if self.SUNDIALS_version >= (3,0,0):
                ext_list[-1].libraries = ["sundials_cvodes", "sundials_nvecserial", "sundials_idas", "sundials_sunlinsoldense", "sundials_sunlinsolspgmr", "sundials_sunmatrixdense", "sundials_sunmatrixsparse"]
            else:
                ext_list[-1].libraries = ["sundials_cvodes", "sundials_nvecserial", "sundials_idas"]
            if self.sundials_with_superlu and self.with_SLU: #If SUNDIALS is compiled with support for SuperLU
                if self.SUNDIALS_version >= (3,0,0):
                    ext_list[-1].libraries.extend(["sundials_sunlinsolsuperlumt"])
                
                ext_list[-1].include_dirs.append(self.SLUincdir)
                ext_list[-1].library_dirs.append(self.SLUlibdir)
                ext_list[-1].libraries.extend(self.superLUFiles)
                
        
            #Kinsol
            ext_list += cythonize(["assimulo"+os.path.sep+"solvers"+os.path.sep+"kinsol.pyx"], 
                        include_path=[".","assimulo","assimulo"+os.sep+"lib"],
                        compile_time_env=compile_time_env, force=True)
            ext_list[-1].include_dirs = [np.get_include(), "assimulo","assimulo"+os.sep+"lib", self.incdirs]
            ext_list[-1].library_dirs = [self.libdirs]
            ext_list[-1].libraries = ["sundials_kinsol", "sundials_nvecserial"]
            
            if self.sundials_with_superlu and self.with_SLU: #If SUNDIALS is compiled with support for SuperLU
                ext_list[-1].include_dirs.append(self.SLUincdir)
                ext_list[-1].library_dirs.append(self.SLUlibdir)
                ext_list[-1].libraries.extend(self.superLUFiles)
        
        for el in ext_list:
            #Debug
            if self.debug_flag:
                el.extra_compile_args = ["-g","-fno-strict-aliasing"]
                el.extra_link_args = ["-g"]
            else:
                el.extra_compile_args = ["-O2", "-fno-strict-aliasing"]
            if self.platform == "mac":
                el.extra_compile_args += ["-Wno-error=return-type"]
            if self.with_openmp:
                el.extra_link_args.append("-fopenmp")
                el.extra_compile_args.append("-fopenmp")
            el.extra_compile_args += self.flag_32bit + self.extra_c_flags

        for el in ext_list:
            if self.is_python3:
                el.cython_directives = {"language_level": 3} 
            el.extra_link_args += extra_link_flags
        return ext_list

    def fortran_extensionlists(self):
        """
        Adds the Fortran extensions using Numpy's distutils extension.
        """
        extra_link_flags = self.static_link_gfortran + self.static_link_gcc + self.flag_32bit + self.extra_fortran_link_flags
        extra_compile_flags = self.flag_32bit + self.extra_c_flags
        
        config = np.distutils.misc_util.Configuration()
        extraargs={'extra_link_args':extra_link_flags[:], 'extra_compile_args':extra_compile_flags[:]}
                  
        extraargs['extra_f77_compile_args'] = extra_compile_flags[:]
        extraargs['extra_f90_compile_args'] = extra_compile_flags[:]
    
        #Hairer
        sources='assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'{0}.f','assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+'{0}.pyf'
        config.add_extension('assimulo.lib.dopri5', sources=[s.format('dopri5') for s in sources], **extraargs)
        config.add_extension('assimulo.lib.rodas', sources=[s.format('rodas_decsol') for s in sources], include_dirs=[np.get_include()],**extraargs)
        config.add_extension('assimulo.lib.radau5', sources=[s.format('radau_decsol') for s in sources], include_dirs=[np.get_include()],**extraargs)
                             
        radar_list=['contr5.f90', 'radar5_int.f90', 'radar5.f90', 'dontr5.f90', 'decsol.f90', 'dc_decdel.f90', 'radar5.pyf']
        src=['assimulo'+os.sep+'thirdparty'+os.sep+'hairer'+os.sep+code for code in radar_list]
        config.add_extension('assimulo.lib.radar5', sources= src, include_dirs=[np.get_include()],**extraargs)
        
        #ODEPACK
        odepack_list = ['opkdmain.f', 'opkda1.f', 'opkda2.f', 'odepack_aux.f90','odepack.pyf']
        src=['assimulo'+os.sep+'thirdparty'+os.sep+'odepack'+os.sep+code for code in odepack_list]
        config.add_extension('assimulo.lib.odepack', sources= src, include_dirs=[np.get_include()],**extraargs)
     
        #ODASSL
        odassl_list=['odassl.pyf','odassl.f','odastp.f','odacor.f','odajac.f','d1mach.f','daxpy.f','ddanrm.f','ddatrp.f','ddot.f',
                      'ddwats.f','dgefa.f','dgesl.f','dscal.f','idamax.f','xerrwv.f']
        src=['assimulo'+os.sep+'thirdparty'+os.sep+'odassl'+os.sep+code for code in odassl_list]
        config.add_extension('assimulo.lib.odassl', sources= src, include_dirs=[np.get_include()],**extraargs)
    
        dasp3_f77_compile_flags = ["-fdefault-double-8","-fdefault-real-8"]
        dasp3_f77_compile_flags += extra_compile_flags
        
        #NOTE, THERE IS A PROBLEM WITH PASSING F77 COMPILER ARGS FOR NUMPY LESS THAN 1.6.1
        dasp3_list = ['dasp3dp.pyf', 'DASP3.f', 'ANORM.f','CTRACT.f','DECOMP.f', 'HMAX.f','INIVAL.f','JACEST.f','PDERIV.f','PREPOL.f','SOLVE.f','SPAPAT.f']
        src=['assimulo'+os.sep+'thirdparty'+os.sep+'dasp3'+os.sep+code for code in dasp3_list]
        config.add_extension('assimulo.lib.dasp3dp',
                              sources= src,
                              include_dirs=[np.get_include()], extra_link_args=extra_link_flags[:],extra_f77_compile_args=dasp3_f77_compile_flags[:],
                              extra_compile_args=extra_compile_flags[:],extra_f90_compile_args=extra_compile_flags[:])
    
        #GLIMDA
        if self.with_BLAS and self.with_LAPACK:
            glimda_list = ['glimda_complete.f','glimda_complete.pyf']
            src=['assimulo'+os.sep+'thirdparty'+os.sep+'glimda'+os.sep+code for code in glimda_list]
            extraargs_glimda={'extra_link_args':extra_link_flags[:], 'extra_compile_args':extra_compile_flags[:], 'library_dirs':[self.BLASdir, self.LAPACKdir], 'libraries':['lapack', self.BLASname]}
            extraargs_glimda["extra_f77_compile_args"] = extra_compile_flags[:]
            config.add_extension('assimulo.lib.glimda', sources= src,include_dirs=[np.get_include()],**extraargs_glimda) 
            extra_link_flags=extra_link_flags[:-2]  # remove LAPACK flags after GLIMDA 
        else:
            L.warning("Could not find Blas or Lapack, disabling support for the solver GLIMDA.")
        
    
        return config.todict()["ext_modules"]
        

prepare=Assimulo_prepare(args, thirdparty_methods)
curr_dir=os.getcwd()
if not os.path.isdir("assimulo"):
    prepare.create_assimulo_dirs_and_populate()
    os.chdir("build") #Change dir
    change_dir = True
else:
    change_dir = False

ext_list = prepare.cython_extensionlists()

#MAJOR HACK DUE TO NUMPY CHANGE IN VERSION 1.6.2 THAT DOES NOT SEEM TO
#HANDLE EXTENSIONS OF BOTH TYPE (DISTUTILS AND NUMPY DISTUTILS) AT THE
#SAME TIME.
for e in ext_list:
    e.extra_f77_compile_args = []
    e.extra_f90_compile_args = []

ext_list += prepare.fortran_extensionlists()

# distutils part


NAME = "Assimulo"
AUTHOR = u"C. Winther (Andersson), C. Führer, J. Åkesson, M. Gäfvert"
AUTHOR_EMAIL = "christian.winther@modelon.com"
VERSION = "trunk" if version_number_arg == "Default" else version_number_arg
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

Assimulo supports Explicit Euler, adaptive Runge-Kutta of 
order 4 and Runge-Kutta of order 4. It also wraps the popular SUNDIALS 
(https://computation.llnl.gov/casc/sundials/main.html) solvers CVode 
(for ODEs) and IDA (for DAEs). Ernst Hairer's 
(http://www.unige.ch/~hairer/software.html) codes Radau5, Rodas and 
Dopri5 are also available. For the full list, see the documentation.

Documentation and installation instructions can be found at: 
http://www.jmodelica.org/assimulo . 

The package requires Numpy, Scipy and Matplotlib and additionally for 
compiling from source, Cython 0.18, Sundials 2.6/2.7/3.1, BLAS and LAPACK 
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

license_info=[place+os.sep+pck+os.sep+'LICENSE_{}'.format(pck.upper()) 
               for pck in  thirdparty_methods for place in ['thirdparty','lib']]
L.debug(license_info)

import numpy.distutils.core as ndc
ndc.setup(name=NAME,
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
      package_data={'assimulo': ['version.txt', 'CHANGELOG', 'README', 'LICENSE']+license_info+['examples'+os.sep+'kinsol_ors_matrix.mtx',
                                'examples'+os.sep+'kinsol_ors_matrix.mtx'] + (['lib'+os.sep+f for f in prepare.extra_fortran_link_files] if prepare.extra_fortran_link_files else [])},
      script_args=prepare.distutil_args)


if change_dir:
    os.chdir(curr_dir) #Change back to original directory

