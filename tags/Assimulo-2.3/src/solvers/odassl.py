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

import numpy as np

from assimulo.ode import *
from assimulo.support import set_type_shape_array
from assimulo.implicit_ode import OverdeterminedDAE

from assimulo.exception import *

from assimulo.lib import odassl

realtype =np.float

class ODASSL_Exception(Exception):
    pass

class ODASSL_Common(object):
    
    def _set_atol(self,atol):
        
        self.options["atol"] = set_type_shape_array(atol) 
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*np.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise ODASSL_Exception("atol must be of length one or same as the dimension of the problem.")

    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        self.options["rtol"] = set_type_shape_array(rtol) 
        if len(self.options["rtol"]) == 1:
            self.options["rtol"] = self.options["rtol"]*np.ones(self._leny)
        elif len(self.options["rtol"]) != self._leny:
            raise ODASSL_Exception("rtol must be of length one or same as the dimension of the problem.")    
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
        """
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)
    
    
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise ODASSL_Exception('The initial step size must be an integer or float.')
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                inith    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    inith = 0.001
            The quantity should be always positive. It will internally be multiplied
            by the sign(tout-t0) to account for the direction of integration.                        
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise ODASSL_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise ODASSL_Exception('Maximal stepsize must be a positiv (scalar) float.')
        
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
          Parameters::
            
                maxh    
                        - Default: maxh=0.0 or None  ignores this option
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["maxh"]
        
    maxh=property(_get_max_h,_set_max_h)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined Jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that Jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined Jacobian
                          False - use an approximation
                          
                        Default:  False  
                    
                        - Should be a Boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)


class ODASSLODE(ODASSL_Common, OverdeterminedDAE):
    """
    Modified version of DASSL for solving overdetermined systems              
    of (singularily) implicit ODEs. The main difference to DASSL         
    is in the corrector iteration part.                                         
                                                                               
    DASSL-Author:  PETZOLD, LINDA                                            
             APPLIED MATHEMATICS DIVISION 8331                                 
             SANDIA NATIONAL LABORATORIES                                      
             LIVERMORE, CA.    94550
    Based on DASSL Version  900103
                                                        
    ODASSL ad-ons : FUEHRER, CLAUS                                         
             DEUTSCHE FORSCHUNGSANSTALT                          
             FUER LUFT- UND RAUMFAHRT (DLR)                                  
             INST. DYNAMIC DER FLUGSYSTEME                                     
             D-8031 WESSLING  (F.R.G)        
    """
    
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
        """
        OverdeterminedDAE.__init__(self, problem) #Calls the base class
            
        #Default values
        self.options["inith"]    = 0.0
        self.options["maxh"]     = 0.0
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6*np.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 5000
        self.options["maxord"]   = 0
        
        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["njac"]        = 0 #Number of jacobian evaluations
        #self.statistics["njacfcn"]     = 0 #Number of function evaluations when evaluating the jacobian
        self.statistics["errfail"]     = 0 #Number of step rejections
        self.statistics["convfail"]         = 0 #Number of LU decompositions
        
        #Internal
        self._leny = len(self.y) #Dimension of the problem
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
        
   
    
    def integrate(self, t, y, yprime, tf, opts):
        ny  = self.problem_info["dim"]
        neq = self.problem_info["neq"]
        lrw = 40+8*ny + neq**2 + 3*neq
        rwork = np.zeros((lrw,))                                                
        liw = 22+neq
        iwork = np.zeros((liw,),np.int)                                                                          
        jac_dummy = lambda t,x,xp: x
        info = np.zeros((15,),np.int) 
        info[1] = 1  # Tolerances are vectors  
        info[2] = normal_mode = 0 if opts["output_list"] != None else 1  # intermediate output mode
        info[6] = 1 if self.options["maxh"] > 0.0 else 0       
        rwork[1] = self.options["maxh"]            
        info[7] = 1 if self.options["inith"] > 0.0 else 0
        rwork[2] = self.options["inith"]
        info[8] =  1 if self.options["maxord"] > 0 else 0  
        iwork[2] = self.options["maxord"]                     
        #info[11] will be set later (see Ticket #230)                        
        #iwork[0] = ML
        #iwork[1] = MU                             
        atol = self.options["atol"]
        rtol = self.options["rtol"]
        for i in range(ny):
            if self.problem.algvar[i] == 0:
                rtol[i]=1.e7
                atol[i]=1.e7
        tlist=[]
        ylist=[]
        ydlist=[]
        
        
        #Store the opts
        self._opts = opts
        
        
        #THIS IS NOT SUPPOSE TO BE NECCESSARY, BUT DUE TO A PROBLEM
        #REPORTED IN TICKET:244 THIS IS HOWEVER NECCESSARY AS A 
        #WORKAROUND FOR NOW...
        def py_residual(t,y,yd):
            return self.problem.res(t,y,yd)
        callback_residual = py_residual
        #----

        if normal_mode == 1: # intermediate output mode
            idid = 1
            while idid==1:
                t,y,yprime,tf,info,idid,rwork,iwork = \
                   odassl.odassl(callback_residual,neq,ny,t,y,yprime,
                         tf,info,rtol,atol,rwork,iwork,jac_dummy)
                tlist.append(t)
                ylist.append(y.copy())
                ydlist.append(yprime.copy())
                if idid==2:
                    flag=ID_PY_COMPLETE
                elif idid < 0:
                    raise ODASSL_Exception("ODASSL failed with flag IDID {IDID}".format(IDID=idid))
        else:   # mode with output_list          
            output_list  = opts["output_list"]
            for tout in output_list: 
                t,y,yprime,tout,info,idid,rwork,iwork = \
                  odassl.odassl(callback_residual,neq,ny,t,y,yprime, \
                         tout,info,rtol,atol,rwork,iwork,jac_dummy)
                tlist.append(t)
                ylist.append(y.copy())
                ydlist.append(yprime.copy())
                if idid > 0 and t >= tf:
                    flag=ID_PY_COMPLETE
                elif idid < 0:
                    raise ODASSL_Exception("ODASSL failed with flag IDID {IDID}".format(IDID=idid))                                      
 
        
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[10]
        self.statistics["nfcn"]        += iwork[11]
        self.statistics["njac"]        += iwork[12]
        self.statistics["errfail"]     += iwork[13]
        self.statistics["convfail"]         += iwork[14]
        
        return flag, tlist, ylist, ydlist
        
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]), verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfcn"]), verbose)
        self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njac"]), verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["errfail"]), verbose)
        self.log_message(' Number of Convergence Test Failures      : '+ str(self.statistics["convfail"]), verbose)
        
        self.log_message('\nSolver options:\n', verbose)
        self.log_message(' Solver                  : ODASSL ',          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)


