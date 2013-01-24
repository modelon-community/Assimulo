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

import numpy as N

from assimulo.exception import *
from assimulo.ode import *

from assimulo.explicit_ode import Explicit_ODE

try:
    from assimulo.lib.odepack import dlsodar
except ImportError:
    print "Could not find ODEPACK"

class LSODAR(Explicit_ODE):
    """
        LOSDAR is a multistep method for solving explicit ordinary 
        differential equations on the form,
        
        .. math::
    
            \dot{y} = f(t,y), \quad y(t_0) = y_0.
            
        LSODAR automatically switches between using an ADAMS method
        or an BDF method and is also able to monitor events.
        
        LSODAR is part of ODEPACK, http://www.netlib.org/odepack/opkd-sum
    """
    
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = False
        self.options["maxsteps"] = 100000

        
        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["njac"]        = 0 #Number of jacobian evaluations
        self.statistics["ng"]          = 0 #Number of root evaluations
        
        self._leny = len(self.y) #Dimension of the problem
        
        # Solver support
        self.supports["state_events"] = True
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
        #self._tlist = []
        #self._ylist = []
    
    def integrate(self, t, y, tf, opts):
        ITOL  = 2 #Both atol and rtol are vectors
        ITASK = 5 #For computation of yout
        ISTATE = 1 #Start of integration
        IOPT = 1 #No optional inputs are used
        RWORK = N.array([0.0]*(22 + self.problem_info["dim"]*max(16,self.problem_info["dim"]+9) + 3*self.problem_info["dimRoot"]))#Real work array
        IWORK = N.array([0]*(20 + self.problem_info["dim"]))
        JT = 1 if self.usejac else 2#Jacobian type indicator
        JROOT = N.array([0]*self.problem_info["dimRoot"])
        
        #Setting work options
        RWORK[0] = tf #Do not integrate past tf
        
        #Setting iwork options
        IWORK[5] = self.maxsteps

        
        #Dummy methods
        g_dummy = (lambda t:x) if not self.problem_info["state_events"] else self.problem.state_events
        jac_dummy = (lambda t,y:N.zeros((len(y),len(y)))) if not self.usejac else self.problem.jac
        
        #Extra args to rhs and state_events
        rhs_extra_args = (self.sw,) if self.problem_info["switches"] else ()
        g_extra_args = (self.sw,) if self.problem_info["switches"] else ()
        
        #Store the opts
        self._opts = opts
        
        #Outputs
        tlist = []
        ylist = []
        
        #Run in normal mode?
        normal_mode = 1 if opts["output_list"] != None else 0
        
        if normal_mode == 0: 
            while (ISTATE == 2 or ISTATE == 1) and t < tf:
            
                y, t, ISTATE, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tf, ITOL, self.rtol*N.ones(self.problem_info["dim"]), self.atol,
                        ITASK, ISTATE, IOPT, RWORK, IWORK, jac_dummy, JT, g_dummy, JROOT,
                        f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                tlist.append(t)
                ylist.append(y.copy())
                self._event_info = roots
            
                #Checking return
                if ISTATE == 2:
                    flag = ID_PY_COMPLETE
                elif ISTATE == 3:
                    flag = ID_PY_EVENT
                else:
                    raise ODEPACK_Exception("LSODAR failed with flag %d"%ISTATE)
            
        else:
            #Change the ITASK
            ITASK = 4 #For computation of yout
            
            output_index = opts["output_index"]
            output_list  = opts["output_list"][output_index:]

            for tout in output_list:
                output_index += 1

                y, t, ISTATE, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tout, ITOL, self.rtol*N.ones(self.problem_info["dim"]), self.atol,
                    ITASK, ISTATE, IOPT, RWORK, IWORK, jac_dummy, JT, g_dummy, JROOT,
                    f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                #Store results
                tlist.append(t)
                ylist.append(y.copy())
                self._event_info = roots
                
                #Checking return
                if ISTATE == 2 and t >= tf:
                    flag = ID_PY_COMPLETE
                    break
                elif ISTATE == 3:
                    flag = ID_PY_EVENT
                    break
                elif ISTATE < 0:
                    raise ODEPACK_Exception("LSODAR failed with flag %d"%ISTATE)
        
            opts["output_index"] = output_index
        
        #Retrieving statistics
        self.statistics["ng"]            += IWORK[9]
        self.statistics["nsteps"]        += IWORK[10]
        self.statistics["nfcn"]          += IWORK[11]
        self.statistics["njac"]          += IWORK[12]
        
        return flag, tlist, ylist
    
    def state_event_info(self):
        """
        Returns the state events indicator as a list where a one (1)
        indicates that the event function have been activated and a (0)
        if not.
        """
        return self._event_info

    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njac"]),    verbose)
        self.log_message(' Number of Root Evaluations               : '+ str(self.statistics["ng"]),       verbose)
        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : LSODAR ',         verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise ODEPACK_Exception("atol must be of length one or same as the dimension of the problem.")

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
        try:
            self.options["rtol"] = float(rtol)
        except (ValueError, TypeError):
            raise ODEPACK_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise ODEPACK_Exception('Relative tolerance must be a positive (scalar) float.')
    
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

    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 10000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
