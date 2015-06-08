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

from assimulo.ode import *
from assimulo.explicit_ode import Explicit_ODE

from assimulo.exception import *
from assimulo.support import set_type_shape_array

from assimulo.lib import rodas

class Rodas_Common(object):
    
    def _set_atol(self,atol):
        
        self.options["atol"] = set_type_shape_array(atol)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise Rodas_Exception("atol must be of length one or same as the dimension of the problem.")

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
            raise Rodas_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise Rodas_Exception('Relative tolerance must be a positive (scalar) float.')
    
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
            raise Rodas_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def _set_fac1(self, fac1):
        try:
            self.options["fac1"] = float(fac1)
        except (ValueError, TypeError):
            raise Rodas_Exception('The fac1 must be an integer or float.')
            
    def _get_fac1(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac1
                        - Default 0.2
                        
                        - Should be a float.
                        
                            Example:
                                fac1 = 0.1
        """
        return self.options["fac1"]
        
    fac1 = property(_get_fac1, _set_fac1)
    
    def _set_fac2(self, fac2):
        try:
            self.options["fac2"] = float(fac2)
        except (ValueError, TypeError):
            raise Rodas_Exception('The fac2 must be an integer or float.')
        
    def _get_fac2(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac2
                        - Default 8.0
                        
                        - Should be a float.
                        
                            Example:
                                fac2 = 10.0
        """
        return self.options["fac2"]
        
    fac2 = property(_get_fac2, _set_fac2)
    
    def _set_safe(self, safe):
        try:
            self.options["safe"] = float(safe)
        except (ValueError, TypeError):
            raise Rodas_Exception('The safe must be an integer or float.')

    def _get_safe(self):
        """
        The safety factor in the step-size prediction.
        
            Parameters::
            
                safe
                        - Default '0.9'
                        
                        - Should be float.
                        
                            Example:
                                safe = 0.8
        """
        return self.options["safe"]
        
    safe = property(_get_safe, _set_safe)
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Rodas_Exception('The initial step must be an integer or float.')
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                inith    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    inith = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise Rodas_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise Rodas_Exception('Maximal stepsize must be a positiv (scalar) float.')
        
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default maxh = final time - current time.
                          
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
                    
                        - Should be a Boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)


class RodasODE(Rodas_Common, Explicit_ODE):
    """
    Rosenbrock method of order (3)4 with step-size control and 
    continuous output.
    
    Based on the FORTRAN code RODAS by E.Hairer and G.Wanner, which can 
    be found here: http://www.unige.ch/~hairer/software.html
    
    Details about the implementation (FORTRAN) can be found in the book,::
    
        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems
        
        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9
    
    """
    
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["inith"]    = 0.01
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 6.0 #Parameters for step-size selection (upper bound)
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000
        
        #Solver support
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        self.supports["state_events"] = True
        
        #Internal
        self._leny = len(self.y) #Dimension of the problem
        
    def initialize(self):
        #Reset statistics
        self.statistics.reset()
            
    def set_problem_data(self):
        if self.problem_info["state_events"]:
            def event_func(t, y):
                return self.problem.state_events(t, y, self.sw)
            def f(t, y):
                return self.problem.rhs(t, y, self.sw)
            self.f = f
            self.event_func = event_func
            self._event_info = [0] * self.problem_info["dimRoot"]
            self.g_old = self.event_func(self.t, self.y)
        else:
            self.f = self.problem.rhs
    
    def interpolate(self, time):
        y = N.empty(self._leny)
        for i in range(self._leny):
            y[i] = rodas.contro(i+1, time, self.cont)
        
        return y
        
    def _solout(self, nrsol, told, t, y, cont, lrc, irtrn):
        """
        This method is called after every successful step taken by Rodas
        """
        self.cont = cont #Saved to be used by the interpolation function.
        
        if self.problem_info["state_events"]:
            flag, t, y = self.event_locator(told, t, y)
            #Convert to Fortram indicator.
            if flag == ID_PY_EVENT: irtrn = -1
        
        if self._opts["report_continuously"]:
            initialize_flag = self.report_solution(t, y, self._opts)
            if initialize_flag: irtrn = -1
        else:
            if self._opts["output_list"] is None:
                self._tlist.append(t)
                self._ylist.append(y.copy())
            else:
                output_list = self._opts["output_list"]
                output_index = self._opts["output_index"]
                try:
                    while output_list[output_index] <= t:
                        self._tlist.append(output_list[output_index])
                        self._ylist.append(self.interpolate(output_list[output_index]))
                        
                        output_index += 1
                except IndexError:
                    pass
                self._opts["output_index"] = output_index
                
                if self.problem_info["state_events"] and flag == ID_PY_EVENT and len(self._tlist) > 0 and self._tlist[-1] != t:
                    self._tlist.append(t)
                    self._ylist.append(y)
        
        return irtrn
    
    def integrate(self, t, y, tf, opts):
        IFCN  = 1 #The function may depend on t
        ITOL  = 1 #Both rtol and atol are vectors
        IJAC  = 1 if self.usejac else 0 #Switch for the jacobian, 0==NO JACOBIAN
        MLJAC = self.problem_info["dim"] #The jacobian is full
        MUJAC = self.problem_info["dim"] #The jacobian is full
        IDFX  = 0 #df/dt is computed internally
        IMAS  = 0 #The mass matrix is the identity
        MLMAS = self.problem_info["dim"] #The mass matrix is full
        MUMAS = self.problem_info["dim"] #The mass matrix is full
        IOUT  = 1 #Solout is called after every accepted step
        WORK  = N.array([0.0]*(2*self.problem_info["dim"]**2+14*self.problem_info["dim"]+20))
        IWORK = N.array([0]*(self.problem_info["dim"]+20))
        
        #Setting work options
        WORK[1] = self.maxh
        WORK[2] = self.fac1
        WORK[3] = self.fac2
        WORK[4] = self.safe
        
        #Setting iwork options
        IWORK[0] = self.maxsteps
        
        #Dummy methods
        mas_dummy = lambda t:x
        jac_dummy = (lambda t:x) if not self.usejac else self.problem.jac
        dfx_dummy = lambda t:x
        
        #Check for initialization
        if opts["initialize"]:
            self.set_problem_data()
            self._tlist = []
            self._ylist = []
        
        #Store the opts
        self._opts = opts
        
        t, y, h, iwork, flag = rodas.rodas(self.f, IFCN, t, y.copy(), tf, self.inith, self.rtol*N.ones(self.problem_info["dim"]), self.atol,
                    ITOL, jac_dummy, IJAC, MLJAC, MUJAC, dfx_dummy, IDFX, mas_dummy, IMAS, MLMAS, MUMAS, self._solout, IOUT, WORK, IWORK)
                    
        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("Rodas failed with flag %d"%flag)
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[16]
        self.statistics["nfcns"]        += iwork[13]
        self.statistics["njacs"]        += iwork[14]
        #self.statistics["nstepstotal"] += iwork[15]
        self.statistics["nfcnjacs"]    += (iwork[14]*self.problem_info["dim"] if not self.usejac else 0)
        self.statistics["nerrfails"]     += iwork[17]
        self.statistics["nlus"]         += iwork[18]
        
        return flag, self._tlist, self._ylist
        
    def state_event_info(self):
        return self._event_info
        
    def set_event_info(self, event_info):
        self._event_info = event_info
        
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        Explicit_ODE.print_statistics(self, verbose) #Calls the base class
        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : Rodas ',          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)


