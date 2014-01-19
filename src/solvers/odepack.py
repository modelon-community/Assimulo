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
        self.statistics["nevents"]     = 0 #Number of events
        
        self._leny = len(self.y) #Dimension of the problem
        self._nordsieck_array = []
        self._nordsieck_order = 0
        self._nordsieck_time  = 0.0
        self._nordsieck_h  = 0.0
        
        # Solver support
        self.supports["state_events"] = True
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
        #self._tlist = []
        #self._ylist = []
        
    def interpolate(self, t):
        """
        Interpolate the solution at time t using the nordsieck history
        array as,
        
        .. math::
        
             interp = \sum_{j=0}^q (t - tn)^(j) * h^(-j) * zn[j] ,
        
        where q is the current order, and zn[j] is the j-th column of 
        the Nordsieck history array.
        """
        interp = N.array([0.0]*self.problem_info["dim"])
        
        for i in range(self._nordsieck_order):
            interp = interp + (t-self._nordsieck_time)**i*self._nordsieck_h**(-i)*self._nordsieck_array[i*self.problem_info["dim"]:(i+1)*self.problem_info["dim"]]
            
        return interp
    
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
        
        #Nordsieck start index
        nordsieck_start_index = 21+3*self.problem_info["dimRoot"] - 1
        
        #Run in normal mode?
        normal_mode = 1 if opts["output_list"] != None else 0
        
        #if normal_mode == 0:
        if opts["report_continuously"] or opts["output_list"] == None:
            while (ISTATE == 2 or ISTATE == 1) and t < tf:
            
                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tf, ITOL, self.rtol*N.ones(self.problem_info["dim"]), self.atol,
                        ITASK, ISTATE, IOPT, RWORK, IWORK, jac_dummy, JT, g_dummy, JROOT,
                        f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                
                self._nordsieck_array = RWORK[nordsieck_start_index:nordsieck_start_index+(IWORK[14]+1)*self.problem_info["dim"]]
                self._nordsieck_order = IWORK[14]
                self._nordsieck_time  = RWORK[12]
                self._nordsieck_h = RWORK[11]
                
                self._event_info = roots
                
                if opts["report_continuously"]:
                    flag_initialize = self.report_solution(t, y, opts)
                    if flag_initialize:
                        #If a step event has occured the integration has to be reinitialized
                        ISTATE = 3
                else:
                    #Store results
                    tlist.append(t)
                    ylist.append(y.copy())
            
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

                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tout, ITOL, self.rtol*N.ones(self.problem_info["dim"]), self.atol,
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
        if flag == ID_PY_EVENT:
            self.statistics["nevents"] += 1
        
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
        self.log_message(' Number of State-Events                   : '+ str(self.statistics["nevents"]), verbose)
        
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
class RKStarterNordsieck(object):
    """
    A family of Runge-Kutta starters producing a 
    Nordsieck array to (re-)start a Nordsieck based multistep
    method with a given order.
    
    See: Mohammadi (2013): https://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=4196026&fileOId=4196027
    """
    def __init__(self,  rhs, H, eval_at=0.,number_of_steps=4):
        """
        Initiates the Runge-Kutta Starter.
        
            Parameters::
                            
                rhs     
                            - The problem's rhs function
                             
                H
                            - starter step size
                            
                eval_at
                            - Where to evaluate the Nordiseck vector
                              evaluation point is (eval_at)*H
                              Possible choices 
                              eval_at=0.
                              eval_at=1.
             
                number_of_steps
                            - the number of steps :math:`(k)` to start the multistep method with
                              This will determine the initial order of the method, which is
                              :math:`k` for BDF methods and :math:`k+1` for Adams Moulton Methods   .
        """
        self.f = rhs
        self.H = H
        if number_of_steps != 4:
            raise RKStarter_Exception('Step number different from 4 not yet implemented')
        else:
            self.number_of_steps = number_of_steps
            
        self.number_of_steps = number_of_steps 
        self.eval_at = float(eval_at)
        if self.eval_at == 0.:
            self.gamma = \
              N.array([[1.,0.,0.,0.,0.],
                       [0.,1.,-5./6.,4./9.,-1./9.],         
                       [0.,0.,0.,0.,0.],
                       [0.,0.,1./2.,-4./9.,1./9.],
                       [0.,0.,7./3.,-19./9.,7./9.],
                       [0.,0.,-3.,10./3.,-4./3.],
                       [0.,0.,1.,-11./9.,5./9.]])
        elif self.eval_at == 1.:
            self.gamma = \
                N.array([[1., 0. ,0., 0., 0.], 
                         [8./9., 1./27.,1./18.,0.,-1./9.],
                         [0.   , 0.,0.,0.,0.],
                         [4./9., -10./27.,-7./18.,0.,1./9.],
                         [40./9.,  20./27.,1./9.,1.,7./9.],
                         [3./3.,-1./3.,-2.,-4./3.],
                         [8./9.,5./9.,1.,5./9.]])
        else:
            raise RKStarter_Exception("Parameter eval_at should be 0. or 1.")                 
    def rk_like4(self, t0, y0, sw0): 
        """
        rk_like computes Runge-Kutta stages
        Note, the currently implementation is **only** correct for
        autonomous systems.
        """
        f = lambda y: self.f(t0, y, sw0)
        h = self.H/4.
        k1 = h*f(y0)
        k2 = h*f(y0 + k1)
        k3 = h*f(y0 + 2. * k2)
        k4 = h*f(y0 + 3./4. * k1 + 9./4. * k3)
        k5 = h*f(y0 + k1/2. + k2 + k3/2. + 2. * k4)
        k6 = h*f(y0+k1/12.+2. * k2 + k3/4. + 2./3. * k4 + 2. * k5)
        return N.array([y0,k1,k2,k3,k4,k5,k6])
    def nordsieck(self,k):
        """
        Nordsieck array computed at initial point
        """
        nord=N.dot(k.T,self.gamma)
        return nord  
    def __call__(self, t0 , y0, sw0=[]):
        """
        Evaluates the Runge-Kutta starting values
        
            Parameters::
            
                y0   
                    - starting value
        """
        k = self.rk_like4(t0, y0, sw0)
        t = t0+self.eval_at*self.H
        return t, self.nordsieck(k).T
