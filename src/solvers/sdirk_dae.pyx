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

from assimulo.explicit_ode import Implicit_ODE

try:
    from assimulo.lib.odepack import dlsodar, dcfode, dintdy
    from assimulo.lib.odepack import set_lsod_common, get_lsod_common
except ImportError:
    print "Could not find ODEPACK functions"

class SDIRK_DAE(Implicit_ODE):
    """
        SDIRK_DAE is a special implicit Runge-Kutta methodstep for solving implicit ordinary 
        differential equations of index index on the form,
        
        .. math::
    
            F(t, y, \dot{y}) = 0, \quad y(t_0) = y_0.
            
        
        
        The core integrator SDIRK_DAE is the original code by Anne Kvaernoe, 
        http://www.math.ntnu.no/~anne/
    """

    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
        """
        Implicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = False
        self.options["hmax"] = 0.
                
        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["njac"]        = 0 #Number of Jacobian evaluations
        
        self._leny = len(self.y) #Dimension of the problem
        
        # Solver support
        self.supports["state_events"] = False
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = False
        
    def initialize(self):
        """
        Initializes the overall simulation process
        (called before _simulate) 
        """ 
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
        #self._tlist = []
        #self._ylist = []
        
        
    def integrate_start(self, t, y):
        
        # first call or classical restart after a discontinuity
        ISTATE=1
        RWORK = N.array([0.0]*(22 + self.problem_info["dim"] * 
                               max(16,self.problem_info["dim"]+9) + 
                               3*self.problem_info["dimRoot"]))
        # Integer work array
        IWORK = N.array([0]*(20 + self.problem_info["dim"]))

 
        return ISTATE, RWORK, IWORK
                                     
    
    def integrate(self, t, y, tf, opts):
        ITOL  = 2 #Only  atol is a  vector
        ITASK = 5 #For one step mode and hitting exactly tcrit, normally tf
        IOPT = 1 #optional inputs are used
        
        # provide work arrays and set common blocks (if needed)
        ISTATE, RWORK, IWORK = self.integrate_start( t, y)
        
        JT = 1 if self.usejac else 2#Jacobian type indicator
        JROOT = N.array([0]*self.problem_info["dimRoot"])
        
        #Setting work options
        RWORK[0] = tf #Do not integrate past tf
        RWORK[5] = self.options["hmax"]
        
        #Setting iwork options
        IWORK[5] = self.maxsteps
        
        #Setting maxord to IWORK
        IWORK[7] = self.maxordn
        IWORK[8] = self.maxords

        
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
            
                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tf, ITOL, 
                        self.rtol*N.ones(self.problem_info["dim"]), self.atol,
                        ITASK, ISTATE, IOPT, RWORK, IWORK, jac_dummy, JT, g_dummy, JROOT,
                        f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
                #print 't= {}, tN={}, y={}, ns={}, hu={}'.format(t , RWORK[12], y, RWORK[nordsieck_start_index],hu)
                self._nordsieck_array = \
                     RWORK[nordsieck_start_index:nordsieck_start_index+(nq+1)*nyh].reshape((nyh,-1),order='F') 
                self._nyh = nyh              
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

                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(self.problem.rhs, y.copy(), t, tout, ITOL, 
                    self.rtol*N.ones(self.problem_info["dim"]), self.atol,
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
        # deciding on restarting options
        self._rkstarter_active = True if ISTATE == 3 and self.rkstarter > 1 else False
        #print 'rkstarter_active set to {} and ISTATE={}'.format(self._rkstarter_active, ISTATE)
        
        #Retrieving statistics
        self.statistics["ng"]            += IWORK[9]
        self.statistics["nsteps"]        += IWORK[10]
        self.statistics["nfcn"]          += IWORK[11]
        self.statistics["njac"]          += IWORK[12]
        self.statistics["nevents"] += 1  if flag == ID_PY_EVENT else 0
        # save RWORK, IWORK for restarting feature
        if self.rkstarter>1:
            self._RWORK=RWORK
            self._IWORK=IWORK        
       
        return flag, tlist, ylist
        
    def get_algorithm_data(self):
        """
        Returns the order and step size used in the last successful step.
        """
        hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
            
        return hu, nqu
    
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
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : LSODAR ',         verbose)
        self.log_message(' Absolute tolerances     : {}'.format(self.options["atol"]),  verbose)
        self.log_message(' Relative tolerances     : {}'.format(self.options["rtol"]),  verbose)
        if self.rkstarter==1:
           self.log_message(' Starter                 : {}'.format('classical'),  verbose)
        else:
           self.log_message(' Starter                 : Runge-Kutta order {}'.format(self.rkstarter),  verbose)
        if self.maxordn < 12 or self.maxords < 5:
            self.log_message(' Maximal order Adams     : {}'.format(self.maxordn),  verbose)
            self.log_message(' Maximal order BDF       : {}'.format(self.maxords),  verbose)
        if self.hmax > 0. :
            self.log_message(' Maximal stepsize hmax   : {}'.format(self.hmax),  verbose)
        self.log_message('',                                                         verbose)

        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of steps                          : {}'.format(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of function evaluations           : {}'.format(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Jacobian evaluations           : {}'.format(self.statistics["njac"]),    verbose)
        self.log_message(' Number of event function evaluations     : {}'.format(self.statistics["ng"]),       verbose)
        self.log_message(' Number of detected state events          : {}'.format(self.statistics["nevents"]), verbose)
        
    
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
    def _get_hmax(self):
        """
        The absolute value of the maximal stepsize for all methods
        
        Parameters::
        
               hmax
                          - Default:  0.  (no maximal step size)
                          
                          - Should be a positive float
        """
        return self.options["hmax"]
    def _set_hmax(self,hmax):
        if not (isinstance(hmax,float) and hmax >= 0.):
           raise ODEPACK_Exception("Maximal step size hmax should be a positive float")
        self.options["hmax"]=hmax
    hmax = property(_get_hmax, _set_hmax)
    def _get_maxordn(self):
        """
        The maximum order used by the Adams-Moulton method (nonstiff case)
        
            Parameters::
            
                maxordn
                            - Default 12
                            
                            - Should be a positive integer
        """
        return self.options["maxordn"]
    
    def _set_maxordn(self, maxordn):
        try:
            maxordn = int(maxordn)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum order must be a positive integer.")
        if maxordn > 12:
            raise ODEPACK_Exception("Maximum order should not exceed 12.")    
        self.options["maxordn"] = maxordn
    
    maxordn = property(_get_maxordn, _set_maxordn)
    def _get_maxords(self):
        """
        The maximum order used by the BDF method (stiff case)
        
            Parameters::
            
                maxords
                            - Default 5
                            
                            - Should be a positive integer
        """
        return self.options["maxords"]
    
    def _set_maxords(self, maxords):
        try:
            maxords = int(maxords)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum order must be a positive integer.")
        if maxords > 5:
            raise ODEPACK_Exception("Maximum order should not exceed 5.")    
        self.options["maxords"] = maxords
    
    maxords = property(_get_maxords, _set_maxords)
    def _get_rkstarter(self):
        """
        This defines how LSODAR is started. 
        (classically or with a fourth order Runge-Kutta starter)
        
            Parameters::
            
                rkstarter
                            - Default False  starts LSODAR in the classical multistep way
                            
                            - Should be a Boolean
        """
        return self.options["rkstarter"]
    
    def _set_rkstarter(self, rkstarter):
        if not rkstarter in {1,2,3,4,5}:
            raise ODEPACK_Exception("Must be a positive integer less than 6")
        self.options["rkstarter"] = rkstarter
    
    rkstarter = property(_get_rkstarter, _set_rkstarter)
class RKStarterNordsieck(object):
    """
    A family of Runge-Kutta starters producing a 
    Nordsieck array to (re-)start a Nordsieck based multistep
    method with a given order.
    
    See: Mohammadi (2013): https://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=4196026&fileOId=4196027
    """
    Gamma_0=[N.array([[1.,0.],                        # 1st order
                      [0.,1.]]),
             N.array([[1.,0.,0.],                     # 2nd order
                      [0.,1.,-1.],
                      [0.,0.,1.]]),
             N.array([[1.,0.,0.,0.],                  # 3rd order
                      [0.,1.,-5./3.,1.],         
                      [0.,0.,3.,-2.],
                      [0.,0.,0.,1.],
                      [0.,0.,-4./3.,0.]]),
             N.array([[1.,0.,0.,0.,0.],               # 4th order
                      [0.,1.,-5./6.,4./9.,-1./9.],         
                      [0.,0.,0.,0.,0.],
                      [0.,0.,1./2.,-4./9.,1./9.],
                      [0.,0.,7./3.,-19./9.,7./9.],
                      [0.,0.,-3.,10./3.,-4./3.],
                      [0.,0.,1.,-11./9.,5./9.]])]  
                      
    scale=N.array([1, 1, 1/2., 1./6., 1./24., 1./120.]).reshape(-1,1)
    
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
        
        if not 1 < number_of_steps < 5:
            raise RKStarter_Exception('Step number larget than 4 not yet implemented')            
        self.number_of_steps = number_of_steps 
        self.eval_at = float(eval_at)
        if not self.eval_at == 0.:
           raise RKStarter_Exception("Parameter eval_at different from 0 not yet implemented.")                    
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
    def rk_like3(self, t0, y0, sw0): 
        """
        rk_like computes Runge-Kutta stages
        Note, the currently implementation is **only** correct for
        autonomous systems.

        """
        f = lambda y: self.f(t0, y, sw0)
        h = self.H/3.
        k1 = h*f(y0)
        k2 = h*f(y0 + k1)
        k3 = h*f(y0 + k1+ k2)
        k4 = h*f(y0 + 3./2. * k1)
        return N.array([y0,k1,k2,k3,k4])
    def rk_like2(self, t0, y0, sw0):
        """
        rk_like2 computes Runge-Kutta 2nd-stages
        Note, the currently implementation is **only** correct for
        autonomous systems.
        """
        f=lambda y: self.f(t0, y, sw0)
        h=self.H/2
        k1=h*f(y0)
        k2=h*f(y0+k1)
        return N.array([y0,k1,k2])        
    def nordsieck(self,k):
        """
        Nordsieck array computed at initial point
        """
        nord=self.scale[:self.number_of_steps+1]*N.dot(self.Gamma_0[self.number_of_steps-1].T,k)
        return nord  
    def __call__(self, t0 , y0, sw0=[]):
        """
        Evaluates the Runge-Kutta starting values
        
            Parameters::
            
                y0   
                    - starting value
        """
        # We construct a call like: rk_like4(self, t0, y0, sw0) 
        k=self.__getattribute__('rk_like{}'.format(self.number_of_steps))(t0, y0, sw0)
        t = t0+self.eval_at*self.H
        return t, self.nordsieck(k)
