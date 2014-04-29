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

cimport numpy as N
import numpy as N
import numpy.linalg as LIN

#from assimulo.ode import *
from assimulo.explicit_ode cimport Explicit_ODE
from assimulo.exception import *

include "constants.pxi" #Includes the constants (textual include)

cdef class ImplicitEuler(Explicit_ODE):
    """
    This solver solves an explicit ordinary differential equation using 
    the implicit Euler method.
    
    We want to approximate the solution to the ordinary differential 
    equation of the form,
    
    .. math::

        \dot{y} = f(t,y), \quad y(t_0) = y_0 .
        
    Using the implicit Euler method, the approximation is defined as 
    follow,
    
    .. math::
    
        y_{n+1} = y_n + hf(t_{n+1},y_{n+1})
        
    with :math:`h` being the step-size and :math:`y_n` the previous 
    solution to the equation.
    """
    cdef N.ndarray yd1
    cdef N.ndarray _old_jac
    cdef object f
    cdef public object event_func
    cdef int _leny
    cdef double _eps
    cdef int _needjac
    cdef int _curjac
    cdef int _steps_since_last_jac
    cdef N.ndarray _yold
    cdef N.ndarray _ynew
    #cdef N.ndarray _event_info
    cdef public N.ndarray g_old
    cdef double _told
    cdef double _h
    cdef double _inith
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01
        self.options["usejac"]   = True if (self.problem_info["jac_fcn"]) else False
        self.options["newt"]     = 7 #Maximum number of newton iterations
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        
        #Internal temporary result vector
        self.yd1 = N.array([0.0]*len(self.y0))
        self._yold = N.array([0.0]*len(self.y0))
        self._ynew = N.array([0.0]*len(self.y0))
        
        #Solver support
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = False
        self.supports["state_events"] = True

        
        self._leny = len(self.y) #Dimension of the problem
        self._eps  = N.finfo('double').eps
        self._needjac = True #Do we need a new jacobian?
        self._curjac = False #Is the current jacobian up to date?
        self._steps_since_last_jac = 0 #Keep track on how long ago we updated the jacobian
        self._inith = 0 #Used for taking an initial step of correct length after an event.
    
    def set_problem_data(self): 
        if self.problem_info["state_events"]: 
            def event_func(t, y): 
                return self.problem.state_events(t, y, self.sw) 
            def f(t, y): 
                return self.problem.rhs(t, y, self.sw)
            self.f = f
            self.event_func = event_func
            self._event_info = N.array([0] * self.problem_info["dimRoot"]) 
            self.g_old = self.event_func(self.t, self.y)
        else: 
            self.f = self.problem.rhs
    
    
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
    
    cpdef step(self,double t,N.ndarray y,double tf,dict opts):
        cdef double h
        h = self.options["h"]
        
        if t+h < tf:
            t, y = self._step(t,y,h)
            return ID_OK, t, y
        else:
            h = min(h, abs(tf-t))
            t, y = self._step(t,y,h)
            return ID_COMPLETE, t, y
    
    cpdef integrate(self, double t,N.ndarray y,double tf, dict opts):
        cdef double h
        cdef list tr,yr
        
        h = self.options["h"]
        h = min(h, abs(tf-t))
        
        if opts["initialize"]:
            self.set_problem_data()
            tr = []
            yr = []
        
        flag = ID_PY_OK
        if self._inith != 0:
            t, y = self._step(t,y,self._inith)
            if self.problem_info["state_events"]: 
                flag, t, y = self.event_locator(t-self._inith , t, y)
            
            if opts["report_continuously"]:
                initialize_flag = self.report_solution(t, y, opts)
                if initialize_flag: flag = ID_PY_EVENT
                
            else:
                tr.append(t)
                yr.append(y)
            
            #If an event was detected, calculate the length of the initial step to take after restarting.
            if flag == ID_PY_EVENT:
                self._inith = self._inith - (t - self._told)
            else:
                self._inith = 0
            h = min(h, abs(tf-t))
        
        while t+h < tf and flag == ID_PY_OK:
            t, y = self._step(t,y,h)
            if self.problem_info["state_events"]: 
                    flag, t, y = self.event_locator(t-h , t, y)
            
            if opts["report_continuously"]:
                initialize_flag = self.report_solution(t, y, opts)
                if initialize_flag: flag = ID_PY_EVENT
            else:
                tr.append(t)
                yr.append(y)

            h = min(h, abs(tf-t))
        else:
            #If no event was detected take the final step to tf.
            if flag == ID_PY_OK:
                t, y = self._step(t, y, h)
                flag = ID_PY_COMPLETE
                if self.problem_info["state_events"]:
                    flag, t, y = self.event_locator(t-h , t, y)
                    if flag == ID_PY_OK: flag = ID_PY_COMPLETE
                
                if opts["report_continuously"]:
                    initialize_flag = self.report_solution(t, y, opts)
                    if initialize_flag: flag = ID_PY_EVENT
                else:
                    tr.append(t)
                    yr.append(y)
                    flag = ID_PY_COMPLETE
        
        #If an event was detected calculate the length of the initial step if not already done.
        if flag == ID_PY_EVENT and self._inith == 0:
            self._inith = h - (t - self._told)
                
        return flag, tr, yr
    
    def _set_newt(self, newt):
        """
        Maximal number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        try:
            self.options["newt"] = int(newt)
        except (ValueError, TypeError):
            raise AssimuloException('The newt must be an integer or float.')
        
    def _get_newt(self):
        """
        Maximal number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        return self.options["newt"]
        
    newt = property(_get_newt,_set_newt)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise AssimuloException("atol must be of length one or same as the dimension of the problem.")

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
            raise AssimuloException('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise AssimuloException('Relative tolerance must be a positive (scalar) float.')
    
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
    
    cpdef initialize(self):
        #Reset statistics
        #for k in self.statistics.keys():
        #    self.statistics[k] = 0
        self.statistics.reset()
    
    def _jacobian(self, t, y):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        self._curjac = True #The jacobian is up to date
        self._needjac = False #A new jacobian is not needed
        self._steps_since_last_jac = 0
        
        if self.usejac: #Retrieve the user-defined jacobian
            jac = self.problem.jac(t,y)
        else:           #Calculate a numeric jacobian
            delt = N.array([(self._eps*max(abs(yi),1.e-5))**0.5 for yi in y])*N.identity(self._leny) #Calculate a disturbance
            Fdelt = N.array([self.f(t,y+e) for e in delt]) #Add the disturbance (row by row) 
            grad = ((Fdelt-self.f(t,y)).T/delt.diagonal()).T
            jac = N.array(grad).T
            
            self.statistics["nfcnjacs"] += 1+self._leny #Add the number of function evaluations
        
        self.statistics["njacs"] += 1 #add the number of jacobian evaluation
        return jac
    
    
    
    cdef double WRMS(self, N.ndarray x, N.ndarray w):
        """
        Calculates the Weighted Root-mean-square.
        """
        cdef double prod
        cdef int N = self._leny
        cdef double sum = 0.0
        
        for i in range(N):
            prod = x[i]*w[i]
            sum += prod*prod
            
        return (sum/N)**0.5
    
    cdef tuple _step(self,double t,N.ndarray y,double h):
        """
        This calculates the next step in the integration.
        """
        cdef double new_norm, old_norm
        cdef double tn1 = t+h
        cdef N.ndarray yn = y.copy() #Old y
        #cdef N.ndarray yn1 = y.copy() #First newton guess
        cdef N.ndarray yn1 = y+h*self.f(t,y) #First newton guess
        cdef N.ndarray I = N.eye(self._leny)
        self.statistics["nfcns"] += 1
        
        FLAG_CONV = False
        
        for j in range(2):
            FLAG_FAIL = False

            if self._needjac: #If the jacobian should be updated or not
                jac = self._jacobian(tn1, yn1)
            else:
                jac = self._old_jac
            
            for i in range(self.newt):
                self.statistics["nniters"] += 1
                
                #jac = self._jacobian(tn1, yn1)
                
                #ynew = yn1 - N.dot(LIN.inv(h*jac-I),(yn-yn1+h*self.problem.rhs(tn1,yn1)))
                ynew = yn1 - LIN.solve(h*jac-I, yn-yn1+h*self.f(tn1,yn1) )
                self.statistics["nfcns"] += 1
                
                #print tn1, self.WRMS(ynew-yn1, 1.0/(self.rtol*N.abs(yn1)+self.atol))
                new_norm = self.WRMS(ynew-yn1, 1.0/(self.rtol*N.abs(yn1)+self.atol))
                
                if new_norm < 0.1: #Newton converged
                    FLAG_CONV = True
                    break
                    
                if i > 0:
                    if new_norm/old_norm > 2: #Newton iterations diverging
                        FLAG_FAIL = True
                        break
                yn1 = ynew
                old_norm = new_norm
            else:
                FLAG_FAIL = True
            
            if FLAG_FAIL:
                self.statistics["nnfails"] += 1
                
                if self._curjac: #The current Jacobian is updated (Newton failed)
                    raise AssimuloException("Newton iteration failed at %f"%t)
                else: #Try with updated jacobian
                    self._needjac = True
            
            if FLAG_CONV:
                self._steps_since_last_jac += 1
                self.statistics["nsteps"] += 1
                
                if self._steps_since_last_jac > 20: #Need to update the jacobian
                    self._needjac = True 
                break
        else:
            raise AssimuloException("Newton iteration failed at %f"%t)
                
        self._curjac = False #The Jacobian is no longer current
        self._old_jac = jac #Store the old jacobian
        
        #Internal values only used for defining the interpolation function.
        self._yold = yn
        self._ynew = ynew.copy()
        self._told = t
        self._h = h
        
        return tn1, ynew
        
    def interpolate(self, time):
        return self._yold + (time - self._told) / self._h * (self._ynew - self._yold)
        
    def _set_h(self,h):
        try:
            self.options["h"] = float(h)
        except:
            raise AssimuloException("Step-size must be a (scalar) float.")
    
    def _get_h(self):
        """
        Defines the step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0.01'.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["h"]
        
    h=property(_get_h,_set_h)
    
    def state_event_info(self):
        return self._event_info
    
    def set_event_info(self, event_info):
        self._event_info = event_info
        
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        Explicit_ODE.print_statistics(self, verbose) #Calls the base class

        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : ImplicitEuler',                     verbose)
        self.log_message(' Solver type       : fixed step size',                      verbose)
        self.log_message(' Step size         : {}\n'.format(self.h), verbose)

cdef class ExplicitEuler(Explicit_ODE):
    """
    This solver solves an explicit ordinary differential equation using 
    the explicit Euler method.
    
    We want to approximate the solution to the ordinary differential 
    equation of the form,
    
    .. math::

        \dot{y} = f(t,y), \quad y(t_0) = y_0 .
        
    Using the explicit Euler method, the approximation is defined as 
    follow,
    
    .. math::
    
        y_{n+1} = y_n + hf(t_n,y_n)
        
    with :math:`h` being the step-size and :math:`y_n` the previous 
    solution to the equation.
    """
    cdef N.ndarray yd1
    cdef object f
    cdef public object event_func
    cdef N.ndarray _yold
    cdef N.ndarray _ynew
    #cdef N.ndarray _event_info
    cdef public N.ndarray g_old
    cdef double _told
    cdef double _h
    cdef double _inith
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01

        
        
        #Internal temporary result vector
        self.yd1 = N.array([0.0]*len(self.y0))
        self._yold = N.array([0.0]*len(self.y0))
        self._ynew = N.array([0.0]*len(self.y0))
        self._inith = 0 #Used for taking an initial step of correct length after an event.
        
        #Solver support
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = False
        self.supports["state_events"] = True
    
    def set_problem_data(self): 
        if self.problem_info["state_events"]: 
            def event_func(t, y): 
                return self.problem.state_events(t, y, self.sw) 
            def f(t, y): 
                return self.problem.rhs(t, y, self.sw)
            self.f = f
            self.event_func = event_func
            self._event_info = N.array([0] * self.problem_info["dimRoot"]) 
            self.g_old = self.event_func(self.t, self.y)
        else: 
            self.f = self.problem.rhs
    
    cpdef step(self,double t,N.ndarray y,double tf,dict opts):
        cdef double h
        h = self.options["h"]
        
        if t+h < tf:
            t, y = self._step(t,y,h)
            return ID_OK, t, y
        else:
            h = min(h, abs(tf-t))
            t, y = self._step(t,y,h)
            return ID_COMPLETE, t, y
    
    cpdef integrate(self, double t,N.ndarray y,double tf, dict opts):
        cdef double h
        cdef list tr,yr
        
        h = self.options["h"]
        h = min(h, abs(tf-t))
        
        if opts["initialize"]:
            self.set_problem_data()
            tr = []
            yr = []
        
        flag = ID_PY_OK
        if self._inith != 0:
            t, y = self._step(t,y,self._inith)
            if self.problem_info["state_events"]: 
                flag, t, y = self.event_locator(t-self._inith , t, y)
            
            if opts["report_continuously"]:
                initialize_flag = self.report_solution(t, y, opts)
                if initialize_flag: flag = ID_PY_EVENT
                
            else:
                tr.append(t)
                yr.append(y.copy())
            
            #If an event was detected calculate the length of the initial step to take after restarting.
            if flag == ID_PY_EVENT:
                self._inith = self._inith - (t - self._told)
            else:
                self._inith = 0
            h = min(h, abs(tf-t))
        
        while t+h < tf and flag == ID_PY_OK:
            t, y = self._step(t,y,h)
            if self.problem_info["state_events"]: 
                    flag, t, y = self.event_locator(t-h , t, y)
            
            if opts["report_continuously"]:
                initialize_flag = self.report_solution(t, y, opts)
                if initialize_flag: flag = ID_PY_EVENT
            else:
                tr.append(t)
                yr.append(y.copy())
            
            h = min(h, abs(tf-t))
        else:
            #If no event was detected take the final step to tf.
            if flag == ID_PY_OK:
                t, y = self._step(t, y, h)
                flag = ID_PY_COMPLETE
                if self.problem_info["state_events"]:
                    flag, t, y = self.event_locator(t-h , t, y)
                    if flag == ID_PY_OK: flag = ID_PY_COMPLETE
                if opts["report_continuously"]:
                    initialize_flag = self.report_solution(t, y, opts)
                    if initialize_flag: flag = ID_PY_EVENT
                else:
                    tr.append(t)
                    yr.append(y.copy())
            
            #If an event was detected calculate the length of the initial step if not already done.
            if flag == ID_PY_EVENT and self._inith == 0:
                self._inith = h - (t - self._told)
                self._inith = self._inith if abs(self._inith) > 1e-15 else 0.0
            
        return flag, tr, yr
    
    cdef tuple _step(self,double t,N.ndarray y,double h):
        """
        This calculates the next step in the integration.
        """
        #self.f(self.yd1,t,y) #The output is stored in yd
        
        #Internal values only used for defining the interpolation function.
        self._yold = y.copy()
        self._ynew = y + h*self.f(t,y)
        self._told = t
        self._h = h
        return t + h, self._ynew
        
    def interpolate(self, time):
        return self._yold + (time - self._told) / self._h * (self._ynew - self._yold)
        
    def _set_h(self,h):
        try:
            self.options["h"] = float(h)
        except:
            raise AssimuloException("Step-size must be a (scalar) float.")
    
    def _get_h(self):
        """
        Defines the step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0.01'.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["h"]
        
    h=property(_get_h,_set_h)
    
    def state_event_info(self):
        return self._event_info
    
    def set_event_info(self, event_info):
        self._event_info = event_info
        
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        Explicit_ODE.print_statistics(self, verbose) #Calls the base class
            
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : ExplicitEuler',                     verbose)
        self.log_message(' Solver type       : fixed step size',                 verbose)
        self.log_message(' Step size         : {}\n'.format(self.h), verbose)
