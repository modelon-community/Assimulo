#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2025 Modelon AB
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

# Testing utilities shared in a variety of solver tests

import numpy as np

from assimulo.problem import Explicit_Problem, Implicit_Problem

class Extended_Problem(Explicit_Problem):
    
    # Sets the initial conditions directly into the problem
    y0 = [0.0, -1.0, 0.0]
    sw0 = [False,True,True]
    event_array = np.array([0.0,0.0,0.0])
    rhs_array   = np.array([0.0,0.0,0.0])
    
    #The right-hand-side function (rhs)
    def rhs(self,t,y,sw):
        """
        This is our function we are trying to simulate. During simulation
        the parameter sw should be fixed so that our function is continuous
        over the interval. The parameters sw should only be changed when the
        integrator has stopped.
        """
        self.rhs_array[0] = (1.0 if sw[0] else -1.0)
        self.rhs_array[1] = 0.0
        self.rhs_array[2] = 0.0

        return self.rhs_array

    #Sets a name to our function
    name = 'ODE with discontinuities and a function with consistency problem'
    
    #The event function
    def state_events(self,t,y,sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        self.event_array[0] = y[1] - 1.0 
        self.event_array[1] = -y[2] + 1.0
        self.event_array[2] = -t + 1.0
        
        return self.event_array    
    
    #Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        event_info = event_info[0] #We only look at the state events information.
        while True: #Event Iteration
            self.event_switch(solver, event_info) #Turns the switches
            
            b_mode = self.state_events(solver.t, solver.y, solver.sw).copy()
            self.init_mode(solver) #Pass in the solver to the problem specified init_mode
            a_mode = self.state_events(solver.t, solver.y, solver.sw).copy()
            
            event_info = self.check_eIter(b_mode, a_mode)
                
            if True not in event_info: #Breaks the iteration loop
                break
    
    #Helper function for handle_event
    def event_switch(self, solver, event_info):
        """
        Turns the switches.
        """
        for i in range(len(event_info)): #Loop across all event functions
            if event_info[i] != 0:
                solver.sw[i] = not solver.sw[i] #Turn the switch
        
    #Helper function for handle_event
    def check_eIter(self, before, after):
        """
        Helper function for handle_event to determine if we have event
        iteration.
        
            Input: Values of the event indicator functions (state_events)
            before and after we have changed mode of operations.
        """
        
        eIter = [False]*len(before)
        
        for i in range(len(before)):
            if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
                eIter[i] = True
                
        return eIter
    
    def init_mode(self, solver):
        """
        Initialize the DAE with the new conditions.
        """
        solver.y[1] = (-1.0 if solver.sw[1] else 3.0)
        solver.y[2] = (0.0 if solver.sw[2] else 2.0)

class Eval_Failure(Explicit_Problem):
    """Problem for testing evaluation failures starting from a given time point and 
    aborting on BaseExceptions."""
    y0 = np.array([1.])
    def __init__(self, t_failure = 0.5, max_evals = 1000):
        self.t_failure = t_failure
        self.max_evals = max_evals
        self.n_evals = 0
    
    def rhs(self, t, y, sw = None):
        self.n_evals += 1
        if t > self.t_failure:
            raise ValueError("passed failure time")
        if (self.max_evals > 0) and (self.n_evals > self.max_evals):
            raise BaseException("Abort")
        return np.array([-1.])


class BaseExceptionAux:
    """Auxiliary class for creating problems (both explicit and implicit) to test
    simulation termination on BaseExceptions.
    
    Set 'fcn', 'jac' or 'event' to True to enable BaseExceptions for the respective functions."""
    def __init__(self, dim, fcn = False, jac = False, event = False, fcn_n = 5, event_n = 5):
        self.dim = dim
        self.fcn_raise = fcn
        self.fcn_n = fcn_n
        self.jac_raise = jac
        self.event_raise = event
        self.event_n = event_n
        self.n_f = 0
        self.n_e = 0

    def f(self, t, y, sw = None):
        if self.fcn_raise:
            self.n_f += 1
            if self.n_f % self.fcn_n == 0:
                raise BaseException('f')
        return -y

    def f_impl(self, t, y, yd, sw = None):
        if self.fcn_raise:
            self.n_f += 1
            if self.n_f % self.fcn_n == 0:
                raise BaseException('f_impl')
        return -y
        
    def jac(self, t, y, sw = None):
        if self.jac_raise:
            raise BaseException('jac')
        else:
            return -np.eye(self.dim)

    def state_events(self,t,y,sw):
        if self.event_raise:
            self.n_e += 1
            if self.n_e % self.event_n == 0:
                raise BaseException('event')
        return np.ones(len(sw))

    def handle_event(self, solver, event_info):
        pass


class ExplicitProbBaseException(Explicit_Problem):
    def __init__(self, dim, fcn = False, jac = False, event = False, fcn_n = 5, event_n = 5):
        self.aux = BaseExceptionAux(dim = dim, fcn = fcn, jac = jac, event = event, fcn_n = fcn_n, event_n = event_n)
        y0 = np.array([1.]*dim)
        sw0 = None
        if event:
            self.state_events = self.aux.state_events
            self.handle_event = self.aux.handle_event
            sw0 = np.array([1.])
        super().__init__(self.aux.f, y0, sw0 = sw0)
        self.jac = self.aux.jac

class ImplicitProbBaseException(Implicit_Problem):
    def __init__(self, dim, fcn = False, jac = False, event = False, fcn_n = 5, event_n = 5):
        self.aux = BaseExceptionAux(dim = dim, fcn = fcn, jac = jac, event = event, fcn_n = fcn_n, event_n = event_n)
        y0 = np.array([1.]*dim)
        yd0 = np.array([0.]*dim)
        sw0 = None
        if event:
            self.state_events = self.aux.state_events
            self.handle_event = self.aux.handle_event
            sw0 = np.array([1.])
        super().__init__(self.aux.f_impl, y0, yd0, sw0 = sw0)
