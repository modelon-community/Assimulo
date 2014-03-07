#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2014 Modelon AB
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
import pylab as P
import nose
from assimulo.solvers import LSODAR
from assimulo.problem import Explicit_Problem
import sys
import os

"""
The bouncing ball example for LSODAR

"""

#Extend Assimulos problem definition
class Extended_Problem(Explicit_Problem):
    
    #Sets the initial conditons directly into the problem
    y0 = [2.0, 0]   # position and (downward) velocity
    sw0=[True,False]
    
    g = -9.81 # gravitational constant
    
    #The right-hand-side function (rhs)
    
    def rhs(self,t,y,sw):
        """
        This is our function we are trying to simulate. 
        """
        yd_0 = y[1]
        yd_1 = self.g
        return N.array([yd_0,yd_1])

    #Sets a name to our function
    name = 'Bouncing Ball Problem'
    
    #The event function
    def state_events(self,t,y,sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        event_0 = y[0] if sw[0] else 5 # hits the ground 
        event_1 = y[1] if sw[1] else 5 # velocity changes sign at topmost point
        
        return N.array([event_0,event_1])
    
    
    #Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        event_info = event_info[0] #We only look at the state events information.
        if event_info[0] !=0:
             solver.sw[0] = False
             solver.sw[1] = True
             solver.y[1] = - 0.88*solver.y[1]
        else:
             solver.sw[0] = True
             solver.sw[1] = False
             
    def initialize(self, solver):
        solver.h_sol=[]
        solver.nq_sol=[]
        
    def handle_result(self, solver, t, y):
        Explicit_Problem.handle_result(self, solver, t, y)
        # Extra output for algorithm analysis
        if solver.report_continuously:
           h, nq = solver.get_algorithm_data()   
           solver.h_sol.extend([h])
           solver.nq_sol.extend([nq])
 
def run_example(with_plots=True):
    """
    Bouncing ball example to demonstrate LSODAR's 
    discontinuity handling.

    Also a way to use :program:`problem.initialize` and :program:`problem.handle_result`
    in order to provide extra information is demonstrated.

    The governing differential equation is

    .. math::

       \\dot y_1 &= y_2\\\\
       \\dot y_2 &= -9.81

    and the switching functions are

    .. math::

       \\mathrm{event}_0 &= y_1 \\;\\;\\;\\text{ if } \\;\\;\\;\\mathrm{sw}_0 = 1\\\\
       \\mathrm{event}_1 &= y_2 \\;\\;\\;\\text{ if }\\;\\;\\; \\mathrm{sw}_1 = 1

    otherwise the events are deactivated by setting the respective value to something different from 0.


    The event handling sets 

    :math:`y_1 = - 0.88 y_1` and :math:`\\mathrm{sw}_1 = 1` if the first event triggers 
    and :math:`\\mathrm{sw}_1 = 0`   if the second event triggers.

    """
    #Create an instance of the problem
    exp_mod = Extended_Problem() #Create the problem

    exp_sim = LSODAR(exp_mod) #Create the solver
    exp_sim.atol=1.e-8
    exp_sim.report_continuously = True
    
    exp_sim.verbosity = 30

    
    #Simulate
    t, y = exp_sim.simulate(10.0) #Simulate 10 seconds 
    
    #Plot
    if with_plots:
        P.subplot(221)
        P.plot(t,y)
        P.title('LSODAR Bouncing ball (one step mode)')
        P.ylabel('States: $y$ and $\dot y$')
        P.subplot(223)
        P.plot(exp_sim.t_sol,exp_sim.h_sol)
        P.title('LSODAR step size plot')
        P.xlabel('Time')
        P.ylabel('Sepsize')
        P.subplot(224)
        P.plot(exp_sim.t_sol,exp_sim.nq_sol)
        P.title('LSODAR order plot')
        P.xlabel('Time')
        P.ylabel('Order')
        P.suptitle(exp_mod.name)
        P.show()
    return exp_mod, exp_sim    
if __name__=="__main__":
    mod,sim = run_example()
    

    
