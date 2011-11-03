#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N
from assimulo.solvers.euler import *
from assimulo.solvers.runge_kutta import *
from assimulo.problem import Explicit_Problem

def run_example():
        
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])

    #Define an Assimulo problem
    exp_mod = Explicit_Problem()
    exp_mod.f = f #Sets the rhs into the problem
    exp_mod.y0 = 4.0 #Sets the initial conditions
    exp_mod.problem_name = 'Simple Explicit Example'
    
    #Explicit Euler
    #==============
    exp_sim = ExplicitEuler(exp_mod) #Create a explicit Euler solver
    exp_sim.options["continuous_output"] = True
    #Simulate
    exp_sim.simulate(3) #Simulate 3 seconds
    
    #exp_sim.reset()
    
    exp_sim.simulate(5,100) #Simulate 2 second more

    #Plot
    exp_sim.plot() #Plot the solution
    
    #==============
    
    
    #RungeKutta4
    #===========
    exp_sim = RungeKutta4(exp_mod) #Create a RungeKutta4 solver
    #exp_sim.options["continuous_output"] = True
    #Simulate
    exp_sim.simulate(5, 100) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()
    
    #===========

    
    #RungeKutta34
    #=============
    exp_sim = RungeKutta34(exp_mod) #Create a RungeKutta34 solver
    
    #Simulate
    exp_sim.simulate(5) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()

    #Reset the solver
    exp_sim.reset()
    
    exp_sim.inith = 0.1 #Sets the initial step, default = 0.01
    
    #Simulate
    exp_sim.simulate(5) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()
    
    #=============
    

if __name__=='__main__':
    run_example()
