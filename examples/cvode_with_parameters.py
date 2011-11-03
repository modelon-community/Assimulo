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
from assimulo.solvers.sundials import CVode
from assimulo.problem import Explicit_Problem

def run_example():
    #This is the same example from the Sundials package (cvsRoberts_FSA_dns.c)
    #----
    #This simple example problem for CVode, due to Robertson, 
    #is from chemical kinetics, and consists of the following three 
    #equations:
    #
    #   dy1/dt = -p1*y1 + p2*y2*y3
    #   dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
    #   dy3/dt = p3*(y2)^2
    
    def f(t, y, p):
        
        yd_0 = -p[0]*y[0]+p[1]*y[1]*y[2]
        yd_1 = p[0]*y[0]-p[1]*y[1]*y[2]-p[2]*y[1]**2
        yd_2 = p[2]*y[1]**2
        
        return N.array([yd_0,yd_1,yd_2])
    
    def handle_result(solver, t ,y):
        solver.t += [t]
        solver.y += [y]
        #if len(solver.t) !=1:
        print solver.interpolate_sensitivity(t, 0, 0)
            #solver.p1 += [solver.interpolate_sensitivity(t, 0, 0)]
            #solver.p2 += [solver.interpolate_sensitivity(t, 0, 1)]
            #solver.p3 += [solver.interpolate_sensitivity(t, 0, 2)]
    
    #The initial conditions
    y0 = [1.0,0.0,0.0]          #Initial conditions for y
    
    #Create an Assimulo explicit problem
    exp_mod = Explicit_Problem(f,y0)
    
    #Sets the options to the problem
    exp_mod.handle_result = handle_result #Change the default handling of the result
    
    exp_mod.p0 = [0.040, 1.0e4, 3.0e7]  #Initial conditions for parameters
    exp_mod.pbar = [0.040, 1.0e4, 3.0e7]

    #Create an Assimulo explicit solver (CVode)
    exp_sim = CVode(exp_mod)
    
    #Sets the paramters
    exp_sim.iter = 'Newton'
    exp_sim.discr = 'BDF'
    exp_sim.rtol = 1.e-4
    exp_sim.atol = N.array([1.0e-8, 1.0e-14, 1.0e-6])
    #exp_sim.pbar = exp_mod.p0 #pbar is used to estimate the tolerances for the parameters
    exp_sim.sensmethod = 'SIMULTANEOUS' #Defines the sensitvity method used
    exp_sim.suppress_sens = False            #Dont suppress the sensitivity variables in the error test.
    exp_sim.options["continuous_output"] = True
    #exp_sim.p1 = [] #Vector for storing the p1 result
    #exp_sim.p2 = [] #Vector for storing the p2 result
    #exp_sim.p3 = [] #Vector for storing the p3 result
    
    #Simulate
    exp_sim.simulate(4,40) #Simulate 4 seconds with 400 communication points
    
    #Plot
    exp_sim.plot() #Plot the solution
    
    

if __name__=='__main__':
    run_example()
