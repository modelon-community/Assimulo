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
import nose
from assimulo.solvers.sundials import IDA
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    #This is the same example from the Sundials package (idasRoberts_FSA_dns.c)
    #----
    #This simple example problem for IDA, due to Robertson, 
    #is from chemical kinetics, and consists of the following three 
    #equations:
    #
    #   dy1/dt = -p1*y1 + p2*y2*y3
    #   dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
    #   0   = y1 + y2 + y3 - 1

    def f(t, y, yd, p):
        
        res1 = -p[0]*y[0]+p[1]*y[1]*y[2]-yd[0]
        res2 = p[0]*y[0]-p[1]*y[1]*y[2]-p[2]*y[1]**2-yd[1]
        res3 = y[0]+y[1]+y[2]-1
        
        return N.array([res1,res2,res3])
        
    def handle_result(solver, t ,y,yd):
        solver.t += [t]
        solver.y += [y]
        solver.p[0] += [solver.interpolate_sensitivity(t, 0, 0)]
        solver.p[1] += [solver.interpolate_sensitivity(t, 0, 1)]
        solver.p[2] += [solver.interpolate_sensitivity(t, 0, 2)]
    
    #The initial conditons
    y0 = [1.0, 0.0, 0.0]        #Initial conditions for y
    yd0 = [0.1, 0.0, 0.0]       #Initial conditions for dy/dt
    p0 = [0.040, 1.0e4, 3.0e7]  #Initial conditions for parameters
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem(f, y0, yd0,p0=p0)
    
    #Sets the options to the problem
    imp_mod.handle_result = handle_result #Change the default handling of the result

    #Create an Assimulo implicit solver (IDA)
    imp_sim = IDA(imp_mod) #Create a IDA solver
    
    #Sets the paramters
    imp_sim.atol = N.array([1.0e-8, 1.0e-14, 1.0e-6])
    imp_sim.algvar = [1.0,1.0,0.0]
    imp_sim.suppress_alg = False #Suppres the algebraic variables on the error test
    imp_sim.continuous_output = True #Store data continuous during the simulation
    imp_sim.pbar = p0
    imp_sim.suppress_sens = False            #Dont suppress the sensitivity variables in the error test.
    imp_sim.p = [[],[],[]] #Vector for storing the p result
    
    #Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
    imp_sim.make_consistent('IDA_YA_YDP_INIT')
    
    #Simulate
    imp_sim.simulate(4,400) #Simulate 4 seconds with 400 communication points
    
    #Basic test
    nose.tools.assert_almost_equal(imp_sim.y[-1][0], 9.05518032e-01, 4)
    nose.tools.assert_almost_equal(imp_sim.y[-1][1], 2.24046805e-05, 4)
    nose.tools.assert_almost_equal(imp_sim.y[-1][2], 9.44595637e-02, 4)
    
    #Plot
    if with_plots:
        imp_sim.plot() #Plot the solution
    
    

if __name__=='__main__':
    run_example()
    
