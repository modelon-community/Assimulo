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
import pylab as P
import nose
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    r"""
    This is the same example from the Sundials package (cvsRoberts_FSA_dns.c)
    Its purpose is to demonstrate the use of parameters in the differential equation.

    This simple example problem for IDA, due to Robertson
    see http://www.dm.uniba.it/~testset/problems/rober.php, 
    is from chemical kinetics, and consists of the system:
    
    .. math:: 
    
       \dot y_1 -( -p_1 y_1 + p_2 y_2 y_3)&=0 \\
       \dot y_2 -(p_1 y_1 - p_2 y_2 y_3 - p_3 y_2^2)&=0  \\
       \dot y_3 -( p_3  y_ 2^2)&=0 
       
    
    on return:
    
       - :dfn:`imp_mod`    problem instance
    
       - :dfn:`imp_sim`    solver instance
    
    """

    def f(t, y, yd, p):
        
        res1 = -p[0]*y[0]+p[1]*y[1]*y[2]-yd[0]
        res2 = p[0]*y[0]-p[1]*y[1]*y[2]-p[2]*y[1]**2-yd[1]
        res3 = y[0]+y[1]+y[2]-1
        
        return N.array([res1,res2,res3])

    #The initial conditons
    y0 = [1.0, 0.0, 0.0]        #Initial conditions for y
    yd0 = [0.1, 0.0, 0.0]       #Initial conditions for dy/dt
    p0 = [0.040, 1.0e4, 3.0e7]  #Initial conditions for parameters
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem(f, y0, yd0,p0=p0)

    #Create an Assimulo implicit solver (IDA)
    imp_sim = IDA(imp_mod) #Create a IDA solver
    
    #Sets the paramters
    imp_sim.atol = N.array([1.0e-8, 1.0e-14, 1.0e-6])
    imp_sim.algvar = [1.0,1.0,0.0]
    imp_sim.suppress_alg = False #Suppres the algebraic variables on the error test
    imp_sim.report_continuously = True #Store data continuous during the simulation
    imp_sim.pbar = p0
    imp_sim.suppress_sens = False            #Dont suppress the sensitivity variables in the error test.

    #Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
    imp_sim.make_consistent('IDA_YA_YDP_INIT')
    
    #Simulate
    t, y, yd = imp_sim.simulate(4,400) #Simulate 4 seconds with 400 communication points
    print(imp_sim.p_sol[0][-1] , imp_sim.p_sol[1][-1], imp_sim.p_sol[0][-1])
    
    #Basic test
    nose.tools.assert_almost_equal(y[-1][0], 9.05518032e-01, 4)
    nose.tools.assert_almost_equal(y[-1][1], 2.24046805e-05, 4)
    nose.tools.assert_almost_equal(y[-1][2], 9.44595637e-02, 4)
    nose.tools.assert_almost_equal(imp_sim.p_sol[0][-1][0], -1.8761, 2) #Values taken from the example in Sundials
    nose.tools.assert_almost_equal(imp_sim.p_sol[1][-1][0], 2.9614e-06, 8)
    nose.tools.assert_almost_equal(imp_sim.p_sol[2][-1][0], -4.9334e-10, 12)
    
    #Plot
    if with_plots:
        P.plot(t, y)
        P.title(imp_mod.name)
        P.xlabel('Time')
        P.ylabel('State')
        P.show()  
        
    return imp_mod, imp_sim  

if __name__=='__main__':
    mod,sim = run_example()
    
