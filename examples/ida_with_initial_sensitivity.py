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

import numpy as np
import nose
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    """
    This example show how to use Assimulo and IDA for simulating sensitivities
    for initial conditions.::
    
        0 = dy1/dt - -(k01+k21+k31)*y1 - k12*y2 - k13*y3 - b1
        0 = dy2/dt - k21*y1 + (k02+k12)*y2
        0 = dy3/dt - k31*y1 + k13*y3
     
        y1(0) = p1, y2(0) = p2, y3(0) = p3
        p1=p2=p3 = 0 
    
    See http://sundials.2283335.n4.nabble.com/Forward-sensitivities-for-initial-conditions-td3239724.html
    
    on return:
    
       - :dfn:`imp_mod`    problem instance
    
       - :dfn:`imp_sim`    solver instance
    
    """
    
    def f(t, y, yd,p):
        y1,y2,y3 = y
        yd1,yd2,yd3 = yd
        k01 = 0.0211
        k02 = 0.0162
        k21 = 0.0111
        k12 = 0.0124
        k31 = 0.0039
        k13 = 0.000035
        b1 = 49.3
        
        res_0 = -yd1 -(k01+k21+k31)*y1+k12*y2+k13*y3+b1
        res_1 = -yd2 + k21*y1-(k02+k12)*y2
        res_2 = -yd3 + k31*y1-k13*y3
        
        return np.array([res_0,res_1,res_2])
    
    #The initial conditions
    y0 = [0.0,0.0,0.0]          #Initial conditions for y
    yd0 = [49.3,0.,0.]
    p0 = [0.0, 0.0, 0.0]  #Initial conditions for parameters
    yS0 = np.array([[1,0,0],[0,1,0],[0,0,1.]])
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem(f,y0,yd0,p0=p0,name='Example: Computing Sensitivities')
    
    #Sets the options to the problem
    imp_mod.yS0=yS0

    #Create an Assimulo explicit solver (IDA)
    imp_sim = IDA(imp_mod)
    
    #Sets the parameters
    imp_sim.rtol = 1e-7
    imp_sim.atol = 1e-6
    imp_sim.pbar = [1,1,1] #pbar is used to estimate the tolerances for the parameters
    imp_sim.report_continuously = True #Need to be able to store the result using the interpolate methods
    imp_sim.sensmethod = 'SIMULTANEOUS' #Defines the sensitivity method used
    imp_sim.suppress_sens = False            #Dont suppress the sensitivity variables in the error test.

    #Simulate
    t, y, yd = imp_sim.simulate(400) #Simulate 400 seconds
    
    #Plot
    if with_plots:
        import pylab as pl
        pl.figure(1)
        pl.subplot(221)
        pl.plot(t, np.array(imp_sim.p_sol[0])[:,0],
               t, np.array(imp_sim.p_sol[0])[:,1],
               t, np.array(imp_sim.p_sol[0])[:,2])
        pl.title("Parameter p1")
        pl.legend(("p1/dy1","p1/dy2","p1/dy3"))
        pl.subplot(222)
        pl.plot(t, np.array(imp_sim.p_sol[1])[:,0],
               t, np.array(imp_sim.p_sol[1])[:,1],
               t, np.array(imp_sim.p_sol[1])[:,2])
        pl.title("Parameter p2")
        pl.legend(("p2/dy1","p2/dy2","p2/dy3"))
        pl.subplot(223)
        pl.plot(t, np.array(imp_sim.p_sol[2])[:,0],
               t, np.array(imp_sim.p_sol[2])[:,1],
               t, np.array(imp_sim.p_sol[2])[:,2])
        pl.title("Parameter p3")
        pl.legend(("p3/dy1","p3/dy2","p3/dy3"))
        pl.subplot(224)
        pl.title('ODE Solution')
        pl.plot(t, y)
        pl.suptitle(imp_mod.name)
        pl.show()
    
    #Basic test
    nose.tools.assert_almost_equal(y[-1][0], 1577.6552477,3)
    nose.tools.assert_almost_equal(y[-1][1], 611.9574565, 3)
    nose.tools.assert_almost_equal(y[-1][2], 2215.88563217, 3)
    nose.tools.assert_almost_equal(imp_sim.p_sol[0][1][0], 1.0)
    
        
    return imp_mod, imp_sim
    
if __name__=='__main__':
    mod,sim = run_example()
