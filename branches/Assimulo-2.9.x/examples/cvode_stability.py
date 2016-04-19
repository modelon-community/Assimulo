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
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Example for the use of the stability limit detection algorithm
    in CVode.
    
    .. math::
       
        \dot y_1 &= y_2 \\
        \dot y_2 &= \mu ((1.-y_1^2) y_2-y_1) \\
        \dot y_3 &= sin(ty_2)

    with :math:`\mu=\frac{1}{5} 10^3`.

    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance

    """
    class Extended_Problem(Explicit_Problem):
        order = []
        
        def handle_result(self, solver, t, y):
            Explicit_Problem.handle_result(self, solver, t, y)
            self.order.append(solver.get_last_order())
            
    eps = 5.e-3
    my = 1./eps
    
    #Define the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        yd_2 = N.sin(t*y[1])
        
        return N.array([yd_0,yd_1, yd_2])
    
    y0 = [2.0,-0.6, 0.1] #Initial conditions
    
    #Define an Assimulo problem
    exp_mod = Extended_Problem(f,y0, 
                          name = "CVode: Stability problem")
    
    #Define an explicit solver
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.stablimdet = True
    exp_sim.report_continuously = True
    
    #Simulate
    t, y = exp_sim.simulate(2.0) #Simulate 2 seconds
    
    #Plot
    if with_plots:
        P.subplot(211)
        P.plot(t,y[:,2])
        P.ylabel("State: $y_1$")
        P.subplot(212)
        P.plot(t,exp_mod.order)
        P.ylabel("Order")
        P.suptitle(exp_mod.name)
        P.xlabel("Time")
        P.show()

    #Basic test
    x1 = y[:,0]
    assert N.abs(x1[-1]-1.8601438) < 1e-1 #For test purpose
    
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()

