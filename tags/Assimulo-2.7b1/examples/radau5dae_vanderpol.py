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
from assimulo.solvers import Radau5DAE
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    r"""
    Example for the use of Radau5DAE to solve
    Van der Pol's equation
    
    .. math::
       
        \dot y_1 &= y_2 \\
        \dot y_2 &= \mu ((1.-y_1^2) y_2-y_1)

    with :math:`\mu= 10^6`.

    on return:
    
       - :dfn:`imp_mod`    problem instance
    
       - :dfn:`imp_sim`    solver instance

    """
    
    #Define the residual
    def f(t,y,yd):
        eps = 1.e-6
        my = 1./eps
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        
        res_0 = yd[0]-yd_0
        res_1 = yd[1]-yd_1
        
        return N.array([res_0,res_1])
    
    y0 = [2.0,-0.6] #Initial conditions
    yd0 = [-.6,-200000.]
    
    #Define an Assimulo problem
    imp_mod = Implicit_Problem(f,y0,yd0)
    imp_mod.name = 'Van der Pol (implicit)'
    
    #Define an explicit solver
    imp_sim = Radau5DAE(imp_mod) #Create a Radau5 solver
    
    #Sets the parameters
    imp_sim.atol = 1e-4 #Default 1e-6
    imp_sim.rtol = 1e-4 #Default 1e-6
    imp_sim.inith = 1.e-4 #Initial step-size

    #Simulate
    t, y, yd = imp_sim.simulate(2.) #Simulate 2 seconds
    
    #Plot
    if with_plots:
        P.subplot(211)
        P.plot(t,y[:,0])#, marker='o')
        P.xlabel('Time')
        P.ylabel('State')
        P.subplot(212)
        P.plot(t,yd[:,0]*1.e-5)#, marker='o')
        P.xlabel('Time')
        P.ylabel('State derivatives scaled with $10^{-5}$')
        P.suptitle(imp_mod.name)
        P.show()
    
    #Basic test
    x1 = y[:,0]
    assert N.abs(x1[-1]-1.706168035) < 1e-3 #For test purpose
    return imp_mod, imp_sim

if __name__=='__main__':
    mod,sim = run_example()

