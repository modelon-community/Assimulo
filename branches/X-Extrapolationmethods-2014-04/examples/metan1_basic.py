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

#sudo python setup.py install

import numpy as N
import pylab as P
import nose
from assimulo.solvers import Metan1
from assimulo.problem import Explicit_Problem
import time


# cd -;sudo python setup.py install;cd -;ipython
 
def run_example(with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    
    
     
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, y0=N.array([4.0]), name = 'Eulex Test Example: $\dot y = - y$')
    #Define an explicit solver
    exp_sim = Metan1(exp_mod) #Create an eulex solver
    
    #Sets the parameters
    #exp_sim.iter  = 'Newton' #Default 'FixedPoint'
    #exp_sim.discr = 'BDF' #Default 'Adams'
    #exp_sim.atol = [1e-4] #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6

    #Simulate
    start=time.clock()
    t1, y1 = exp_sim.simulate(5,90) #Simulate 5 seconds
    #t2, y2 = exp_sim.simulate(7,50) #Simulate 2 seconds more
    stop=time.clock()
    print'tottttal time is {} '.format(stop-start)
    #Basic test
    #nose.tools.assert_almost_equal(y2[-1], 0.00347746, 5)
    #nose.tools.assert_almost_equal(exp_sim.get_last_step(), 0.0222169642893, 3)
    
    #Plot
    if with_plots:
        P.plot(t1, y1, color="b")
        #P.plot(t2, y2, color="r")
        P.title(exp_mod.name)
        P.ylabel('y')
        P.xlabel('Time')
        P.show()
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()
