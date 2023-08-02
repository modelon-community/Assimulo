#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010-2023 Modelon AB
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
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode

def run_example(with_plots=True):
    """
    Simulations for the Gyro (Heavy Top) example in Celledoni/Safstrom: 
        Journal of Physics A, Vol 39, 5463-5478, 2006
        
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
    
    """

    def curl(v):
        return np.array([[0, v[2], -v[1]], [-v[2], 0, v[0]], [v[1], -v[0], 0]])

    # Defines the rhs
    def f(t, u):
        """
        Simulations for the Gyro (Heavy Top) example in Celledoni/Safstrom: 
        Journal of Physics A, Vol 39, 5463-5478, 2006
        """
        I1 = 1000.
        I2 = 5000.
        I3 = 6000.
        u0 = [0, 0, 1.]
        pi = u[0:3]
        Q = (u[3:12]).reshape((3, 3))
        Qu0 = np.dot(Q, u0)
        omega = np.array([pi[0]/I1, pi[1]/I2, pi[2]/I3])
        pid = np.dot(curl(omega), pi)
        Qd = np.dot(curl(omega),Q)
        return np.hstack([pid, Qd.reshape((9,))])

    # Initial conditions
    y0 = np.hstack([[1000.*10, 5000.*10, 6000*10], np.eye(3).reshape((9,))])
    
    # Create an Assimulo explicit problem
    exp_mod = Explicit_Problem(f, y0, name="Gyroscope Example")
    
    # Create an Assimulo explicit solver (CVode)
    exp_sim = CVode(exp_mod)
    
    # Sets the parameters
    exp_sim.discr = 'BDF'
    exp_sim.iter = 'Newton'
    exp_sim.maxord = 2 # Sets the maxorder
    exp_sim.atol = 1.e-10
    exp_sim.rtol = 1.e-10
    
    # Simulate
    t, y = exp_sim.simulate(0.1)
    
    # Plot
    if with_plots:
        import pylab as P
        P.plot(t, y/10000.)
        P.xlabel('Time')
        P.ylabel('States, scaled by $10^4$')
        P.title(exp_mod.name)
        P.show()
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0], 692.800241862)
    nose.tools.assert_almost_equal(y[-1][8], 7.08468221e-1)
    
    return exp_mod, exp_sim    


if __name__=='__main__':
    mod, sim = run_example()
