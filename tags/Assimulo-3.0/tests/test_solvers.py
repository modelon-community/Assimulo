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

import nose
import numpy as np
from assimulo import testattr
from assimulo.exception import *
from assimulo.problem import *
from assimulo.solvers import Radau5DAE, Radau5ODE, Dopri5, RodasODE

def res(t,y,yd,sw):
    return np.array([yd+y])
def state_events(t, y, yd, sw):
    return y - 0.5
def rhs(t,y,sw):
    return np.array([-y])
def estate_events(t, y, sw):
    return y - 0.5

problem = Implicit_Problem(res, [1.0], [-1.0])
problem.state_events = state_events

eproblem = Explicit_Problem(rhs, [1.0])
eproblem.state_events = estate_events

class Test_Solvers:
    
    @testattr(stddist = True)
    def test_radau5dae_state_events(self):
        solver = Radau5DAE(problem)
        
        t,y,yd = solver.simulate(2,33)
        
        nose.tools.assert_almost_equal(float(y[-1]), 0.135, 3)
        
    @testattr(stddist = True)
    def test_radau5ode_state_events(self):
        solver = Radau5ODE(eproblem)
        
        t,y = solver.simulate(2,33)
        
        nose.tools.assert_almost_equal(float(y[-1]), 0.135, 3)

    @testattr(stddist = True)
    def test_dopri5_state_events(self):
        solver = Dopri5(eproblem)
        
        t,y = solver.simulate(2,33)
        
        nose.tools.assert_almost_equal(float(y[-1]), 0.135, 3)
        
    @testattr(stddist = True)
    def test_rodasode_state_events(self):
        solver = RodasODE(eproblem)
        
        t,y = solver.simulate(2,33)
        
        nose.tools.assert_almost_equal(float(y[-1]), 0.135, 3)

    
