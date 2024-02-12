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

import pytest
import numpy as np
from assimulo.problem import Explicit_Problem, Implicit_Problem
from assimulo.solvers import Radau5DAE, Dopri5, RodasODE

def res(t,y,yd,sw):
    return np.array([yd+y])
def state_events(t, y, yd, sw):
    return y - 0.5
def rhs(t,y,sw):
    return np.array([-y])
def estate_events(t, y, sw):
    return y - 0.5
def handle_event(solver, event_info):
    pass

class Test_Solvers:
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        cls.problem = Implicit_Problem(res, [1.0], [-1.0])
        cls.problem.state_events = state_events
        cls.problem.handle_event = handle_event

        cls.eproblem = Explicit_Problem(rhs, [1.0])
        cls.eproblem.state_events = estate_events
        cls.eproblem.handle_event = handle_event
    
    def test_radau5dae_state_events(self):
        solver = Radau5DAE(self.problem)
        
        t,y,yd = solver.simulate(2,33)
        
        assert float(y[-1]) == pytest.approx(0.135, abs = 1e-3)
        
    def test_dopri5_state_events(self):
        solver = Dopri5(self.eproblem)
        
        t,y = solver.simulate(2,33)
        
        assert float(y[-1]) == pytest.approx(0.135, abs = 1e-3)
        
    def test_rodasode_state_events(self):
        solver = RodasODE(self.eproblem)
        
        t,y = solver.simulate(2,33)
        
        assert float(y[-1]) == pytest.approx(0.135, abs = 1e-3)
