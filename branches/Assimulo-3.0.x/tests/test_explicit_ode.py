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
from assimulo import testattr
from assimulo.explicit_ode import *
from assimulo.problem import Explicit_Problem
from assimulo.exception import *

class Test_Explicit_ODE:
    pass
    
    @testattr(stddist = True)
    def test_elapsed_step_time(self):
        rhs = lambda t,y: y
        
        prob = Explicit_Problem(rhs, 0.0)
        solv = Explicit_ODE(prob)
        
        assert solv.get_elapsed_step_time() == -1.0
        
    @testattr(stddist = True)
    def test_problem_name_attribute(self):
        rhs = lambda t,y: y
        
        prob = Explicit_Problem(rhs, 0.0)
        assert prob.name == "---"
        prob = Explicit_Problem(rhs, 0.0, name="Test")
        assert prob.name == "Test"
    
    @testattr(stddist = True)
    def test_re_init(self):
        
        rhs = lambda t,y: y
        
        prob = Explicit_Problem(rhs, 0.0)
        solv = Explicit_ODE(prob)
        
        assert solv.t == 0.0
        assert solv.y[0] == 0.0
        
        solv.re_init(1.0, 2.0)
        
        assert solv.t == 1.0
        assert solv.y[0] == 2.0
