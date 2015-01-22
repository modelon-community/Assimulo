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
from assimulo.solvers.kinsol import *
from assimulo.problem import Algebraic_Problem
from assimulo.exception import *

class Test_KINSOL:
    
    @testattr(stddist = True)
    def test_problem_name_attribute(self):
        res = lambda y: y
        model  = Algebraic_Problem(res, 1)
        assert model.name == "---"
        model  = Algebraic_Problem(res, 1, name="Test")
        assert model.name == "Test"
        
    @testattr(stddist = True)
    def test_properties_simple(self):
        res = lambda y: y
        model  = Algebraic_Problem(res, 1)
        solver = KINSOL(model)
        
        solver.max_iter = 150
        assert solver.max_iter == 150
        
        solver.no_initial_setup = True
        assert solver.no_initial_setup == True
        
        solver.max_solves_between_setup_calls = 15
        assert solver.max_solves_between_setup_calls == 15
        
        solver.max_newton_step = 1.0
        assert solver.max_newton_step == 1.0
        
        solver.no_min_epsilon = True
        assert solver.no_min_epsilon == True
        
        solver.max_beta_fails = 15
        assert solver.max_beta_fails == 15
