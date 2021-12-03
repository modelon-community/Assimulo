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
        nose.tools.assert_equal(model.name, "---")
        model  = Algebraic_Problem(res, 1, name="Test")
        nose.tools.assert_equal(model.name, "Test")
        
    @testattr(stddist = True)
    def test_properties_simple(self):
        res = lambda y: y
        model  = Algebraic_Problem(res, 1)
        solver = KINSOL(model)
        
        solver.max_iter = 150
        nose.tools.assert_equal(solver.max_iter, 150)
        
        solver.no_initial_setup = True
        nose.tools.assert_true(solver.no_initial_setup)
        
        solver.max_solves_between_setup_calls = 15
        nose.tools.assert_equal(solver.max_solves_between_setup_calls, 15)
        
        solver.max_newton_step = 1.0
        nose.tools.assert_equal(solver.max_newton_step, 1.0)
        
        solver.no_min_epsilon = True
        nose.tools.assert_true(solver.no_min_epsilon)
        
        solver.max_beta_fails = 15
        nose.tools.assert_equal(solver.max_beta_fails, 15)
