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
from assimulo.solvers.rosenbrock import *
from assimulo.problem import Explicit_Problem
from assimulo.exception import *

class Test_RodasODE:
    def setUp(self):
        #Define the rhs
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return N.array([yd_0,yd_1])
        
        #Define the jacobian
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            j = N.array([[0.0,1.0],[my*((-2.0*y[0])*y[1]-1.0), my*(1.0-y[0]**2)]])
            return j
        
        y0 = [2.0,-0.6] #Initial conditions
        
        #Define an Assimulo problem
        exp_mod = Explicit_Problem(f,y0, name = 'Van der Pol (explicit)')
        exp_mod.jac = jac
        self.mod = exp_mod
    
    @testattr(stddist = True)
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = RodasODE(self.mod)
        
        sim.usejac = False
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] > 0
        
        sim = RodasODE(self.mod)
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] == 0
    
