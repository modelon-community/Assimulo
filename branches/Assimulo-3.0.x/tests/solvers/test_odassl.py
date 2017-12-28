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
from assimulo.solvers.odassl import *
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.problem import Overdetermined_Problem
from assimulo.exception import *
import numpy as np

class Test_ODASSL:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y,yd: yd + 1
        y0 = [1.0, 1.0, 1.0]
        yd0 = [-1.0, -1.0, -1.0]
        
        self.problem = Overdetermined_Problem(f,y0, yd0)
        self.simulator = ODASSL(self.problem)
    
    @testattr(stddist = True)
    def test_overdetermined(self):
        f = lambda t,y,yd: np.hstack((yd + 1, yd +1))
        y0 = [1.0, 1.0, 1.0]
        yd0 = [-1.0, -1.0, -1.0]
        
        self.problem = Overdetermined_Problem(f,y0, yd0)
        self.simulator = ODASSL(self.problem)
        
        self.simulator.simulate(1)
        
    @testattr(stddist = True)
    def test_implicit_problem(self):
        f = lambda t,y,yd: yd + 1
        y0 = [1.0, 1.0, 1.0]
        yd0 = [-1.0, -1.0, -1.0]
        
        self.problem = Implicit_Problem(f,y0, yd0)
        self.simulator = ODASSL(self.problem)
        
        self.simulator.simulate(1)

    @testattr(stddist = True)
    def test_atol(self):
        
        #Test a simulation
        self.simulator.simulate(1)
        
        self.simulator.reset()
        self.simulator.atol = 1e-6
        
        #Test a simulation
        self.simulator.simulate(1)
        
        self.problem.algvar = [1,1,0]
        simulator = ODASSL(self.problem)
        simulator.atol = 1e-6
        
        #Test a simulation
        self.simulator.simulate(1)
    
    @testattr(stddist = True)
    def test_rtol(self):
        
        #Test a simulation
        self.simulator.simulate(1)
        
        self.simulator.reset()
        self.simulator.rtol = 1e-6
        
        #Test a simulation
        self.simulator.simulate(1)
        
        self.problem.algvar = [1,1,0]
        simulator = ODASSL(self.problem)
        simulator.rtol = 1e-6
        
        #Test a simulation
        self.simulator.simulate(1)    
