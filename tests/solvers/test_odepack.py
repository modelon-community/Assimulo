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
from assimulo.solvers import LSODAR
from assimulo.problem import Explicit_Problem
from assimulo.exception import *

import numpy as N

class Test_LSODAR:
    """
    Tests the LSODAR solver.
    """
    def setUp(self):
        """
        This sets up the test case.
        """
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return N.array([yd_0,yd_1])
        
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            J = N.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return J
        
        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod_t0 = Explicit_Problem(f,y0,1.0)
        
        exp_mod.jac = jac
        self.mod = exp_mod
            
        #Define an explicit solver
        self.sim = LSODAR(exp_mod) #Create a LSODAR solve
        
        #Sets the parameters
        self.sim.atol = 1e-6 #Default 1e-6
        self.sim.rtol = 1e-6 #Default 1e-6
        self.sim.usejac = False

    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.simulate(1.) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
        
    @testattr(stddist = True)
    def test_simulation_with_jac(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.usejac = True
        self.sim.simulate(1.) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.simulate(1.,100) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
        
    @testattr(stddist = True)    
    def test_simulation_ncp_with_jac(self):
        """
        Test a simulation with ncp.
        """
        self.sim.usejac= True
        self.sim.simulate(1.,100) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
