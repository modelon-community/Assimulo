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
from assimulo.solvers import GLIMDA
from assimulo.problem import Implicit_Problem, Explicit_Problem
from assimulo.exception import GLIMDA_Exception

import numpy as np

class Test_GLIMDA:
    """
    Tests the GLIMDA solver.
    """
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This sets up the test case.
        """
        #Define the residual
        def f(t,y,yd):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            res_0 = yd[0]-yd_0
            res_1 = yd[1]-yd_1
            
            return np.array([res_0,res_1])
        
        y0 = [2.0,-0.6] #Initial conditions
        yd0 = [-.6,-200000.]
        
        #Define an Assimulo problem
        cls.mod = Implicit_Problem(f,y0,yd0)
        cls.mod_t0 = Implicit_Problem(f,y0,yd0,1.0)
            
        #Define an explicit solver
        cls.sim = GLIMDA(cls.mod) #Create a Radau5 solve
        cls.sim_t0 = GLIMDA(cls.mod_t0)
        
        #Sets the parameters
        cls.sim.atol = 1e-4 #Default 1e-6
        cls.sim.rtol = 1e-4 #Default 1e-6
        cls.sim.inith = 1.e-4 #Initial step-size
    
    def test_simulate_explicit(self):
        """
        Test a simulation of an explicit problem using GLIMDA.
        """
        f = lambda t,y:np.array(-y)
        y0 = [1.0]
        
        problem = Explicit_Problem(f,y0)
        simulator = GLIMDA(problem)
        
        assert simulator.yd0[0] == -simulator.y0[0]
        
        t,y = simulator.simulate(1.0)
        
        assert float(y[-1]) == pytest.approx(float(np.exp(-1.0)),4)
    
    def test_maxord(self):
        """
        Tests the maximum order of GLIMDA.
        """
        assert self.sim.maxord == 3#Default
        assert self.sim.options["maxord"] == 3
        
        self.sim.maxord = 2
        
        assert self.sim.maxord == 2#Default
        assert self.sim.options["maxord"] == 2
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_maxord(4)
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_maxord(0)
    
    def test_minord(self):
        """
        Tests the minimum order of GLIMDA.
        """
        assert self.sim.minord == 1#Default
        assert self.sim.options["minord"] == 1
        
        self.sim.minord = 2
        
        assert self.sim.minord == 2#Default
        assert self.sim.options["minord"] == 2
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_minord(4)
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_minord(0)
        
    def test_maxsteps(self):
        """
        Tests the maximum allowed steps of GLIMDA
        """
        assert self.sim.maxsteps == 100000
        assert self.sim.options["maxsteps"] == 100000
        
        self.sim.maxsteps = 100
        
        assert self.sim.maxsteps == 100
        assert self.sim.options["maxsteps"] == 100
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_maxsteps(-1)
    
    def test_newt(self):
        """
        Tests the maximum allowed number of Newton iterations GLIMDA
        """
        assert self.sim.newt == 5
        assert self.sim.options["newt"] == 5
        
        self.sim.newt = 3
        
        assert self.sim.newt == 3
        assert self.sim.options["newt"] == 3
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_newt(-1)
        
    def test_minh(self):
        """
        Tests the minimum stepsize of GLIMDA.
        """
        assert self.sim.minh == np.finfo(np.double).eps
        assert self.sim.options["minh"] == np.finfo(np.double).eps
        
        self.sim.minh = 1e-5
        
        assert self.sim.minh == 1e-5
        assert self.sim.options["minh"] == 1e-5
        
    def test_order(self):
        """
        Tests the order of GLIMDA.
        """
        assert self.sim.order == 0
        assert self.sim.options["order"] == 0
        
        self.sim.order = 1
        
        assert self.sim.order == 1
        assert self.sim.options["order"] == 1
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_order(-1)
    
    def test_maxh(self):
        """
        Tests the maximum stepsize of GLIMDA.
        """
        assert self.sim.maxh == np.inf
        assert self.sim.options["maxh"] == np.inf
        
        self.sim.maxh = 1e5
        
        assert self.sim.maxh == 1e5
        assert self.sim.options["maxh"] == 1e5
        
    def test_maxretry(self):
        """
        Tests the maximum number of retries of GLIMDA.
        """
        assert self.sim.maxretry == 15
        assert self.sim.options["maxretry"] == 15
        
        self.sim.maxretry = 10
        
        assert self.sim.maxretry == 10
        assert self.sim.options["maxretry"] == 10
        
        with pytest.raises(GLIMDA_Exception):
            self.sim._set_maxretry(-1)
