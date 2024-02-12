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
from assimulo import testattr
from assimulo.solvers.rosenbrock import RodasODE
from assimulo.problem import Explicit_Problem
from assimulo.exception import TimeLimitExceeded
import numpy as np
import scipy.sparse as sps

float_regex = "\s*[+-]?\d*.\d*((e|E)[+-]?\d*)?"

class Test_RodasODE:
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        #Define the rhs
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return np.array([yd_0,yd_1])
        
        #Define the jacobian
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            j = np.array([[0.0,1.0],[my*((-2.0*y[0])*y[1]-1.0), my*(1.0-y[0]**2)]])
            return j
        
        def jac_sparse(t,y):
            eps = 1.e-6
            my = 1./eps
            J = np.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return sps.csc_matrix(J)
        
        y0 = [2.0,-0.6] #Initial conditions
        
        #Define an Assimulo problem
        exp_mod = Explicit_Problem(f,y0, name = 'Van der Pol (explicit)')
        exp_mod_sp = Explicit_Problem(f,y0, name = 'Van der Pol (explicit)')
        exp_mod.jac = jac
        exp_mod_sp.jac = jac_sparse
        cls.mod = exp_mod
        cls.mod_sp = exp_mod_sp
    
    @testattr(stddist = True)
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = RodasODE(self.mod)
        
        sim.usejac = False
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] > 0
        
        sim = RodasODE(self.mod)
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] == 0
    
    @testattr(stddist = True)
    def test_usejac_csc_matrix(self):
        sim = RodasODE(self.mod_sp)
        
        sim.usejac = True
        
        sim.simulate(2.) #Simulate 2 seconds
    
        assert sim.statistics["nfcnjacs"] == 0
        
        assert sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)

    @testattr(stddist = True)
    def test_time_limit(self):
        """ Test that simulation is canceled when a set time limited is exceeded. """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = RodasODE(prob)
        
        sim.h = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = f'The time limit was exceeded at integration time {float_regex}.'
        with pytest.raises(TimeLimitExceeded, match = err_msg):
            sim.simulate(1.)
