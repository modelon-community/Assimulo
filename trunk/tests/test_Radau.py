#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import nose
import numpy as N
from assimulo.lib.radau_core import Radau_Common, Radau_Exception
from assimulo.explicit_ode import Radau5
from assimulo.implicit_ode import Radau5 as Radau5Imp
from assimulo.problem import Explicit_Problem, Implicit_Problem

class Test_Radau_Common:
    """
    Tests the common attributes of the Radau solvers.
    """
    def setUp(self):
        """
        This sets up the test case.
        """
        self.sim = Radau_Common()
        
    def test_fac1(self):
        """
        This tests the functionality of the property fac1.
        """
        self.sim.fac1 = 0.01
        assert self.sim.fac1 == 0.01
        self.sim.fac1 = 0.001
        assert self.sim.fac1 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac1, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac1, [-1.0])
    
    
    def test_fac2(self):
        """
        This tests the functionality of the property fac2.
        """
        self.sim.fac2 = 0.01
        assert self.sim.fac2 == 0.01
        self.sim.fac2 = 0.001
        assert self.sim.fac2 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac2, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac2, [-1.0])
    
    def test_fnewt(self):
        """
        This tests the functionality of the property fnewt.
        """
        self.sim.fnewt = 0.01
        assert self.sim.fnewt == 0.01
        self.sim.fnewt = 0.001
        assert self.sim.fnewt == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fnewt, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fnewt, [-1.0])
    
    def test_h(self):
        """
        This tests the functionality of the property h.
        """
        self.sim.h = 0.01
        assert self.sim.h == 0.01
        self.sim.h = 0.001
        assert self.sim.h == 0.001
    
    def test_initial_step(self):
        """
        This tests the functionality of the property initial step.
        """
        self.sim.initstep = 0.01
        assert self.sim.initstep == 0.01
        self.sim.initstep = 0.001
        assert self.sim.initstep == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_initial_step, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_initial_step, [-1.0])
    
    def test_newt(self):
        """
        This tests the functionality of the property newt.
        """
        self.sim.newt = 1
        assert self.sim.newt == 1
        self.sim.newt = 10
        assert self.sim.newt == 10
        self.sim.newt = 9.8
        assert self.sim.newt == 9
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_newt, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_newt, [-1.0])
    
    def test_quot1(self):
        """
        This tests the functionality of the property quot1.
        """
        self.sim.quot1 = 0.01
        assert self.sim.quot1 == 0.01
        self.sim.quot1 = 0.001
        assert self.sim.quot1 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot1, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot1, [-1.0])
        
    def test_quot2(self):
        """
        This tests the functionality of the property quot2.
        """
        self.sim.quot2 = 0.01
        assert self.sim.quot2 == 0.01
        self.sim.quot2 = 0.001
        assert self.sim.quot2 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot2, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot2, [-1.0])
    
    def test_safe(self):
        """
        This tests the functionality of the property safe.
        """
        self.sim.safe = 0.01
        assert self.sim.safe == 0.01
        self.sim.safe = 0.001
        assert self.sim.safe == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_safe, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_safe, [-1.0])
    
    def test_thet(self):
        """
        This tests the functionality of the property thet.
        """
        self.sim.thet = 0.01
        assert self.sim.thet == 0.01
        self.sim.thet = 0.001
        assert self.sim.thet == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_thet, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_thet, [-1.0])
    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim.usejac = True
        assert self.sim.usejac == True
        self.sim.usejac = False
        assert self.sim.usejac == False
        self.sim.usejac = 1
        assert self.sim.usejac == True
        self.sim.usejac = []
        assert self.sim.usejac == False
        
    def test_print_statistics(self):
        """
        This tests the functionality of the method print statistics.
        """
        pass
    
class Test_Explicit_Radau:
    """
    Tests the explicit Radau solver.
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
        exp_mod = Explicit_Problem()
        exp_mod.f = f
        exp_mod.jac = jac
            
        #Define an explicit solver
        y0 = [2.0,-0.6] #Initial conditions
        self.sim = Radau5(exp_mod,y0) #Create a Radau5 solve
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.initstep = 1.e-4 #Initial step-size
        self.sim.usejac = False
    
    def test_col(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim._nsteps == 285
        assert self.sim._nfcn   == 2494
        assert self.sim._njac   == 188
        assert self.sim._njacfcn == 564
        assert len(self.sim.y) == 201
        
        nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y[-1][0], 1.7061680350, 4)
        
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        assert self.sim._nsteps == 285
        assert self.sim._nfcn   == 2494
        assert self.sim._njac   == 188
        assert self.sim._njacfcn == 564
        assert len(self.sim.y)  == 286
        
        nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.8753606, 4)
        nose.tools.assert_almost_equal(self.sim.y[-1][0], 1.7061680350, 4)

    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds
        
        assert self.sim._nsteps == 285
        assert self.sim._nfcn   == 2494
        assert self.sim._njac   == 188
        assert self.sim._njacfcn == 0
        assert len(self.sim.y)  == 286
        
        nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.8753606, 4)
        nose.tools.assert_almost_equal(self.sim.y[-1][0], 1.7061680350, 4)

    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim._nsteps == 284
        assert self.sim._nsteps == self.sim._njac
        
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(N.diff(self.sim.t))-N.finfo('double').eps <= 0.01
        
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        self.sim.newt = 10
        self.sim.simulate(1.0)
        
        assert self.sim._nniterfail == 1
    
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        assert self.sim._nsteps == 142
    

class Test_Implicit_Radau:
    """
    Tests the implicit Radau solver.
    """
    def setUp(self):
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
            
            return N.array([res_0,res_1])
        
        #Define an Assimulo problem
        imp_mod = Implicit_Problem()
        imp_mod.f = f
            
        #Define an explicit solver
        y0 = [2.0,-0.6] #Initial conditions
        yd0 = [-.6,-200000.]
        self.sim = Radau5Imp(imp_mod,y0,yd0) #Create a Radau5 solve
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.initstep = 1.e-4 #Initial step-size
        
    def test_simulation(self):
        """
        Test a simulation of the vanderpol equations.
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim._nsteps == 272
    
    
    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(2.)
        assert self.sim._nsteps == 437
        
        
