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
from assimulo.solvers.euler import *
from assimulo.problem import Explicit_Problem
from assimulo.exception import *


class Test_Explicit_Euler:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y: 1.0
        y0 = 1.0
        
        self.problem = Explicit_Problem(f, y0)
        self.simulator = ExplicitEuler(self.problem)
    
    @testattr(stddist = True)
    def test_h(self):
        
        nose.tools.assert_almost_equal(self.simulator.h, 0.01)
        self.simulator.h = 1.0
        nose.tools.assert_almost_equal(self.simulator.h, 1.0)
        nose.tools.assert_raises(AssimuloException, self.simulator._set_h, [1])
        
    @testattr(stddist = True)
    def test_time_event(self):
        f = lambda t,y: N.array(1.0)
        global tnext
        global nevent
        tnext = 0.0
        nevent = 0
        def time_events(t,y,sw):
            global tnext,nevent
            events = [1.0, 2.0, 2.5, 3.0]
            for ev in events:
                if t < ev:
                    tnext = ev
                    break
                else:
                    tnext = None
            nevent += 1
            return tnext
            
        def handle_event(solver, event_info):
            solver.y+= 1.0
            global tnext
            nose.tools.assert_almost_equal(solver.t, tnext)
            assert event_info[0] == []
            assert event_info[1] == True
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = ExplicitEuler(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    @testattr(stddist = True)
    def test_integrator(self):
        """
        This tests the functionality of using the normal mode.
        """
        values = self.simulator.simulate(1)
        
        nose.tools.assert_almost_equal(self.simulator.t_sol[-1], 1.0)
        nose.tools.assert_almost_equal(float(self.simulator.y_sol[-1]), 2.0)
    
    @testattr(stddist = True)
    def test_step(self):
        """
        This tests the functionality of using one step mode.
        """
        self.simulator.report_continuously = True
        
        self.simulator.h = 0.1
        self.simulator.simulate(1)
        
        nose.tools.assert_almost_equal(self.simulator.t_sol[-1], 1.0)
        nose.tools.assert_almost_equal(float(self.simulator.y_sol[-1]), 2.0)
        
    @testattr(stddist = True)
    def test_exception(self):
        """
        This tests that exceptions are no caught when evaluating the RHS in ExpEuler.
        """
        def f(t,y):
            raise NotImplementedError
        
        prob = Explicit_Problem(f,0.0)
        sim = ExplicitEuler(prob)
        
        nose.tools.assert_raises(NotImplementedError, sim.simulate, 1.0)

    @testattr(stddist = True)
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,sw: N.array([1.0])
        state_events = lambda t,x,sw: N.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Explicit_Problem(f,[0.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = ExplicitEuler(mod)
        assert sim.sw[0] == True
        sim.simulate(3)
        assert sim.sw[0] == False
    
class Test_Implicit_Euler:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y: 1.0
        y0 = 1.0
        
        self.problem = Explicit_Problem(f, y0)
        self.simulator = ImplicitEuler(self.problem)
    
    @testattr(stddist = True)
    def test_reset_statistics(self):
        assert self.simulator.statistics["nsteps"] == 0
        
        self.simulator.simulate(5)
        nsteps = self.simulator.statistics["nsteps"]
        self.simulator.simulate(6)
        
        assert self.simulator.statistics["nsteps"] < nsteps
    
    @testattr(stddist = True)
    def test_h(self):
        
        nose.tools.assert_almost_equal(self.simulator.h, 0.01)
        self.simulator.h = 1.0
        nose.tools.assert_almost_equal(self.simulator.h, 1.0)
        nose.tools.assert_raises(AssimuloException, self.simulator._set_h, [1])
        
    @testattr(stddist = True)
    def test_time_event(self):
        f = lambda t,y: N.array(1.0)
        global tnext
        global nevent
        tnext = 0.0
        nevent = 0
        def time_events(t,y,sw):
            global tnext,nevent
            events = [1.0, 2.0, 2.5, 3.0]
            for ev in events:
                if t < ev:
                    tnext = ev
                    break
                else:
                    tnext = None
            nevent += 1
            return tnext
            
        def handle_event(solver, event_info):
            solver.y+= 1.0
            global tnext
            nose.tools.assert_almost_equal(solver.t, tnext)
            assert event_info[0] == []
            assert event_info[1] == True
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = ImplicitEuler(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    @testattr(stddist = True)
    def test_integrator(self):
        """
        This tests the functionality of using the normal mode.
        """
        values = self.simulator.simulate(1)
        
        nose.tools.assert_almost_equal(self.simulator.t_sol[-1], 1.0)
        nose.tools.assert_almost_equal(float(self.simulator.y_sol[-1]), 2.0)
    
    @testattr(stddist = True)
    def test_step(self):
        """
        This tests the functionality of using one step mode.
        """
        self.simulator.report_continuously = True
        
        self.simulator.h = 0.1
        self.simulator.simulate(1)
        
        nose.tools.assert_almost_equal(self.simulator.t_sol[-1], 1.0)
        nose.tools.assert_almost_equal(float(self.simulator.y_sol[-1]), 2.0)
    
    @testattr(stddist = True)
    def test_stiff_problem(self):
        f = lambda t,y: -15.0*y
        y0 = 1.0
        
        problem = Explicit_Problem(f, y0)
        simulator = ImplicitEuler(problem)

        t,y = simulator.simulate(1)
        
        y_correct = lambda t: N.exp(-15*t)
        
        abs_err = N.abs(y[:,0]-y_correct(N.array(t)))
        assert N.max(abs_err) < 0.1
        
    @testattr(stddist = True)
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,sw: N.array([1.0])
        state_events = lambda t,x,sw: N.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Explicit_Problem(f,[0.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = ImplicitEuler(mod)
        assert sim.sw[0] == True
        sim.simulate(3)
        assert sim.sw[0] == False
