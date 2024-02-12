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
from assimulo.solvers.euler import ExplicitEuler, ImplicitEuler
from assimulo.problem import Explicit_Problem
from assimulo.exception import AssimuloException, TimeLimitExceeded
import numpy as np
import scipy.sparse as sps

float_regex = r"\s*[+-]?\d*.\d*((e|E)[+-]?\d*)?"

class Extended_Problem(Explicit_Problem):
    
    #Sets the initial conditions directly into the problem
    y0 = [0.0, -1.0, 0.0]
    sw0 = [False,True,True]
    event_array = np.array([0.0,0.0,0.0])
    rhs_array   = np.array([0.0,0.0,0.0])
    
    #The right-hand-side function (rhs)
    def rhs(self,t,y,sw):
        """
        This is our function we are trying to simulate. During simulation
        the parameter sw should be fixed so that our function is continuous
        over the interval. The parameters sw should only be changed when the
        integrator has stopped.
        """
        self.rhs_array[0] = (1.0 if sw[0] else -1.0)
        self.rhs_array[1] = 0.0
        self.rhs_array[2] = 0.0

        return self.rhs_array

    #Sets a name to our function
    name = 'ODE with discontinuities and a function with consistency problem'
    
    #The event function
    def state_events(self,t,y,sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        self.event_array[0] = y[1] - 1.0 
        self.event_array[1] = -y[2] + 1.0
        self.event_array[2] = -t + 1.0
        
        return self.event_array    
    
    #Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        event_info = event_info[0] #We only look at the state events information.
        while True: #Event Iteration
            self.event_switch(solver, event_info) #Turns the switches
            
            b_mode = self.state_events(solver.t, solver.y, solver.sw).copy()
            self.init_mode(solver) #Pass in the solver to the problem specified init_mode
            a_mode = self.state_events(solver.t, solver.y, solver.sw).copy()
            
            event_info = self.check_eIter(b_mode, a_mode)
                
            if True not in event_info: #Breaks the iteration loop
                break
    
    #Helper function for handle_event
    def event_switch(self, solver, event_info):
        """
        Turns the switches.
        """
        for i in range(len(event_info)): #Loop across all event functions
            if event_info[i] != 0:
                solver.sw[i] = not solver.sw[i] #Turn the switch
        
    #Helper function for handle_event
    def check_eIter(self, before, after):
        """
        Helper function for handle_event to determine if we have event
        iteration.
        
            Input: Values of the event indicator functions (state_events)
            before and after we have changed mode of operations.
        """
        
        eIter = [False]*len(before)
        
        for i in range(len(before)):
            if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
                eIter[i] = True
                
        return eIter
    
    def init_mode(self, solver):
        """
        Initialize the DAE with the new conditions.
        """
        solver.y[1] = (-1.0 if solver.sw[1] else 3.0)
        solver.y[2] = (0.0 if solver.sw[2] else 2.0)


class Test_Explicit_Euler:
    
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This function sets up the test case.
        """
        f = lambda t,y: 1.0
        y0 = 1.0
        
        cls.problem = Explicit_Problem(f, y0)
        cls.simulator = ExplicitEuler(cls.problem)
    
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = ExplicitEuler(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)
    
    def test_h(self):
        
        assert self.simulator.h == pytest.approx(0.01)
        self.simulator.h = 1.0
        assert self.simulator.h == pytest.approx(1.0)
        with pytest.raises(AssimuloException):
            self.simulator._set_h([1])
        
    def test_time_event(self):
        f = lambda t,y: np.array(1.0)
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
            assert solver.t == pytest.approx(tnext)
            assert event_info[0] == []
            assert event_info[1]
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = ExplicitEuler(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_integrator(self):
        """
        This tests the functionality of using the normal mode.
        """
        values = self.simulator.simulate(1)
        
        assert self.simulator.t_sol[-1] == pytest.approx(1.0)
        assert float(self.simulator.y_sol[-1]) == pytest.approx(2.0)
    
    def test_step(self):
        """
        This tests the functionality of using one step mode.
        """
        self.simulator.report_continuously = True
        
        self.simulator.h = 0.1
        self.simulator.simulate(1)
        
        assert self.simulator.t_sol[-1] == pytest.approx(1.0)
        assert float(self.simulator.y_sol[-1]) == pytest.approx(2.0)
        
    def test_exception(self):
        """
        This tests that exceptions are no caught when evaluating the RHS in ExpEuler.
        """
        def f(t,y):
            raise NotImplementedError
        
        prob = Explicit_Problem(f,0.0)
        sim = ExplicitEuler(prob)
        
        with pytest.raises(NotImplementedError):
            sim.simulate(1.0)

    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,sw: np.array([1.0])
        state_events = lambda t,x,sw: np.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Explicit_Problem(f,[0.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = ExplicitEuler(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]

    def test_time_limit(self):
        """ Test that simulation is canceled when a set time limited is exceeded. """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = ExplicitEuler(prob)
        
        sim.h = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = f'The time limit was exceeded at integration time {float_regex}.'
        with pytest.raises(TimeLimitExceeded, match = err_msg):
            sim.simulate(1.)
    
class Test_Implicit_Euler:
    
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This function sets up the test case.
        """
        f = lambda t,y: 1.0
        y0 = 1.0
        
        cls.problem = Explicit_Problem(f, y0)
        cls.simulator = ImplicitEuler(cls.problem)
    
    def test_reset_statistics(self):
        assert self.simulator.statistics["nsteps"] == 0
        
        self.simulator.simulate(5)
        nsteps = self.simulator.statistics["nsteps"]
        self.simulator.simulate(6)
        
        assert self.simulator.statistics["nsteps"] < nsteps
    
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: np.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: sps.csc_matrix(np.array([[0.,1.],[0.,0.]])) #Defines the jacobian
        
        exp_mod = Explicit_Problem(f, [1.0,0.0])
        exp_mod.jac = jac
        
        exp_sim = ImplicitEuler(exp_mod)
        exp_sim.simulate(5.,100)
        
        assert exp_sim.statistics["nfcnjacs"] == 0
        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.995500, abs = 1e-4)
        
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.995500, abs = 1e-4)
        assert exp_sim.statistics["nfcnjacs"] > 0
    
    def test_h(self):
        
        assert self.simulator.h == pytest.approx(0.01)
        self.simulator.h = 1.0
        assert self.simulator.h == pytest.approx(1.0)
        with pytest.raises(AssimuloException):
            self.simulator._set_h([1])
        
    def test_time_event(self):
        f = lambda t,y: np.array(1.0)
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
            assert solver.t == pytest.approx(tnext)
            assert event_info[0] == []
            assert event_info[1]
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = ImplicitEuler(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_integrator(self):
        """
        This tests the functionality of using the normal mode.
        """
        values = self.simulator.simulate(1)
        
        assert self.simulator.t_sol[-1] == pytest.approx(1.0)
        assert float(self.simulator.y_sol[-1]) == pytest.approx(2.0)
    
    def test_step(self):
        """
        This tests the functionality of using one step mode.
        """
        self.simulator.report_continuously = True
        
        self.simulator.h = 0.1
        self.simulator.simulate(1)
        
        assert self.simulator.t_sol[-1] == pytest.approx(1.0)
        assert float(self.simulator.y_sol[-1]) == pytest.approx(2.0)
    
    def test_stiff_problem(self):
        f = lambda t,y: -15.0*y
        y0 = 1.0
        
        problem = Explicit_Problem(f, y0)
        simulator = ImplicitEuler(problem)

        t,y = simulator.simulate(1)
        
        y_correct = lambda t: np.exp(-15*t)
        
        abs_err = np.abs(y[:,0]-y_correct(np.array(t)))
        assert np.max(abs_err) < 0.1
        
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,sw: np.array([1.0])
        state_events = lambda t,x,sw: np.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Explicit_Problem(f,[0.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = ImplicitEuler(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]

    def test_time_limit(self):
        """ Test that simulation is canceled when a set time limited is exceeded. """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = ImplicitEuler(prob)
        
        sim.h = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = f'The time limit was exceeded at integration time {float_regex}.'
        with pytest.raises(TimeLimitExceeded, match = err_msg):
            sim.simulate(1.)
