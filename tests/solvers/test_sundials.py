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

import re
import pytest
from assimulo.solvers.sundials import CVode, IDA, CVodeError, IDAError, get_sundials_version
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.exception import AssimuloException, TimeLimitExceeded, TerminateSimulation
import numpy as np
import scipy.sparse as sps
from .utils import (
    Extended_Problem,
    Eval_Failure,
    ExplicitProbBaseException,
    ImplicitProbBaseException
)


class Test_CVode:
    
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This function sets up the test case.
        """
        f = lambda t,y:np.array(y)
        y0 = [1.0]
        
        cls.problem = Explicit_Problem(f,y0)
        cls.simulator = CVode(cls.problem)
        cls.simulator.verbosity = 0
    
    def test_backward_integration(self):
        def f(t, y):
            x, v = y
            return [x, -x - 0.1*v]
            
        mod = Explicit_Problem(f, y0=[1, 0], t0=10)
        sim = CVode(mod)
        sim.backward = True
        t, y = sim.simulate(0, ncp_list=np.arange(1, 10))
        
        assert np.all(t == np.arange(0,11)[::-1])
        
        mod = Explicit_Problem(f, y0=[1, 0], t0=10)
        sim = CVode(mod)
        sim.backward = True
        t, y = sim.simulate(0, ncp_list=np.arange(1, 10)[::-1])
        
        assert np.all(t == np.arange(0,11)[::-1])

    def test_event_localizer(self):
        """ Test that CVode internal event localization works correctly."""
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = CVode(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        exp_sim.external_event_detection = False # default; CVode internal event detection
        
        #Simulate
        t, y = exp_sim.simulate(10.0, 1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)
    
    def test_event_localizer_external(self):
        """ Test that CVode with Assimulo event localization works correctly."""
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = CVode(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        exp_sim.external_event_detection = True # Assimulo event detection
        
        #Simulate
        t, y = exp_sim.simulate(10.0, 1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)
    
    def test_get_error_weights(self):
        with pytest.raises(CVodeError):
            self.simulator.get_error_weights()
        
        self.simulator.simulate(1.0)
        weights = self.simulator.get_error_weights()
        assert weights[0] < 1e6
        
    def test_get_used_initial_step(self):
        self.simulator.simulate(1.0)
        
        step = self.simulator.get_used_initial_step()
        assert step == pytest.approx(0.001, abs = 1e-3)
        
        self.simulator.reset()
        
        self.simulator.inith = 1e-8
        self.simulator.simulate(1.0)
        
        step = self.simulator.get_used_initial_step()
        assert np.abs(step-1e-8) < 1e-2
        
    
    def test_get_local_errors(self):
        with pytest.raises(CVodeError):
            self.simulator.get_local_errors()
    
        self.simulator.simulate(1.0)
        
        err = self.simulator.get_local_errors()
        assert err[0] < 1e-5
    
    def test_get_last_order(self):
        with pytest.raises(CVodeError):
            self.simulator.get_last_order()
        
        self.simulator.simulate(1.0)
        
        qlast = self.simulator.get_last_order()
        assert qlast == 4
        
    def test_max_convergence_failures(self):
        assert self.simulator.maxncf == self.simulator.options["maxncf"]
        self.simulator.maxncf = 15
        assert self.simulator.maxncf == 15
        
        with pytest.raises(AssimuloException):
            self.simulator._set_max_conv_fails(-1)
        
    def test_max_error_tests_failures(self):
        assert self.simulator.maxnef == self.simulator.options["maxnef"]
        self.simulator.maxnef = 15
        assert self.simulator.maxnef == 15
        assert self.simulator.options["maxnef"] == 15
        
        with pytest.raises(AssimuloException):
            self.simulator._set_max_err_fails(-1)
        
    def test_max_nonlinear_iterations(self):
        assert self.simulator.maxcor == self.simulator.options["maxcor"]
        self.simulator.maxcor = 15
        assert self.simulator.maxcor == 15
        assert self.simulator.options["maxcor"] == 15
        
        # with pytest.raises(AssimuloException):
        #     self.simulator._set_max_err_fails(-1)
        
    def test_get_current_order(self):  
        
        with pytest.raises(CVodeError):
            self.simulator.get_current_order()

        self.simulator.simulate(1.0)
        
        qcur = self.simulator.get_current_order()
        assert qcur == 4
        
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        assert self.simulator.y == 1.0
        assert self.simulator.discr == 'BDF'
        assert self.simulator.iter == 'Newton'
        assert self.simulator.maxord == 5
        
        self.simulator.discr = 'Adams'
        assert self.simulator.discr == 'Adams'
        assert self.simulator.maxord == 12
    
    def test_time_event(self):
        f = lambda t,y: [1.0]
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
        exp_sim = CVode(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
        
    def test_clear_event_log(self):
        f = lambda t,y: [1.0]
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
        exp_sim = CVode(exp_mod)
        exp_sim.verbosity = 10
        exp_sim(5.,100)
        
        assert len(exp_sim.event_data) == 4
        
        tnext = 0.0
        nevent = 0
        
        exp_sim.reset()
        assert len(exp_sim.event_data) == 0
        
        exp_sim(5.,100)
        assert len(exp_sim.event_data) == 4
    
    def test_time_limit(self):
        f = lambda t,y: -y
        
        exp_mod = Explicit_Problem(f,1.0)
        exp_sim = CVode(exp_mod)
        
        exp_sim.maxh = 1e-8
        exp_sim.time_limit = 1 #One second
        exp_sim.report_continuously = True
        
        with pytest.raises(TimeLimitExceeded):
            exp_sim.simulate(1)
    
    def test_statistics_stored(self):
        """
        Test that the statistics is stored even if there is a TimeLimit exception
        """
        f = lambda t,y: -y
        
        exp_mod = Explicit_Problem(f,1.0)
        exp_sim = CVode(exp_mod)
        
        exp_sim.maxh = 1e-8
        exp_sim.time_limit = 1 #One second
        exp_sim.report_continuously = True
        
        err_msg = "The time limit was exceeded at integration time"
        with pytest.raises(TimeLimitExceeded, match = re.escape(err_msg)):
            exp_sim.simulate(1.0)
        assert any(exp_sim.statistics[k] > 0 for k in exp_sim.statistics.keys()), "No statistics was found to be stored"
    
    def test_discr_method(self):
        """
        This tests the functionality of the property discr.
        """
        
        with pytest.raises(Exception):
            self.simulator._set_discr_method('Test')
        with pytest.raises(Exception):
            self.simulator._set_discr_method(1)
        with pytest.raises(Exception):
            self.simulator._set_discr_method([1.0, 1])
        with pytest.raises(Exception):
            self.simulator._set_discr_method({'Test':'case'})
        with pytest.raises(Exception):
            self.simulator._set_discr_method(5.1)
        with pytest.raises(Exception):
            self.simulator._set_discr_method(['Test'])
        
        self.simulator.discr = 'BDF'
        assert self.simulator.discr == 'BDF'
        self.simulator.discr = 'Adams'
        assert self.simulator.discr == 'Adams'
    
    def test_change_discr(self):
        """
        This tests that the change from Functional to Newton works
        """
        f = lambda t,y: np.array([1.0])
        y0 = 4.0 #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.iter = "FixedPoint"
        exp_sim.simulate(1)
        assert exp_sim.statistics["njacs"] == 0
        exp_sim.iter = "Newton"
        exp_sim.simulate(2)
        assert exp_sim.statistics["njacs"] > 0
        
    def test_change_norm(self):
        
        assert self.simulator.options["norm"] == "WRMS"
        self.simulator.norm = 'WRMS'
        assert self.simulator.norm == 'WRMS'
        self.simulator.norm = 'EUCLIDEAN'
        assert self.simulator.options["norm"] == "EUCLIDEAN"
        assert self.simulator.norm == 'EUCLIDEAN'
        
        f = lambda t,y: np.array([1.0])
        y0 = 4.0 #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.norm = "WRMS"
        exp_sim.simulate(1)
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.norm = "EUCLIDEAN"
        exp_sim.simulate(1)
    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: np.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: np.array([[0.,1.],[0.,0.]]) #Defines the jacobian
        
        exp_mod = Explicit_Problem(f, [1.0,0.0])
        exp_mod.jac = jac
        
        exp_sim = CVode(exp_mod)
        exp_sim.discr='BDF'
        exp_sim.iter='Newton'
        exp_sim.simulate(5.,100)
        
        assert exp_sim.statistics["nfcnjacs"] == 0
        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.75000143, abs = 1e-4)
        
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.75000143, abs = 1e-4)
        assert exp_sim.statistics["nfcnjacs"] > 0
    
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: np.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: sps.csc_matrix(np.array([[0.,1.],[0.,0.]])) #Defines the jacobian
        
        exp_mod = Explicit_Problem(f, [1.0,0.0])
        exp_mod.jac = jac
        
        exp_sim = CVode(exp_mod)
        exp_sim.discr='BDF'
        exp_sim.iter='Newton'
        exp_sim.simulate(5.,100)
        
        assert exp_sim.statistics["nfcnjacs"] == 0
        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.75000143, abs = 1e-4)
        
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        assert exp_sim.y_sol[-1][0] == pytest.approx(-121.75000143, abs = 1e-4)
        assert exp_sim.statistics["nfcnjacs"] > 0
    
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
        
        sim = CVode(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]
    
    def test_iter_method(self):
        """
        This tests the functionality of the property iter.
        """
        
        with pytest.raises(Exception):
            self.simulator._set_iter_method('Test')
        with pytest.raises(Exception):
            self.simulator._set_iter_method(1)
        with pytest.raises(Exception):
            self.simulator._set_iter_method(0)
        with pytest.raises(Exception):
            self.simulator._set_iter_method(['Test'])
        with pytest.raises(Exception):
            self.simulator._set_iter_method([1.0, 1])
        with pytest.raises(Exception):
            self.simulator._set_iter_method(11.1)
        
        self.simulator.iter = 'Newton'
        assert self.simulator.iter == 'Newton'
        self.simulator.iter = 'FixedPoint'
        assert self.simulator.iter == 'FixedPoint'
    
    def test_initial_step(self):
        """
        This tests the functionality of the property initstep.
        """
        
        with pytest.raises(Exception):
            self.simulator._set_initial_step('Test')
        with pytest.raises(Exception):
            self.simulator._set_initial_step(['Test'])
        
        assert self.simulator.inith == 0.0
        self.simulator.inith = 10.0
        assert self.simulator.inith == 10.0
        self.simulator.inith = 1
        assert self.simulator.inith == 1.0
    
    def test_interpolate(self):
        """
        This tests the functionality of the method interpolate.
        """
        f = lambda t,x: np.array(x**0.25)
        
        prob = Explicit_Problem(f,[1.0])

        sim = CVode(prob)
        sim.simulate(10., 100)
        y100 = sim.y_sol
        t100 = sim.t_sol
        sim.reset()
        sim.simulate(10.)
        assert y100[-2][0] == pytest.approx(sim.interpolate(9.9, 0)[0], abs = 1e-5)
    
    def test_ncp_list(self):
        f = lambda t,y:np.array(-y)
        y0 = [4.0]
        
        prob = Explicit_Problem(f,y0)
        sim = CVode(prob)
        
        t, y = sim.simulate(7, ncp_list=np.arange(0, 7, 0.1)) #Simulate 5 seconds
        
        assert y[-1][0] == pytest.approx(0.00364832, abs = 1e-4)
        
    def test_handle_result(self):
        """
        This function tests the handle result.
        """
        f = lambda t,x: x**0.25
        def handle_result(solver,t,y):
            assert solver.t == t
        
        prob = Explicit_Problem(f, [1.0])
        prob.handle_result = handle_result
        
        sim = CVode(prob)
        sim.report_continuously = True
        sim.simulate(10.)
    
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        self.simulator.discr='Adams'
        
        with pytest.raises(Exception):
            self.simulator._set_max_ord("Test")
        with pytest.raises(Exception):
            self.simulator._set_max_ord([1,1])
        
        self.simulator.maxord = -1
        assert self.simulator.maxord == 1
        self.simulator.maxord = 2
        assert self.simulator.maxord == 2
        self.simulator.maxord = 13
        assert self.simulator.maxord == 12
        
        self.simulator.discr='BDF'
        
        with pytest.raises(Exception):
            self.simulator._set_max_ord("Test")
        with pytest.raises(Exception):
            self.simulator._set_max_ord([1,1])
        
        self.simulator.maxord = -1
        assert self.simulator.maxord == 1
        self.simulator.maxord = 2
        assert self.simulator.maxord == 2
        self.simulator.maxord = 6
        assert self.simulator.maxord == 5
    
    def test_spgmr(self):
        f = lambda t,y: np.array([y[1], -9.82])
        fsw = lambda t,y,sw: np.array([y[1], -9.82])
        fp = lambda t,y,p: np.array([y[1], -9.82])
        fswp = lambda t,y,sw,p: np.array([y[1], -9.82])
        jacv = lambda t,y,fy,v: np.dot(np.array([[0,1.],[0,0]]),v)
        jacvsw = lambda t,y,fy,v,sw: np.dot(np.array([[0,1.],[0,0]]),v)
        jacvp = lambda t,y,fy,v,p: np.dot(np.array([[0,1.],[0,0]]),v)
        jacvswp = lambda t,y,fy,v,sw,p: np.dot(np.array([[0,1.],[0,0]]),v)
        y0 = [1.0,0.0] #Initial conditions
        
        def run_sim(exp_mod):
            exp_sim = CVode(exp_mod) #Create a CVode solver
            exp_sim.linear_solver = 'SPGMR' #Change linear solver

            #Simulate
            t, y = exp_sim.simulate(5, 1000) #Simulate 5 seconds with 1000 communication points
        
            #Basic tests
            assert y[-1][0] == pytest.approx(-121.75000000, abs = 1e-4)
            assert y[-1][1] == pytest.approx(-49.100000000)
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod.jacv = jacv #Sets the jacobian
        run_sim(exp_mod)
        
        #Need someway of suppressing error messages from deep down in the Cython wrapper
        #See http://stackoverflow.com/questions/1218933/can-i-redirect-the-stdout-in-python-into-some-sort-of-string-buffer
        try:
            from cStringIO import StringIO
        except ImportError:
            from io import StringIO
        import sys
        stderr = sys.stderr
        sys.stderr = StringIO()
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod.jacv = jacvsw #Sets the jacobian
        with pytest.raises(CVodeError):
            run_sim(exp_mod)
        
        exp_mod = Explicit_Problem(fswp,y0,sw0=[True],p0=1.0)
        exp_mod.jacv = jacvsw #Sets the jacobian
        with pytest.raises(CVodeError):
            run_sim(exp_mod)
        
        #Restore standard error
        sys.stderr = stderr
        
        exp_mod = Explicit_Problem(fp,y0,p0=1.0)
        exp_mod.jacv = jacvp #Sets the jacobian
        run_sim(exp_mod)
        
        exp_mod = Explicit_Problem(fsw,y0,sw0=[True])
        exp_mod.jacv = jacvsw #Sets the jacobian
        run_sim(exp_mod)
        
        exp_mod = Explicit_Problem(fswp,y0,sw0=[True],p0=1.0)
        exp_mod.jacv = jacvswp #Sets the jacobian
        run_sim(exp_mod)
    
    def test_max_order_discr(self):
        """
        This tests the maximum order when the discretization is changed.
        """
        self.simulator.discr = "Adams"
        self.simulator.maxord = 7
        assert self.simulator.maxord == 7
        
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 12
        self.simulator.discr = 'BDF'
        assert self.simulator.maxord == 5
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 12
        self.simulator.maxord = 4
        self.simulator.discr = 'BDF'
        assert self.simulator.maxord == 5
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 12
   
    def test_pretype(self):
        """
        This tests the precondition option.
        """
        assert self.simulator.precond == 'PREC_NONE'
        self.simulator.precond = 'prec_none'
        assert self.simulator.precond == 'PREC_NONE'
        
        with pytest.raises(Exception):
            self.simulator._set_pre_cond(-1.0)
        with pytest.raises(Exception):
            self.simulator._set_pre_cond('PREC_BOTH1')
    
    def test_maxkrylov(self):
        """
        This test the maximum number of krylov subspaces.
        """
        assert self.simulator.maxkrylov == 5
        self.simulator.maxkrylov = 3
        assert self.simulator.maxkrylov == 3
        self.simulator.maxkrylov = 4.5
        assert self.simulator.maxkrylov == 4
        
        with pytest.raises(Exception):
            self.simulator._set_max_krylov('Test')
        
    def test_stablimit(self):
        assert not self.simulator.stablimit
        self.simulator.stablimit = True
        assert self.simulator.stablimit
        assert self.simulator.options["stablimit"]
    
    def test_linearsolver(self):
        """
        This test the choice of the linear solver.
        """
        assert self.simulator.linear_solver == 'DENSE'
        self.simulator.linear_solver = 'dense'
        assert self.simulator.linear_solver == 'DENSE'
        self.simulator.linear_solver = 'spgmr'
        assert self.simulator.linear_solver == 'SPGMR'
        
        with pytest.raises(Exception):
            self.simulator._set_linear_solver('Test')
    
    def test_terminate_simulation(self):
        """
        This tests the functionality of raising TerminateSimulation exception in handle_result.
        """
        class Extended_Problem(Explicit_Problem):
            def handle_event(self, solver, event_info):
                if solver.t > 1.5:
                    raise TerminateSimulation
            rhs = lambda self,t,y,sw: np.array([1.0])
            y0 = [1.0]
            sw0 = [False,True]
            state_events = lambda self,t,y,sw:  np.array([t-1.0, t-2.0])

        exp_mod = Extended_Problem()
        simulator = CVode(exp_mod)
        simulator(3.)
        
        assert simulator.t == pytest.approx(2.000000, abs = 1e-4)
    
    def test_completed_step(self):
        """
        This tests the functionality of the method completed_step.
        """
        global nsteps
        nsteps = 0
        f = lambda t,x: x**0.25
        def completed_step(solver):
            global nsteps
            nsteps += 1
        mod = Explicit_Problem(f, 1.0)
        mod.step_events = completed_step

        sim = CVode(mod)
        
        sim.simulate(2., 100)
        assert len(sim.t_sol) == 101
        assert nsteps == sim.statistics["nsteps"]
        
        sim = CVode(mod)
        nsteps = 0
        sim.simulate(2.)
        assert len(sim.t_sol) == sim.statistics["nsteps"]+1
        assert nsteps == sim.statistics["nsteps"]

    def test_rtol_vector(self):
        """This tests the functionality of using an rtol vector, if supported."""
        f = lambda t, y: y
        prob = Explicit_Problem(f, np.array([1, 1]))
        sim = CVode(prob)

        sim.rtol = [1e-2, 1e-2] # reduces to scalar
        assert sim.rtol == 1e-2

        if sim.supports['rtol_as_vector']:
            sim.rtol = [1e-2, 1e-3]
            assert sim.rtol[0] == 1e-2
            assert sim.rtol[1] == 1e-3

            sim.simulate(1.)
        else:
            with pytest.raises(AssimuloException):
                sim._set_rtol([1e-2, 1e-3])
            with pytest.raises(AssimuloException):
                sim._set_rtol(np.array([1e-2, 1e-3]))

    def test_rtol_zero(self):
        """ Test CVode with rtol = 0. """
        f = lambda t, y: y
        prob = Explicit_Problem(f, np.array([1]))
        sim = CVode(prob)

        sim.rtol = 0.
        assert sim.rtol == 0.

    def test_rtol_vector_with_zeroes(self):
        """ Test CVode with rtol vector containing zeroes. """
        f = lambda t, y: y
        prob = Explicit_Problem(f, np.array([1, 1]))
        sim = CVode(prob)

        if sim.supports['rtol_as_vector']:
            sim.rtol = [1., 0.]
            assert sim.rtol[0] == 1.
            assert sim.rtol[1] == 0.

            sim.simulate(1.)
        else:
            with pytest.raises(AssimuloException):
                sim._set_rtol([1., 0.])

    def test_rtol_vector_sense(self):
        """ Test CVode with rtol vector and sensitivity analysis. """
        n = 2
        f = lambda t, y, p: p*y
        prob = Explicit_Problem(f, np.ones(n), p0 = np.ones(n))
        prob.yS0 = np.zeros((2, 2))

        sim = CVode(prob)

        # not supported
        with pytest.raises(AssimuloException):
            sim._set_rtol([1., 0.])
        # Ok
        sim.rtol = 1e-6
        sim.rtol = [1e-6]
        sim.rtol = np.array([1e-6])

    def test_no_progress_force_min_h(self):
        """ Test example where CVode fails to make progress and minimal stepsize 
        will be forced."""
        prob = Eval_Failure(t_failure = 0.5, max_evals = -1)
        sim = CVode(prob)
        sim.maxstepshnil = 10 # default
        assert sim.minh == 0
        err_msg = "The right-hand side function had repeated recoverable errors."
        with pytest.raises(CVodeError, match = re.escape(err_msg)):
            sim.simulate(1.)
        assert sim.minh > 0

    def test_no_progress_maxstepshnil_not_active(self):
        """ Test example where CVode fails to make progress, but forcing minh 
        is deactivated."""
        prob = Eval_Failure(t_failure = 0.5, max_evals = 1000)
        sim = CVode(prob)
        sim.maxstepshnil = 0
        assert sim.minh == 0
        err_msg = "failed in an unrecoverable manner"
        with pytest.raises(CVodeError, match = re.escape(err_msg)):
            sim.simulate(1.)
        assert sim.minh == 0

    @pytest.mark.parametrize("ncp, ncp_list",
        [
            (10, None),
            (0, [0., 0.7, 1.])
        ]
        )
    def test_maxstepshnil_not_enforced(self, ncp, ncp_list):
        """ Test example where CVode fails to make progress, but minh will not be forced.
        CVode fails in ordinary ways."""
        prob = Eval_Failure(t_failure = 0.5, max_evals = -1)
        sim = CVode(prob)
        assert sim.minh == 0
        sim.report_continuously = False
        sim.maxsteps = 200
        err_msg = "The solver took max internal steps but could not reach tout"
        with pytest.raises(CVodeError, match = re.escape(err_msg)):
            sim.simulate(1., ncp = ncp, ncp_list = ncp_list)
        assert sim.minh == 0

    @pytest.mark.parametrize("val", [-1, 0, 20])
    def test_maxstepshnil_set_valid(self, val):
        """Test setting valid values for the 'maxstepshnil' option."""
        prob = Explicit_Problem(lambda t, x: -x, np.array([1.]))
        sim = CVode(prob)
        sim.maxstepshnil = val
        assert val == sim.maxstepshnil

    @pytest.mark.parametrize("val", [1.23, "no", [1]])
    def test_maxstepshnil_set_invalid(self, val):
        """Test setting invalid values for the 'maxstepshnil' option."""
        prob = Explicit_Problem(lambda t, x: -x, np.array([1.]))
        sim = CVode(prob)
        msg = "'maxstepshnil' must be an integer."
        with pytest.raises(TypeError, match = re.escape(msg)):
            sim.maxstepshnil = val

    def test_base_exception_interrupt_fcn(self):
        """Test that BaseExceptions in right-hand side terminate the simulation."""
        prob = ExplicitProbBaseException(dim = 2, fcn = True)
        sim = CVode(prob)
        msg = "The user-provided rhs function failed in an unrecoverable manner"
        with pytest.raises(CVodeError, match = re.escape(msg)):
            sim.simulate(1.)

    def test_base_exception_interrupt_jac(self):
        """Test that BaseExceptions in jacobian terminate the simulation."""
        prob = ExplicitProbBaseException(dim = 2, jac = True)
        sim = CVode(prob)
        sim.usejac = True

        msg = "The linear solvers setup function failed in an unrecoverable manner"
        with pytest.raises(CVodeError, match = re.escape(msg)):
            sim.simulate(1.)

    def test_base_exception_interrupt_jac_sparse(self):
        """Test that BaseExceptions in jacobian terminate the simulation."""
        prob = ExplicitProbBaseException(dim = 2, jac = True)
        prob.jac_nnz = 1
        sim = CVode(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        msg = "The linear solvers setup function failed in an unrecoverable manner"
        with pytest.raises(CVodeError, match = re.escape(msg)):
            sim.simulate(1.)

    def test_base_exception_interrupt_event_indicator(self):
        """Test that BaseExceptions in event indicator function resp. solout callback correctly terminate solution."""
        prob = ExplicitProbBaseException(dim = 1, event = True, event_n = 3)
        sim = CVode(prob)

        msg = "The rootfinding function failed in an unrecoverable manner"
        with pytest.raises(CVodeError, match = re.escape(msg)):
            sim.simulate(1.)

class Test_IDA:
    
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This function sets up the test case.
        """
        f = lambda t,y,yd: y
        y0 = [1.0]
        yd0 = [1.0]
        
        cls.problem = Implicit_Problem(f,y0,yd0)
        cls.simulator = IDA(cls.problem)
    
    def test_time_limit(self):
        f = lambda t,y,yd: yd-y
        
        exp_mod = Implicit_Problem(f,1.0,1.0)
        exp_sim = IDA(exp_mod)
        
        exp_sim.maxh = 1e-8
        exp_sim.time_limit = 1 #One second
        exp_sim.report_continuously = True
        
        with pytest.raises(TimeLimitExceeded):
            exp_sim.simulate(1)
    
    def test_simulate_explicit(self):
        """
        Test a simulation of an explicit problem using IDA.
        """
        f = lambda t,y:np.array(-y)
        y0 = [1.0]
        
        problem = Explicit_Problem(f,y0)
        simulator = IDA(problem)
        
        assert simulator.yd0[0] == -simulator.y0[0]
        
        t,y = simulator.simulate(1.0)
        
        assert y[-1][0] == pytest.approx(np.exp(-1.0),4)
    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        assert not self.simulator.suppress_alg
        assert self.simulator.algvar[0] == 1.0
        assert self.simulator.sw is None
        assert self.simulator.maxsteps == 10000
        assert self.simulator.y[0] == 1.0
    
    def test_interpolate(self):
        """
        This tests the functionality of the method interpolate.
        """
        f = lambda t,y,yd: y**0.25-yd
        
        prob = Implicit_Problem(f,[1.0],[1.0])

        sim = IDA(prob)
        sim.simulate(10., 100)
        y100 = sim.y_sol
        t100 = sim.t_sol
        sim.reset()
        sim.simulate(10.)
        assert y100[-2] == pytest.approx(sim.interpolate(9.9, 0), abs = 1e-5)
    
    def test_handle_result(self):
        """
        This function tests the handle result.
        """
        f = lambda t,x,xd: x**0.25-xd
        def handle_result(solver, t ,y, yd):
            assert solver.t == t
        
        prob = Implicit_Problem(f, [1.0],[1.0])
        prob.handle_result = handle_result
        
        sim = IDA(prob)

        sim.report_continuously = True
        
        sim.simulate(10.)
    
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        with pytest.raises(Exception):
            self.simulator._set_max_ord("Test")
        with pytest.raises(Exception):
            self.simulator._set_max_ord([1,1])
        
        
        self.simulator.maxord = -1
        assert self.simulator.maxord == 1
        self.simulator.maxord = 2
        assert self.simulator.maxord == 2
        self.simulator.maxord = 6
        assert self.simulator.maxord == 5
    
    def test_tout1(self):
        """
        This tests the functionality of the property tout1.
        """
        with pytest.raises(Exception):
            self.simulator._set_tout1('Test')
        with pytest.raises(Exception):
            self.simulator._set_tout1([1,1])
        with pytest.raises(Exception):
            self.simulator._set_tout1('Test')
        
        assert self.simulator.tout1 == 0.0001
        self.simulator.tout1 = -0.001
        assert self.simulator.tout1 == -0.001
        self.simulator.tout1 = 1
        assert self.simulator.tout1 == 1.0
        
    def test_lsoff(self):
        """
        This tests the functionality of the property lsoff.
        """
        assert not self.simulator.lsoff
        self.simulator.lsoff = True
        assert self.simulator.lsoff
        self.simulator.lsoff = False
        assert not self.simulator.lsoff
    
    def test_initstep(self):
        """
        This tests the funtionality of the property initstep.
        """
        
        def f(t,y,yd):
            res_0 = yd[0] - y[1]
            res_1 = yd[1] +9.82-0.01*y[1]**2
            return np.array([res_0,res_1])
            
        mod = Implicit_Problem(f,y0=[5.0,0.0], yd0=[0.0,9.82])
        
        
        sim = IDA(mod)
        sim.simulate(2.0)

        assert sim.y_sol[-1][0] == pytest.approx(-13.4746473811, abs = 1e-7)
        
        sim.reset()
        sim.inith = 1e-10
        sim.simulate(2.0)

        assert sim.y_sol[-1][0] == pytest.approx(-13.4746596311, abs = 1e-7)
        
    def test_time_event(self):
        """
        This tests the functionality of the time event function.
        """
        f = lambda t,x,xd,sw: xd-x
        
        def time_events(t, y, yd, sw):
            if sw[0]:
                return 1.0
            if sw[1]:
                return 3.0
            return None
        
        def handle_event(solver, event_info):
            
            if event_info[1]:
                solver.y  = np.array([1.0])
                solver.yd = np.array([1.0])
                
                if not solver.sw[0]:
                    solver.sw[1] = False
                
                if solver.sw[0]:
                    solver.sw[0] = False
        
        mod = Implicit_Problem(f,[1.0],[1.0])
        mod.time_events = time_events
        mod.handle_event = handle_event
        mod.sw0 = [True, True]
        
        sim = IDA(mod)
        
        sim.simulate(5.0)

        assert sim.y_sol[38] == pytest.approx(1.0000000, abs = 1e-5)
        assert sim.y_sol[87] == pytest.approx(1.0000000, abs = 1e-5)
        assert sim.t_sol[-1] == pytest.approx(5.0000000, abs = 1e-5)
    
    def test_clear_event_log(self):
        """
        This tests the functionality of the time event function.
        """
        f = lambda t,x,xd,sw: xd-x
        
        def time_events(t, y, yd, sw):
            if sw[0]:
                return 1.0
            if sw[1]:
                return 3.0
            return None
        
        def handle_event(solver, event_info):
            
            if event_info[1]:
                solver.y  = np.array([1.0])
                solver.yd = np.array([1.0])
                
                if not solver.sw[0]:
                    solver.sw[1] = False
                
                if solver.sw[0]:
                    solver.sw[0] = False
        
        mod = Implicit_Problem(f,[1.0],[1.0], sw0=[True, True])
        mod.time_events = time_events
        mod.handle_event = handle_event
        
        sim = IDA(mod)
        sim.verbosity = 10
        assert len(sim.event_data) == 0
        sim.simulate(5.0)
        assert len(sim.event_data) > 0
        
        sim.reset()
        assert len(sim.event_data) == 0
        sim.simulate(5.0)
        assert len(sim.event_data) > 0
        
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x,xd: np.array([xd[0]-x[1], xd[1]-9.82])       #Defines the rhs
        jac = lambda c,t,x,xd: np.array([[c,-1.],[0.,c]]) #Defines the jacobian

        imp_mod = Implicit_Problem(f,[1.0,0.0],[0.,-9.82])
        imp_mod.jac = jac
        
        imp_sim = IDA(imp_mod)
        
        imp_sim.simulate(3,100)

        assert imp_sim.statistics["nfcnjacs"] == 0
        assert imp_sim.y_sol[-1][0] == pytest.approx(45.1900000, abs = 1e-4)
        
        imp_sim.reset()
        imp_sim.usejac=False
        imp_sim.simulate(3.,100)

        assert imp_sim.y_sol[-1][0] == pytest.approx(45.1900000, abs = 1e-4)
        assert imp_sim.statistics["nfcnjacs"] > 0
    
    def test_terminate_simulation(self):
        """
        This tests the functionality of raising TerminateSimulation exception in handle_result.
        """
        class Extended_Problem(Implicit_Problem):
            def handle_event(self,solver, event_info):
                if solver.t > 1.5:
                    raise TerminateSimulation
            res = lambda self,t,y,yd,sw: np.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: np.array([t-1.0, t-2.0])
            y0 = [1.0]
            yd0 = [1.0]
            sw0 = [False]

        prob = Extended_Problem()
    
        sim = IDA(prob)
        sim.simulate(2.5)
        
        assert sim.t == pytest.approx(2.000000, abs = 1e-4)

    def test_terminate_simulation_external_event(self):
        """
        This tests the functionality of raising TerminateSimulation exception in handle_result. External event detection.
        """
        class Extended_Problem(Implicit_Problem):
            def handle_event(self,solver, event_info):
                if solver.t > 1.5:
                    raise TerminateSimulation
            res = lambda self, t, y, yd, sw: np.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: np.array([t-1.0, t-2.0])
            y0 = [1.0]
            yd0 = [1.0]
            sw0 = [False]

        prob = Extended_Problem()
    
        sim = IDA(prob)
        sim.external_event_detection = True
        sim.simulate(2.5)
        
        assert sim.t == pytest.approx(2.000000, abs = 1e-4)
    
    def test_algvar(self):
        """
        This tests the functionality of the property algvar.
        """
        with pytest.raises(Exception):
            self.simulator._set_algvar([1,'hej',1])
        with pytest.raises(Exception):
            self.simulator._set_algvar({'Test':'case'})
        with pytest.raises(Exception):
            self.simulator._set_algvar([-1,0,1])
        with pytest.raises(Exception):
            self.simulator._set_algvar([1.0,1.0])
        with pytest.raises(Exception):
            self.simulator._set_algvar([3.0, 1.0, 1.0])
    
    def test_time_event_2(self):
        f = lambda t,y,yd: y-yd
        global tnext
        global nevent
        tnext = 0.0
        nevent = 0
        def time_events(t,y,yd,sw):
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
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = IDA(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_suppress_alg(self):
        """
        This tests the functionality of the property suppress_alg.
        """
        self.simulator.suppress_alg = True
        assert self.simulator.suppress_alg
        self.simulator.suppress_alg = False
        assert not self.simulator.suppress_alg
        
    def test_make_consistency(self):
        """
        This tests the functionality of the method make_consistency.
        """
        def f(t,y,yd):
            res_1 = y[0] + y[1]+1.0
            res_2 = y[1]
            return np.array([res_1, res_2])
        y0 = [2.0, 2.0]
        yd0 = [1.0 , 0.0]
        
        
        my_Prob = Implicit_Problem(f, y0, yd0)
        
        simulator = IDA(my_Prob)
        
        [flag, y, yd] = simulator.make_consistent('IDA_Y_INIT')
        
        assert y[1] == pytest.approx(0.00000)
        assert y[0] == pytest.approx(-1.0000)
        assert yd[0] == pytest.approx(1.0000)
        assert yd[1] == pytest.approx(0.0000)

    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,xd,sw: np.array([xd[0]- 1.0])
        state_events = lambda t,x,xd,sw: np.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Implicit_Problem(f, [0.0],[1.0])
        mod.f = f
        mod.sw0 = [True]
        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = IDA(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]
    
    def test_completed_step(self):
        """
        This tests the functionality of the method completed_step.
        """
        global nsteps
        nsteps = 0
        def f(t,y,yd):
            res_1 = y[0] + y[1]+1.0
            res_2 = y[1]
            return np.array([res_1, res_2])
        def completed_step(solver):
            global nsteps
            nsteps += 1
        
        y0 = [-1.0, 0.0]
        yd0 = [1.0 , 0.0]    
        
        mod = Implicit_Problem(f, y0, yd0)
        mod.step_events = completed_step
        
        sim = IDA(mod)
        
        sim.simulate(2., 100)
        assert len(sim.t_sol) == 101
        assert nsteps == sim.statistics["nsteps"]
        
        sim = IDA(mod)
        nsteps = 0
        sim.simulate(2.)
        assert len(sim.t_sol) == sim.statistics["nsteps"] + 1
        assert nsteps == sim.statistics["nsteps"]

    def test_base_exception_interrupt_fcn(self):
        """Test that BaseExceptions in right-hand side terminate the simulation. Radau5 + C + implicit problem."""
        prob = ImplicitProbBaseException(dim = 2, fcn = True)
        sim = IDA(prob)

        err_msg = "The user-provided residual function failed in an unrecoverable manner."
        with pytest.raises(IDAError, match = re.escape(err_msg)):
            sim.simulate(1.)


class Test_Sundials:
    
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This sets up the test case.
        """
        class Prob_IDA(Implicit_Problem):
            res = lambda self,t,y,yd,sw: np.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: np.array([t-1.0, t])
            y0 = [1.0]
            yd0 = [1.0]
            sw0 = [False, True]
            
        res = Prob_IDA()
        
        class Prob_CVode(Explicit_Problem):
            rhs = lambda self,t,y,sw: np.array([1.0])
            state_events = lambda self,t,y,sw: np.array([t-1.0, t])
            y0 = [1.0]
            sw0 = [False, True]

        f = Prob_CVode()
        
        cls.simulators = [IDA(res), CVode(f)]
        
        
        f = lambda t,y,yd,p: np.array([0.0])
        y0 = [1.0]
        yd0 = [1.0]
        p0 = [1.0]
        
        mod = Implicit_Problem(f, y0,yd0,p0=p0)
        cls.sim = IDA(mod)
    
    def test_atol(self):
        """
        This tests the functionality of the property atol.
        """
        assert self.simulators[1].atol == 1.0e-6
        assert self.simulators[0].atol == 1.0e-6
        
        for i in range(len(self.simulators)):
            with pytest.raises(Exception):
                self.simulators([i]._set_atol, -1.0)
            with pytest.raises(Exception):
                self.simulators([i]._set_atol, [1.0, 1.0])
            with pytest.raises(Exception):
                self.simulators([i]._set_atol, "Test")
            
            self.simulators[i].atol = 1.0e-5
            assert self.simulators[i].atol == 1.0e-5
            self.simulators[i].atol = 1.0
            assert self.simulators[i].atol == 1.0
            self.simulators[i].atol = 1
            assert self.simulators[i].atol == 1.0
            self.simulators[i].atol = 1001.0
            assert self.simulators[i].atol == 1001.0
            self.simulators[i].atol = [np.array([1e-5])]
            assert len(self.simulators[i].atol.shape) == 1
            assert self.simulators[i].atol == 1e-5
    
    def test_rtol(self):
        """
        This tests the functionality of the property rtol.
        """
        for sim in self.simulators:
            with pytest.raises(Exception):
                sim._set_rtol(-1.0)
            with pytest.raises(Exception):
                sim._set_rtol([1.0, 2.0]) ## size mismatch
            with pytest.raises(Exception):
                sim._set_rtol("Test")
            
            sim.rtol = 1.0e-5
            assert sim.rtol == 1.0e-5
            sim.rtol = 1.0
            assert sim.rtol == 1.0
            sim.rtol = 1001.0
            assert sim.rtol == 1001.0
            sim.rtol = 1001
            assert sim.rtol == 1001.0

    def test_maxh(self):
        """
        This tests the functionality of the property maxh.
        """
        for i in range(len(self.simulators)):
            with pytest.raises(Exception):
                self.simulators([i]._set_max_h, [1.0, 1.0])
            with pytest.raises(Exception):
                self.simulators([i]._set_max_h, "Test")
            
            self.simulators[i].maxh = 1.0e-5
            assert self.simulators[i].maxh == 1.0e-5
            self.simulators[i].maxh = 1.0
            assert self.simulators[i].maxh == 1.0
            self.simulators[i].maxh = 1001.0
            assert self.simulators[i].maxh == 1001.0

    def test_dqtype(self):
        """
        Tests the property of dqtype.
        """
        
        assert self.sim.dqtype == 'CENTERED' #Test the default value.
        
        self.sim.dqtype = 'FORWARD'
        assert self.sim.dqtype == 'FORWARD'
        self.sim.dqtype = 'CENTERED'
        assert self.sim.dqtype == 'CENTERED'
        
        self.sim.dqtype = 'forward'
        assert self.sim.dqtype == 'FORWARD'
        self.sim.dqtype = 'centered'
        assert self.sim.dqtype == 'CENTERED'
        
        with pytest.raises(Exception):
            self.sim._set_dqtype(1)
        with pytest.raises(Exception):
            self.sim._set_dqtype('IDA_CE')
        with pytest.raises(Exception):
            self.sim._set_dqtype([1])
        with pytest.raises(Exception):
            self.sim._set_dqtype(-1)
    
    def test_dqrhomax(self):
        """
        Tests the property of DQrhomax.
        """
        assert self.sim.dqrhomax == 0.0 #Test the default value.
        
        self.sim.dqrhomax = 1.0
        assert self.sim.dqrhomax == 1.0
        self.sim.dqrhomax = 10
        assert self.sim.dqrhomax == 10
        
        with pytest.raises(Exception):
            self.sim._set_dqrhomax(-1)
        with pytest.raises(Exception):
            self.sim._set_dqrhomax('str')
        with pytest.raises(Exception):
            self.sim._set_dqrhomax([])
        with pytest.raises(Exception):
            self.sim._set_dqrhomax(-10)
    
    def test_usesens(self):
        """
        Tests the property of usesens.
        """
        assert self.sim.usesens#Test the default value.
        self.sim.usesens = False
        assert not self.sim.usesens
        self.sim.usesens = 0
        assert not self.sim.usesens
        self.sim.usesens = 1
        assert self.sim.usesens

    def test_sensmethod(self):
        """
        Tests the property of sensmethod.
        """
        assert self.sim.sensmethod == 'STAGGERED' #Test the default value
        
        self.sim.sensmethod = 'SIMULTANEOUS'
        assert self.sim.sensmethod == 'SIMULTANEOUS'
        self.sim.sensmethod = 'STAGGERED'
        assert self.sim.sensmethod == 'STAGGERED'
        
        self.sim.sensmethod = 'simultaneous'
        assert self.sim.sensmethod == 'SIMULTANEOUS'
        self.sim.sensmethod = 'staggered'
        assert self.sim.sensmethod == 'STAGGERED'
        
        with pytest.raises(Exception):
            self.sim._set_sensitivity_method(1)
        with pytest.raises(Exception):
            self.sim._set_sensitivity_method('IDA_CE')
        with pytest.raises(Exception):
            self.sim._set_sensitivity_method([1])
        with pytest.raises(Exception):
            self.sim._set_sensitivity_method(-1)
    
    def test_suppress_sens(self):
        """
        Tests the property of suppress_sens.
        """
        assert not self.sim.suppress_sens
        self.sim.suppress_sens = False
        assert not self.sim.suppress_sens
        self.sim.suppress_sens = 0
        assert not self.sim.suppress_sens
        self.sim.suppress_sens = 1
        assert self.sim.suppress_sens
    
    def test_maxsensiter(self):
        """
        Tests the property of maxsensiter.
        """
        assert self.sim.maxcorS == 3 #Test the default value
        self.sim.maxcorS = 1
        assert self.sim.maxcorS == 1
        self.sim.maxcorS = 10.5
        assert self.sim.maxcorS == 10
        
        with pytest.raises(Exception):
            self.sim._set_max_cor_S('str')
        with pytest.raises(Exception):
            self.sim._set_max_cor_S([])
    
    def test_pbar(self):
        """
        Tests the property of pbar.
        """
        f = lambda t,y,p:np.array([0.0]*len(y))
        y0 = [1.0]*2
        p0 = [1000.0, -100.0]
        exp_mod = Explicit_Problem(f,y0,p0=p0)
        
        exp_sim = CVode(exp_mod)

        assert exp_sim.pbar[0] == pytest.approx(1000.00000, abs = 1e-4)
        assert exp_sim.pbar[1] == pytest.approx(100.000000, abs = 1e-4)
        
        f = lambda t,y,yd,p: np.array([0.0]*len(y))
        yd0 = [0.0]*2
        imp_mod = Implicit_Problem(f,y0,yd0,p0=p0)
        
        imp_sim = IDA(imp_mod)
        
        assert imp_sim.pbar[0] == pytest.approx(1000.00000, abs = 1e-4)
        assert imp_sim.pbar[1] == pytest.approx(100.000000, abs = 1e-4)

    def test_get_sundials_version(self):
        """Test fetching the sundials version."""
        version = get_sundials_version()
        assert isinstance(version, tuple), "Expected version to be a tuple"
