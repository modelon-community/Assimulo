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
from assimulo.solvers.radau5 import Radau5DAE, _Radau5DAE
from assimulo.solvers.radau5 import Radau5ODE, _Radau5ODE
from assimulo.solvers.radau5 import Radau5Error
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.lib.radau_core import Radau_Exception
from assimulo.exception import TimeLimitExceeded
import scipy.sparse as sps
import numpy as np

from .utils import (
    Extended_Problem,
    Eval_Failure,
    ExplicitProbBaseException,
    ImplicitProbBaseException,
    ExplicitTimeEventCloseToFinalTime,
    ImplicitTimeEventCloseToFinalTime,
)


import re
float_regex = r"[\s]*[\d]*.[\d]*((e|E)(\+|\-)\d\d|)"


class Test_Explicit_Radau5_Py:
    """
    Tests the explicit Radau solver (Python implementation).
    """
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This sets up the test case.
        """
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return np.array([yd_0,yd_1])
        
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            J = np.zeros([2,2])
            
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
        cls.mod = exp_mod
            
        #Define an explicit solver
        cls.sim = _Radau5ODE(exp_mod) #Create a Radau5 solve
        cls.sim_t0 = _Radau5ODE(exp_mod_t0)
        
        #Sets the parameters
        cls.sim.atol = 1e-4 #Default 1e-6
        cls.sim.rtol = 1e-4 #Default 1e-6
        cls.sim.inith = 1.e-4 #Initial step-size
        cls.sim.usejac = False
    
    @pytest.mark.skip("Does not support state events")
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = _Radau5ODE(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)
    
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
        exp_sim = _Radau5ODE(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = _Radau5ODE(self.mod)
        
        assert sim._leny == 2
    
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.report_continuously = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300
        
        #assert self.sim.y[-2][0] == pytest.approx(1.71505001, abs = 1e-4)
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
        self.sim.report_continuously = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        #assert self.sim.y[-2][0] == pytest.approx(1.71505001, abs = 1e-4)
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
        self.sim_t0.simulate(3.)
        assert self.sim_t0.t_sol[0] == pytest.approx(1.0000000, abs = 1e-4)
        assert self.sim_t0.t_sol[-1] == pytest.approx(3.0000000, abs = 1e-4)
        assert self.sim_t0.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der Pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
    
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["nfcnjacs"] == 0
        
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)

    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["nsteps"] == self.sim.statistics["njacs"]
    
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(np.diff(self.sim.t_sol))-np.finfo('double').eps <= 0.01
        
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        self.sim.newt = 10
        self.sim.simulate(1.0)
        
        assert self.sim.statistics["nnfails"] == 1
    
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        assert self.sim.statistics["nsteps"] < 150
        
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        assert self.sim.statistics["nsteps"] < steps*1.5
        
    def test_atol(self):
        """
        This test the absolute tolerance.
        """
        self.sim.simulate(1.0)
        
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        
        self.sim.rtol = 1e-8
        self.sim.atol = 1e-8
        
        self.sim.simulate(1.0)
        steps2 = self.sim.statistics["nsteps"]
        
        assert steps2 > steps
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        assert steps3 == steps2

        err_msg = "atol must be of length one or same as the dimension of the problem."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.atol = [1e-6,1e-6,1e-6]
        
class Test_Explicit_Radau5:
    """
    Tests the explicit Radau solver.
    """
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This sets up the test case.
        """
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return np.array([yd_0,yd_1])
        
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            J = np.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return J
        
        def jac_sparse(t,y):
            eps = 1.e-6
            my = 1./eps
            J = np.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return sps.csc_matrix(J)
        
        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod_t0 = Explicit_Problem(f,y0,1.0)
        exp_mod_sp = Explicit_Problem(f,y0)
        
        exp_mod.jac = jac
        exp_mod_sp.jac = jac_sparse
        cls.mod = exp_mod
            
        #Define an explicit solver
        cls.sim = Radau5ODE(exp_mod) #Create a Radau5 solve
        cls.sim_t0 = Radau5ODE(exp_mod_t0)
        cls.sim_sp = Radau5ODE(exp_mod_sp)
        
        #Sets the parameters
        cls.sim.atol = 1e-4 #Default 1e-6
        cls.sim.rtol = 1e-4 #Default 1e-6
        cls.sim.inith = 1.e-4 #Initial step-size
        cls.sim.usejac = False

    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = Radau5ODE(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)
        assert exp_sim.get_statistics()['nstateevents'] == 1, "Incorrect number of state events"
    
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = Radau5ODE(self.mod)
        
        sim.usejac = False
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] > 0
        
        sim = Radau5ODE(self.mod)
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] == 0
    
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
        exp_sim = Radau5ODE(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = Radau5ODE(self.mod)
        
        assert sim._leny == 2
    
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.report_continuously = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300
        
        #assert self.sim.y[-2][0] == pytest.approx(1.71505001, abs = 1e-4)
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
        self.sim.report_continuously = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        #assert self.sim.y[-2][0] == pytest.approx(1.71505001, abs = 1e-4)
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
        self.sim_t0.simulate(3.)
        assert self.sim_t0.t_sol[0] == pytest.approx(1.0000000, abs = 1e-4)
        assert self.sim_t0.t_sol[-1] == pytest.approx(3.0000000, abs = 1e-4)
        assert self.sim_t0.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
    
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["nfcnjacs"] == 0
        
        assert self.sim.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
    
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim_sp.usejac = True
        
        self.sim_sp.simulate(2.) #Simulate 2 seconds
    
        assert self.sim_sp.statistics["nfcnjacs"] == 0
        
        assert self.sim_sp.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
    
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["nsteps"] == self.sim.statistics["njacs"]
    
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(np.diff(self.sim.t_sol))-np.finfo('double').eps <= 0.01
    
    @pytest.mark.skip("Statistic not recorded by Radau5")
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        self.sim.simulate(1.0)
        self.sim.reset()
        self.sim.newt = 10
        self.sim.simulate(1.0)
        assert self.sim.statistics["nniterfail"] == 1
    
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        assert self.sim.statistics["nsteps"] < 150
        
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        assert self.sim.statistics["nsteps"] < steps*1.5
        
    def test_weighted_error(self):
        
        def handle_result(solver, t, y):
            err = solver.get_weighted_local_errors()
            assert len(err) == len(y)
        
        self.mod.handle_result = handle_result
            
        #Define an explicit solver
        sim = Radau5ODE(self.mod) #Create a Radau5 solve
        sim.get_weighted_local_errors()
        sim.simulate(1)
        
    def test_atol(self):
        """
        This test the absolute tolerance.
        """
        self.sim.simulate(1.0)
        
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        
        self.sim.rtol = 1e-8
        self.sim.atol = 1e-8
        
        self.sim.simulate(1.0)
        steps2 = self.sim.statistics["nsteps"]
        
        assert steps2 > steps
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        assert steps3 == steps2
        
        err_msg = "atol must be of length one or same as the dimension of the problem."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.atol = [1e-6,1e-6,1e-6]
        
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
        
        sim = Radau5ODE(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]

    def test_nmax_steps(self):
        """
        This tests the error upon exceeding a set maximum number of steps
        """
        sim = Radau5ODE(self.mod)

        sim.maxh = 1.e-1
        sim.maxsteps = 9

        err_msg = f'Radau5 failed with flag -5. At time {float_regex}. Message: Maximal number of steps = 9 exceeded.'
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1.)

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_step_size_too_small(self):
        """
        This tests the error for too small step-sizes
        """
        sim = Radau5ODE(self.mod)

        sim.atol = 1.e10
        sim.rtol = 1.e10

        sim.inith = 1.e-1
        sim.maxh = 1.e-1

        err_msg = f"Radau5 failed with flag -6. At time {float_regex}. Message: Stepsize too small with h = {float_regex}."
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1. + 1.e-16)

    def test_repeated_unexpected_step_rejections(self):
        """
        This tests the error for repeated unexpected step rejections
        """
        def f(t, y):
            raise np.linalg.LinAlgError()
        y0 = np.array([1.])
        prob = Explicit_Problem(f, y0)
        sim = Radau5ODE(prob)

        err_msg = 'Repeated unexpected step rejections.'
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1.)

    def test_sparse_solver_jac_disabled(self):
        """
        This tests the trying to simulate using the sparse linear solver, with no analytical jacobian provided.
        """
        f = lambda t, y: [y]
        y0 = np.array([1.])
        prob = Explicit_Problem(f, y0)

        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = False

        sim.simulate(1.)
        assert sim.linear_solver == 'DENSE'

    def test_solver_no_jac(self):
        """
        This tests the error when trying to simulate using an analytical jacobian, with none provided
        """
        f = lambda t, y: [y]
        y0 = np.array([1.])
        prob = Explicit_Problem(f, y0)

        sim = Radau5ODE(prob)
        sim.usejac = True

        err_msg = "Use of an analytical Jacobian is enabled, but problem does contain a 'jac' function."
        with pytest.raises(Radau_Exception, match = err_msg):
            sim.simulate(1.)

    def test_solver_sparse_jac_wrong_format(self):
        """
        This tests the error when using a sparse jacobian of the wrong format
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sps.spdiags([1], 0, 1, 1, format = 'csr')
        y0 = np.array([1.])
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 1

        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = f'Radau5 failed with flag -11. At time {float_regex}. Message: Jacobian given in wrong format, required sparsity format: CSC.'
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1.)

    def test_solver_sparse_jac_nnz_too_small(self):
        """
        This tests the error when using a sparse jacobian with nnz larger than specified
        """
        n = 5
        f = lambda t, y: y
        jac = lambda t, y: sps.eye(n, n, dtype = np.double, format = 'csc')
        y0 = np.array([1.]*n)
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 1

        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = f'Radau5 failed with flag -9. At time {float_regex}. Message: Number of nonzero elements in provided jacobian is too small, specified = 1, actual = 5.'
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1.)

    def test_solver_sparse_jac_nnz_zero(self):
        """
        This tests that using a sparse jacobian with nnz = 0 is valid.
        """
        n = 5
        f = lambda t, y: [0.]*n
        jac = lambda t, y: sps.csc_matrix((n, n), dtype = np.double)
        y0 = np.array([1.]*n)
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 0

        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        sim.simulate(1.)

    def test_sparse_solver_no_nnz(self):
        """
        This tests the error when trying to simulate using the sparse linear solver, without specifying the number of non-zero elements
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sps.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = np.array([1.])
        prob = Explicit_Problem(f, y0)
        prob.jac = jac

        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = "Number of non-zero elements of sparse Jacobian must be non-negative. Detected default value of '-1', has 'problem.jac_fcn_nnz' been set?"
        with pytest.raises(Radau_Exception, match = err_msg):
            sim.simulate(1.)

    def test_sparse_solver_invalid_nnz_type(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; wrong type.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sps.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = np.array([1.])

        for nnz in [None, "test"]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian must be an integer, received: {}."
            with pytest.raises(Radau_Exception, match = err_msg.format(nnz)):
                sim.simulate(1.)

    def test_sparse_solver_invalid_nnz_negative(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; negative.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sps.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = np.array([1.])

        for nnz in [-2, -10]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian must be non-negative, given value = {}."
            with pytest.raises(Radau_Exception, match = err_msg.format(nnz)):
                sim.simulate(1.)

    def test_sparse_solver_invalid_nnz_too_large(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; too_large.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sps.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = np.array([1.])

        for nnz in [5, 100]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian infeasible, must be smaller than the problem dimension squared."
            with pytest.raises(Radau_Exception, match = err_msg):
                sim.simulate(1.)

    def test_sparse_solver_jacobian(self):
        """Testing sparse solver to not produce segmentation faults for Jacobian."""
        ## Take trivial problem with somewhat arbitrary jacobians
        ## Test that functions for internal processing of jacobian do not produces segfaults
        jacobians = [
            (lambda t, y: sps.csc_matrix(np.array([[1., 1., 1.], [1., 1., 1.], [1., 1., 1.]])), 9), 
            (lambda t, y: sps.csc_matrix(np.array([[0., 1., 1.], [1., 0., 1.], [1., 1., 0.]])), 6),
            (lambda t, y: sps.csc_matrix(np.array([[0., 1., 1.], [1., 1., 1.], [1., 1., 1.]])), 8),
            (lambda t, y: sps.csc_matrix(np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sps.csc_matrix(np.array([[0., 0., 0.], [1., 0., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sps.csc_matrix(np.array([[0., 0., 0.], [0., 0., 0.], [0., 1., 0.]])), 1),
            (lambda t, y: sps.csc_matrix(np.array([[0., 0., 1.], [0., 0., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sps.csc_matrix(np.array([[1., 0., 0.], [0., 0., 0.], [0., 0., 0.]])), 1),
        ]

        for i, (jac, nnz) in enumerate(jacobians):
            f = lambda t, y: y
            y0 = 1.*np.ones(3)
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            assert sim.simulate(1.), f"Jacobian #{i} failed: {jac(0, 0)}"

    def test_linear_solver(self):
        """
        This tests the functionality of the property linear_solver.
        """
        self.sim.linear_solver = 'dense'
        assert self.sim.linear_solver == 'DENSE'
        self.sim.linear_solver = 'sparse'
        assert self.sim.linear_solver == 'SPARSE'
        self.sim.linear_solver = 'DENSE'
        assert self.sim.linear_solver == 'DENSE'
        self.sim.linear_solver = 'SPARSE'
        assert self.sim.linear_solver == 'SPARSE'

        err_msg = "'linear_solver' parameter needs to be either 'DENSE' or 'SPARSE'. Set value: {}"
        with pytest.raises(Radau_Exception, match = err_msg.format('default')):
            self.sim.linear_solver = 'default'
        with pytest.raises(Radau_Exception, match = err_msg.format('GMRES')):
            self.sim.linear_solver = 'GMRES'

        err_msg = "'linear_solver' parameter needs to be the STRING 'DENSE' or 'SPARSE'. Set value: {}, type: {}"
        with pytest.raises(Radau_Exception, match = err_msg.format('0', "<class 'int'>")):
            self.sim.linear_solver = 0

    def test_base_exception_interrupt_fcn(self):
        """Test that BaseExceptions in right-hand side terminate the simulation. Radau5 + C + explicit problem."""
        prob = ExplicitProbBaseException(dim = 2, fcn = True)
        sim = Radau5ODE(prob)
        with pytest.raises(BaseException, match = "f"):
            sim.simulate(1.)

    def test_base_exception_interrupt_jac(self):
        """Test that BaseExceptions in jacobian terminate the simulation. Radau5 + C + explicit problem."""
        prob = ExplicitProbBaseException(dim = 2, jac = True)
        sim = Radau5ODE(prob)
        sim.usejac = True

        with pytest.raises(BaseException, match = "jac"):
            sim.simulate(1.)

    def test_base_exception_interrupt_jac_sparse(self):
        """Test that BaseExceptions in jacobian terminate the simulation. Radau5 + C + explicit problem + sparse jac."""
        prob = ExplicitProbBaseException(dim = 2, jac = True)
        prob.jac_nnz = 2
        sim = Radau5ODE(prob)
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        with pytest.raises(BaseException, match = "jac"):
            sim.simulate(1.)

    def test_base_exception_interrupt_event_indicator(self):
        """Test that BaseExceptions in event indicator function resp. solout callback correctly terminate solution."""
        prob = ExplicitProbBaseException(dim = 1, event = True, event_n = 3)
        sim = Radau5ODE(prob)

        with pytest.raises(BaseException, match = "event"):
            sim.simulate(1.)

    def test_time_limit(self):
        """ Test that simulation is canceled when a set time limited is exceeded. """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = Radau5ODE(prob)
        
        sim.maxh = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = f'The time limit was exceeded at integration time {float_regex}.'
        with pytest.raises(TimeLimitExceeded, match = err_msg):
            sim.simulate(1.)
    
    def test_statistics_stored(self):
        """
        Test that the statistics is stored even if there is a TimeLimit exception
        """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = Radau5ODE(prob)
        
        sim.maxh = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = "The time limit was exceeded at integration time"
        with pytest.raises(TimeLimitExceeded, match = re.escape(err_msg)):
            sim.simulate(1.0)
        assert any(sim.statistics[k] > 0 for k in sim.statistics.keys()), "No statistics was found to be stored"

    def test_no_progress(self):
        """Test example where solver cannot make progress past a given time."""
        prob = Eval_Failure(t_failure = 0.5, max_evals = -1)
        sim = Radau5ODE(prob)

        err_msg = "passed failure time"
        with pytest.raises(ValueError, match = re.escape(err_msg)):
            sim.simulate(1.0)

    @pytest.mark.parametrize("tfinal", [1e-6, 1, 1e6])
    @pytest.mark.parametrize("report_continuously", [True, False])
    def test_event_close_to_final_time(self, tfinal, report_continuously):
        """Test event close to final time."""
        exp_sim = Radau5ODE(ExplicitTimeEventCloseToFinalTime(tfinal))
        exp_sim.report_continuously = report_continuously
        tt, _ = exp_sim.simulate(tfinal = tfinal, ncp = 2)
        assert tt[-1] < tfinal # check final interval is skipped

    def test_backwards_report_continuously(self):
        """Test that backward simulation functions with report_continuously = True."""
        mod = Explicit_Problem(lambda t, y: [0], y0 = [1], t0 = 1)
        sim = Radau5ODE(mod)
        sim.backward = True
        sim.report_continuously = True
        sim.simulate(0, ncp = 10)


class Test_Implicit_Radau5:
    """
    Tests the implicit Radau solver.
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
            
        #Define an implicit solver
        cls.sim = Radau5DAE(cls.mod) #Create a Radau5 solve
        cls.sim_t0 = Radau5DAE(cls.mod_t0)
        
        #Sets the parameters
        cls.sim.atol = 1e-4 #Default 1e-6
        cls.sim.rtol = 1e-4 #Default 1e-6
        cls.sim.inith = 1.e-4 #Initial step-size

    def test_implementation_get(self):
        """
            Test getting of implementation property of Radau5DAE.
        """
        assert self.sim.implementation == 'f'

    def test_implementation_set(self):
        """
            Test setting of implementation property of Radau5DAE.
        """
        err_msg = "Radau5DAE does not support setting the 'implementation' attribute, since it only supports the Fortran implementation of Radau5."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.implementation = 'c'

    def test_linear_solver_get(self):
        """
            Test getting of linear_solver property of Radau5DAE.
        """
        assert self.sim.linear_solver == 'DENSE'

    def test_linear_solver_set(self):
        """
            Test setting of linear_solver property of Radau5DAE.
        """
        err_msg = "Radau5DAE does not support setting the 'linear_solver' attribute, since it only supports the DENSE linear solver in Fortran implementation of Radau5."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.linear_solver = 'SPARSE'
    
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = Radau5DAE(self.mod)
        
        sim.usejac = False
        sim.simulate(1)
        
        assert sim.statistics["nfcnjacs"] > 0
    
    def test_simulate_explicit(self):
        """
        Test a simulation of an explicit problem using Radau5DAE.
        """
        f = lambda t,y:np.array(-y)
        y0 = [1.0]
        
        problem = Explicit_Problem(f,y0)
        simulator = Radau5DAE(problem)
        
        assert simulator.yd0[0] == -simulator.y0[0]
        
        t,y = simulator.simulate(1.0)
        
        assert y[-1][0] == pytest.approx(np.exp(-1.0),4)
    
    def test_time_event(self):
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
            #solver.y+= 1.0
            global tnext
            assert solver.t == pytest.approx(tnext)
            assert event_info[0] == []
            assert event_info[1]
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event

        #CVode
        exp_sim = Radau5DAE(exp_mod)
        exp_sim.verbosity = 0
        exp_sim(5.,100)

        assert nevent == 5

    def test_init(self):
        """
        This tests the functionality of Radau5 Implicit Init.
        """
        #Test both y0 in problem and not.

        sim = Radau5DAE(self.mod)
        
        assert sim._leny == 2
    
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(.5) #Simulate 2 seconds
        assert self.sim.statistics["nsteps"] == self.sim.statistics["njacs"]

    def test_simulation(self):
        """
        Test a simulation of the van der Pol equations (1).
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        assert self.sim.y_sol[-1][0] == pytest.approx(1.706272, abs = 1e-3)

        self.sim.reset()

        self.sim.report_continuously = True

        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        assert self.sim.y_sol[-1][0] == pytest.approx(1.706166, abs = 1e-3)

        self.sim_t0.simulate(3.)
        assert self.sim_t0.t_sol[0] == pytest.approx(1.0000000, abs = 1e-4)
        assert self.sim_t0.t_sol[-1] == pytest.approx(3.0000000, abs = 1e-4)
        assert self.sim_t0.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)

    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True

        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201

        self.sim.reset()
        self.sim.report_continuously = False

        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201

    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(np.diff(self.sim.t_sol))-np.finfo('double').eps <= 0.01
        
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        res = lambda t,x,xd,sw: np.array([1.0 - xd])
        state_events = lambda t,x,xd,sw: np.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Implicit_Problem(res,[0.0], [1.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = Radau5DAE(mod)
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]

    def test_nmax_steps(self):
        """
        This tests the error upon exceeding a set maximum number of steps
        """
        sim = Radau5DAE(self.mod)

        sim.maxh = 1.e-1
        sim.maxsteps = 9

        err_msg = "The solver took max internal steps but could not reach the next output time."
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1.)

    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    def test_step_size_too_small(self):
        """
        This tests the error for too small step-sizes
        """
        f = lambda t, y, yd: -y
        y0 = np.array([1.])
        yd0 = np.array([0.])

        prob = Implicit_Problem(f, y0, yd0)

        sim = Radau5DAE(prob)

        sim.atol = 1.e10
        sim.rtol = 1.e10

        sim.inith = 1.e-1
        sim.maxh = 1.e-1
        
        err_msg = f"The step size became too small. At time {float_regex}."
        with pytest.raises(Radau5Error, match = err_msg):
            sim.simulate(1. + 1.e-16)

    def test_repeated_unexpected_step_rejections(self):
        """
        This tests the error for repeated unexpected step rejections in Radau5DAE.
        """
        def f(t, y, yd):
            raise np.linalg.LinAlgError()
        prob = Implicit_Problem(f, np.array([1.]), np.array([1.]))
        sim = Radau5DAE(prob)

        # XXX: Error is raised, but may be due to singular jacobians
        # err_msg = 'Repeated unexpected step rejections.'
        # with pytest.raises(Radau5Error, match = err_msg):
        with pytest.raises(Radau5Error):
            sim.simulate(1.)

    @pytest.mark.parametrize("tfinal", [1e-6, 1, 1e6])
    @pytest.mark.parametrize("report_continuously", [True, False])
    def test_event_close_to_final_time(self, tfinal, report_continuously):
        """Test event close to final time."""
        exp_sim = Radau5DAE(ImplicitTimeEventCloseToFinalTime(tfinal))
        exp_sim.report_continuously = report_continuously
        tt, _, _ = exp_sim.simulate(tfinal = tfinal, ncp = 2)
        assert tt[-1] < tfinal # check final interval is skipped

    def test_backwards_report_continuously(self):
        """Test that backward simulation functions with report_continuously = True."""
        mod = Implicit_Problem(lambda t, y, yd: yd, y0 = [1], yd0 = [0], t0 = 1)
        sim = Radau5DAE(mod)
        sim.backward = True
        sim.report_continuously = True
        sim.simulate(0, ncp = 10)


class Test_Implicit_Radau5_Py:
    """
    Tests the implicit Radau solver (Python implementation).
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
        cls.sim = _Radau5DAE(cls.mod) #Create a Radau5 solve
        cls.sim_t0 = _Radau5DAE(cls.mod_t0)
        
        #Sets the parameters
        cls.sim.atol = 1e-4 #Default 1e-6
        cls.sim.rtol = 1e-4 #Default 1e-6
        cls.sim.inith = 1.e-4 #Initial step-size
    
    def test_time_event(self):
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
            #solver.y+= 1.0
            global tnext
            assert solver.t == pytest.approx(tnext)
            assert event_info[0] == []
            assert event_info[1]
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = _Radau5DAE(exp_mod)
        exp_sim.verbosity = 0
        exp_sim(5.,100)
        
        assert nevent == 5
    
    def test_init(self):
        """
        This tests the functionality of Radau5 Implicit Init.
        """
        #Test both y0 in problem and not.

        sim = _Radau5DAE(self.mod)
        
        assert sim._leny == 2
    
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(.5) #Simulate 2 seconds

        assert self.sim.statistics["nsteps"] == self.sim.statistics["njacs"]
        
    def test_simulation(self):
        """
        Test a simulation of the van der Pol equations (2).
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        assert self.sim.y_sol[-1][0] == pytest.approx(1.706272, abs = 1e-3)
        
        self.sim.reset()
        
        self.sim.report_continuously = True
        
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        assert self.sim.y_sol[-1][0] == pytest.approx(1.706947, abs = 1e-2)
        
        self.sim_t0.simulate(3.)
        assert self.sim_t0.t_sol[0] == pytest.approx(1.0000000, abs = 1e-4)
        assert self.sim_t0.t_sol[-1] == pytest.approx(3.0000000, abs = 1e-4)
        assert self.sim_t0.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
    
    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(np.diff(self.sim.t_sol))-np.finfo('double').eps <= 0.01

    def test_base_exception_interrupt_fcn(self):
        """Test that BaseExceptions in right-hand side terminate the simulation. Radau5 + C + implicit problem."""
        prob = ImplicitProbBaseException(dim = 2, fcn = True)
        sim = Radau5DAE(prob)

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with pytest.raises(Radau5Error, match = re.escape(err_msg)):
            sim.simulate(1.)

class Test_Radau_Common:
    """
    Tests the common attributes of the Radau solvers.
    """
    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls):
        """
        This sets up the test case.
        """

        f = lambda t,y:[1.0,2.0]

        #Define an explicit Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        exp_mod = Explicit_Problem(f,y0)
        cls.sim = Radau5ODE(exp_mod)

    def test_fac1(self):
        """
        This tests the functionality of the property fac1.
        """
        self.sim.fac1 = 0.01
        assert self.sim.fac1 == 0.01
        self.sim.fac1 = 0.001
        assert self.sim.fac1 == 0.001

        err_msg = "The attribute 'fac1' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fac1 = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fac1 = [-1.0]
    
    def test_fac2(self):
        """
        This tests the functionality of the property fac2.
        """
        self.sim.fac2 = 0.01
        assert self.sim.fac2 == 0.01
        self.sim.fac2 = 0.001
        assert self.sim.fac2 == 0.001
        
        err_msg = "The attribute 'fac2' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fac2 = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fac2 = [-1.0]
    
    def test_fnewt(self):
        """
        This tests the functionality of the property fnewt.
        """
        self.sim.fnewt = 0.01
        assert self.sim.fnewt == 0.01
        self.sim.fnewt = 0.001
        assert self.sim.fnewt == 0.001
        
        err_msg = "The attribute 'fnewt' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fnewt = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.fnewt = [-1.0]
    
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
        self.sim.inith = 0.01
        assert self.sim.inith == 0.01
        self.sim.inith = 0.001
        assert self.sim.inith == 0.001

        err_msg = 'The initial step must be an integer or float.'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.inith = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.inith = [-1.0]
        
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

        err_msg = "The attribute 'newt' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.newt = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.newt = [-1.0]
        
    def test_quot1(self):
        """
        This tests the functionality of the property quot1.
        """
        self.sim.quot1 = 0.01
        assert self.sim.quot1 == 0.01
        self.sim.quot1 = 0.001
        assert self.sim.quot1 == 0.001

        err_msg = "The attribute 'quot1' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.quot1 = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.quot1 = [-1.0]
    
    def test_quot2(self):
        """
        This tests the functionality of the property quot2.
        """
        self.sim.quot2 = 0.01
        assert self.sim.quot2 == 0.01
        self.sim.quot2 = 0.001
        assert self.sim.quot2 == 0.001
        
        err_msg = "The attribute 'quot2' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.quot2 = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.quot2 = [-1.0]
    
    def test_safe(self):
        """
        This tests the functionality of the property safe.
        """
        self.sim.safe = 0.01
        assert self.sim.safe == 0.01
        self.sim.safe = 0.001
        assert self.sim.safe == 0.001

        err_msg = "The attribute 'safe' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.safe = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.safe = [-1.0]
        
    def test_thet(self):
        """
        This tests the functionality of the property thet.
        """
        self.sim.thet = 0.01
        assert self.sim.thet == 0.01
        self.sim.thet = 0.001
        assert self.sim.thet == 0.001
        
        err_msg = "The attribute 'thet' must be an integer or float."
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.thet = 'Test'
        with pytest.raises(Radau_Exception, match = err_msg):
            self.sim.thet = [-1.0]
    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim.usejac = True
        assert self.sim.usejac
        self.sim.usejac = False
        assert not self.sim.usejac
        self.sim.usejac = 1
        assert self.sim.usejac
        self.sim.usejac = []
        assert not self.sim.usejac
