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
from assimulo.solvers.radau5 import Radau5DAE, _Radau5DAE
from assimulo.solvers.radau5 import Radau5ODE, _Radau5ODE
from assimulo.solvers.radau5 import Radau5Error
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.exception import *
from assimulo.lib.radau_core import Radau_Exception, Radau_Common
import scipy.sparse as sp
import numpy as N

import re
float_regex = "[\s]*[\d]*.[\d]*((e|E)(\+|\-)\d\d|)"

class KeyboardInterruptAux:
    """Auxiliary class for creating problems (both explicit and implicit) 
    that simulate a Keyboardinterrupt (Ctrl + c) in the rhs f or jacobian.
    
    Set 'fcn' or 'jac' to True to enable KeyboardInterrupt Exceptions in f or jac."""
    def __init__(self, dim, fcn = False, jac = False, fcn_n = 5):
        self.dim = dim
        self.fcn_raise = fcn
        self.fcn_n = fcn_n
        self.jac_raise = jac
        self.n = 0

    def f(self, t, y):
        self.n += 1
        if self.fcn_raise and (self.n % self.fcn_n == 0):
            raise KeyboardInterrupt()
        else:
            return -y

    def f_impl(self, t, y, yd):
        self.n += 1
        if self.fcn_raise and (self.n % self.fcn_n == 0):
            raise KeyboardInterrupt()
        else:
            return -y
        
    def jac(self, t, y):
        if self.jac_raise:
            raise KeyboardInterrupt()
        else:
            return -N.eye(self.dim)

class Extended_Problem(Explicit_Problem):
    
    #Sets the initial conditions directly into the problem
    y0 = [0.0, -1.0, 0.0]
    sw0 = [False,True,True]
    event_array = N.array([0.0,0.0,0.0])
    rhs_array   = N.array([0.0,0.0,0.0])
    
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
                
            if not True in event_info: #Breaks the iteration loop
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

class Test_Explicit_Radau5_Py:
    """
    Tests the explicit Radau solver (Python implementation).
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
        self.sim = _Radau5ODE(exp_mod) #Create a Radau5 solve
        self.sim_t0 = _Radau5ODE(exp_mod_t0)
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
        self.sim.usejac = False
    
    @testattr(stddist = True)
    def test_event_localizer(self):
        pass
    #     exp_mod = Extended_Problem() #Create the problem

    #     exp_sim = _Radau5ODE(exp_mod) #Create the solver
        
    #     exp_sim.verbosity = 0
    #     exp_sim.report_continuously = True
        
    #     #Simulate
    #     t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
    #     #Basic test
    #     nose.tools.assert_almost_equal(y[-1][0],8.0)
    #     nose.tools.assert_almost_equal(y[-1][1],3.0)
    #     nose.tools.assert_almost_equal(y[-1][2],2.0)
    
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = _Radau5ODE(exp_mod)
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
    
    @testattr(stddist = True)
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = _Radau5ODE(self.mod)
        
        nose.tools.assert_equal(sim._leny, 2)
    
    @testattr(stddist = True)
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.report_continuously = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)
        
        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        print
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim.report_continuously = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
        
    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der Pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)

    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nsteps"], self.sim.statistics["njacs"])
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        nose.tools.assert_less_equal(max(N.diff(self.sim.t_sol))-N.finfo('double').eps, 0.01)
        
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        self.sim.newt = 10
        self.sim.simulate(1.0)
        
        nose.tools.assert_equal(self.sim.statistics["nnfails"], 1)
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        nose.tools.assert_less(self.sim.statistics["nsteps"], 150)
        
    @testattr(stddist = True)
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], steps*1.5)
        
    @testattr(stddist = True)
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
        
        nose.tools.assert_greater(steps2, steps)
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        nose.tools.assert_equal(steps3, steps2)

        err_msg = "atol must be of length one or same as the dimension of the problem."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.atol = [1e-6,1e-6,1e-6]
        
class Test_Explicit_Fortran_Radau5:
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
        
        def jac_sparse(t,y):
            eps = 1.e-6
            my = 1./eps
            J = N.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return sp.csc_matrix(J)
        
        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod_t0 = Explicit_Problem(f,y0,1.0)
        exp_mod_sp = Explicit_Problem(f,y0)
        
        exp_mod.jac = jac
        exp_mod_sp.jac = jac_sparse
        self.mod = exp_mod
            
        #Define an explicit solver
        self.sim = Radau5ODE(exp_mod) #Create a Radau5 solve
        self.sim.implementation = 'f'
        self.sim_t0 = Radau5ODE(exp_mod_t0)
        self.sim_t0.implementation = 'f'
        self.sim_sp = Radau5ODE(exp_mod_sp)
        self.sim_sp.implementation = 'f'
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
        self.sim.usejac = False

    @testattr(stddist = True)
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = Radau5ODE(exp_mod) #Create the solver
        exp_sim.implementation = 'f'
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        nose.tools.assert_almost_equal(y[-1][0],8.0)
        nose.tools.assert_almost_equal(y[-1][1],3.0)
        nose.tools.assert_almost_equal(y[-1][2],2.0)
    
    @testattr(stddist = True)
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'
        
        sim.usejac = False
        sim.simulate(1)
        
        nose.tools.assert_greater(sim.statistics["nfcnjacs"], 0)
        
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'
        sim.simulate(1)
        
        nose.tools.assert_equal(sim.statistics["nfcnjacs"], 0)
    
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = Radau5ODE(exp_mod)
        exp_sim.implementation = 'f'
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
    
    @testattr(stddist = True)
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'
        
        nose.tools.assert_equal(sim._leny, 2)
    
    @testattr(stddist = True)
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.report_continuously = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)
        
        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim.report_continuously = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
        
    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim_sp.usejac = True
        
        self.sim_sp.simulate(2.) #Simulate 2 seconds
    
        nose.tools.assert_equal(self.sim_sp.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim_sp.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nsteps"], self.sim.statistics["njacs"])
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        nose.tools.assert_less_equal(max(N.diff(self.sim.t_sol))-N.finfo('double').eps, 0.01)
        
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        pass
        # self.sim.simulate(1.0)
        # self.sim.reset()
        # self.sim.newt = 10
        # self.sim.simulate(1.0)
        
        # nose.tools.assert_equal(self.sim.statistics["nniterfail"], 1)
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        nose.tools.assert_less(self.sim.statistics["nsteps"], 150)
        
    @testattr(stddist = True)
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], steps*1.5)
        
    @testattr(stddist = True)
    def test_weighted_error(self):
        
        def handle_result(solver, t, y):
            err = solver.get_weighted_local_errors()
            nose.tools.assert_equal(len(err), len(y))
        
        self.mod.handle_result = handle_result
            
        #Define an explicit solver
        sim = Radau5ODE(self.mod) #Create a Radau5 solve
        sim.implementation = 'f'
        
        sim.get_weighted_local_errors()
        
        sim.simulate(1)
        
        
    @testattr(stddist = True)
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
        
        nose.tools.assert_greater(steps2, steps)
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        nose.tools.assert_equal(steps3, steps2)
        
        err_msg = "atol must be of length one or same as the dimension of the problem."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.atol = [1e-6,1e-6,1e-6]
        
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
        
        sim = Radau5ODE(mod)
        sim.implementation = 'f'
        nose.tools.assert_true(sim.sw[0])
        sim.simulate(3)
        nose.tools.assert_false(sim.sw[0])

    @testattr(stddist = True)
    def test_nmax_steps(self):
        """
        This tests the error upon exceeding a set maximum number of steps
        """
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'

        sim.maxh = 1.e-1
        sim.maxsteps = 9

        err_msg = "The solver took max internal steps but could not reach the next output time."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_step_size_too_small(self):
        """
        This tests the error for too small step-sizes
        """
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'

        sim.atol = 1.e10
        sim.rtol = 1.e10

        sim.inith = 1.e-1
        sim.maxh = 1.e-1
        
        err_msg = f"The step size became too small. At time {float_regex}."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1. + 1.e-16)

    @testattr(stddist = True)
    def test_repeated_unexpected_step_rejections(self):
        """
        This tests the error for repeated unexpected step rejections
        """
        def f(t, y):
            raise N.linalg.LinAlgError()
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)
        sim = Radau5ODE(prob)
        sim.implementation = 'f'

        err_msg = f'Repeated unexpected step rejections.'
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_attribute(self):
        """
        This tests the error when trying to simulate using the Fortran based solver with sparse linear solver setting enabled
        """
        sim = Radau5ODE(self.mod)
        sim.implementation = 'f'
        sim.linear_solver = 'SPARSE'

        err_msg = "Sparse Linear solver not supported for Fortran based implementation, instead use 'implementation' = 'c' or 'linear_solver' = 'DENSE'."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_linear_solver(self):
        """
        This tests the functionality of the property linear_solver.
        """
        self.sim.linear_solver = 'dense'
        nose.tools.assert_equal(self.sim.linear_solver, 'DENSE')
        self.sim.linear_solver = 'sparse'
        nose.tools.assert_equal(self.sim.linear_solver, 'SPARSE')
        self.sim.linear_solver = 'DENSE'
        nose.tools.assert_equal(self.sim.linear_solver, 'DENSE')
        self.sim.linear_solver = 'SPARSE'
        nose.tools.assert_equal(self.sim.linear_solver, 'SPARSE')

        err_msg = "'linear_solver' parameter needs to be either 'DENSE' or 'SPARSE'. Set value: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('default')):
            self.sim.linear_solver = 'default'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('GMRES')):
            self.sim.linear_solver = 'GMRES'

        err_msg = "'linear_solver' parameter needs to be the STRING 'DENSE' or 'SPARSE'. Set value: {}, type: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('0', "<class 'int'>")):
            self.sim.linear_solver = 0

    @testattr(stddist = True)
    def test_implementation(self):
        """
        This tests the functionality of the property solver.
        """
        self.sim.implementation = 'f'
        nose.tools.assert_equal(self.sim.implementation, 'f')
        self.sim.implementation = 'c'
        nose.tools.assert_equal(self.sim.implementation, 'c')
        self.sim.implementation = 'f'
        nose.tools.assert_equal(self.sim.implementation, 'f')
        self.sim.implementation = 'c'
        nose.tools.assert_equal(self.sim.implementation, 'c')

        err_msg = "'implementation' parameter needs to be either 'f' or 'c'. Set value: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('Python')):
            self.sim.implementation = 'Python'
        err_msg = "'implementation' parameter needs to be the STRING 'c' or 'f'. Set value: {}, type: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('True', "<class 'bool'>")):
            self.sim.implementation = True
        
    def test_keyboard_interrupt_fcn(self):
        """Test that KeyboardInterrupts in right-hand side terminate the simulation. Radau5 + Fortran + explicit problem."""

        y0 = N.array([1., 1.])
        aux = KeyboardInterruptAux(dim = len(y0), fcn = True)
        prob = Explicit_Problem(aux.f, y0)
        sim = Radau5ODE(prob)
        sim.implementation = 'f'

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_keyboard_interrupt_jac(self):
        """Test that KeyboardInterrupts in jacobian terminate the simulation. Radau5 + Fortran + explicit problem."""

        y0 = N.array([1., 1.])
        aux = KeyboardInterruptAux(dim = len(y0), jac = True)
        prob = Explicit_Problem(aux.f, y0)
        prob.jac = aux.jac
        sim = Radau5ODE(prob)
        sim.usejac = True
        sim.implementation = 'f'

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)


    @testattr(stddist = True)
    def test_sparse_solver_attribute(self):
        """
        This tests the error when trying to simulate using the Fotran based solver with sparse linear solver setting enabled
        """
        sim = Radau5ODE(self.mod)
        sim.solver = 'f'
        sim.linear_solver = 'SPARSE'
        nose.tools.assert_raises(Radau_Exception, sim.simulate, 1.)

class Test_Explicit_C_Radau5:
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
        
        def jac_sparse(t,y):
            eps = 1.e-6
            my = 1./eps
            J = N.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return sp.csc_matrix(J)
        
        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod_t0 = Explicit_Problem(f,y0,1.0)
        exp_mod_sp = Explicit_Problem(f,y0)
        
        exp_mod.jac = jac
        exp_mod_sp.jac = jac_sparse
        self.mod = exp_mod
            
        #Define an explicit solver
        self.sim = Radau5ODE(exp_mod) #Create a Radau5 solve
        self.sim.implementation = 'c'
        self.sim_t0 = Radau5ODE(exp_mod_t0)
        self.sim_t0.implementation = 'c'
        self.sim_sp = Radau5ODE(exp_mod_sp)
        self.sim_sp.implementation = 'c'
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
        self.sim.usejac = False

    @testattr(stddist = True)
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = Radau5ODE(exp_mod) #Create the solver
        exp_sim.implementation = 'c'
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        nose.tools.assert_almost_equal(y[-1][0],8.0)
        nose.tools.assert_almost_equal(y[-1][1],3.0)
        nose.tools.assert_almost_equal(y[-1][2],2.0)
    
    @testattr(stddist = True)
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = Radau5ODE(self.mod)
        sim.implementation = 'c'
        
        sim.usejac = False
        sim.simulate(1)
        
        nose.tools.assert_greater(sim.statistics["nfcnjacs"], 0)
        
        sim = Radau5ODE(self.mod)
        sim.implementation = 'c'
        sim.simulate(1)
        
        nose.tools.assert_equal(sim.statistics["nfcnjacs"], 0)
    
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = Radau5ODE(exp_mod)
        exp_sim.implementation = 'c'
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
    
    @testattr(stddist = True)
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = Radau5ODE(self.mod)
        sim.implementation = 'c'
        
        nose.tools.assert_equal(sim._leny, 2)
    
    @testattr(stddist = True)
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.report_continuously = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)
        
        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim.report_continuously = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
        
    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], 300)

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim_sp.usejac = True
        
        self.sim_sp.simulate(2.) #Simulate 2 seconds
    
        nose.tools.assert_equal(self.sim_sp.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim_sp.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nsteps"], self.sim.statistics["njacs"])
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        nose.tools.assert_less_equal(max(N.diff(self.sim.t_sol))-N.finfo('double').eps, 0.01)
        
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        pass
        # self.sim.simulate(1.0)
        # self.sim.reset()
        # self.sim.newt = 10
        # self.sim.simulate(1.0)
        # nose.tools.assert_equal(self.sim.statistics["nniterfail"], 1)
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        nose.tools.assert_less(self.sim.statistics["nsteps"], 150)
        
    @testattr(stddist = True)
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        nose.tools.assert_less(self.sim.statistics["nsteps"], steps*1.5)
        
    @testattr(stddist = True)
    def test_weighted_error(self):
        
        def handle_result(solver, t, y):
            err = solver.get_weighted_local_errors()
            nose.tools.assert_equal(len(err), len(y))
        
        self.mod.handle_result = handle_result
            
        #Define an explicit solver
        sim = Radau5ODE(self.mod) #Create a Radau5 solve
        sim.implementation = 'c'
        
        sim.get_weighted_local_errors()
        
        sim.simulate(1)
        
        
    @testattr(stddist = True)
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
        
        nose.tools.assert_greater(steps2, steps)
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        nose.tools.assert_equal(steps3, steps2)
        
        err_msg = "atol must be of length one or same as the dimension of the problem."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.atol = [1e-6,1e-6,1e-6]
        
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
        
        sim = Radau5ODE(mod)
        sim.implementation = 'c'
        nose.tools.assert_true(sim.sw[0])
        sim.simulate(3)
        nose.tools.assert_false(sim.sw[0])

    @testattr(stddist = True)
    def test_nmax_steps(self):
        """
        This tests the error upon exceeding a set maximum number of steps
        """
        sim = Radau5ODE(self.mod)
        sim.implementation = 'c'

        sim.maxh = 1.e-1
        sim.maxsteps = 9

        err_msg = "The solver took max internal steps but could not reach the next output time."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_step_size_too_small(self):
        """
        This tests the error for too small step-sizes
        """
        sim = Radau5ODE(self.mod)
        sim.implementation = 'c'

        sim.atol = 1.e10
        sim.rtol = 1.e10

        sim.inith = 1.e-1
        sim.maxh = 1.e-1

        err_msg = f"The step size became too small. At time {float_regex}."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1. + 1.e-16)

    @testattr(stddist = True)
    def test_repeated_unexpected_step_rejections(self):
        """
        This tests the error for repeated unexpected step rejections
        """
        def f(t, y):
            raise N.linalg.LinAlgError()
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)
        sim = Radau5ODE(prob)
        sim.implementation = 'c'

        err_msg = f'Repeated unexpected step rejections.'
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_jac_disabled(self):
        """
        This tests the trying to simulate using the sparse linear solver, with no analytical jacobian provided.
        """
        f = lambda t, y: [y]
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = False

        sim.simulate(1.)
        nose.tools.assert_equal(sim.linear_solver, 'DENSE')

    @testattr(stddist = True)
    def test_solver_no_jac(self):
        """
        This tests the error when trying to simulate using an analytical jacobian, with none provided
        """
        f = lambda t, y: [y]
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.usejac = True

        err_msg = "Use of an analytical Jacobian is enabled, but problem does contain a 'jac' function."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_solver_sparse_jac_wrong_format(self):
        """
        This tests the error when using a sparse jacobian of the wrong format
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sp.spdiags([1], 0, 1, 1, format = 'csr')
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 1

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = 'Sparse Jacobian given in wrong format, expects CSC format.'
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_solver_sparse_jac_nnz_too_small(self):
        """
        This tests the error when using a sparse jacobian with nnz larger than specified
        """
        n = 5
        f = lambda t, y: y
        jac = lambda t, y: sp.eye(n, n, dtype = N.double, format = 'csc')
        y0 = N.array([1.]*n)
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 1

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = 'Failure in sparse Jacobian evaluation, specified number of nonzero elements too small.'
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_solver_sparse_jac_nnz_zero(self):
        """
        This tests that using a sparse jacobian with nnz = 0 is valid.
        """
        n = 5
        f = lambda t, y: [0.]*n
        jac = lambda t, y: sp.csc_matrix((n, n), dtype = N.double)
        y0 = N.array([1.]*n)
        prob = Explicit_Problem(f, y0)
        prob.jac = jac
        prob.jac_nnz = 0

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_no_nnz(self):
        """
        This tests the error when trying to simulate using the sparse linear solver, without specifying the number of non-zero elements
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sp.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = N.array([1.])
        prob = Explicit_Problem(f, y0)
        prob.jac = jac

        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = "Number of non-zero elements of sparse Jacobian must be non-negative. Detected default value of '-1', has 'problem.jac_fcn_nnz' been set?"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_invalid_nnz_type(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; wrong type.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sp.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = N.array([1.])

        for nnz in [None, "test"]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.implementation = 'c'
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian must be an integer, received: {}."
            with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format(nnz)):
                sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_invalid_nnz_negative(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; negative.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sp.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = N.array([1.])

        for nnz in [-2, -10]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.implementation = 'c'
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian must be non-negative, given value = {}."
            with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format(nnz)):
                sim.simulate(1.)

    @testattr(stddist = True)
    def test_sparse_solver_invalid_nnz_too_large(self):
        """
        This tests the error when trying to simulate using the sparse linear solver with invalid inputs for nnz; too_large.
        """
        f = lambda t, y: [y]
        jac = lambda t, y: sp.spdiags([1], 0, 1, 1, format = 'csc')
        y0 = N.array([1.])

        for nnz in [5, 100]:
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.implementation = 'c'
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            err_msg = "Number of non-zero elements of sparse Jacobian infeasible, must be smaller than the problem dimension squared."
            with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
                sim.simulate(1.)

    def test_sparse_solver_jacobian(self):
        """Testing sparse solver to not produce segmentation faults for Jacobian."""
        ## Take trivial problem with somewhat arbitrary jacobians
        ## Test that functions for internal processing of jacobian do not produces segfaults
        jacobians = [
            (lambda t, y: sp.csc_matrix(N.array([[1., 1., 1.], [1., 1., 1.], [1., 1., 1.]])), 9), 
            (lambda t, y: sp.csc_matrix(N.array([[0., 1., 1.], [1., 0., 1.], [1., 1., 0.]])), 6),
            (lambda t, y: sp.csc_matrix(N.array([[0., 1., 1.], [1., 1., 1.], [1., 1., 1.]])), 8),
            (lambda t, y: sp.csc_matrix(N.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sp.csc_matrix(N.array([[0., 0., 0.], [1., 0., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sp.csc_matrix(N.array([[0., 0., 0.], [0., 0., 0.], [0., 1., 0.]])), 1),
            (lambda t, y: sp.csc_matrix(N.array([[0., 0., 1.], [0., 0., 0.], [0., 0., 0.]])), 1),
            (lambda t, y: sp.csc_matrix(N.array([[1., 0., 0.], [0., 0., 0.], [0., 0., 0.]])), 1),
        ]

        for i, (jac, nnz) in enumerate(jacobians):
            f = lambda t, y: y
            y0 = 1.*N.ones(3)
            prob = Explicit_Problem(f, y0)
            prob.jac = jac
            prob.jac_nnz = nnz

            sim = Radau5ODE(prob)
            sim.implementation = 'c'
            sim.linear_solver = 'SPARSE'
            sim.usejac = True

            nose.tools.ok_(sim.simulate(1.), msg = f"Jacobian #{i} failed: {jac(0, 0)}")

    @testattr(stddist = True)
    def test_linear_solver(self):
        """
        This tests the functionality of the property linear_solver.
        """
        self.sim.linear_solver = 'dense'
        nose.tools.assert_equal(self.sim.linear_solver, 'DENSE')
        self.sim.linear_solver = 'sparse'
        nose.tools.assert_equal(self.sim.linear_solver, 'SPARSE')
        self.sim.linear_solver = 'DENSE'
        nose.tools.assert_equal(self.sim.linear_solver, 'DENSE')
        self.sim.linear_solver = 'SPARSE'
        nose.tools.assert_equal(self.sim.linear_solver, 'SPARSE')

        err_msg = "'linear_solver' parameter needs to be either 'DENSE' or 'SPARSE'. Set value: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('default')):
            self.sim.linear_solver = 'default'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('GMRES')):
            self.sim.linear_solver = 'GMRES'

        err_msg = "'linear_solver' parameter needs to be the STRING 'DENSE' or 'SPARSE'. Set value: {}, type: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('0', "<class 'int'>")):
            self.sim.linear_solver = 0

    @testattr(stddist = True)
    def test_implementation(self):
        """
        This tests the functionality of the property solver.
        """
        self.sim.implementation = 'f'
        nose.tools.assert_equal(self.sim.implementation, 'f')
        self.sim.implementation = 'c'
        nose.tools.assert_equal(self.sim.implementation, 'c')
        self.sim.implementation = 'f'
        nose.tools.assert_equal(self.sim.implementation, 'f')
        self.sim.implementation = 'c'
        nose.tools.assert_equal(self.sim.implementation, 'c')

        err_msg = "'implementation' parameter needs to be either 'f' or 'c'. Set value: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('Python')):
            self.sim.implementation = 'Python'
        err_msg = "'implementation' parameter needs to be the STRING 'c' or 'f'. Set value: {}, type: {}"
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg.format('True', "<class 'bool'>")):
            self.sim.implementation = True

    @testattr(stddist = True)
    def test_keyboard_interrupt_fcn(self):
        """Test that KeyboardInterrupts in right-hand side terminate the simulation. Radau5 + C + explicit problem."""

        y0 = N.array([1., 1.])
        aux = KeyboardInterruptAux(dim = len(y0), fcn = True)
        prob = Explicit_Problem(aux.f, y0)
        sim = Radau5ODE(prob)
        sim.implementation = 'c'

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_keyboard_interrupt_jac(self):
        """Test that KeyboardInterrupts in jacobian terminate the simulation. Radau5 + C + explicit problem."""

        y0 = N.array([1., 1.])
        aux = KeyboardInterruptAux(dim = len(y0), jac = True)
        prob = Explicit_Problem(aux.f, y0)
        prob.jac = aux.jac
        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.usejac = True

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_keyboard_interrupt_jac_sparse(self):
        """Test that KeyboardInterrupts in jacobian terminate the simulation. Radau5 + C + explicit problem + sparse jac."""

        y0 = N.array([1., 1.])
        aux = KeyboardInterruptAux(dim = len(y0), jac = True)
        prob = Explicit_Problem(aux.f, y0)
        prob.jac = aux.jac
        prob.jac_nnz = 1
        sim = Radau5ODE(prob)
        sim.implementation = 'c'
        sim.linear_solver = 'SPARSE'
        sim.usejac = True

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)


class Test_Implicit_Radau5:
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
        
        y0 = [2.0,-0.6] #Initial conditions
        yd0 = [-.6,-200000.]
        
        #Define an Assimulo problem
        self.mod = Implicit_Problem(f,y0,yd0)
        self.mod_t0 = Implicit_Problem(f,y0,yd0,1.0)
            
        #Define an implicit solver
        self.sim = Radau5DAE(self.mod) #Create a Radau5 solve
        self.sim_t0 = Radau5DAE(self.mod_t0)
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size

    @testattr(stddist = True)
    def test_implementation_get(self):
        """
            Test getting of implementation property of Radau5DAE.
        """
        nose.tools.assert_equal(self.sim.implementation, 'f')

    @testattr(stddist = True)
    def test_implementation_set(self):
        """
            Test setting of implementation property of Radau5DAE.
        """
        err_msg = "Radau5DAE does not support setting the 'implementation' attribute, since it only supports the Fortran implementation of Radau5."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.implementation = 'c'

    @testattr(stddist = True)
    def test_linear_solver_get(self):
        """
            Test getting of linear_solver property of Radau5DAE.
        """
        nose.tools.assert_equal(self.sim.linear_solver, 'DENSE')

    @testattr(stddist = True)
    def test_linear_solver_set(self):
        """
            Test setting of linear_solver property of Radau5DAE.
        """
        err_msg = "Radau5DAE does not support setting the 'linear_solver' attribute, since it only supports the DENSE linear solver in Fortran implementation of Radau5."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.linear_solver = 'SPARSE'
    
    @testattr(stddist = True)
    def test_nbr_fcn_evals_due_to_jac(self):
        sim = Radau5DAE(self.mod)
        
        sim.usejac = False
        sim.simulate(1)
        
        nose.tools.assert_greater(sim.statistics["nfcnjacs"], 0)
    
    @testattr(stddist = True)
    def test_simulate_explicit(self):
        """
        Test a simulation of an explicit problem using Radau5DAE.
        """
        f = lambda t,y:N.array(-y)
        y0 = [1.0]
        
        problem = Explicit_Problem(f,y0)
        simulator = Radau5DAE(problem)
        
        nose.tools.assert_equal(simulator.yd0[0], -simulator.y0[0])
        
        t,y = simulator.simulate(1.0)
        
        nose.tools.assert_almost_equal(float(y[-1]), float(N.exp(-1.0)),4)
    
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event

        #CVode
        exp_sim = Radau5DAE(exp_mod)
        exp_sim.verbosity = 0
        exp_sim(5.,100)

        nose.tools.assert_equal(nevent, 5)

    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of Radau5 Implicit Init.
        """
        #Test both y0 in problem and not.

        sim = Radau5DAE(self.mod)
        
        nose.tools.assert_equal(sim._leny, 2)
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(.5) #Simulate 2 seconds
        nose.tools.assert_equal(self.sim.statistics["nsteps"], self.sim.statistics["njacs"])

    @testattr(stddist = True)    
    def test_simulation(self):
        """
        Test a simulation of the van der Pol equations (1).
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706272, 3)

        self.sim.reset()

        self.sim.report_continuously = True

        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706166, 3)

        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)

    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True

        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)

        self.sim.reset()
        self.sim.report_continuously = False

        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)


    @testattr(stddist = True)
    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        nose.tools.assert_less_equal(max(N.diff(self.sim.t_sol))-N.finfo('double').eps, 0.01)
        
    @testattr(stddist = True)
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        res = lambda t,x,xd,sw: N.array([1.0 - xd])
        state_events = lambda t,x,xd,sw: N.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Implicit_Problem(res,[0.0], [1.0])
        mod.sw0 = [True]

        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = Radau5DAE(mod)
        nose.tools.assert_true(sim.sw[0])
        sim.simulate(3)
        nose.tools.assert_false(sim.sw[0])

    @testattr(stddist = True)
    def test_nmax_steps(self):
        """
        This tests the error upon exceeding a set maximum number of steps
        """
        sim = Radau5DAE(self.mod)

        sim.maxh = 1.e-1
        sim.maxsteps = 9

        err_msg = "The solver took max internal steps but could not reach the next output time."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)

    @testattr(stddist = True)
    def test_step_size_too_small(self):
        """
        This tests the error for too small step-sizes
        """
        f = lambda t, y, yd: -y
        y0 = N.array([1.])
        yd0 = N.array([0.])

        prob = Implicit_Problem(f, y0, yd0)

        sim = Radau5DAE(prob)

        sim.atol = 1.e10
        sim.rtol = 1.e10

        sim.inith = 1.e-1
        sim.maxh = 1.e-1
        
        err_msg = f"The step size became too small. At time {float_regex}."
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1. + 1.e-16)

    @testattr(stddist = True)
    def test_repeated_unexpected_step_rejections(self):
        """
        This tests the error for repeated unexpected step rejections in Radau5DAE.
        """
        def f(t, y, yd):
            raise N.linalg.LinAlgError()
        prob = Implicit_Problem(f, N.array([1.]), N.array([1.]))
        sim = Radau5DAE(prob)

        err_msg = f'Repeated unexpected step rejections.'
        with nose.tools.assert_raises_regex(Radau5Error, err_msg):
            sim.simulate(1.)


class Test_Implicit_Radau5_Py:
    """
    Tests the implicit Radau solver (Python implementation).
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
        
        y0 = [2.0,-0.6] #Initial conditions
        yd0 = [-.6,-200000.]
        
        #Define an Assimulo problem
        self.mod = Implicit_Problem(f,y0,yd0)
        self.mod_t0 = Implicit_Problem(f,y0,yd0,1.0)
            
        #Define an explicit solver
        self.sim = _Radau5DAE(self.mod) #Create a Radau5 solve
        self.sim_t0 = _Radau5DAE(self.mod_t0)
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
    
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = _Radau5DAE(exp_mod)
        exp_sim.verbosity = 0
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
    
    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of Radau5 Implicit Init.
        """
        #Test both y0 in problem and not.

        sim = _Radau5DAE(self.mod)
        
        nose.tools.assert_equal(sim._leny, 2)
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(.5) #Simulate 2 seconds

        nose.tools.assert_equal(self.sim.statistics["nsteps"], self.sim.statistics["njacs"])
        
    @testattr(stddist = True)    
    def test_simulation(self):
        """
        Test a simulation of the van der Pol equations (2).
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706272, 3)
        
        self.sim.reset()
        
        self.sim.report_continuously = True
        
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706947, 2)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.report_continuously = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
        
        self.sim.reset()
        self.sim.report_continuously = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        nose.tools.assert_equal(len(self.sim.t_sol), 201)
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        nose.tools.assert_less_equal(max(N.diff(self.sim.t_sol))-N.finfo('double').eps, 0.01)

    @testattr(stddist = True)
    def test_keyboard_interrupt_fcn(self):
        """Test that KeyboardInterrupts in right-hand side terminate the simulation. Radau5 + C + implicit problem."""

        y0 = N.array([1., 1.])
        yd = N.array([0., 0.])
        aux = KeyboardInterruptAux(dim = len(y0), fcn = True)
        prob = Implicit_Problem(aux.f_impl, y0, yd)
        sim = Radau5DAE(prob)

        err_msg = "Unrecoverable exception encountered during callback to problem (right-hand side/jacobian)."
        with nose.tools.assert_raises_regex(Radau5Error, re.escape(err_msg)):
            sim.simulate(1.)

class Test_Radau_Common:
    """
    Tests the common attributes of the Radau solvers.
    """
    def setUp(self):
        """
        This sets up the test case.
        """

        f = lambda t,y:[1.0,2.0]

        #Define an explicit Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        exp_mod = Explicit_Problem(f,y0)
        self.sim = Radau5ODE(exp_mod)

    @testattr(stddist = True)    
    def test_fac1(self):
        """
        This tests the functionality of the property fac1.
        """
        self.sim.fac1 = 0.01
        nose.tools.assert_equal(self.sim.fac1, 0.01)
        self.sim.fac1 = 0.001
        nose.tools.assert_equal(self.sim.fac1, 0.001)

        err_msg = "The attribute 'fac1' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fac1 = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fac1 = [-1.0]
    
    @testattr(stddist = True)
    def test_fac2(self):
        """
        This tests the functionality of the property fac2.
        """
        self.sim.fac2 = 0.01
        nose.tools.assert_equal(self.sim.fac2, 0.01)
        self.sim.fac2 = 0.001
        nose.tools.assert_equal(self.sim.fac2, 0.001)
        
        err_msg = "The attribute 'fac2' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fac2 = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fac2 = [-1.0]
    
    @testattr(stddist = True)
    def test_fnewt(self):
        """
        This tests the functionality of the property fnewt.
        """
        self.sim.fnewt = 0.01
        nose.tools.assert_equal(self.sim.fnewt, 0.01)
        self.sim.fnewt = 0.001
        nose.tools.assert_equal(self.sim.fnewt, 0.001)
        
        err_msg = "The attribute 'fnewt' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fnewt = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.fnewt = [-1.0]
    
    @testattr(stddist = True)
    def test_h(self):
        """
        This tests the functionality of the property h.
        """
        self.sim.h = 0.01
        nose.tools.assert_equal(self.sim.h, 0.01)
        self.sim.h = 0.001
        nose.tools.assert_equal(self.sim.h, 0.001)
    
    @testattr(stddist = True)
    def test_initial_step(self):
        """
        This tests the functionality of the property initial step.
        """
        self.sim.inith = 0.01
        nose.tools.assert_equal(self.sim.inith, 0.01)
        self.sim.inith = 0.001
        nose.tools.assert_equal(self.sim.inith, 0.001)

        err_msg = 'The initial step must be an integer or float.'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.inith = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.inith = [-1.0]
        
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the functionality of the property newt.
        """
        self.sim.newt = 1
        nose.tools.assert_equal(self.sim.newt, 1)
        self.sim.newt = 10
        nose.tools.assert_equal(self.sim.newt, 10)
        self.sim.newt = 9.8
        nose.tools.assert_equal(self.sim.newt, 9)

        err_msg = "The attribute 'newt' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.newt = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.newt = [-1.0]
        
    @testattr(stddist = True)
    def test_quot1(self):
        """
        This tests the functionality of the property quot1.
        """
        self.sim.quot1 = 0.01
        nose.tools.assert_equal(self.sim.quot1, 0.01)
        self.sim.quot1 = 0.001
        nose.tools.assert_equal(self.sim.quot1, 0.001)

        err_msg = "The attribute 'quot1' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.quot1 = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.quot1 = [-1.0]
    
    @testattr(stddist = True)    
    def test_quot2(self):
        """
        This tests the functionality of the property quot2.
        """
        self.sim.quot2 = 0.01
        nose.tools.assert_equal(self.sim.quot2, 0.01)
        self.sim.quot2 = 0.001
        nose.tools.assert_equal(self.sim.quot2, 0.001)
        
        err_msg = "The attribute 'quot2' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.quot2 = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.quot2 = [-1.0]
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the functionality of the property safe.
        """
        self.sim.safe = 0.01
        nose.tools.assert_equal(self.sim.safe, 0.01)
        self.sim.safe = 0.001
        nose.tools.assert_equal(self.sim.safe, 0.001)

        err_msg = "The attribute 'safe' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.safe = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.safe = [-1.0]
        
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests the functionality of the property thet.
        """
        self.sim.thet = 0.01
        nose.tools.assert_equal(self.sim.thet, 0.01)
        self.sim.thet = 0.001
        nose.tools.assert_equal(self.sim.thet, 0.001)
        
        err_msg = "The attribute 'thet' must be an integer or float."
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.thet = 'Test'
        with nose.tools.assert_raises_regex(Radau_Exception, err_msg):
            self.sim.thet = [-1.0]
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim.usejac = True
        nose.tools.assert_true(self.sim.usejac)
        self.sim.usejac = False
        nose.tools.assert_false(self.sim.usejac)
        self.sim.usejac = 1
        nose.tools.assert_true(self.sim.usejac)
        self.sim.usejac = []
        nose.tools.assert_false(self.sim.usejac)
