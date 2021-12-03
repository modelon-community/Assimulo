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
from assimulo.solvers.sundials import *
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.exception import *
import numpy as np
import scipy.sparse as sp

class Extended_Problem(Explicit_Problem):
    
    #Sets the initial conditons directly into the problem
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

class Test_CVode:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y:N.array(y)
        y0 = [1.0]
        
        self.problem = Explicit_Problem(f,y0)
        self.simulator = CVode(self.problem)
        self.simulator.verbosity = 0
    
    @testattr(stddist = True)
    def test_backward_integration(self):
        def f(t, y):
            x, v = y
            return [x, -x - 0.1*v]
            
        mod = Explicit_Problem(f, y0=[1, 0], t0=10)
        sim = CVode(mod)
        sim.backward = True
        t, y = sim.simulate(0, ncp_list=np.arange(1, 10))
        
        nose.tools.assert_true(np.all(t == np.arange(0,11)[::-1]))
        
        mod = Explicit_Problem(f, y0=[1, 0], t0=10)
        sim = CVode(mod)
        sim.backward = True
        t, y = sim.simulate(0, ncp_list=np.arange(1, 10)[::-1])
        
        nose.tools.assert_true(np.all(t == np.arange(0,11)[::-1]))

    @testattr(stddist = True)
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = CVode(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        nose.tools.assert_almost_equal(y[-1][0],8.0)
        nose.tools.assert_almost_equal(y[-1][1],3.0)
        nose.tools.assert_almost_equal(y[-1][2],2.0)
    
    @testattr(stddist = True)
    def test_get_error_weights(self):
        nose.tools.assert_raises(CVodeError, self.simulator.get_error_weights)
        
        self.simulator.simulate(1.0)
        
        weights = self.simulator.get_error_weights()
        nose.tools.assert_less(weights[0], 1e6)
        
    @testattr(stddist = True)
    def test_get_used_initial_step(self):
        self.simulator.simulate(1.0)
        
        step = self.simulator.get_used_initial_step()
        nose.tools.assert_almost_equal(step, 0.001, 3)
        
        self.simulator.reset()
        
        self.simulator.inith = 1e-8
        self.simulator.simulate(1.0)
        
        step = self.simulator.get_used_initial_step()
        nose.tools.assert_less(N.abs(step-1e-8), 1e-2)
        
    
    @testattr(stddist = True)
    def test_get_local_errors(self):
        nose.tools.assert_raises(CVodeError, self.simulator.get_local_errors)
    
        self.simulator.simulate(1.0)
        
        err = self.simulator.get_local_errors()
        nose.tools.assert_less(err[0], 1e-5)
    
    @testattr(stddist = True)
    def test_get_last_order(self):
        nose.tools.assert_raises(CVodeError, self.simulator.get_last_order)
        
        self.simulator.simulate(1.0)
        
        qlast = self.simulator.get_last_order()
        nose.tools.assert_equal(qlast, 4)
        
    @testattr(stddist = True)
    def test_max_convergence_failures(self):
        nose.tools.assert_equal(self.simulator.maxncf, self.simulator.options["maxncf"])
        self.simulator.maxncf = 15
        nose.tools.assert_equal(self.simulator.maxncf, 15)
        
        nose.tools.assert_raises(AssimuloException, self.simulator._set_max_conv_fails, -1)
        
    @testattr(stddist = True)
    def test_max_error_tests_failures(self):
        nose.tools.assert_equal(self.simulator.maxnef, self.simulator.options["maxnef"])
        self.simulator.maxnef = 15
        nose.tools.assert_equal(self.simulator.maxnef, 15)
        nose.tools.assert_equal(self.simulator.options["maxnef"], 15)
        
        nose.tools.assert_raises(AssimuloException, self.simulator._set_max_err_fails, -1)
        
    @testattr(stddist = True)
    def test_max_nonlinear_iterations(self):
        nose.tools.assert_equal(self.simulator.maxcor, self.simulator.options["maxcor"])
        self.simulator.maxcor = 15
        nose.tools.assert_equal(self.simulator.maxcor, 15)
        nose.tools.assert_equal(self.simulator.options["maxcor"], 15)
        
        #nose.tools.assert_raises(AssimuloException, self.simulator._set_max_err_fails, -1)
        
    @testattr(stddist = True)
    def test_get_current_order(self):  
        
        nose.tools.assert_raises(CVodeError, self.simulator.get_current_order)

        self.simulator.simulate(1.0)
        
        qcur = self.simulator.get_current_order()
        nose.tools.assert_equal(qcur, 4)


        
    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        # nose.tools.assert_equal(self.simulator.f, 'Test function')
        nose.tools.assert_equal(self.simulator.y, 1.0)
        nose.tools.assert_equal(self.simulator.discr, 'BDF')
        nose.tools.assert_equal(self.simulator.iter, 'Newton')
        nose.tools.assert_equal(self.simulator.maxord, 5)
        
        self.simulator.discr = 'Adams'
        nose.tools.assert_equal(self.simulator.discr, 'Adams')
        nose.tools.assert_equal(self.simulator.maxord, 12)
    
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
        exp_sim = CVode(exp_mod)
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
        
    @testattr(stddist = True)
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
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = CVode(exp_mod)
        exp_sim.verbosity = 10
        exp_sim(5.,100)
        
        nose.tools.assert_equal(len(exp_sim.event_data), 4)
        
        tnext = 0.0
        nevent = 0
        
        exp_sim.reset()
        nose.tools.assert_equal(len(exp_sim.event_data), 0)
        
        exp_sim(5.,100)
        nose.tools.assert_equal(len(exp_sim.event_data), 4)
    
    @testattr(stddist = True)
    def test_time_limit(self):
        f = lambda t,y: -y
        
        exp_mod = Explicit_Problem(f,1.0)
        exp_sim = CVode(exp_mod)
        
        exp_sim.maxh = 1e-8
        exp_sim.time_limit = 1 #One second
        exp_sim.report_continuously = True
        
        nose.tools.assert_raises(TimeLimitExceeded, exp_sim.simulate, 1)
    
    @testattr(stddist = True)    
    def test_discr_method(self):
        """
        This tests the functionality of the property discr.
        """
        
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, 'Test')
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, 1)
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, [1.0, 1])
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, {'Test':'case'})
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, 5.1)
        nose.tools.assert_raises(Exception, self.simulator._set_discr_method, ['Test'])
        
        self.simulator.discr = 'BDF'
        nose.tools.assert_equal(self.simulator.discr, 'BDF')
        self.simulator.discr = 'Adams'
        nose.tools.assert_equal(self.simulator.discr, 'Adams')
    
    @testattr(stddist = True)
    def test_change_discr(self):
        """
        This tests that the change from Functional to Newton works
        """
        f = lambda t,y: N.array([1.0])
        y0 = 4.0 #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.iter = "FixedPoint"
        exp_sim.simulate(1)
        nose.tools.assert_equal(exp_sim.statistics["njacs"], 0)
        exp_sim.iter = "Newton"
        exp_sim.simulate(2)
        nose.tools.assert_greater(exp_sim.statistics["njacs"], 0)
        
    @testattr(stddist = True)
    def test_change_norm(self):
        
        nose.tools.assert_equal(self.simulator.options["norm"], "WRMS")
        self.simulator.norm = 'WRMS'
        nose.tools.assert_equal(self.simulator.norm, 'WRMS')
        self.simulator.norm = 'EUCLIDEAN'
        nose.tools.assert_equal(self.simulator.options["norm"], "EUCLIDEAN")
        nose.tools.assert_equal(self.simulator.norm, 'EUCLIDEAN')
        
        f = lambda t,y: N.array([1.0])
        y0 = 4.0 #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.norm = "WRMS"
        exp_sim.simulate(1)
        
        exp_mod = Explicit_Problem(f,y0)
        exp_sim = CVode(exp_mod) #Create a CVode solver
        
        exp_sim.norm = "EUCLIDEAN"
        exp_sim.simulate(1)
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: N.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: N.array([[0.,1.],[0.,0.]]) #Defines the jacobian
        
        exp_mod = Explicit_Problem(f, [1.0,0.0])
        exp_mod.jac = jac
        
        exp_sim = CVode(exp_mod)
        exp_sim.discr='BDF'
        exp_sim.iter='Newton'
        exp_sim.simulate(5.,100)
        
        nose.tools.assert_equal(exp_sim.statistics["nfcnjacs"], 0)
        nose.tools.assert_almost_equal(exp_sim.y_sol[-1][0], -121.75000143, 4)
        
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        nose.tools.assert_almost_equal(exp_sim.y_sol[-1][0], -121.75000143, 4)
        nose.tools.assert_greater(exp_sim.statistics["nfcnjacs"], 0)
    
    @testattr(stddist = True)
    def test_usejac_csc_matrix(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: N.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: sp.csc_matrix(N.array([[0.,1.],[0.,0.]])) #Defines the jacobian
        
        exp_mod = Explicit_Problem(f, [1.0,0.0])
        exp_mod.jac = jac
        
        exp_sim = CVode(exp_mod)
        exp_sim.discr='BDF'
        exp_sim.iter='Newton'
        exp_sim.simulate(5.,100)
        
        nose.tools.assert_equal(exp_sim.statistics["nfcnjacs"], 0)
        nose.tools.assert_almost_equal(exp_sim.y_sol[-1][0], -121.75000143, 4)
        
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        nose.tools.assert_almost_equal(exp_sim.y_sol[-1][0], -121.75000143, 4)
        nose.tools.assert_greater(exp_sim.statistics["nfcnjacs"], 0)
    
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
        
        sim = CVode(mod)
        nose.tools.assert_true(sim.sw[0])
        sim.simulate(3)
        nose.tools.assert_false(sim.sw[0])
    
    @testattr(stddist = True)
    def test_iter_method(self):
        """
        This tests the functionality of the property iter.
        """
        
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, 'Test')
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, 1)
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, 0)
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, ['Test'])
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, [1.0, 1])
        nose.tools.assert_raises(Exception, self.simulator._set_iter_method, 11.1)
        
        self.simulator.iter = 'Newton'
        nose.tools.assert_equal(self.simulator.iter, 'Newton')
        self.simulator.iter = 'FixedPoint'
        nose.tools.assert_equal(self.simulator.iter, 'FixedPoint')
    
    @testattr(stddist = True)
    def test_initial_step(self):
        """
        This tests the functionality of the property initstep.
        """
        
        nose.tools.assert_raises(Exception, self.simulator._set_initial_step, 'Test')
        nose.tools.assert_raises(Exception, self.simulator._set_initial_step, ['Test'])
        
        nose.tools.assert_equal(self.simulator.inith, 0.0)
        self.simulator.inith = 10.0
        nose.tools.assert_equal(self.simulator.inith, 10.0)
        self.simulator.inith = 1
        nose.tools.assert_equal(self.simulator.inith, 1.0)
    
    @testattr(stddist = True)
    def test_interpolate(self):
        """
        This tests the functionality of the method interpolate.
        """
        f = lambda t,x: N.array(x**0.25)
        
        prob = Explicit_Problem(f,[1.0])

        sim = CVode(prob)
        sim.simulate(10., 100)
        y100 = sim.y_sol
        t100 = sim.t_sol
        sim.reset()
        sim.simulate(10.)
        nose.tools.assert_almost_equal(float(y100[-2]), float(sim.interpolate(9.9,0)),5)
    
    @testattr(stddist = True)
    def test_ncp_list(self):
        f = lambda t,y:N.array(-y)
        y0 = [4.0]
        
        prob = Explicit_Problem(f,y0)
        sim = CVode(prob)
        
        t, y = sim.simulate(7, ncp_list=N.arange(0, 7, 0.1)) #Simulate 5 seconds
        
        nose.tools.assert_almost_equal(float(y[-1]), 0.00364832, 4)
        
    @testattr(stddist = True)
    def test_handle_result(self):
        """
        This function tests the handle result.
        """
        f = lambda t,x: x**0.25
        def handle_result(solver,t,y):
            nose.tools.assert_equal(solver.t, t)
        
        prob = Explicit_Problem(f, [1.0])
        prob.handle_result = handle_result
        
        sim = CVode(prob)
        sim.report_continuously = True
        sim.simulate(10.)
    
    @testattr(stddist = True)    
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        self.simulator.discr='Adams'
        
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, "Test")
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, [1,1])
        
        self.simulator.maxord = -1
        nose.tools.assert_equal(self.simulator.maxord, 1)
        self.simulator.maxord = 2
        nose.tools.assert_equal(self.simulator.maxord, 2)
        self.simulator.maxord = 13
        nose.tools.assert_equal(self.simulator.maxord, 12)
        
        self.simulator.discr='BDF'
        
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, "Test")
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, [1,1])
        
        self.simulator.maxord = -1
        nose.tools.assert_equal(self.simulator.maxord, 1)
        self.simulator.maxord = 2
        nose.tools.assert_equal(self.simulator.maxord, 2)
        self.simulator.maxord = 6
        nose.tools.assert_equal(self.simulator.maxord, 5)
    
    @testattr(stddist = True)
    def test_spgmr(self):
        f = lambda t,y: N.array([y[1], -9.82])
        fsw = lambda t,y,sw: N.array([y[1], -9.82])
        fp = lambda t,y,p: N.array([y[1], -9.82])
        fswp = lambda t,y,sw,p: N.array([y[1], -9.82])
        jacv = lambda t,y,fy,v: N.dot(N.array([[0,1.],[0,0]]),v)
        jacvsw = lambda t,y,fy,v,sw: N.dot(N.array([[0,1.],[0,0]]),v)
        jacvp = lambda t,y,fy,v,p: N.dot(N.array([[0,1.],[0,0]]),v)
        jacvswp = lambda t,y,fy,v,sw,p: N.dot(N.array([[0,1.],[0,0]]),v)
        y0 = [1.0,0.0] #Initial conditions
        
        def run_sim(exp_mod):
            exp_sim = CVode(exp_mod) #Create a CVode solver
            exp_sim.linear_solver = 'SPGMR' #Change linear solver

            #Simulate
            t, y = exp_sim.simulate(5, 1000) #Simulate 5 seconds with 1000 communication points
        
            #Basic tests
            nose.tools.assert_almost_equal(y[-1][0],-121.75000000,4)
            nose.tools.assert_almost_equal(y[-1][1],-49.100000000)
        
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
        nose.tools.assert_raises(CVodeError,run_sim,exp_mod)
        
        exp_mod = Explicit_Problem(fswp,y0,sw0=[True],p0=1.0)
        exp_mod.jacv = jacvsw #Sets the jacobian
        nose.tools.assert_raises(CVodeError,run_sim,exp_mod)
        
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
    
    @testattr(stddist = True)
    def test_max_order_discr(self):
        """
        This tests the maximum order when the discretization is changed.
        """
        self.simulator.discr = "Adams"
        self.simulator.maxord = 7
        nose.tools.assert_equal(self.simulator.maxord, 7)
        
        self.simulator.discr = 'Adams'
        nose.tools.assert_equal(self.simulator.maxord, 12)
        self.simulator.discr = 'BDF'
        nose.tools.assert_equal(self.simulator.maxord, 5)
        self.simulator.discr = 'Adams'
        nose.tools.assert_equal(self.simulator.maxord, 12)
        self.simulator.maxord = 4
        self.simulator.discr = 'BDF'
        nose.tools.assert_equal(self.simulator.maxord, 5)
        self.simulator.discr = 'Adams'
        nose.tools.assert_equal(self.simulator.maxord, 12)
   
    @testattr(stddist = True)
    def test_pretype(self):
        """
        This tests the precondition option.
        """
        nose.tools.assert_equal(self.simulator.precond, 'PREC_NONE')
        self.simulator.precond = 'prec_none'
        nose.tools.assert_equal(self.simulator.precond, 'PREC_NONE')
        
        nose.tools.assert_raises(Exception, self.simulator._set_pre_cond, -1.0)
        nose.tools.assert_raises(Exception, self.simulator._set_pre_cond, 'PREC_BOTH1')
    
    @testattr(stddist = True)
    def test_maxkrylov(self):
        """
        This test the maximum number of krylov subspaces.
        """
        nose.tools.assert_equal(self.simulator.maxkrylov, 5)
        self.simulator.maxkrylov = 3
        nose.tools.assert_equal(self.simulator.maxkrylov, 3)
        self.simulator.maxkrylov = 4.5
        nose.tools.assert_equal(self.simulator.maxkrylov, 4)
        
        nose.tools.assert_raises(Exception, self.simulator._set_max_krylov, 'Test')
        
    @testattr(stddist = True)
    def test_stablimit(self):
        nose.tools.assert_false(self.simulator.stablimit)
        self.simulator.stablimit = True
        nose.tools.assert_true(self.simulator.stablimit)
        nose.tools.assert_true(self.simulator.options["stablimit"])
    
    @testattr(stddist = True)
    def test_linearsolver(self):
        """
        This test the choice of the linear solver.
        """
        nose.tools.assert_equal(self.simulator.linear_solver, 'DENSE')
        self.simulator.linear_solver = 'dense'
        nose.tools.assert_equal(self.simulator.linear_solver, 'DENSE')
        self.simulator.linear_solver = 'spgmr'
        nose.tools.assert_equal(self.simulator.linear_solver, 'SPGMR')
        
        nose.tools.assert_raises(Exception, self.simulator._set_linear_solver, 'Test')
    
    @testattr(stddist = True)
    def test_terminate_simulation(self):
        """
        This tests the functionality of raising TerminateSimulation exception in handle_result.
        """
        class Extended_Problem(Explicit_Problem):
            def __init__(self):
                pass
            def handle_event(self, solver, event_info):
                if solver.t > 1.5:
                    raise TerminateSimulation
            rhs = lambda self,t,y,sw: N.array([1.0])
            y0 = [1.0]
            sw0 = [False,True]
            state_events = lambda self,t,y,sw:  N.array([t-1.0, t-2.0])

        exp_mod = Extended_Problem()
        simulator = CVode(exp_mod)
        simulator(3.)
        
        nose.tools.assert_almost_equal(simulator.t, 2.000000, 4)
    
    @testattr(stddist = True)
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
        nose.tools.assert_equal(len(sim.t_sol), 101)
        nose.tools.assert_equal(nsteps, sim.statistics["nsteps"])
        
        sim = CVode(mod)
        nsteps = 0
        sim.simulate(2.)
        nose.tools.assert_equal(len(sim.t_sol), sim.statistics["nsteps"]+1)
        nose.tools.assert_equal(nsteps, sim.statistics["nsteps"])
        
class Test_IDA:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y,yd: y
        y0 = [1.0]
        yd0 = [1.0]
        
        self.problem = Implicit_Problem(f,y0,yd0)
        self.simulator = IDA(self.problem)
    
    @testattr(stddist = True)
    def test_time_limit(self):
        f = lambda t,y,yd: yd-y
        
        exp_mod = Implicit_Problem(f,1.0,1.0)
        exp_sim = IDA(exp_mod)
        
        exp_sim.maxh = 1e-8
        exp_sim.time_limit = 1 #One second
        exp_sim.report_continuously = True
        
        nose.tools.assert_raises(TimeLimitExceeded, exp_sim.simulate, 1)
    
    @testattr(stddist = True)
    def test_simulate_explicit(self):
        """
        Test a simulation of an explicit problem using IDA.
        """
        f = lambda t,y:N.array(-y)
        y0 = [1.0]
        
        problem = Explicit_Problem(f,y0)
        simulator = IDA(problem)
        
        nose.tools.assert_equal(simulator.yd0[0], -simulator.y0[0])
        
        t,y = simulator.simulate(1.0)
        
        nose.tools.assert_almost_equal(float(y[-1]), float(N.exp(-1.0)),4)
    
    @testattr(stddist = True)    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        nose.tools.assert_false(self.simulator.suppress_alg)
        nose.tools.assert_equal(self.simulator.algvar[0], 1.0)
        nose.tools.assert_equal(self.simulator.sw, None)
        nose.tools.assert_equal(self.simulator.maxsteps, 10000)
        nose.tools.assert_equal(self.simulator.y[0], 1.0)
    
    @testattr(stddist = True)
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
        nose.tools.assert_almost_equal(y100[-2], sim.interpolate(9.9,0),5)
    
    @testattr(stddist = True)
    def test_handle_result(self):
        """
        This function tests the handle result.
        """
        f = lambda t,x,xd: x**0.25-xd
        def handle_result(solver, t ,y, yd):
            nose.tools.assert_equal(solver.t, t)
        
        prob = Implicit_Problem(f, [1.0],[1.0])
        prob.handle_result = handle_result
        
        sim = IDA(prob)

        sim.report_continuously = True
        
        sim.simulate(10.)
    
    @testattr(stddist = True)    
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, "Test")
        nose.tools.assert_raises(Exception, self.simulator._set_max_ord, [1,1])
        
        
        self.simulator.maxord = -1
        nose.tools.assert_equal(self.simulator.maxord, 1)
        self.simulator.maxord = 2
        nose.tools.assert_equal(self.simulator.maxord, 2)
        self.simulator.maxord = 6
        nose.tools.assert_equal(self.simulator.maxord, 5)
    
    @testattr(stddist = True)    
    def test_tout1(self):
        """
        This tests the functionality of the property tout1.
        """
        nose.tools.assert_raises(Exception, self.simulator._set_tout1, 'Test')
        nose.tools.assert_raises(Exception, self.simulator._set_tout1, [1,1])
        nose.tools.assert_raises(Exception, self.simulator._set_tout1, 'Test')
        
        nose.tools.assert_equal(self.simulator.tout1, 0.0001)
        self.simulator.tout1 = -0.001
        nose.tools.assert_equal(self.simulator.tout1, -0.001)
        self.simulator.tout1 = 1
        nose.tools.assert_equal(self.simulator.tout1, 1.0)
        
    @testattr(stddist = True)
    def test_lsoff(self):
        """
        This tests the functionality of the property lsoff.
        """
        nose.tools.assert_false(self.simulator.lsoff)
        self.simulator.lsoff = True
        nose.tools.assert_true(self.simulator.lsoff)
        self.simulator.lsoff = False
        nose.tools.assert_false(self.simulator.lsoff)
    
    @testattr(stddist = True)
    def test_initstep(self):
        """
        This tests the funtionality of the property initstep.
        """
        
        def f(t,y,yd):
            res_0 = yd[0] - y[1]
            res_1 = yd[1] +9.82-0.01*y[1]**2
            return N.array([res_0,res_1])
            
        mod = Implicit_Problem(f,y0=[5.0,0.0], yd0=[0.0,9.82])
        
        
        sim = IDA(mod)
        sim.simulate(2.0)

        nose.tools.assert_almost_equal(sim.y_sol[-1][0], -13.4746473811, places=7)
        
        sim.reset()
        sim.inith = 1e-10
        sim.simulate(2.0)

        nose.tools.assert_almost_equal(sim.y_sol[-1][0], -13.4746596311, places=7)
        
    @testattr(stddist = True)
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
                solver.y  = N.array([1.0])
                solver.yd = N.array([1.0])
                
                if not solver.sw[0]:
                    solver.sw[1] = False
                
                if solver.sw[0]:
                    solver.sw[0] = False
        
        mod = Implicit_Problem(f,[1.0],[1.0])
        mod.time_events = time_events
        mod.handle_event = handle_event
        mod.switches0 = [True, True]
        
        sim = IDA(mod)
        
        sim.simulate(5.0)

        nose.tools.assert_almost_equal(sim.y_sol[38], 1.0000000, 5)
        nose.tools.assert_almost_equal(sim.y_sol[87], 1.0000000, 5)
        
        sim = IDA(mod, [1.0],[1.0])
        sim.simulate(2.0)
        
        nose.tools.assert_almost_equal(sim.t_sol[-1], 2.0000000, 5)
    
    @testattr(stddist = True)
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
                solver.y  = N.array([1.0])
                solver.yd = N.array([1.0])
                
                if not solver.sw[0]:
                    solver.sw[1] = False
                
                if solver.sw[0]:
                    solver.sw[0] = False
        
        mod = Implicit_Problem(f,[1.0],[1.0], sw0=[True, True])
        mod.time_events = time_events
        mod.handle_event = handle_event
        
        sim = IDA(mod)
        sim.verbosity = 10
        nose.tools.assert_equal(len(sim.event_data), 0)
        sim.simulate(5.0)
        nose.tools.assert_greater(len(sim.event_data), 0)
        
        sim.reset()
        nose.tools.assert_equal(len(sim.event_data), 0)
        sim.simulate(5.0)
        nose.tools.assert_greater(len(sim.event_data), 0)
        
    @testattr(stddist = True)    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x,xd: N.array([xd[0]-x[1], xd[1]-9.82])       #Defines the rhs
        jac = lambda c,t,x,xd: N.array([[c,-1.],[0.,c]]) #Defines the jacobian

        imp_mod = Implicit_Problem(f,[1.0,0.0],[0.,-9.82])
        imp_mod.jac = jac
        
        imp_sim = IDA(imp_mod)
        
        imp_sim.simulate(3,100)

        nose.tools.assert_equal(imp_sim.statistics["nfcnjacs"], 0)
        nose.tools.assert_almost_equal(imp_sim.y_sol[-1][0], 45.1900000, 4)
        
        imp_sim.reset()
        imp_sim.usejac=False
        imp_sim.simulate(3.,100)

        nose.tools.assert_almost_equal(imp_sim.y_sol[-1][0], 45.1900000, 4)
        nose.tools.assert_greater(imp_sim.statistics["nfcnjacs"], 0)
    
    @testattr(stddist = True)
    def test_terminate_simulation(self):
        """
        This tests the functionality of raising TerminateSimulation exception in handle_result.
        """
        class Extended_Problem(Implicit_Problem):
            def __init__(self):
                pass
            def handle_event(self,solver, event_info):
                if solver.t > 1.5:
                    raise TerminateSimulation
            res = lambda self,t,y,yd,sw: N.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: N.array([t-1.0, t-2.0])
            y0 = [1.0]
            yd0 = [1.0]
            sw0 = [False]

        prob = Extended_Problem()
    
        sim = IDA(prob)
        sim.simulate(2.5)
        
        nose.tools.assert_almost_equal(sim.t, 2.000000, 4)
    
    @testattr(stddist = True)    
    def test_algvar(self):
        """
        This tests the functionality of the property algvar.
        """
        #self.simulator.Integrator.dim = 3
        
        #nose.tools.assert_raises(Exception, self.simulator._set_algvar, 1)
        #nose.tools.assert_raises(Exception, self.simulator._set_algvar, 1.0)
        nose.tools.assert_raises(Exception, self.simulator._set_algvar, [1,'hej',1])
        nose.tools.assert_raises(Exception, self.simulator._set_algvar, {'Test':'case'})
        nose.tools.assert_raises(Exception, self.simulator._set_algvar, [-1,0,1])
        nose.tools.assert_raises(Exception, self.simulator._set_algvar, [1.0,1.0])
        nose.tools.assert_raises(Exception, self.simulator._set_algvar, [3.0, 1.0, 1.0])
        
        #vector = [1.0,0.0,1.0]
        #vectorb = [True,False,True]
        #vectori = [1,0,1]
        
        #self.simulator.algvar = vectorb
        #self.simulator.algvar = vectori
        #self.simulator.algvar = vector
        #nose.tools.assert_equal(self.simulator.algvar[0], vector[0])
        #nose.tools.assert_equal(self.simulator.algvar[1], vector[1])
        #nose.tools.assert_equal(self.simulator.algvar[2], vector[2])
    
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
            solver.y+= 1.0
            global tnext
            nose.tools.assert_almost_equal(solver.t, tnext)
            nose.tools.assert_equal(event_info[0], [])
            nose.tools.assert_true(event_info[1])
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = IDA(exp_mod)
        exp_sim(5.,100)
        
        nose.tools.assert_equal(nevent, 5)
    
    @testattr(stddist = True)    
    def test_suppress_alg(self):
        """
        This tests the functionality of the property suppress_alg.
        """
        self.simulator.suppress_alg = True
        nose.tools.assert_true(self.simulator.suppress_alg)
        self.simulator.suppress_alg = False
        nose.tools.assert_false(self.simulator.suppress_alg)
        
    @testattr(stddist = True)
    def test_make_consistency(self):
        """
        This tests the functionality of the method make_consistency.
        """
        def f(t,y,yd):
            res_1 = y[0] + y[1]+1.0
            res_2 = y[1]
            return N.array([res_1, res_2])
        y0 = [2.0, 2.0]
        yd0 = [1.0 , 0.0]
        
        
        my_Prob = Implicit_Problem(f, y0, yd0)
        
        simulator = IDA(my_Prob)
        
        [flag, y, yd] = simulator.make_consistent('IDA_Y_INIT')
        
        nose.tools.assert_almost_equal(y[1], 0.00000)
        nose.tools.assert_almost_equal(y[0], -1.0000)
        nose.tools.assert_almost_equal(yd[0], 1.0000)
        nose.tools.assert_almost_equal(yd[1], 0.0000)

    @testattr(stddist = True)
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,xd,sw: N.array([xd[0]- 1.0])
        state_events = lambda t,x,xd,sw: N.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.sw = [False] #Override the switches to point to another instance
        
        mod = Implicit_Problem(f, [0.0],[1.0])
        mod.f = f
        mod.sw0 = [True]
        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = IDA(mod)
        nose.tools.assert_true(sim.sw[0])
        sim.simulate(3)
        nose.tools.assert_false(sim.sw[0])
    
    @testattr(stddist = True)
    def test_completed_step(self):
        """
        This tests the functionality of the method completed_step.
        """
        global nsteps
        nsteps = 0
        def f(t,y,yd):
            res_1 = y[0] + y[1]+1.0
            res_2 = y[1]
            return N.array([res_1, res_2])
        def completed_step(solver):
            global nsteps
            nsteps += 1
        
        y0 = [-1.0, 0.0]
        yd0 = [1.0 , 0.0]    
        
        mod = Implicit_Problem(f, y0, yd0)
        mod.step_events = completed_step
        
        sim = IDA(mod)
        
        sim.simulate(2., 100)
        nose.tools.assert_equal(len(sim.t_sol), 101)
        nose.tools.assert_equal(nsteps, sim.statistics["nsteps"])
        
        sim = IDA(mod)
        nsteps = 0
        sim.simulate(2.)
        nose.tools.assert_equal(len(sim.t_sol), sim.statistics["nsteps"] + 1)
        nose.tools.assert_equal(nsteps, sim.statistics["nsteps"])



class Test_Sundials:
    
    def setUp(self):
        """
        This sets up the test case.
        """
        class Prob_IDA(Implicit_Problem):
            def __init__(self):
                pass
            res = lambda self,t,y,yd,sw: N.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: N.array([t-1.0, t])
            y0 = [1.0]
            yd0 = [1.0]
            sw0 = [False, True]
            
        res = Prob_IDA()
        
        class Prob_CVode(Explicit_Problem):
            def __init__(self):
                pass
            rhs = lambda self,t,y,sw: N.array([1.0])
            state_events = lambda self,t,y,sw: N.array([t-1.0, t])
            y0 = [1.0]
            sw0 = [False, True]

        f = Prob_CVode()
        
        self.simulators = [IDA(res), CVode(f)]
        
        
        f = lambda t,y,yd,p: N.array([0.0])
        y0 = [1.0]
        yd0 = [1.0]
        p0 = [1.0]
        
        mod = Implicit_Problem(f, y0,yd0,p0=p0)
        self.sim = IDA(mod)
    
    @testattr(stddist = True)
    def test_atol(self):
        """
        This tests the functionality of the property atol.
        """
        nose.tools.assert_equal(self.simulators[1].atol, 1.0e-6)
        nose.tools.assert_equal(self.simulators[0].atol, 1.0e-6)
        
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Exception, self.simulators[i]._set_atol, -1.0)
            nose.tools.assert_raises(Exception, self.simulators[i]._set_atol, [1.0, 1.0])
            nose.tools.assert_raises(Exception, self.simulators[i]._set_atol, "Test")
            
            self.simulators[i].atol = 1.0e-5
            nose.tools.assert_equal(self.simulators[i].atol, 1.0e-5)
            self.simulators[i].atol = 1.0
            nose.tools.assert_equal(self.simulators[i].atol, 1.0)
            self.simulators[i].atol = 1
            nose.tools.assert_equal(self.simulators[i].atol, 1.0)
            self.simulators[i].atol = 1001.0
            nose.tools.assert_equal(self.simulators[i].atol, 1001.0)
            self.simulators[i].atol = [N.array([1e-5])]
            nose.tools.assert_equal(len(self.simulators[i].atol.shape), 1)
            nose.tools.assert_equal(self.simulators[i].atol, 1e-5)
            """
            self.simulators[i].Integrator.dim = 3
            nose.tools.assert_raises(Exception, self.simulators[i]._set_atol, [1.0, 1.0])
            nose.tools.assert_raises(Exception, self.simulators[i]._set_atol, [1.0, 1.0, -1.0])
            self.simulators[i].atol = [1.0, 1.0, 1.0]
            nose.tools.assert_equal(self.simulators[i].atol, [1.0, 1.0, 1.0])
            self.simulators[i].atol = N.array([1.0, 1.0, 1.0])
            nose.tools.assert_equal(self.simulators[i].atol[0], 1.0)
            self.simulators[i].atol = N.array([1, 5, 1.0])
            nose.tools.assert_equal(self.simulators[i].atol[0], 1.0)
            """
    
    
    @testattr(stddist = True)
    def test_rtol(self):
        """
        This tests the functionality of the property rtol.
        """
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Exception, self.simulators[i]._set_rtol, -1.0)
            nose.tools.assert_raises(Exception, self.simulators[i]._set_rtol, [1.0, 1.0])
            nose.tools.assert_raises(Exception, self.simulators[i]._set_rtol, "Test")
            
            self.simulators[i].rtol = 1.0e-5
            nose.tools.assert_equal(self.simulators[i].rtol, 1.0e-5)
            self.simulators[i].rtol = 1.0
            nose.tools.assert_equal(self.simulators[i].rtol, 1.0)
            self.simulators[i].rtol = 1001.0
            nose.tools.assert_equal(self.simulators[i].rtol, 1001.0)
            self.simulators[i].rtol = 1001
            nose.tools.assert_equal(self.simulators[i].rtol, 1001.0)
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the functionality of the property maxh.
        """
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Exception, self.simulators[i]._set_max_h, [1.0, 1.0])
            nose.tools.assert_raises(Exception, self.simulators[i]._set_max_h, "Test")
            
            self.simulators[i].maxh = 1.0e-5
            nose.tools.assert_equal(self.simulators[i].maxh, 1.0e-5)
            self.simulators[i].maxh = 1.0
            nose.tools.assert_equal(self.simulators[i].maxh, 1.0)
            self.simulators[i].maxh = 1001.0
            nose.tools.assert_equal(self.simulators[i].maxh, 1001.0)

    @testattr(stddist = True)    
    def test_dqtype(self):
        """
        Tests the property of dqtype.
        """
        
        nose.tools.assert_equal(self.sim.dqtype, 'CENTERED') #Test the default value.
        
        self.sim.dqtype = 'FORWARD'
        nose.tools.assert_equal(self.sim.dqtype, 'FORWARD')
        self.sim.dqtype = 'CENTERED'
        nose.tools.assert_equal(self.sim.dqtype, 'CENTERED')
        
        self.sim.dqtype = 'forward'
        nose.tools.assert_equal(self.sim.dqtype, 'FORWARD')
        self.sim.dqtype = 'centered'
        nose.tools.assert_equal(self.sim.dqtype, 'CENTERED')
        
        nose.tools.assert_raises(Exception,self.sim._set_dqtype, 1)
        nose.tools.assert_raises(Exception,self.sim._set_dqtype, 'IDA_CE')
        nose.tools.assert_raises(Exception,self.sim._set_dqtype, [1])
        nose.tools.assert_raises(Exception,self.sim._set_dqtype, -1)
    
    @testattr(stddist = True)
    def test_dqrhomax(self):
        """
        Tests the property of DQrhomax.
        """
        nose.tools.assert_equal(self.sim.dqrhomax, 0.0) #Test the default value.
        
        self.sim.dqrhomax = 1.0
        nose.tools.assert_equal(self.sim.dqrhomax, 1.0)
        self.sim.dqrhomax = 10
        nose.tools.assert_equal(self.sim.dqrhomax, 10)
        
        nose.tools.assert_raises(Exception,self.sim._set_dqrhomax, -1)
        nose.tools.assert_raises(Exception,self.sim._set_dqrhomax, 'str')
        nose.tools.assert_raises(Exception,self.sim._set_dqrhomax, [])
        nose.tools.assert_raises(Exception,self.sim._set_dqrhomax, -10)
    
    @testattr(stddist = True)    
    def test_usesens(self):
        """
        Tests the property of usesens.
        """
        nose.tools.assert_true(self.sim.usesens)#Test the default value.
        self.sim.usesens = False
        nose.tools.assert_false(self.sim.usesens)
        self.sim.usesens = 0
        nose.tools.assert_false(self.sim.usesens)
        self.sim.usesens = 1
        nose.tools.assert_true(self.sim.usesens)

    @testattr(stddist = True)
    def test_sensmethod(self):
        """
        Tests the property of sensmethod.
        """
        nose.tools.assert_equal(self.sim.sensmethod, 'STAGGERED') #Test the default value
        
        self.sim.sensmethod = 'SIMULTANEOUS'
        nose.tools.assert_equal(self.sim.sensmethod, 'SIMULTANEOUS')
        self.sim.sensmethod = 'STAGGERED'
        nose.tools.assert_equal(self.sim.sensmethod, 'STAGGERED')
        
        self.sim.sensmethod = 'simultaneous'
        nose.tools.assert_equal(self.sim.sensmethod, 'SIMULTANEOUS')
        self.sim.sensmethod = 'staggered'
        nose.tools.assert_equal(self.sim.sensmethod, 'STAGGERED')
        
        nose.tools.assert_raises(Exception,self.sim._set_sensitivity_method, 1)
        nose.tools.assert_raises(Exception,self.sim._set_sensitivity_method, 'IDA_CE')
        nose.tools.assert_raises(Exception,self.sim._set_sensitivity_method, [1])
        nose.tools.assert_raises(Exception,self.sim._set_sensitivity_method, -1)
    
    @testattr(stddist = True)    
    def test_suppress_sens(self):
        """
        Tests the property of suppress_sens.
        """
        nose.tools.assert_false(self.sim.suppress_sens)
        self.sim.suppress_sens = False
        nose.tools.assert_false(self.sim.suppress_sens)
        self.sim.suppress_sens = 0
        nose.tools.assert_false(self.sim.suppress_sens)
        self.sim.suppress_sens = 1
        nose.tools.assert_true(self.sim.suppress_sens)
    
    @testattr(stddist = True)
    def test_maxsensiter(self):
        """
        Tests the property of maxsensiter.
        """
        nose.tools.assert_equal(self.sim.maxcorS, 3) #Test the default value
        self.sim.maxcorS = 1
        nose.tools.assert_equal(self.sim.maxcorS, 1)
        self.sim.maxcorS = 10.5
        nose.tools.assert_equal(self.sim.maxcorS, 10)
        
        #nose.tools.assert_raises(Exception, self.sim._set_max_cor_S, 0)
        nose.tools.assert_raises(Exception, self.sim._set_max_cor_S, 'str')
        nose.tools.assert_raises(Exception, self.sim._set_max_cor_S, [])
        #nose.tools.assert_raises(Exception, self.sim._set_max_cor_S, -10)
    
    @testattr(stddist = True)
    def test_pbar(self):
        """
        Tests the property of pbar.
        """
        f = lambda t,y,p:N.array([0.0]*len(y))
        y0 = [1.0]*2
        p0 = [1000.0, -100.0]
        exp_mod = Explicit_Problem(f,y0,p0=p0)
        
        exp_sim = CVode(exp_mod)

        nose.tools.assert_almost_equal(exp_sim.pbar[0], 1000.00000,4)
        nose.tools.assert_almost_equal(exp_sim.pbar[1], 100.000000,4)
        
        f = lambda t,y,yd,p: N.array([0.0]*len(y))
        yd0 = [0.0]*2
        imp_mod = Implicit_Problem(f,y0,yd0,p0=p0)
        
        imp_sim = IDA(imp_mod)
        
        nose.tools.assert_almost_equal(imp_sim.pbar[0], 1000.00000,4)
        nose.tools.assert_almost_equal(imp_sim.pbar[1], 100.000000,4)
