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
from assimulo.lib.odepack import dsrcar, dcfode
from assimulo.solvers.odepack import RKStarterNordsieck
from assimulo.solvers import LSODAR
from assimulo.problem import Explicit_Problem
from assimulo.exception import TimeLimitExceeded

import numpy as np
import scipy.sparse as sps

float_regex = r"[\s]*[\d]*.[\d]*((e|E)(\+|\-)\d\d|)"

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


class Test_LSODAR:
    """
    Tests the LSODAR solver.
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
        cls.sim = LSODAR(exp_mod) #Create a LSODAR solve
        cls.sim_sp = LSODAR(exp_mod_sp)
        
        #Sets the parameters
        cls.sim.atol = 1e-6 #Default 1e-6
        cls.sim.rtol = 1e-6 #Default 1e-6
        cls.sim.usejac = False

    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = LSODAR(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)

    def test_simulation(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.simulate(1.) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)
    
    def test_setcoefficients(self):
        elco,tesco=dcfode(1)
        assert elco[0,2] == pytest.approx(5./12., abs = 1e-9) # AM-2
        assert tesco[0,2] == pytest.approx(2., abs = 1e-9) # AM-2 error coeff  
    
    def test_readcommon(self):
        """
        This tests the LSODAR's subroutine dsrcar  (read functionality).
        """
        self.sim.simulate(1.) #Simulate 2 seconds
        r=np.ones((245,),'d')
        i=np.ones((55,),'i')
        dsrcar(r,i,1)
        assert r[217] == pytest.approx(2.22044605e-16, abs = 1e-20)
        assert i[36] == 3
        
    def test_writereadcommon(self):
        """
        This tests the LSODAR's subroutine dsrcar  (write and read functionality).
        """
        r=np.ones((245,),'d')
        i=np.ones((55,),'i')
        dsrcar(r,i,2)
        r[0]=100.
        i[0]=10
        dsrcar(r,i,1)
        assert r[0] == pytest.approx(1., abs = 1e-4)
        assert i[0] == 1  
    
    @pytest.mark.skip("Test broken")
    def test_rkstarter(self):
        """
        This test checks the correctness of the Nordsieck array generated 
        from a RK starter
        """
        A=np.array([[0.,1.],[-4.,0.]])
        def f(t,x,sw0):
            return np.dot(A,np.array(x))
        H = 1.e-8
        # Compute the exact solution at h=0,H/4,H/2,3H/4,H
        T=np.array([0.,H/4.,H/2.,3./4.*H,H])
        y0=np.array([1.,0.])
        from scipy.linalg import expm
        exact=np.array([np.dot(expm(A*t),y0) for t in T])
        #polynomial interpolation
        from scipy import polyfit
        coeff = polyfit(T,exact,4)
        d1coeff=np.array([4,3,2,1]).reshape(-1,1)*coeff[:-1,:]
        d2coeff=np.array([3,2,1]).reshape(-1,1)*d1coeff[:-1,:]
        d3coeff=np.array([2,1]).reshape(-1,1)*d2coeff[:-1,:]
        d4coeff=np.array([1]).reshape(-1,1)*d3coeff[:-1,:]
        h=H/4.
        nordsieck_at_0=np.array([coeff[-1,:],h*d1coeff[-1,:],h**2/2.*d2coeff[-1,:],
                                     h**3/6.*d3coeff[-1,:],h**4/24.*d4coeff[-1,:]])
        rkNordsieck=RKStarterNordsieck(f,H)
        computed=rkNordsieck(0,y0)       
        np.testing.assert_allclose(computed[1], nordsieck_at_0, atol=H/100., verbose=True)      
        
    def test_interpol(self):
        # a with interpolation and report_continuously
        self.sim.report_continuously=True
        t_sol,y_sol=self.sim.simulate(1.,ncp_list=[0.5])
        self.sim.reset()
        t_sol1,y_sol1=self.sim.simulate(0.5)
        ind05=np.nonzero(np.array(t_sol)==0.5)[0][0]
        assert y_sol[ind05,0] == pytest.approx(y_sol1[-1,0], abs = 1e-6)
        
    def test_simulation_with_jac(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.usejac = True
        self.sim.simulate(1.) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)
    
    def test_simulation_ncp(self):
        self.sim.simulate(1.,100) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)
        
    def test_usejac_csc_matrix(self):
        self.sim_sp.usejac = True
        
        self.sim_sp.simulate(2.) #Simulate 2 seconds
    
        assert self.sim_sp.statistics["nfcnjacs"] == 0
        assert self.sim_sp.y_sol[-1][0] == pytest.approx(1.7061680350, abs = 1e-4)
        
    def test_simulation_ncp_list(self):
        self.sim.simulate(1.,ncp_list=[0.5]) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)
        
    def test_maxh(self):
        
        self.sim.hmax = 1.0
        assert self.sim.options["maxh"] == 1.0
        assert self.sim.maxh == 1.0
        
        self.sim.maxh = 2.0
        assert self.sim.options["maxh"] == 2.0
        assert self.sim.maxh == 2.0
        
    def test_simulation_ncp_list_2(self):
        self.sim.simulate(1.,ncp_list=[0.5,4]) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)
        
    def test_simulation_ncp_with_jac(self):
        """
        Test a simulation with ncp.
        """
        self.sim.usejac= True
        self.sim.simulate(1.,100) #Simulate 2 seconds

        assert self.sim.y_sol[-1][0] == pytest.approx(-1.863646028, abs = 1e-4)

    def test_time_limit(self):
        """ Test that simulation is canceled when a set time limited is exceeded. """
        import time
        def f(t, y):
            time.sleep(.1)
            return -y
        
        prob = Explicit_Problem(f,1.0)
        sim = LSODAR(prob)
        
        sim.maxh = 1e-5
        sim.time_limit = 1
        sim.report_continuously = True

        err_msg = f'The time limit was exceeded at integration time {float_regex}.'
        with pytest.raises(TimeLimitExceeded, match = err_msg):
            sim.simulate(1.)
