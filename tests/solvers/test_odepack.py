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
import numpy.testing
from assimulo import testattr
from assimulo.lib.odepack import dsrcar, dcfode
from assimulo.solvers import LSODAR, odepack
from assimulo.problem import Explicit_Problem
from assimulo.exception import *

import numpy as N
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


class Test_LSODAR:
    """
    Tests the LSODAR solver.
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
        self.sim = LSODAR(exp_mod) #Create a LSODAR solve
        self.sim_sp = LSODAR(exp_mod_sp)
        
        #Sets the parameters
        self.sim.atol = 1e-6 #Default 1e-6
        self.sim.rtol = 1e-6 #Default 1e-6
        self.sim.usejac = False

    @testattr(stddist = True)
    def test_event_localizer(self):
        exp_mod = Extended_Problem() #Create the problem

        exp_sim = LSODAR(exp_mod) #Create the solver
        
        exp_sim.verbosity = 0
        exp_sim.report_continuously = True
        
        #Simulate
        t, y = exp_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
        
        #Basic test
        nose.tools.assert_almost_equal(y[-1][0],8.0)
        nose.tools.assert_almost_equal(y[-1][1],3.0)
        nose.tools.assert_almost_equal(y[-1][2],2.0)

    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.simulate(1.) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
    @testattr(stddist = True)
    def test_setcoefficients(self):
        elco,tesco=dcfode(1)
        nose.tools.assert_almost_equal(elco[0,2],5./12.,9) # AM-2
        nose.tools.assert_almost_equal(tesco[0,2],2.,9) # AM-2 error coeff  
    @testattr(stddist = True)
    def test_readcommon(self):
        """
        This tests the LSODAR's subroutine dsrcar  (read functionality).
        """
        self.sim.simulate(1.) #Simulate 2 seconds
        r=N.ones((245,),'d')
        i=N.ones((55,),'i')
        dsrcar(r,i,1)
        nose.tools.assert_almost_equal(r[217], 2.22044605e-16, 20)
        nose.tools.assert_equal(i[36], 3)
        
    @testattr(stddist = True)
    def test_writereadcommon(self):
        """
        This tests the LSODAR's subroutine dsrcar  (write and read functionality).
        """
        r=N.ones((245,),'d')
        i=N.ones((55,),'i')
        dsrcar(r,i,2)
        r[0]=100.
        i[0]=10
        dsrcar(r,i,1)
        nose.tools.assert_almost_equal(r[0], 1., 4)
        nose.tools.assert_equal(i[0], 1)  
    def test_rkstarter(self):
        """
        This test checks the correctness of the Nordsieck array generated 
        from a RK starter
        """
        pass
        """
        A=N.array([[0.,1.],[-4.,0.]])
        def f(t,x,sw0):
            return N.dot(A,N.array(x))
        H = 1.e-8
        # Compute the exact solution at h=0,H/4,H/2,3H/4,H
        T=N.array([0.,H/4.,H/2.,3./4.*H,H])
        y0=N.array([1.,0.])
        from scipy.linalg import expm
        exact=N.array([N.dot(expm(A*t),y0) for t in T])
        #polynomial interpolation
        from scipy import polyfit
        coeff = polyfit(T,exact,4)
        d1coeff=N.array([4,3,2,1]).reshape(-1,1)*coeff[:-1,:]
        d2coeff=N.array([3,2,1]).reshape(-1,1)*d1coeff[:-1,:]
        d3coeff=N.array([2,1]).reshape(-1,1)*d2coeff[:-1,:]
        d4coeff=N.array([1]).reshape(-1,1)*d3coeff[:-1,:]
        h=H/4.
        nordsieck_at_0=N.array([coeff[-1,:],h*d1coeff[-1,:],h**2/2.*d2coeff[-1,:],
                                     h**3/6.*d3coeff[-1,:],h**4/24.*d4coeff[-1,:]])
        rkNordsieck=odepack.RKStarterNordsieck(f,H)
        computed=rkNordsieck(0,y0)       
        numpy.testing.assert_allclose(computed[1], nordsieck_at_0, atol=H/100., verbose=True)      
        """              
        
    @testattr(stddist = True)
    def test_interpol(self):
        # a with interpolation and report_continuously
        self.sim.report_continuously=True
        t_sol,y_sol=self.sim.simulate(1.,ncp_list=[0.5])
        self.sim.reset()
        t_sol1,y_sol1=self.sim.simulate(0.5)
        ind05=N.nonzero(N.array(t_sol)==0.5)[0][0]
        #print y_sol[ind05],y_sol1[-1]
        nose.tools.assert_almost_equal(y_sol[ind05,0],y_sol1[-1,0],6)
        
        
    def test_simulation_with_jac(self):
        """
        This tests the LSODAR with a simulation of the van der pol problem.
        """
        self.sim.usejac = True
        self.sim.simulate(1.) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        self.sim.simulate(1.,100) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
        
    @testattr(stddist = True)
    def test_usejac_csc_matrix(self):
        self.sim_sp.usejac = True
        
        self.sim_sp.simulate(2.) #Simulate 2 seconds
    
        nose.tools.assert_equal(self.sim_sp.statistics["nfcnjacs"], 0)
        
        nose.tools.assert_almost_equal(self.sim_sp.y_sol[-1][0], 1.7061680350, 4)
        
    @testattr(stddist = True)    
    def test_simulation_ncp_list(self):
        self.sim.simulate(1.,ncp_list=[0.5]) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
        
    @testattr(stddist = True)    
    def test_maxh(self):
        
        self.sim.hmax = 1.0
        nose.tools.assert_equal(self.sim.options["maxh"], 1.0)
        nose.tools.assert_equal(self.sim.maxh, 1.0)
        
        self.sim.maxh = 2.0
        nose.tools.assert_equal(self.sim.options["maxh"], 2.0)
        nose.tools.assert_equal(self.sim.maxh, 2.0)
        
    @testattr(stddist = True)    
    def test_simulation_ncp_list_2(self):
        self.sim.simulate(1.,ncp_list=[0.5,4]) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
        
    @testattr(stddist = True)    
    def test_simulation_ncp_with_jac(self):
        """
        Test a simulation with ncp.
        """
        self.sim.usejac= True
        self.sim.simulate(1.,100) #Simulate 2 seconds

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], -1.863646028, 4)
