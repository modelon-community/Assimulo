import nose
from assimulo.explicit_ode import *
from assimulo.problem import Explicit_Problem

class Test_Explicit_ODE:
    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        def f(t, y):
            pass
        Test = Explicit_Problem()
        Test.f = f
        
        nose.tools.assert_raises(Explicit_ODE_Exception, Explicit_ODE, Test, 'test')
        nose.tools.assert_raises(Explicit_ODE_Exception, Explicit_ODE, Test, 1, 'test')
        nose.tools.assert_raises(Explicit_ODE_Exception, Explicit_ODE, Test, 'test', 'test')
        nose.tools.assert_raises(Explicit_ODE_Exception, Explicit_ODE, Test, [1.0 , 1.0, 'test'])
        nose.tools.assert_raises(Explicit_ODE_Exception, Explicit_ODE, None, [1.0 , 1.0, 1.0])
        
        
        simulator = Explicit_ODE(Test, [1.0 , 1.0], 1)
        #assert simulator.t[0] == 1.0
        #assert simulator.y[0][0] == 1.0
        assert simulator.t_cur == 1.0
        assert simulator.y_cur[0] == 1.0
        
        
    def test_call(self):
        """
        This tests the functionality of the method __call__.
        """
        y0 = [0.0]
        
        def f(t, x):
            return x

        problem = Explicit_Problem()
        problem.f = f
        
        simulator = CVode(problem,y0)
        nose.tools.assert_raises(Explicit_ODE_Exception, simulator, -1.0)
        nose.tools.assert_raises(Explicit_ODE_Exception, simulator, 'test')
        simulator.reset()
        #[t,y] = simulator(1.0,10)
        simulator(1.0,10)
        t = simulator.t
        y = simulator.y
        
        assert len(t) == 11 #11 Due to t0 is counted as well
        #simulator = Explicit_Euler(problem,y0)
        #[t,y] = simulator(1.0,10)
        #print len(t)
        #assert len(t) == 11 #11 Due to t0 is counted as well
        #simulator = RungeKutta34(problem,y0)
        #[t,y] = simulator(1.0,10)
        #assert len(t) == 11 #11 Due to t0 is counted as well
        #simulator = RungeKutta4(problem,y0)
        #[t,y] = simulator(1.0,10)
        #assert len(t) == 11 #11 Due to t0 is counted as well
    
    def test_reset(self):
        """
        This tests the functionality of the method reset.
        """
        y0 = [0.0]
        
        f = lambda t, x: x

        problem = Explicit_Problem()
        problem.f = f
        
        simulator = Explicit_ODE(problem, y0)
        
        simulator.atol = 0.1
        simulator.rtol = 0.01
        simulator.maxsteps = 1000
        simulator.post_process = True
        simulator.verbosity = 2
        
        simulator.reset()
        
        assert simulator.atol == 0.1
        assert simulator.rtol == 0.01
        assert simulator.maxsteps == 1000
        assert simulator.post_process == True
        assert simulator.verbosity == 2
    
class Test_Explicit_Euler:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        def f(t, y):
            return 1.0
        problem = Explicit_Problem()
        problem.f = f
        y0 = 1.0
        self.simulator = Explicit_Euler(problem,y0)
        
    def test_integrator(self):
        """
        This tests the functionality of the method integrate.
        """
        t = 0.0
        tfinal = 1.0
        y = 1.0
        dt = 0.1
        
        values = self.simulator._integrator(t,y,tfinal,dt)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 0.100000)
        nose.tools.assert_almost_equal(y, 1.100000)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 0.200000)
        nose.tools.assert_almost_equal(y, 1.200000)
        
    def test_step(self):
        """
        This tests the functionality of the method step.
        """
        self.simulator.h = 0.1
        [tplus, yplus] = self.simulator.step(t=0.0, y=1.0)
        
        assert tplus == 0.1
        assert yplus == 1.1
    
class Test_RungeKutta34:
    
    def setUp(self):
        """
        This function sets up the test case.
        """ 
        def f(t, y):
            return 1.0
        problem = Explicit_Problem()
        problem.f = f
        y0 = 1
        self.simulator = RungeKutta34(problem,y0)
    
    def test_integrator(self):
        """
        This tests the functionality of the method integrate.
        """
        t = 0.0
        tfinal = 1.0
        y = 1.0
        nt = 0.0
        self.simulator.initstep = 0.1
        
        values = self.simulator._integrator(t,y,tfinal,nt)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 0.100000)
        nose.tools.assert_almost_equal(y, 1.100000)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 1.0)
        nose.tools.assert_almost_equal(y, 2.0)
        
    def test_adjust_stepsize(self):
        """
        This tests the functionality of the method adjust_stesize.
        """
        self.simulator.h = 0.1
        self.simulator.atol = 1.0
        self.simulator.error = 0.5
        
        self.simulator.adjust_stepsize()

        nose.tools.assert_almost_equal(self.simulator.h, 0.1000)
        
        self.simulator.atol = 0.1
        self.simulator.adjust_stepsize()
        
        nose.tools.assert_almost_equal(self.simulator.h, 0.05623413)
    
    def test_step(self):
        """
        This tests the functionality of the method step.
        """
        self.simulator.h = 0.1
        [tplus, yplus] = self.simulator.step(t=0.0, y=1.0)
        
        assert tplus == 0.1
        assert yplus == 1.1
    
class Test_RungeKutta4:
    
    def setUp(self):
        """
        This function sets up the test case.
        """ 
        def f(t, y):
            return 1.0
        problem = Explicit_Problem()
        problem.f = f
        y0 = 1
        self.simulator = RungeKutta4(problem,y0)
        
    def test_integrate(self):
        """
        This tests the functionality of the method integrate.
        """
        t = 0.0
        tfinal = 0.10
        y = 1.0
        dt = 0.05
        
        values = self.simulator._integrator(t,y,tfinal,dt)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 0.050000)
        nose.tools.assert_almost_equal(y, 1.050000)
        [t, y] = values.next()
        nose.tools.assert_almost_equal(t, 0.10)
        nose.tools.assert_almost_equal(y, 1.10)
        
    def test_step(self):
        """
        This tests the functionality of the method step.
        """
        self.simulator.h = 0.1
        [tplus, yplus] = self.simulator.step(t=0.0, y=1.0)
        
        assert tplus == 0.1
        assert yplus == 1.1
    
class Test_CVode:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        f = lambda t,y:N.array(y)
        problem = Explicit_Problem()
        problem.f = f
        y0 = [1.0]
        
        self.simulator = CVode(problem,y0)
        
    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        
        #assert self.simulator.f == 'Test function'
        assert self.simulator.y_cur == 1.0
        assert self.simulator.discr == 'Adams'
        assert self.simulator.iter == 'FixedPoint'
        assert self.simulator.maxord == 12
        
        self.simulator.discr = 'BDF'
        assert self.simulator.discr == 'BDF'
        assert self.simulator.maxord == 5
        
        nose.tools.assert_raises(Explicit_ODE_Exception, CVode, 'Test function', 'test', [1.0])
        nose.tools.assert_raises(Explicit_ODE_Exception, CVode, 'Test function', [1.0], switches0='Error')
        
        
        f = lambda t,y,sw:N.array(y)
        state_events = lambda t,x,sw: N.array(x)
        jac = lambda t,x,sw: N.zeros([len(x),len(x)])
        problem = Explicit_Problem()
        problem.f = f
        problem.state_events = state_events
        problem.jac = jac
        y0 = [1.0]

        switches = [True]
        
        simulator = CVode(problem ,y0, switches0=switches)
        
        assert simulator.f == problem.f
        assert simulator.switches == switches
        assert simulator.y_cur == 1.0
        assert simulator.problem_data['RHS'] == simulator.f
        assert simulator.problem_data['JAC'] == simulator.jac
        assert simulator.problem_data['ROOT'] == simulator.state_events
         
        
    def test_discr_method(self):
        """
        This tests the functionality of the property discr.
        """
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, 'Test')
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, 1)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, [1.0, 1])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, {'Test':'case'})
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, 5.1)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_discr_method, ['Test'])
        
        self.simulator.discr = 'BDF'
        assert self.simulator.discr == 'BDF'
        self.simulator.discr = 'Adams'
        assert self.simulator.discr == 'Adams'
    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x: N.array([x[1], -9.82])       #Defines the rhs
        jac = lambda t,x: N.array([[0.,1.],[0.,0.]]) #Defines the jacobian
        
        exp_mod = Explicit_Problem()
        exp_mod.f = f
        exp_mod.jac = jac
        
        exp_sim = CVode(exp_mod, [1.0,0.0])
        exp_sim.discr='BDF'
        exp_sim.iter='Newton'
        exp_sim.simulate(5.,100)
        
        assert exp_sim.stats[3] == 0
        nose.tools.assert_almost_equal(exp_sim.y[-1][0], -121.75000143, 4)
        exp_sim.reset()
        exp_sim.usejac=False
        exp_sim.simulate(5.,100)

        nose.tools.assert_almost_equal(exp_sim.y[-1][0], -121.75000143, 4)
        assert exp_sim.stats[3] > 0
        
    def test_switches(self):
        """
        This tests that the switches are actually turned when override.
        """
        f = lambda t,x,sw: N.array([1.0])
        state_events = lambda t,x,sw: N.array([x[0]-1.])
        def handle_event(solver, event_info):
            solver.switches = [False] #Override the switches to point to another instance
        
        mod = Explicit_Problem()
        mod.f = f
        mod.state_events = state_events
        mod.handle_event = handle_event
        
        sim = CVode(mod, [0.0], switches0=[True])
        assert sim.switches[0] == True
        sim.simulate(3)
        assert sim.switches[0] == False
    
    def test_iter_method(self):
        """
        This tests the functionality of the property iter.
        """
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, 'Test')
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, 1)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, 0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, ['Test'])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, [1.0, 1])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_iter_method, 11.1)
        
        self.simulator.iter = 'Newton'
        assert self.simulator.iter == 'Newton'
        self.simulator.iter = 'FixedPoint'
        assert self.simulator.iter == 'FixedPoint'
    
    def test_initial_step(self):
        """
        This tests the functionality of the property initstep.
        """
        
        nose.tools.assert_raises(Explicit_ODE_Exception, self.simulator._set_initial_step, 'Test')
        nose.tools.assert_raises(Explicit_ODE_Exception, self.simulator._set_initial_step, ['Test'])
        
        assert self.simulator.initstep == 0.0
        self.simulator.initstep = 10.0
        assert self.simulator.initstep == 10.0
        self.simulator.initstep = 1
        assert self.simulator.initstep == 1.0
    
    def test_interpolate(self):
        """
        This tests the functionality of the method interpolate.
        """
        f = lambda t,x: N.array(x**0.25)
        
        prob = Explicit_Problem()
        prob.f = f
        
        sim = CVode(prob, [1.0])
        sim.simulate(10., 100)
        y100 = sim.y
        t100 = sim.t
        sim.reset()
        sim.simulate(10.)
        
        nose.tools.assert_almost_equal(y100[-2], sim.interpolate(9.9,0),5)
        
    def test_handle_result(self):
        """
        This function tests the handle result.
        """
        f = lambda t,x: x**0.25
        def handle_result(solver,t,y):
            solver.temp+=1
        
        
        prob = Explicit_Problem()
        prob.f = f
        prob.handle_result = handle_result
        
        sim = CVode(prob, [1.0])
        sim.temp = 0
        sim.store_cont = True
        sim.simulate(10., 100)
        print sim.temp
        assert sim.temp == 101
        sim.simulate(20.)
        print sim.temp
        assert sim.temp == 118
        
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        self.simulator.discr='Adams'
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, "Test")
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, 1.0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, -1.0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, [1,1])
        
        self.simulator.maxord = -1
        assert self.simulator.maxord == 1
        self.simulator.maxord = 2
        assert self.simulator.maxord == 2
        self.simulator.maxord = 13
        assert self.simulator.maxord == 12
        
        self.simulator.discr='BDF'
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, "Test")
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, 1.0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, -1.0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_max_ord, [1,1])
        
        self.simulator.maxord = -1
        assert self.simulator.maxord == 1
        self.simulator.maxord = 2
        assert self.simulator.maxord == 2
        self.simulator.maxord = 6
        assert self.simulator.maxord == 5

    def test_max_order_discr(self):
        """
        This tests the maximum order when the discretization is changed.
        """
        self.simulator.maxord = 7
        assert self.simulator.maxord == 7
        
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 7
        self.simulator.discr = 'BDF'
        assert self.simulator.maxord == 5
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 12
        self.simulator.maxord = 4
        self.simulator.discr = 'BDF'
        assert self.simulator.maxord == 4
        self.simulator.discr = 'Adams'
        assert self.simulator.maxord == 4

    def test_is_disc(self):
        """
        This tests the functionality of the property is_disc.
        """
        class Prob_CVode(Explicit_Problem):
            f = lambda self,t,y,sw: N.array([1.0])
            y0 = [1.0]
            state_events = lambda self,t,y,sw: N.array([t-1.0, t])
        
        switches=[False,True]
        f = Prob_CVode()
        simulator = CVode(f, switches0=switches)
        simulator(2.)
        
        #assert simulator.t[-1] == 1.0 #For now, this error serves as prof of discontinuities
        #assert simulator.is_disc == True
    
    def test_completed_step(self):
        """
        This tests the functionality of the method completed_step.
        """
        f = lambda t,x: x**0.25
        def completed_step(solver):
            solver._nstepevents += 1
        mod = Explicit_Problem()
        mod.f = f
        mod.completed_step = completed_step
        
        y0 = [1.0]
        sim = CVode(mod, y0)
        sim._nstepevents = 0
        
        sim.simulate(2., 100)
        assert len(sim.t) == 101
        assert sim._nstepevents == 19
        
        sim = CVode(mod, y0)
        sim._nstepevents = 0
        sim.simulate(2.)
        assert len(sim.t) == 20
        assert sim._nstepevents == 19
