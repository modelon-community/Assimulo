import nose
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem

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
        assert simulator.t[0] == 1.0
        assert simulator.y[0][0] == 1.0
        
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
        [t,y] = simulator(1.0,10)
        
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
        
    def test_integrate(self):
        """
        This tests the functionality of the method integrate.
        """
        t = 0.0
        tfinal = 1.0
        y = 1.0
        nt = 10
        
        values = self.simulator.integrate(t,y,tfinal,nt)
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
    
    def test_integrate(self):
        """
        This tests the functionality of the method integrate.
        """
        t = 0.0
        tfinal = 1.0
        y = 1.0
        nt = 0.0
        self.simulator.initstep = 0.1
        
        values = self.simulator.integrate(t,y,tfinal,nt)
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
        nt = 2
        
        values = self.simulator.integrate(t,y,tfinal,nt)
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
        f = 'Test function'
        problem = Explicit_Problem()
        problem.f = f
        y0 = [1.0]
        
        self.simulator = CVode(problem,y0)
        
    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        
        assert self.simulator.f == 'Test function'
        assert self.simulator.y[0][0] == 1.0
        assert self.simulator.discr == 'Adams'
        assert self.simulator.iter == 'FixedPoint'
        assert self.simulator.maxord == 12
        
        self.simulator.discr = 'BDF'
        assert self.simulator.discr == 'BDF'
        assert self.simulator.maxord == 5
        
        nose.tools.assert_raises(Explicit_ODE_Exception, CVode, 'Test function', 'test', [1.0])
        nose.tools.assert_raises(Explicit_ODE_Exception, CVode, 'Test function', [1.0], switches0='Error')
        
        
        f = 'Test function'
        event_fcn = lambda t,x,sw: x
        jac = lambda t,x,sw: N.zeros([len(x),len(x)])
        problem = Explicit_Problem()
        problem.f = f
        problem.event_fcn = event_fcn
        problem.jac = jac
        y0 = [1.0]

        switches = [True]
        
        simulator = CVode(problem ,y0, switches0=switches)
        
        assert simulator.f == problem.f
        assert simulator.switches == switches
        assert simulator.y[0][0] == 1.0
        assert simulator.problem_spec[0][0] == simulator.f
        assert simulator.problem_spec[0][1] == simulator.jac
        assert simulator.problem_spec[1][0] == simulator.event_fcn
        assert simulator.problem_spec[1][1] == switches
         
        
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
        f = lambda t,x: x
        jac = lambda t,x: N.array([x*x])
        
        prob = Explicit_Problem()
        prob.f = f
        prob.jac = jac
        
        sim = CVode(prob, [0])
        
        assert sim._RHS[0] == f
        assert sim._RHS[1] == jac
        assert sim.problem_spec[0][0] == f
        assert sim.problem_spec[0][1] == jac
        sim.usejac = False
        assert sim._RHS[0] == f
        assert len(sim._RHS) == 1
        assert sim.problem_spec[0][0] == f
        assert len(sim.problem_spec[0]) == 1
        
    
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
        f = lambda t,x: x**0.25
        
        prob = Explicit_Problem()
        prob.f = f
        
        sim = CVode(prob, [1.0])
        sim.simulate(10., 100)
        y100 = sim.y
        t100 = sim.t
        sim.reset()
        sim.simulate(10.)
        
        nose.tools.assert_almost_equal(y100[-2], sim.interpolate(9.9,0),5)
        
    def test_post_process(self):
        """
        This function tests the post processing.
        """
        f = lambda t,x: x**0.25
        def post_process(solver):
            solver.temp+=1
        
        
        prob = Explicit_Problem()
        prob.f = f
        prob.post_process = post_process
        
        sim = CVode(prob, [1.0])
        sim.temp = 0
        sim.post_process = True
        sim.simulate(10., 100)
        assert sim.temp == 100
        sim.simulate(20.)
        assert sim.temp == 110
        
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

    def test_is_disc(self):
        """
        This tests the functionality of the property is_disc.
        """
        class Prob_CVode(Explicit_Problem):
            f = lambda self,t,y,sw: [1.0]
            y0 = [1.0]
            event_fcn = lambda self,t,y,sw: [t-1.0, t]
        
        switches=[False,True]
        f = Prob_CVode()
        simulator = CVode(f, switches0=switches)
        simulator(2.)
        
        #assert simulator.t[-1] == 1.0 #For now, this error serves as prof of discontinuities
        #assert simulator.is_disc == True
