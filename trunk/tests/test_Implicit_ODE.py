import nose
import numpy as N
from Assimulo.Implicit_ODE import *
from Assimulo.Problem import Implicit_Problem


class Test_Implicit_ODE:
    
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        def f(self, t, y, yd):
            pass
        res = Implicit_Problem()
        res.f = f
        
        nose.tools.assert_raises(Implicit_ODE_Exception, Implicit_ODE, res, 1, 'test')
        nose.tools.assert_raises(Implicit_ODE_Exception, Implicit_ODE, res, 1, 1, 'test')
        nose.tools.assert_raises(Implicit_ODE_Exception, Implicit_ODE, res, 'test', 'test', 'test')
        nose.tools.assert_raises(Implicit_ODE_Exception, Implicit_ODE, res, [1.0 , 1.0], [1.0, 'test'])
        nose.tools.assert_raises(Implicit_ODE_Exception, Implicit_ODE, None, [1.0 , 1.0, 1.0], [1.0, 1, 1])
        
        
        simulator = Implicit_ODE(res, [1 , 1.0], [2, 2.0], 1)
        assert simulator.t[0] == 1.0
        assert simulator.y[0][0] == 1.0
        assert simulator.yd[0][0] == 2.0
        
        
    def test_call(self):
        """
        This tests the functionality of the method __call__.
        """
        y0 = [0.0]
        yd0 = [1.0]
        my_Prob = Implicit_Problem()
        my_Prob.f = lambda t,x,xd: x

        simulator = IDA(my_Prob,y0,yd0)
        nose.tools.assert_raises(Implicit_ODE_Exception, simulator, -1.0)
        nose.tools.assert_raises(Implicit_ODE_Exception, simulator, 'test')
        print simulator.problem_spec[0][0]
        [t,y,yd] = simulator(1.0,10)
        
        assert len(t) == 11 #11 Due to t0 is counted as well
        
    
class Test_IDA:
    
    def setUp(self):
        """
        This function sets up the test case.
        """
        my_test_prob = Implicit_Problem()
        my_test_prob.f = 'Test function'
        y0 = [1.0]
        yd0 = [1.0]

        self.simulator = IDA(my_test_prob, y0, yd0)
        
        
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        
        assert self.simulator.res_fcn == 'Test function'
        assert self.simulator.suppress_alg == False
        assert self.simulator.algvar == [1.0]
        assert self.simulator.switches == None
        assert self.simulator.maxsteps == 10000
        assert self.simulator.verbosity == self.simulator.NORMAL
        assert self.simulator.y[0][0] == 1.0
        
        nose.tools.assert_raises(Implicit_ODE_Exception, IDA, 'Test function', 'test', [1.0])
        nose.tools.assert_raises(Implicit_ODE_Exception, IDA, 'Test function', [1.0], [1.0], switches0='Error')
        
        my_Prob = Implicit_Problem()
        my_Prob.f = 'Test function'
        my_Prob.event_fcn = lambda t,x,xd,sw: x
        
        def jac(c,t,y,yd,sw):
            re = N.zeros([len(y),len(y)])
            return re
        
        my_Prob.jac = jac
        y0 = [1.0, 1.0, 1]
        yd0 = [1, 1, 1]
        
        switches = [True, False]

        simulator = IDA(my_Prob,y0,yd0, switches0=switches)
        
        assert simulator.res_fcn == 'Test function'
        assert simulator.switches == switches
        assert simulator.yd[0][0] == 1.0
        assert simulator.problem_spec[0][0] == simulator.res_fcn
        assert simulator.problem_spec[0][1] == simulator.jac
        assert simulator.problem_spec[1][0] == simulator.event_fcn
        assert simulator.problem_spec[1][1] == switches
    
    def test_max_order(self):
        """
        This tests the functionality of the property maxord.
        """
        
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

    def test_tout1(self):
        """
        This tests the functionality of the property tout1.
        """
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_calcIC_tout1, 'Test')
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_calcIC_tout1, [1,1])
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_calcIC_tout1, 'Test')
        
        assert self.simulator.tout1 == 0.001
        self.simulator.tout1 = -0.001
        assert self.simulator.tout1 == -0.001
        self.simulator.tout1 = 1
        assert self.simulator.tout1 == 1.0
        
    def test_lsoff(self):
        """
        This tests the functionality of the property lsoff.
        """
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_lsoff, 'Test')
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_lsoff, 1.0)
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_lsoff, 0.0)
        nose.tools.assert_raises(Implicit_ODE_Exception, self.simulator._set_lsoff, [1,1])
        
        assert self.simulator.lsoff == False
        self.simulator.lsoff = True
        assert self.simulator.lsoff == True
        self.simulator.lsoff = False
        assert self.simulator.lsoff == False
    
    def test_initstep(self):
        """
        This tests the funtionality of the property initstep.
        """
        
        def f(t,y,yd):
            res_0 = yd[0] - y[1]
            res_1 = yd[1] +9.82-0.01*y[1]**2
            return N.array([res_0,res_1])
            
        mod = Implicit_Problem()
        mod.f=f
        sim = IDA(mod, y0=[5.0,0.0], yd0=[0.0,9.82])
        
        sim.simulate(2.0)

        nose.tools.assert_almost_equal(sim.y[-1][0], -13.4746473811, places=7)
        
        sim = IDA(mod, y0=[5.0,0.0], yd0=[0.0,9.82])
        sim.initstep = 1e-10
        
        sim.simulate(2.0)

        nose.tools.assert_almost_equal(sim.y[-1][0], -13.4746596311, places=7)
        
    
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        f = lambda t,x,xd: x
        jac = lambda c,t,x,xd: N.array([x*x])
        
        prob = Implicit_Problem()
        prob.f = f
        prob.jac = jac
        
        sim = IDA(prob, [0],[0])
        
        assert sim._RES[0] == f
        assert sim._RES[1] == jac
        assert sim.problem_spec[0][0] == f
        assert sim.problem_spec[0][1] == jac
        sim.usejac = False
        assert sim._RES[0] == f
        assert len(sim._RES) == 1
        assert sim.problem_spec[0][0] == f
        assert len(sim.problem_spec[0]) == 1
    
    def test_run(self):
        """
        This tests the functionality of the property run. (With jacobian)
        """
        pass
        
        
    def test_algvar(self):
        """
        This tests the functionality of the property algvar.
        """
        self.simulator.Integrator.dim = 3
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, 1)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, 1.0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, [1,'hej',1])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, {'Test':'case'})
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, [-1,0,1])
        

        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, [1.0,1.0])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_algvar, [3.0, 1.0, 1.0])
        
        vector = [1.0,0.0,1.0]
        vectorb = [True,False,True]
        vectori = [1,0,1]
        
        self.simulator.algvar = vectorb
        self.simulator.algvar = vectori
        self.simulator.algvar = vector
        nose.tools.assert_equal(self.simulator.algvar[0], vector[0])
        nose.tools.assert_equal(self.simulator.algvar[1], vector[1])
        nose.tools.assert_equal(self.simulator.algvar[2], vector[2])
        
        
    def test_suppress_alg(self):
        """
        This tests the functionality of the property suppress_alg.
        """
        
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, "Test")
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, [1,2])
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, {'Test':'case'})
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, 3)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, 0)
        nose.tools.assert_raises(Sundials_Exception, self.simulator._set_suppress_alg, 0.1)
        
        self.simulator.suppress_alg = True
        assert self.simulator.suppress_alg == True
        self.simulator.suppress_alg = False
        assert self.simulator.suppress_alg == False
        
    def test_make_consistency(self):
        """
        This tests the functionality of the method make_consistency.
        """
        def f(t,y,yd):
            res_1 = y[0] + y[1]+1.0
            res_2 = y[1]
            print res_1
            return [res_1, res_2]
        my_Prob = Implicit_Problem()
        my_Prob.f = f
        
        y0 = [2.0, 2.0]
        yd0 = [1.0 , 0.0]
        simulator = IDA(my_Prob, y0, yd0)
        print simulator.Integrator.jacobian
        [y, yd] = simulator.make_consistency('IDA_Y_INIT')
        
        nose.tools.assert_almost_equal(y[1], 0.00000)
        nose.tools.assert_almost_equal(y[0], -1.0000)
        nose.tools.assert_almost_equal(yd[0], 1.0000)
        nose.tools.assert_almost_equal(yd[1], 0.0000)
     
    def test_is_disc(self):
        """
        This tests the functionality of the property is_disc.
        """
        class Prob_IDA(Implicit_Problem):
            f = lambda self,t,y,yd,sw: [y[0]-1.0]
            event_fcn = lambda self,t,y,yd,sw: [t-1.0, t]
            y0 = [1.0]
            yd0 = [1.0]
        switches = [False,True]
        res = Prob_IDA()
        simulator = IDA(res, switches0=switches)
        simulator(2.)
        
        #assert simulator.t[-1] == 1.0 #For now, this error serves as prof of discontinuities
        #assert simulator.is_disc == True
