import nose
import numpy as N
from assimulo import testattr
from assimulo.sundials import Sundials_Exception, Sundials
from assimulo.implicit_ode import IDA
from assimulo.explicit_ode import CVode
from assimulo.problem import Implicit_Problem
from assimulo.problem import Explicit_Problem

class Test_Sundials:
    
    def setUp(self):
        """
        This sets up the test case.
        """
        class Prob_IDA(Implicit_Problem):
            f = lambda self,t,y,yd,sw: N.array([y[0]-1.0])
            state_events = lambda self,t,y,yd,sw: N.array([t-1.0, t])
            y0 = [1.0]
            yd0 = [1.0]
            
        res = Prob_IDA()
        
        class Prob_CVode(Explicit_Problem):
            f = lambda self,t,y,sw: N.array([1.0])
            y0 = [1.0]
            state_events = lambda self,t,y,sw: N.array([t-1.0, t])
        
        switches=[False, True]
        f = Prob_CVode()
        
        self.simulators = [IDA(res, switches0=switches), CVode(f, switches0=switches)]
        
        mod = Implicit_Problem()
        mod.f = lambda t,y,yd,p: N.array([0.0])
        y0 = [1.0]
        yd0 = [1.0]
        p0 = [1.0]
    
        self.sim = IDA(mod, y0, yd0,p0=p0)
    
    @testattr(stddist = True)
    def test_atol(self):
        """
        This tests the functionality of the property atol.
        """
        assert self.simulators[1].atol == 1.0e-6
        assert self.simulators[0].atol == 1.0e-6
        
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_atol, 5)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_atol, -1.0)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_atol, [1.0, 1.0])
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_atol, "Test")
            
            self.simulators[i].atol = 1.0e-5
            assert self.simulators[i].atol == 1.0e-5
            self.simulators[i].atol = 1.0
            assert self.simulators[i].atol == 1.0
            self.simulators[i].atol = 1001.0
            assert self.simulators[i].atol == 1001.0
            
            self.simulators[i].Integrator.dim = 3
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_atol, [1.0, 1.0])
            self.simulators[i].atol = [1.0, 1.0, 1.0]
            assert self.simulators[i].atol == [1.0, 1.0, 1.0]
            self.simulators[i].atol = N.array([1.0, 1.0, 1.0])
            assert self.simulators[i].atol[0] == 1.0
    
    
    @testattr(stddist = True)
    def test_rtol(self):
        """
        This tests the functionality of the property rtol.
        """
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_rtol, 5)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_rtol, -1.0)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_rtol, [1.0, 1.0])
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_rtol, "Test")
            
            self.simulators[i].rtol = 1.0e-5
            assert self.simulators[i].rtol == 1.0e-5
            self.simulators[i].rtol = 1.0
            assert self.simulators[i].rtol == 1.0
            self.simulators[i].rtol = 1001.0
            assert self.simulators[i].rtol == 1001.0
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the functionality of the property maxh.
        """
        for i in range(len(self.simulators)):
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_max_h, 5)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_max_h, -1.0)
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_max_h, [1.0, 1.0])
            nose.tools.assert_raises(Sundials_Exception, self.simulators[i]._set_max_h, "Test")
            
            self.simulators[i].maxh = 1.0e-5
            assert self.simulators[i].maxh == 1.0e-5
            self.simulators[i].maxh = 1.0
            assert self.simulators[i].maxh == 1.0
            self.simulators[i].maxh = 1001.0
            assert self.simulators[i].maxh == 1001.0
    
    @testattr(stddist = True)    
    def test_stats_print(self):
        """
        This tests the functionality of the method stats_print.
        """
        pass
    
    @testattr(stddist = True)
    def test_stats(self):
        """
        This tests the functionality of the property stats.
        """
        pass
    
    @testattr(stddist = True)    
    def test_plot_stepsize_order(self):
        """
        This tests the functionality of the method plot_stepsize_order.
        """
        pass
    
    @testattr(stddist = True)    
    def test_disc_info(self):
        """
        This tests the functionality of the property disc_info.
        """
        
        self.simulators[1](2)
        
        [t, info] = self.simulators[1].disc_info

        #assert t == 1.0
        assert info[0] != 0
        assert info[1] == 0
        
        self.simulators[0](2)
        
        [t, info] = self.simulators[0].disc_info

        #assert t == 1.0
        assert info[0] != 0
        assert info[1] == 0
    
    @testattr(stddist = True)    
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
        
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqtype, 1)
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqtype, 'IDA_CE')
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqtype, [1])
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqtype, -1)
    
    @testattr(stddist = True)
    def test_dqrhomax(self):
        """
        Tests the property of DQrhomax.
        """
        assert self.sim.dqrhomax == 0.0 #Test the default value.
        
        self.sim.dqrhomax = 1.0
        assert self.sim.dqrhomax == 1.0
        self.sim.dqrhomax = 10
        assert self.sim.dqrhomax == 10
        
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqrhomax, -1)
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqrhomax, 'str')
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqrhomax, [])
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_dqrhomax, -10)
    
    @testattr(stddist = True)    
    def test_usesens(self):
        """
        Tests the property of usesens.
        """
        assert self.sim.usesens == True #Test the default value.
        
        self.sim.usesens = False
        assert self.sim.usesens == False
        self.sim.usesens = 0
        assert self.sim.usesens == False
        self.sim.usesens = 1
        assert self.sim.usesens == True

    @testattr(stddist = True)
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
        
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_sensitivity_method, 1)
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_sensitivity_method, 'IDA_CE')
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_sensitivity_method, [1])
        nose.tools.assert_raises(Sundials_Exception,self.sim._set_sensitivity_method, -1)
    
    @testattr(stddist = True)    
    def test_suppress_sens(self):
        """
        Tests the property of suppress_sens.
        """
        assert self.sim.suppress_sens == True
        self.sim.suppress_sens = False
        assert self.sim.suppress_sens == False
        self.sim.suppress_sens = 0
        assert self.sim.suppress_sens == False
        self.sim.suppress_sens = 1
        assert self.sim.suppress_sens == True
    
    @testattr(stddist = True)
    def test_maxsensiter(self):
        """
        Tests the property of maxsensiter.
        """
        assert self.sim.maxsensiter == 3 #Test the default value
        self.sim.maxsensiter = 1
        assert self.sim.maxsensiter == 1
        self.sim.maxsensiter = 10.5
        assert self.sim.maxsensiter == 10
        
        nose.tools.assert_raises(Sundials_Exception, self.sim._set_max_nonlin, 0)
        nose.tools.assert_raises(Sundials_Exception, self.sim._set_max_nonlin, 'str')
        nose.tools.assert_raises(Sundials_Exception, self.sim._set_max_nonlin, [])
        nose.tools.assert_raises(Sundials_Exception, self.sim._set_max_nonlin, -10)
    
    @testattr(stddist = True)
    def test_pbar(self):
        """
        Tests the property of pbar.
        """
        f = lambda t,y,p:N.array([0.0]*len(y))
        exp_mod = Explicit_Problem()
        exp_mod.f = f
        
        y0 = [1.0]*2
        p0 = [1000.0, -100.0]
        
        exp_sim = CVode(exp_mod,y0,p0=p0)
        
        nose.tools.assert_almost_equal(exp_sim.pbar[0], 1000.00000,4)
        nose.tools.assert_almost_equal(exp_sim.pbar[1], 100.000000,4)
        
        f = lambda t,y,yd,p: N.array([0.0]*len(y))
        imp_mod = Implicit_Problem()
        imp_mod.f = f
        
        yd0 = [0.0]*2
        
        imp_sim = IDA(imp_mod, y0,yd0,p0=p0)
        
        nose.tools.assert_almost_equal(exp_sim.pbar[0], 1000.00000,4)
        nose.tools.assert_almost_equal(exp_sim.pbar[1], 100.000000,4)
