import nose
import numpy as N
from Assimulo.Sundials import Sundials_Exception, Sundials
from Assimulo.Implicit_ODE import IDA
from Assimulo.Explicit_ODE import CVode
from Assimulo.Problem import Implicit_Problem
from Assimulo.Problem import Explicit_Problem

class Test_Sundials:
    
    def setUp(self):
        """
        This sets up the test case.
        """
        class Prob_IDA(Implicit_Problem):
            f = lambda self,t,y,yd,sw: [y[0]-1.0]
            event_fcn = lambda self,t,y,yd,sw: [t-1.0, t]
            y0 = [1.0]
            yd0 = [1.0]
            
        res = Prob_IDA()
        
        class Prob_CVode(Explicit_Problem):
            f = lambda self,t,y,sw: [1.0]
            y0 = [1.0]
            event_fcn = lambda self,t,y,sw: [t-1.0, t]
        
        switches=[False, True]
        f = Prob_CVode()
        
        self.simulators = [IDA(res, switches0=switches), CVode(f, switches0=switches)]
    
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
        
    def test_stats_print(self):
        """
        This tests the functionality of the method stats_print.
        """
        pass
        
    def test_stats(self):
        """
        This tests the functionality of the property stats.
        """
        pass
        
    def test_plot_stepsize_order(self):
        """
        This tests the functionality of the method plot_stepsize_order.
        """
        pass
        
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
