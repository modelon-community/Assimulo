import nose
from Assimulo.ODE import *

class Test_ODE:
    
    def setUp(self):
        self.simulator = ODE()
        
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        simulator = ODE()
        
        assert simulator.verbosity == simulator.NORMAL
        assert simulator.maxsteps == 10000
        
    def test_verbosity(self):
        """
        This tests the functionality of the property verbosity.
        """
        
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, 'Test')
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, 1.3)
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, [1, 31])
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, [1])
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, 5)
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_verbosity, -1)
        
        self.simulator.verbosity=1
        assert self.simulator.verbosity==1
        self.simulator.verbosity=4
        assert self.simulator.verbosity==4
        
    def test_maxsteps(self):
        """
        This tests the functionality of the property maxsteps.
        """
        
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_steps, 'Test')
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_steps, 5.0)
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_steps, -1)
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_steps, [10, 1])
        nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_steps, {'Test':'case'})
        
        self.simulator.maxsteps=100
        assert self.simulator.maxsteps==100
        self.simulator.maxsteps=1000000
        assert self.simulator.maxsteps==1000000
        self.simulator.maxsteps=1
        assert self.simulator.maxsteps==1
        
    #def test_max_eIter(self):
    #    """
    #    This tests the functionality of the property max_eIter.
    #    """
    #    
    #    nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_eIteration, 'Test')
    #    nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_eIteration, 1.0)
    #    nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_eIteration, -1)
    #    nose.tools.assert_raises(ODE_Exception, self.simulator._set_max_eIteration, -1.0)
    #    
    #    self.simulator.max_eIter = 100
    #    assert self.simulator.max_eIter == 100
    #    self.simulator.max_eIter = 1
    #    assert self.simulator.max_eIter == 1
        
    #def test_check_eIter(self):
    #    """
    #    This tests the functionality of the method check_eIter.
    #    """
    #    
    #    b = [1, 1, -1]
    #    a = [1, 1, 1]
    #    
    #    eIter = self.simulator.check_eIter(b,a)
    #    
    #    assert eIter[0] == False
    #    assert eIter[1] == False
    #    assert eIter[2] == True
    #    
    #    b = [1, 1, 1]
    #    a = [1, 1, 1]
    #    
    #    eIter = self.simulator.check_eIter(b,a)
    #    
    #    assert eIter[0] == False
    #    assert eIter[1] == False
    #    assert eIter[2] == False
