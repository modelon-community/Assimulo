import nose
import numpy as N
from assimulo.non_linear_problem import *
from assimulo.kinsol import *


class Test_KINSOL:
    
    def setUp(self):
        
        """
        sets up the test case
        """
        
        self.solver = KINSOL()
        
        class Prob1(NL_Problem):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
        self.p1 = Prob1()
        
        class Prob_no_f(NL_Problem):
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,_x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
        self.no_f = Prob_no_f()
        
        class Prob_no_x0(NL_Problem):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            
            def set_x0(self,x0):
                self._x0 =x0
                
            def get_x0(self):
                return self._x0
            
        self.no_x0 = Prob_no_f()
        
        class Prob_no_getx0(NL_Problem):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,x0):
                self._x0 = x0
                
        self.no_getx0 = Prob_no_getx0()
        
        class Prob_Jac(NL_Problem):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            jac_called = False
            
            def jac(self,x):
                self.jac_called = True
                print " Jacobian called"
                return N.zeros(3)
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
        self.pb_jac = Prob_Jac()    
        
        
            
    def test_solve(self):
        # tests for problems without all member attributes/methods
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.no_f)
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.no_x0)
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.no_getx0)
        
        # test solver for simple problem
        res = self.solver.solve(self.p1)
        nose.tools.assert_almost_equal(res[0],1.0,5)
        nose.tools.assert_almost_equal(res[1],2.0,5)
        nose.tools.assert_almost_equal(res[2],3.0,5)
        
        self.p1.set_x0('a string')
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.p1)
        self.p1.set_x0(5)
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.p1)
        self.p1.set_x0(0.6)
        nose.tools.assert_raises(KINSOL_Exception,self.solver.solve,self.p1)
        
    def test_jac_usage(self):
        # test is supplied jacobian is called
        try:
            self.solver.solve(self.pb_jac,use_jac=True)
        except :
            pass
        
        nose.tools.assert_true(self.pb_jac.jac_called)
        
        self.pb_jac.jac_called = False
        try:
            self.solver.solve(self.pb_jac,use_jac=False)
        except :
            pass
        
        nose.tools.assert_false(self.pb_jac.jac_called,'Jacobian used although use_jac = false')
        

        
        
        