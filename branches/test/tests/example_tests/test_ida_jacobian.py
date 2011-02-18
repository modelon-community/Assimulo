import nose
import numpy as N
from assimulo import testattr
from assimulo.implicit_ode import IDA
from assimulo.problem import Implicit_Problem

def run_example():
    
    def f(t,y,yd):
        
        res_0 = yd[0]-y[2]
        res_1 = yd[1]-y[3]
        res_2 = yd[2]+y[4]*y[0]
        res_3 = yd[3]+y[4]*y[1]+9.82
        #res_4 = y[0]**2+y[1]**2-1
        res_4 = y[2]**2+y[3]**2-y[4]*(y[0]**2+y[1]**2)-y[1]*9.82

        return N.array([res_0,res_1,res_2,res_3,res_4])
    
    def jac(c,t,y,yd):
        jacobian = N.zeros([len(y),len(y)])
        
        #Derivative
        jacobian[0,0] = 1*c
        jacobian[1,1] = 1*c
        jacobian[2,2] = 1*c
        jacobian[3,3] = 1*c
        
        #Differentiated
        jacobian[0,2] = -1
        jacobian[1,3] = -1
        jacobian[2,0] = y[4]
        jacobian[3,1] = y[4]
        jacobian[4,0] = y[0]*2*y[4]*-1
        jacobian[4,1] = y[1]*2*y[4]*-1-9.82
        jacobian[4,2] = y[2]*2
        jacobian[4,3] = y[3]*2
        
        #Algebraic
        jacobian[2,4] = y[0]
        jacobian[3,4] = y[1]
        jacobian[4,4] = -(y[0]**2+y[1]**2)
        
        return jacobian
        
    
    global imp_mod
    global imp_sim
    
    imp_mod = Implicit_Problem()
    imp_mod.f = f #Sets the residual function
    imp_mod.jac = jac #Sets the jacobian
    imp_mod.algvar = [1.0,1.0,1.0,1.0,0.0] #Set the algebraic components
    imp_mod.problem_name = 'Test Jacobian'

    y0 = [1.0,0.0,0.0,0.0,5] #Initial conditions
    yd0 = [0.0,0.0,0.0,-9.82,0.0] #Initial conditions
    imp_sim = IDA(imp_mod,y0,yd0) #Create a IDA solver
    imp_sim.atol = 1e-6 #Default 1e-6
    imp_sim.rtol = 1e-6 #Default 1e-6
    imp_sim.suppress_alg = True 
    imp_sim.make_consistent('IDA_YA_YDP_INIT')
    imp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points
    #imp_sim.plot() #Plot the solution
    
    #assert exp_sim.problemname == exp_mod.Problem_Name
    #assert exp_sim.y[0] == y0
    nose.tools.assert_almost_equal(imp_sim.y[-1][0], 0.94019949, places=7)
    nose.tools.assert_almost_equal(imp_sim.y[-1][1], -0.34095123, places=7)
    nose.tools.assert_almost_equal(imp_sim.y[-1][2], -0.88198309, places=7)
    nose.tools.assert_almost_equal(imp_sim.yd[-1][0], -0.88198927, places=7)
    nose.tools.assert_almost_equal(imp_sim.yd[-1][1], -2.43227069, places=7)
    nose.tools.assert_almost_equal(imp_sim.yd[-1][2], -9.43939383, places=7)


if __name__=='__main__':
    run_example()

@testattr(stddist = True)
def test_jac():
    """
    Runs the test.
    """
    run_example()
