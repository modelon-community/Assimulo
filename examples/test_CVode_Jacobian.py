import nose
import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem


def run_example():
    
    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82

        return N.array([yd_0,yd_1])
        
    def jac(t,y):
        j = N.array([[0,1],[0,0]])
        return j
        
    global exp_mod
    global exp_sim
    
    exp_mod = Explicit_Problem()
    exp_mod.f = f #Sets the rhs
    exp_mod.jac = jac #Sets the jacobian
    exp_mod.Problem_Name = 'Example using Jacobian'

    y0 = [1.0,0.0] #Initial conditions
    exp_sim = CVode(exp_mod,y0) #Create a CVode solver
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim(5,1000) #Simulate 5 seconds
    exp_sim.plot() #Plot the solution
    
    assert exp_sim.problemname == exp_mod.Problem_Name
    assert exp_sim.y[0][0] == y0[0]
    nose.tools.assert_almost_equal(exp_sim.y[-1][0], -121.75011017723, places=8)
    nose.tools.assert_almost_equal(exp_sim.y[-1][1], -49.1000000, places=4)


if __name__=='__main__':
    run_example()
