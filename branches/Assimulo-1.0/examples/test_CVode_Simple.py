import nose
import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem

def run_example():
    
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    global exp_mod
    global exp_sim
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.Problem_Name = 'Simple CVode Example'

    y0 = 4.0 #Initial conditions
    exp_sim = CVode(exp_mod,y0) #Create a CVode solver
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.simulate(5) #Simulate 3 seconds
    #exp_sim.plot() #Plot the solution
    
    assert exp_sim.problemname == exp_mod.Problem_Name
    assert exp_sim.y[0] == y0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=8)

def test_CVode():
    """
    Runs the tests
    """
    run_example()

if __name__=='__main__':
    run_example()
