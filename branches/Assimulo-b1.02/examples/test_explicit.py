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
    exp_mod.y0 = 4.0
    exp_mod.Problem_Name = 'Simple Explicit Example'

    #Explicit Euler
    exp_sim = Explicit_Euler(exp_mod) #Create a explicit Euler solver
    exp_sim.simulate(3,100) #Simulate 3 seconds
    exp_sim.simulate(5,100) #Simulate 2 second more
    #exp_sim.plot() #Plot the solution
    
    assert exp_sim.problemname == exp_mod.Problem_Name
    assert exp_sim.y[0] == 4.0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.0252255, places=4)
    
    #RungeKutta4
    exp_sim = RungeKutta4(exp_mod)
    exp_sim.simulate(5, 100) #Simulate 5 seconds
    #exp_sim.plot()
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)
    
    #RungeKutta34
    exp_sim = RungeKutta34(exp_mod)
    exp_sim.simulate(5) #Simulate 5 seconds
    #exp_sim.plot()
    print exp_sim.t[-1]
    assert len(exp_sim.t) == 62
    assert exp_sim.t[-1] == 5.0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)
    exp_sim.reset()
    exp_sim.initstep = 0.1
    exp_sim.simulate(5) #Simulate 5 seconds
    #exp_sim.plot()
    print exp_sim.t[-1]
    assert exp_sim.t[-1] == 5.0
    assert len(exp_sim.t) == 60
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)

def test_explicit():
    run_example()

if __name__=='__main__':
    run_example()
