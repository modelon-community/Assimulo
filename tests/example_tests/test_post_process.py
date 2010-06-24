import nose
import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Implicit_ODE import *
from Assimulo.Problem import Explicit_Problem
from Assimulo.Problem import Implicit_Problem

def run_example():
    
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    def post_process(solver, t, y):
        solver.temp += 1
    global exp_mod
    global exp_sim
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.post_process = post_process
    exp_mod.problem_name = 'Simple CVode Example'

    y0 = 4.0 #Initial conditions
    exp_sim = CVode(exp_mod,y0) #Create a CVode solver
    exp_sim.temp = 0
    exp_sim.post_process = True
    exp_sim.verbosity = exp_sim.SCREAM
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.simulate(5.,10) #Simulate 5 seconds

    assert len(exp_sim.y) == 11
    assert exp_sim.stats != None
    assert exp_sim.problem_name == exp_mod.problem_name
    assert exp_sim.y[0] == y0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=8)
    
    exp_sim = CVode(exp_mod,y0)
    exp_sim.temp = 0
    exp_sim.post_process = True
    assert exp_sim.stats == None
    exp_sim.simulate(5.) #Simulate 5 seconds
    assert exp_sim.stats != None
    

def run_example2():
    def f(t,y):
        ydot = 1.0
        return N.array([ydot])
    
    def time_event_fcn(t,y,sw):
        events = [1.0, 2.0, 2.5, 3.0]
        for ev in events:
            if t < ev:
                tnext = ev
                break
            else:
                tnext = None
        #print tnext
        return tnext
        
    def handle_event(solver, event_info):
        solver.y[-1]+= 1.0
    
    global exp_mod
    global exp_sim
    
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.time_event_fcn = time_event_fcn
    exp_mod.handle_event = handle_event
    exp_mod.y0 = 0.0
    exp_mod.problem_name = 'Test explicit time event with post process.'
    
    #CVode
    exp_sim = CVode(exp_mod)
    exp_sim.post_process = True
    exp_sim(5.,10)
    #exp_sim.plot()
    

def test_CVode():
    """
    Runs the tests
    """
    run_example()
    
def test_post_and_time_event():
    """
    Tests an explicit solver with time event and with post process.
    """
    run_example2()

if __name__=='__main__':
    run_example2()
