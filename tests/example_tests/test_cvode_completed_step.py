import nose
import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem


def run_example():

    def f(t,y):
        ydot = 1.0
        return N.array([ydot])

    
    def completed_step(solver):
        solver.numb = solver.numb+1
        return 0
        
    global exp_sim
    
        
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.completed_step = completed_step
    exp_mod.y0 = 0.0
    exp_mod.problem_name = 'Test CVode completed step.'
    
    #CVode
    exp_sim = CVode(exp_mod)
    
    #Testing
    exp_sim.numb = 0
    
    exp_sim(5.)

    assert len(exp_sim.y)-1 == exp_sim.numb

    exp_sim.reset()
    
    exp_sim.numb = 0
    
    exp_sim(5.,100)

    assert 7 == exp_sim.numb
    
    
    def post_process(solver, t, y):
        solver.tt += [t]
        solver.yy += y
        solver.test_post =solver.test_post+1
        
    exp_mod.post_process = post_process
    
    exp_sim = CVode(exp_mod)
    exp_sim.post_process = True
    exp_sim.test_post = 0
    exp_sim.numb = 0
    
    exp_sim.tt = []
    exp_sim.yy = []
    
    exp_sim(5.)
    assert len(exp_sim.y)-1 == exp_sim.numb
    assert exp_sim.test_post == exp_sim.numb
    
    exp_sim.reset()
    exp_sim.tt = []
    exp_sim.yy = []
    
    exp_sim.test_post = 0
    exp_sim.numb = 0

    exp_sim(5.,10)
    assert 7 == exp_sim.numb
    assert 10 == exp_sim.test_post
    assert exp_sim.tt[0] == 0.5
    assert exp_sim.tt[1] == 1.0
    assert exp_sim.tt[-1] == 5.0
    
def run_example2():
    global exp_mod
    global exp_sim
    
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])

    def completed_step(solver):
        if solver.t[-1] > 10 and solver.reini:
            solver.reini = False
            return 1
        else:
            return 0
    
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.completed_step = completed_step
    exp_mod.y0 = 4.0
    exp_mod.problem_name = 'Example Completed Steps'
    
    exp_sim = CVode(exp_mod)
    exp_sim.reini = True
    
    exp_sim.simulate(20.)
    
    x = exp_sim.stats.keys()[4]
    assert exp_sim.stats[x] == 77

def test_cvode_completed_step():
    """
    Runs the test_cvode_completed_step.py
    """
    run_example()
    
def test_cvode_completed_step_statistics():
    """
    Runs the test_cvode_completed_step_statistics.
    """
    run_example2()

if __name__=='__main__':
    run_example()
    run_example2()
