import nose
import numpy as N
from assimulo.explicit_ode import *
from assimulo.problem import Explicit_Problem


def run_example():
    
    def f(t,y):
        ydot = 1.0
        return N.array([ydot])
    
    def time_events(t,y,sw):
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
        solver.y_cur+= 1.0
        assert event_info[0] == []
        assert event_info[1] == True
    
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.time_events = time_events
    exp_mod.handle_event = handle_event
    exp_mod.y0 = 0.0
    exp_mod.problem_name = 'Test explicit time event.'
    
    #CVode
    exp_sim = CVode(exp_mod)
    exp_sim(5.,100)
    nose.tools.assert_almost_equal(exp_sim.y[19], 0.950000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[20], 1.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[21], 2.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[62], 5.950000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[64], 7.00000, 4)
    assert len(exp_sim.t) == 105

    #Explicit_Euler
    exp_sim = Explicit_Euler(exp_mod)
    exp_sim(5.,100)
    nose.tools.assert_almost_equal(exp_sim.y[19], 0.950000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[20], 1.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[21], 2.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[62], 5.900000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[63], 5.950000, 4)
    #print len(exp_sim.t)
    assert len(exp_sim.t) == 108
    
    #RungeKutta4
    exp_sim = RungeKutta4(exp_mod)
    exp_sim(5.,100)
    nose.tools.assert_almost_equal(exp_sim.y[19], 0.950000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[20], 1.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[21], 2.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[62], 5.900000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[63], 5.950000, 4)
    assert len(exp_sim.t) == 108
    
    #RungeKutta34
    exp_sim = RungeKutta34(exp_mod)
    exp_sim(5.)
    nose.tools.assert_almost_equal(exp_sim.y[3], 1.0000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[4], 2.0000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[8], 4.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[16], 7.00000, 4)
    nose.tools.assert_almost_equal(exp_sim.y[-1], 9.00000, 4)
    assert len(exp_sim.t) == 20
        
def test_explicit_time_event():
    """
    Runs the test_explicit_time_event.py
    """
    run_example()

if __name__=='__main__':
    run_example()
