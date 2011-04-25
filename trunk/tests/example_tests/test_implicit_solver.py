import nose
import numpy as N
from assimulo import testattr
from assimulo.implicit_ode import *
from assimulo.problem import Implicit_Problem


def run_example():
    
    def f(t,y,yd):
        ydot = 1.0
        return N.array([yd[0]-y[0]])
    
    def time_events(t,y,yd,sw):
        events = [1.0, 2.0, 2.5, 3.0]
        for ev in events:
            if t < ev:
                tnext = ev
                break
            else:
                tnext = None
        print tnext
        return tnext
        
    def handle_event(solver, event_info):
        print solver.y_cur, solver.yd_cur
        solver.y_cur+= 1.0
        solver.yd_cur+= 1.0
        assert event_info[0] == []
        assert event_info[1] == True
        #solver.make_consistent("IDA_YA_YDP_INIT")
        
    def completed_step(solver):
        #solver.numb = solver.numb+1
        return 0
    
    imp_mod = Implicit_Problem()
    imp_mod.f = f
    imp_mod.time_events = time_events
    imp_mod.handle_event = handle_event
    imp_mod.completed_step = completed_step
    imp_mod.y0 = 1.0
    imp_mod.yd0 = 1.0
    imp_mod.problem_name = 'Test implicit time event.'
    
    #IDA
    imp_sim = IDA(imp_mod)
    imp_sim(3.5)
    #nose.tools.assert_almost_equal(imp_sim.y[19], 0.950000, 4)
    #nose.tools.assert_almost_equal(imp_sim.y[20], 1.00000, 4)
    #nose.tools.assert_almost_equal(imp_sim.y[21], 2.00000, 4)
    #nose.tools.assert_almost_equal(imp_sim.y[62], 5.950000, 4)
    nose.tools.assert_almost_equal(imp_sim.y[-1], 54.14745721, 4)
    nose.tools.assert_almost_equal(imp_sim.t[-1], 3.50000, 4)


@testattr(stddist = True)
def test_implicit_time_event():
    """
    Runs the test_implicit_time_event.py
    """
    run_example()

if __name__=='__main__':
    run_example()
