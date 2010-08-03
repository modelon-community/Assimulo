import numpy as N
import nose
from Assimulo.Implicit_ODE import IDA
from Assimulo.Problem import Implicit_Problem

"""
An example with event iteration and with three switches.

t=0     , [False, True, True]   (Start of simulation)
t=1 (1) , [False, True, False]  (Found a root at t=1)
t=1 (2) , [False, False, False] (Second iteration at t=1)
t=1 (3) , [True, False, False]  (Third iteration at t=1)
t=10    , [True, False, False]  (End of simulation)

"""

#Extend Assimulos problem definition
class Extended_Problem(Implicit_Problem):
    
    #Sets the initial conditons directly into the problem
    y0 = [0.0, -1.0, 0.0]
    yd0 = [-1.0, 0.0, 0.0]
    switches0 = [False,True,True]
    algvar = [1.0, 0.0, 0.0] #Determine which variables are differential and algebraic
    
    
    #The residual
    def f(self,t,y,yd,sw):
        """
        This is our function we are trying to simulate. During simulation
        the parameter sw should be fixed so that our function is continuous
        over the interval. The parameters sw should only be changed when the
        integrator has stopped.
        """
        res_0 = -yd[0] + (1.0 if sw[0] else -1.0)
        res_1 = -y[1] + (-1.0 if sw[1] else 3.0)
        res_2 = -y[2] + (0.0 if sw[2] else 2.0)
        
        return N.array([res_0,res_1,res_2])

    #Sets a name to our function
    problem_name = 'Function with consistency problem'
    
    #The event function
    def state_events(self,t,y,yd,sw):
        """
        This is our function that keep track of our events, when the sign
        of any of the events has changed, we have an event.
        """
        event_0 = y[1] - 1.0 
        event_1 = -y[2] + 1.0
        event_2 = -t + 1.0
        
        return N.array([event_0,event_1,event_2])
    
    
    #Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        while True: #Event Iteration
            self.event_switch(solver, event_info) #Turns the switches
            
            b_mode = self.state_events(solver.t_cur, solver.y_cur, solver.yd_cur, solver.switches)
            self.init_mode(solver) #Pass in the solver to the problem specified init_mode
            a_mode = self.state_events(solver.t_cur, solver.y_cur, solver.yd_cur, solver.switches)
            
            event_info = self.check_eIter(b_mode, a_mode)
            
            if solver.verbosity >= solver.SCREAM:
                print 'Event iteration?: ', event_info
                
            if not True in event_info: #Breaks the iteration loop
                break
    
    #Helper function for handle_event
    def event_switch(self, solver, event_info):
        """
        Turns the switches.
        """
        for i in range(len(event_info)): #Loop across all event functions
            if event_info[i] != 0:
                solver.switches[i] = not solver.switches[i] #Turn the switch
        
        if solver.verbosity >= solver.LOUD:
            print 'New switches: ', solver.switches
    
    #Helper function for handle_event
    def check_eIter(self, before, after):
        """
        Helper function for handle_event to determine if we have event
        iteration.
        
            Input: Values of the event indicator functions (state_events)
            before and after we have changed mode of operations.
        """
        
        eIter = [False]*len(before)
        
        for i in range(len(before)):
            if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
                eIter[i] = True
                
        return eIter
    
    def init_mode(self, solver):
        """
        Initialize the DAE with the new conditions.
        """
        solver.make_consistency('IDA_YA_YDP_INIT') #Calculate new initial conditions.
                                                   #see SUNDIALS IDA documentation
                                                   #on the option 'IDA_YA_YDP_INIT'



def test_disc():
    """
    Runs the tests
    """
    iter_mod = Extended_Problem() #Create the problem
    iter_sim = IDA(iter_mod) #Create the solver
    iter_sim.verbosity = iter_sim.SCREAM #Set the verbosity
    iter_sim.simulate(10.0,100) #Simulate 10 seconds with 100 communications points
    
    assert iter_sim.yd[-1][0] == 1.0
    assert iter_sim.yd[-1][1] == 0.0
    assert iter_sim.yd[-1][2] == 0.0

    nose.tools.assert_almost_equal(iter_sim.y[-1][0], 8.0000000, 5)
    nose.tools.assert_almost_equal(iter_sim.y[-1][1], 3.0000000, 5)
    nose.tools.assert_almost_equal(iter_sim.y[-1][2], 2.0000000, 5)
    

if __name__=="__main__":
    test_disc()
    

    
