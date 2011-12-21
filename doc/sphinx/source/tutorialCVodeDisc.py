#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Tutorial example showing how to use the explicit solver CVode for a discontinious problem. 
To run the example simply type,

    run tutorialCVodeDisc.py (in IPython)
    
or,

    python tutorialCVodeDisc.py (in a command prompt)
"""
import numpy as N
import pylab as P
from assimulo.problem import Explicit_Problem
from assimulo.solvers.sundials import CVode


def run_example():

    def pendulum(t,y,sw):
        """
        The ODE to be simulated. The parameter sw should be fixed during 
        the simulation and only be changed during the event handling.
        """
        l=1.0
        g=9.81
        yd_0 = y[1]
        yd_1 = -g/l*N.sin(y[0])
            
        return N.array([yd_0, yd_1])

    def state_events(t,y,sw):
        """
        This is our function that keep track of our events, when the sign
        of any of the events has changed, we have an event.
        """
        if sw[0]:
            e_0 = y[0]+N.pi/4.
        else:
            e_0 = y[0]

        return N.array([e_0])


    def handle_event(solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        state_info = event_info[0] #We are only interested in state events info

        if state_info[0] != 0: #Check if the first event function have been triggered
            
            if solver.sw[0]: #If the switch is True the pendulum bounces
                solver.y[1] = -0.9*solver.y[1] #Change the velocity and lose energy
                
            solver.sw[0] = not solver.sw[0] #Change event function

    #Initial values
    y0 = [N.pi/2.0, 0.0] #Initial states
    t0 = 0.0             #Initial time
    switches0 = [True]   #Initial switches

    #Create an Assimulo Problem
    mod = Explicit_Problem(pendulum, y0, t0, sw0=switches0)
    
    mod.state_events = state_events #Sets the state events to the problem
    mod.handle_event = handle_event #Sets the event handling to the problem
    mod.name = 'Pendulum with events'   #Sets the name of the problem

    #Create an Assimulo solver (CVode)
    sim = CVode(mod)
    
    #Specifies options 
    sim.discr = 'Adams'     #Sets the discretization method
    sim.iter = 'FixedPoint' #Sets the iteration method
    sim.rtol = 1.e-8        #Sets the relative tolerance
    sim.atol = 1.e-6        #Sets the absolute tolerance
    
    #Simulation
    ncp = 200     #Number of communication points
    tfinal = 10.0 #Final time
    
    t,y = sim.simulate(tfinal, ncp) #Simulate
    
    #Print event information
    sim.print_event_data()
    
    #Plot
    P.plot(t,y)
    P.show()
    
if __name__=='__main__':
    run_example()
