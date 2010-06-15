import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem

def run_example():
        
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])

    #Define an Assimulo problem
    exp_mod = Explicit_Problem()
    exp_mod.f = f #Sets the rhs into the problem
    exp_mod.y0 = 4.0 #Sets the initial conditions
    exp_mod.problem_name = 'Simple Explicit Example'

    #Explicit Euler
    #==============
    exp_sim = Explicit_Euler(exp_mod) #Create a explicit Euler solver
    
    #Simulate
    exp_sim.simulate(3,100) #Simulate 3 seconds
    exp_sim.simulate(5,100) #Simulate 2 second more
    
    #Plot
    exp_sim.plot() #Plot the solution
    
    #==============

    
    #RungeKutta4
    #===========
    exp_sim = RungeKutta4(exp_mod) #Create a RungeKutta4 solver
    
    #Simulate
    exp_sim.simulate(5, 100) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()
    
    #===========
    
    
    #RungeKutta34
    #=============
    exp_sim = RungeKutta34(exp_mod) #Create a RungeKutta34 solver
    
    #Simulate
    exp_sim.simulate(5) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()

    #Reset the solver
    exp_sim.reset()
    
    exp_sim.initstep = 0.1 #Sets the initial step, default = 0.01
    
    #Simulate
    exp_sim.simulate(5) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot()
    
    #=============


if __name__=='__main__':
    run_example()
