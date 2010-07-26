import numpy as N
from Assimulo.Explicit_ODE import Radau5
from Assimulo.Problem import Explicit_Problem

def run_example():
    
    #Define the rhs
    def f(t,y):
        eps = 1.e-6
        my = 1./eps
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        
        return N.array([yd_0,yd_1])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.problem_name = 'Van der Pol'
    
    #Define an explicit solver
    y0 = [2.0,-0.6] #Initial conditions
    exp_sim = Radau5(exp_mod,y0) #Create a Radau5 solver
    
    #Sets the parameters
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.initstep = 1.e-4 #Initial step-size
    
    #Simulate
    exp_sim.simulate(2.) #Simulate 2 seconds
    
    #Plot
    exp_sim.plot(mask=[1,0],marker='o') #Plot the solution
    
    assert exp_sim._nsteps == 285 #For test purpose

if __name__=='__main__':
    run_example()
