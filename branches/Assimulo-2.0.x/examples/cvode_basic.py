import numpy as N
from assimulo.solvers.cvode import CVode
from assimulo.problem import Explicit_Problem

def run_example():
    
    #Define the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, y0=4)
    exp_mod.problem_name = 'Simple CVode Example'
    
    #Define an explicit solver
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    
    #Simulate
    exp_sim.simulate(5) #Simulate 5 seconds
    
    #Plot
    exp_sim.plot() #Plot the solution
    

if __name__=='__main__':
    run_example()
