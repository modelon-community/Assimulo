import numpy as N
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem


def run_example():
    
    #Defines the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82

        return N.array([yd_0,yd_1])
    
    #Defines the jacobian
    def jac(t,y):
        j = N.array([[0,1],[0,0]])
        return j
    
    #Defines an Assimulo explicit problem
    exp_mod = Explicit_Problem()
    exp_mod.f = f #Sets the rhs
    exp_mod.jac = jac #Sets the jacobian
    exp_mod.problem_name = 'Example using Jacobian'

    
    y0 = [1.0,0.0] #Initial conditions
    exp_sim = CVode(exp_mod,y0) #Create a CVode solver
    
    #Set the parameters
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    
    #Simulate
    exp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points
        
    #Plot
    exp_sim.plot() #Plot the solution


if __name__=='__main__':
    run_example()
