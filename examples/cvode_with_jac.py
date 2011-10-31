import numpy as N
from assimulo.solvers.cvode import CVode
from assimulo.problem import Explicit_Problem


def run_example():
    
    #Defines the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82
        #print y, yd_0, yd_1
        return N.array([yd_0,yd_1])
    
    #Defines the jacobian
    def jac(t,y):
        j = N.array([[0,1.],[0,0]])
        return j
    
    #Defines an Assimulo explicit problem
    y0 = [1.0,0.0] #Initial conditions

    exp_mod = Explicit_Problem(f,y0)
    
    exp_mod.jac = jac #Sets the jacobian
    exp_mod.name = 'Example using Jacobian'

    
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Set the parameters
    exp_sim.options["iter"] = 'Newton' #Default 'FixedPoint'
    exp_sim.options["discr"] = 'BDF' #Default 'Adams'
    exp_sim.options["atol"] = [1e-4, 1e-4] #Default 1e-6
    exp_sim.options["rtol"] = 1e-4 #Default 1e-6
    #exp_sim.options["continuous_output"] = True
    exp_sim.options["usejac"] = False
    #Simulate
    exp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points
        
    #Plot
    exp_sim.plot(linestyle="dashed",marker="o") #Plot the solution


if __name__=='__main__':
    run_example()
