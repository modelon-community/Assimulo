import numpy as N
from assimulo.explicit_ode import Radau5 as Radau5Exp
from assimulo.problem import Explicit_Problem, Implicit_Problem
from assimulo.implicit_ode import Radau5 as Radau5Imp

def run_example_explicit():
    
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
    exp_mod.problem_name = 'Van der Pol (explicit)'
    
    #Define an explicit solver
    y0 = [2.0,-0.6] #Initial conditions
    exp_sim = Radau5Exp(exp_mod,y0) #Create a Radau5 solver
    
    #Sets the parameters
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.initstep = 1.e-4 #Initial step-size
    
    #Simulate
    exp_sim.simulate(2.) #Simulate 2 seconds
    
    #Plot
    exp_sim.plot(mask=[1,0],marker='o') #Plot the solution
    
    assert exp_sim._nsteps == 285 #For test purpose

def run_example_implicit():
    
    #Define the residual
    def f(t,y,yd):
        eps = 1.e-6
        my = 1./eps
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        
        res_0 = yd[0]-yd_0
        res_1 = yd[1]-yd_1
        
        return N.array([res_0,res_1])
        
    #Define an Assimulo problem
    imp_mod = Implicit_Problem()
    imp_mod.f = f
    imp_mod.problem_name = 'Van der Pol (implicit)'
    
    #Define an explicit solver
    y0 = [2.0,-0.6] #Initial conditions
    yd0 = [-.6,-200000.]
    imp_sim = Radau5Imp(imp_mod,y0,yd0) #Create a Radau5 solver
    
    #Sets the parameters
    imp_sim.atol = 1e-4 #Default 1e-6
    imp_sim.rtol = 1e-4 #Default 1e-6
    imp_sim.initstep = 1.e-4 #Initial step-size

    #Simulate
    imp_sim.simulate(2.) #Simulate 2 seconds
    
    #Plot
    imp_sim.plot(mask=[1,0],marker='o') #Plot the solution
    
    assert imp_sim._nsteps == 272

if __name__=='__main__':
    run_example_explicit()
    run_example_implicit()
