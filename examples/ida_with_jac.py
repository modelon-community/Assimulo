import numpy as N
from Assimulo.Implicit_ODE import IDA
from Assimulo.Problem import Implicit_Problem

def run_example():
    
    #Defines the residual
    def f(t,y,yd):
        
        res_0 = yd[0]-y[2]
        res_1 = yd[1]-y[3]
        res_2 = yd[2]+y[4]*y[0]
        res_3 = yd[3]+y[4]*y[1]+9.82
        #res_4 = y[0]**2+y[1]**2-1
        res_4 = y[2]**2+y[3]**2-y[4]*(y[0]**2+y[1]**2)-y[1]*9.82

        return N.array([res_0,res_1,res_2,res_3,res_4])
    
    #Defines the jacobian
    def jac(c,t,y,yd):
        jacobian = N.zeros([len(y),len(y)])
        
        #Derivative
        jacobian[0,0] = 1*c
        jacobian[1,1] = 1*c
        jacobian[2,2] = 1*c
        jacobian[3,3] = 1*c
        
        #Differentiated
        jacobian[0,2] = -1
        jacobian[1,3] = -1
        jacobian[2,0] = y[4]
        jacobian[3,1] = y[4]
        jacobian[4,0] = y[0]*2*y[4]*-1
        jacobian[4,1] = y[1]*2*y[4]*-1-9.82
        jacobian[4,2] = y[2]*2
        jacobian[4,3] = y[3]*2
        
        #Algebraic
        jacobian[2,4] = y[0]
        jacobian[3,4] = y[1]
        jacobian[4,4] = -(y[0]**2+y[1]**2)
        
        return jacobian
        
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem()
    
    #Sets the options to the problem
    imp_mod.f = f #Sets the residual function
    imp_mod.jac = jac #Sets the jacobian
    imp_mod.algvar = [1.0,1.0,1.0,1.0,0.0] #Set the algebraic components
    imp_mod.problem_name = 'Test Jacobian'

    #The initial conditons
    y0 = [1.0,0.0,0.0,0.0,5] #Initial conditions
    yd0 = [0.0,0.0,0.0,-9.82,0.0] #Initial conditions
    
    #Create an Assimulo implicit solver (IDA)
    imp_sim = IDA(imp_mod,y0,yd0) #Create a IDA solver
    
    #Sets the paramters
    imp_sim.atol = 1e-6 #Default 1e-6
    imp_sim.rtol = 1e-6 #Default 1e-6
    imp_sim.suppress_alg = True #Suppres the algebraic variables on the error test
    
    #Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
    imp_sim.make_consistency('IDA_YA_YDP_INIT')
    
    #Simulate
    imp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points
    
    #Plot
    imp_sim.plot() #Plot the solution


if __name__=='__main__':
    run_example()
    
