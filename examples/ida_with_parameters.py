##
##     The sensitivity calculations are not fully implemented!
##
import numpy as N
from Assimulo.Implicit_ODE import IDA
from Assimulo.Problem import Implicit_Problem

def run_example():
    #This is the same example from the Sundials package (idasRoberts_FSA_dns.c)
    #----
    #This simple example problem for IDA, due to Robertson, 
    #is from chemical kinetics, and consists of the following three 
    #equations:
    #
    #   dy1/dt = -p1*y1 + p2*y2*y3
    #   dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
    #   0   = y1 + y2 + y3 - 1

    def f(t, y, yd, p):
        
        res1 = -p[0]*y[0]+p[1]*y[1]*y[2]-yd[0]
        res2 = p[0]*y[0]-p[1]*y[1]*y[2]-p[2]*y[1]**2-yd[1]
        res3 = y[0]+y[1]+y[2]-1
        
        return N.array([res1,res2,res3])
        
    def handle_result(solver, t ,y,yd):
        solver.t += [t]
        solver.y += [y]
        solver.p1 += [solver.interpolate_sensitivity(t, 0, 0)]
        solver.p2 += [solver.interpolate_sensitivity(t, 0, 1)]
        solver.p3 += [solver.interpolate_sensitivity(t, 0, 2)]
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem()
    
    #Sets the options to the problem
    imp_mod.f = f #Sets the residual function
    imp_mod.handle_result = handle_result #Change the default handling of the result
    
    #The initial conditons
    y0 = [1.0, 0.0, 0.0]        #Initial conditions for y
    yd0 = [0.1, 0.0, 0.0]       #Initial conditions for dy/dt
    p0 = [0.040, 1.0e4, 3.0e7]  #Initial conditions for parameters
    global imp_sim
    #Create an Assimulo implicit solver (IDA)
    imp_sim = IDA(imp_mod,y0,yd0,p0=p0) #Create a IDA solver
    
    #Sets the paramters
    imp_sim.atol = N.array([1.0e-8, 1.0e-14, 1.0e6])
    imp_sim.algvar = [1.0,1.0,0.0]
    imp_sim.suppress_alg = True #Suppres the algebraic variables on the error test
    imp_sim.store_cont = True #Store data continuous during the simulation
    imp_sim.p1 = [] #Vector for storing the p1 result
    imp_sim.p2 = [] #Vector for storing the p2 result
    imp_sim.p3 = [] #Vector for storing the p3 result
    
    #Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
    imp_sim.make_consistency('IDA_YA_YDP_INIT')
    
    #Simulate
    imp_sim.simulate(4,400) #Simulate 4 seconds with 400 communication points
    
    #Plot
    imp_sim.plot() #Plot the solution
    
    

if __name__=='__main__':
    run_example()
    
