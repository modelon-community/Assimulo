##
##     The sensitivity calculations are not fully implemented!
##
import numpy as N
import pylab as P
from assimulo.implicit_ode import IDA
from assimulo.problem import Implicit_Problem

def run_example():
    """
    This example show how to use Assimulo and IDA for simulating sensitivities
    for initial conditions.
    
    0 = dy1/dt - -(k01+k21+k31)*y1 - k12*y2 - k13*y3 - b1
    0 = dy2/dt - k21*y1 + (k02+k12)*y2
    0 = dy3/dt - k31*y1 + k13*y3
 
    y1(0) = p1, y2(0) = p2, y3(0) = p3
    p1=p2=p3 = 0 
    
    See http://sundials.2283335.n4.nabble.com/Forward-sensitivities-for-initial-conditions-td3239724.html
    """
    
    def f(t, y, yd,p):
        y1,y2,y3 = y
        yd1,yd2,yd3 = yd
        k01 = 0.0211
        k02 = 0.0162
        k21 = 0.0111
        k12 = 0.0124
        k31 = 0.0039
        k13 = 0.000035
        b1 = 49.3
        
        res_0 = -yd1 -(k01+k21+k31)*y1+k12*y2+k13*y3+b1
        res_1 = -yd2 + k21*y1-(k02+k12)*y2
        res_2 = -yd3 + k31*y1-k13*y3
        
        return N.array([res_0,res_1,res_2])
    
    def handle_result(solver, t ,y,yd):
        solver.t += [t]
        solver.y += [y]
        if len(solver.t) !=1:
            solver.p1 += [solver.interpolate_sensitivity(t, 0, 0)]
            solver.p2 += [solver.interpolate_sensitivity(t, 0, 1)]
            solver.p3 += [solver.interpolate_sensitivity(t, 0, 2)]
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem()
    
    #Sets the options to the problem
    imp_mod.f = f #Sets the residual function
    imp_mod.handle_result = handle_result #Change the default handling of the result
    
    #The initial conditions
    y0 = [0.0,0.0,0.0]          #Initial conditions for y
    yd0 = [49.3,0.,0.]
    p0 = [0.0, 0.0, 0.0]  #Initial conditions for parameters
    yS0 = N.array([[1,0,0],[0,1,0],[0,0,1.]])

    #Create an Assimulo explicit solver (IDA)
    imp_sim = IDA(imp_mod, y0, yd0=yd0, p0=p0)
    
    #Sets the paramters
    imp_sim.rtol = 1e-7
    imp_sim.atol = 1e-6
    imp_sim.yS0 = yS0 #Specify the initial condition for the sensitivities
    imp_sim.pbar = [1,1,1] #pbar is used to estimate the tolerances for the parameters
    imp_sim.store_cont = True #Need to be able to store the result using the interpolate methods
    imp_sim.sensmethod = 'SIMULTANEOUS' #Defines the sensitvity method used
    imp_sim.suppress_sens = False            #Dont suppress the sensitivity variables in the error test.
    
    imp_sim.p1 = [] #Vector for storing the p1 result
    imp_sim.p2 = [] #Vector for storing the p2 result
    imp_sim.p3 = [] #Vector for storing the p3 result
    
    #Simulate
    imp_sim.simulate(400) #Simulate 400 seconds
    
    #Plot
    P.figure(2)
    P.plot(imp_sim.t[1:], N.array(imp_sim.p1)[:,0],
           imp_sim.t[1:], N.array(imp_sim.p1)[:,1],
           imp_sim.t[1:], N.array(imp_sim.p1)[:,2])
    P.title("Parameter p1")
    P.legend(("p1/dy1","p1/dy2","p1/dy3"))
    P.figure(3)
    P.plot(imp_sim.t[1:], N.array(imp_sim.p2)[:,0],
           imp_sim.t[1:], N.array(imp_sim.p2)[:,1],
           imp_sim.t[1:], N.array(imp_sim.p2)[:,2])
    P.title("Parameter p2")
    P.legend(("p2/dy1","p2/dy2","p2/dy3"))
    P.figure(4)
    P.plot(imp_sim.t[1:], N.array(imp_sim.p3)[:,0],
           imp_sim.t[1:], N.array(imp_sim.p3)[:,1],
           imp_sim.t[1:], N.array(imp_sim.p3)[:,2])
    P.title("Parameter p3")
    P.legend(("p3/dy1","p3/dy2","p3/dy3"))
    imp_sim.plot() #Plot the solution

if __name__=='__main__':
    run_example()
