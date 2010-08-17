#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Tutorial example showing how to use the explicit solver CVode. To run the example
simply type,

    run tutorialCVode.py (in IPython)
    
or,

    python tutorialCVode.py (in a command prompt)
"""
import numpy as N
from Assimulo.Explicit_ODE import CVode #Imports the solver CVode from Assimulo
from Assimulo.Problem import Explicit_Problem #Imports the problem formulation from Assimulo   

def run_example():

    def rhs(t,y):
        A =N.array([[0,1],[-2,-1]])
        yd=N.dot(A,y)
        
        return yd
        
    y0=N.array([1.0,1.0])
    t0=0.0
        
    model = Explicit_Problem() #Create an Assimulo problem
    model.f = rhs #This is how the rhs connects to the Assimulo problem
    model.problem_name = 'Linear Test ODE'
         

    sim = CVode(model, y0, t0) #Create the solver CVode

    tfinal = 10.0 #Specify the final time
        
    sim.simulate(tfinal) #Use the .simulate method to simulate and provide the final time
        
    sim.plot()

if __name__=='__main__':
    run_example()
