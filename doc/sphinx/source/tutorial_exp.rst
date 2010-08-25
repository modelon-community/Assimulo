Explicit Problems (Commonly ODEs)
=================================

In the next few sections we show how to use the solver CVode for solving an explicit ordinary differential equation (commonly ODE) on the form,

.. math::

    \dot{y} = f(t,y), \quad y(t_0) = y_0

Problem formulation
-----------------------

The problem consists of a 'right-hand-side function' (in the explicit ODE case) together with initial conditions for the time and the states. 
This has to be packed into a problem class:

The 'right-hand-side function' takes as input the time (t) and the states (y) and returns the calculated state derivatives (yd).

An example of a 'right-hand-side function' rhs is shown below (Python)::

    import numpy as N

    def rhs(t,y):
        A =N.array([[0,1],[-2,-1]])
        yd=N.dot(A,y)
        
        return yd

The initial conditions to the rhs needs to also to be specified::

    y0=N.array([1.0,1.0])
    t0=0.0


Both, the rhs-function and the initial conditions are packed now into the problem class, 
being the Python equivalent to an explicit ODE::
    
    from assimulo.problem import Explicit_Problem  #Imports the problem formulation from Assimulo
    
    model = Explicit_Problem()             #Create an Assimulo problem
    model.f = rhs                          #This is how the rhs connects to the Assimulo problem
    model.y0 = y0                          # Here we provide the initial conditions
    model.problem_name = 'Linear Test ODE' #Specifies the name of problem (optional)

Creating an Assimulo solver
------------------------------    
And now we create the actual solver object::

    from assimulo.explicit_ode import CVode #Imports the solver CVode from Assimulo

    sim = CVode(model, t0)

Simulate
----------

To simulate the problem using the default values, simply specify the final time of the simulation and simulate::

    tfinal = 10.0        #Specify the final time
    
    sim.simulate(tfinal) #Use the .simulate method to simulate and provide the final time
    
This returns all sorts of information in the prompt, the statistics of the solver, how many function calls were needed, among others. 
Also information about the solver, which options the problem was solved with. 
The *simulate* method can also take the number of communication points for which the solution should be returned. 
This is specified by a second argument, e.g. *sim.simulate(tfinal,200)* means that the result vector should contain 200 equally spaced points.

To plot the simulation result, use the plot method::

    sim.plot() #Plots the result
    
The plot is given below,

.. image:: tutorialCVodePlot.svg
   :align: center
   :scale: 50 %

together with the statistics. ::

    Final Run Statistics: Linear Test ODE 

     Number of Error Test Failures             = 4
     Number of F-Eval During Jac-Eval          = 0
     Number of Function Evaluations            = 153
     Number of Jacobian Evaluations            = 0
     Number of Nonlinear Convergence Failures  = 0
     Number of Nonlinear Iterations            = 149
     Number of Root Evaluations                = 0
     Number of Steps                           = 84

    Solver options:

     Solver                  :  CVode
     Linear Multistep Method :  Adams
     Nonlinear Solver        :  FixedPoint
     Maxord                  :  12
     Tolerances (absolute)   :  1e-06
     Tolerances (relative)   :  1e-06

    Elapsed simulation time: 0.0 seconds.

For the complete example, :download:`tutorialCVode.py`

Setting options and parameters
-------------------------------------

To control the integration, SUNDIALS provides a number of parameters and options which of a few have been lifted up to Python.

Here are some:

    - **atol** The absolute tolerance. This controls the global error increment in every step. It can be set as a scalar or (preferably) as a vector, which defines the absolute tolerance for every solution component.
    
    - **rtol** The relative tolerance. It is a scalar.
    
    - **maxord** The maximal order. It cannot exceed 12 in case of Adams methods or 5 in case of BDF.
    
    - **discr** The discretization method, Adams or BDF. (Only for CVode)
    
    - **iter** The type of corrector iteration, FixedPoint or Newton (Only for CVode)

Example.::

    sim.atol=N.array([1.0,0.1])*1.e-5
    sim.rtol=1.e-8
    sim.maxord=3
    sim.discr='BDF'
    sim.iter='Newton'

For the full range of available options see each solver, for example `CVode <solver_CVode.html>`_ or `IDA <solver_IDA.html>`_ .
