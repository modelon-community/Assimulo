

=============
Examples
=============


Explicit ODE
============

This example shows how to use Assimulo to solve an explicit ODE together with the solver CVode. 
In the example an object dropped in a gravitational-field is simulated.
 
(The full example code can be found in the examples catalog, ``test_CVode_Jacobian.py``.) 


We first start by importing the necessary classes from Assimulo::

    import numpy as N
    from Assimulo.Explicit_ODE import CVode
    from Assimulo.Problem import Explicit_Problem
    
Here we imported an explicit solver, CVode together with the base class for an explicit problem. Also we imported the numpy package for use of the numpy array.

Next, the 'right-hand-side' has to be defined::

    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82

        return N.array([yd_0,yd_1])
        
In this example, we also provide the Jacobian. If the Jacobian is not provided it will be approximated (when using CVode or IDA). 
To define the Jacobian we type::

    def jac(t,y):
        j = N.array([[0,1],[0,0]])
        return j
        
The Jacobian must return a matrix of appropriate size.

Now that we have defined the problem, we have to set-up the problem class for Assimulo. ::

    exp_mod = Explicit_Problem()
    exp_mod.f = f #Sets the rhs
    exp_mod.jac = jac #Sets the jacobian
    exp_mod.Problem_Name = 'Example using Jacobian'

We are now almost ready to create our solver. But first we need the initial conditions, they can either be set directly in the problem class, ::

    exp_mod.y0 = [1.0, 0.0]
    
Or, they can be provided when the solver object is created, ::

    y0 = [1.0,0.0] #Initial conditions
    exp_sim = CVode(exp_mod,y0) #Create a CVode solver
    

The problem is now ready to be simulated ::

    exp_sim.simulate(5,1000) #Simulate 5 seconds
    
Which simulates the problem five seconds using 1000 communication points. To change the default settings of the solver, see below ::

    exp_sim.iter = 'Newton'   #Default 'FixedPoint'
    exp_sim.discr = 'BDF'     #Default 'Adams'
    exp_sim.atol = 1e-4       #Default 1e-6
    exp_sim.rtol = 1e-4       #Default 1e-6
    
Here we changed the method from using the Adams-Moulton to the BDF together with Newton iteration. We also changed the tolerances. This is just a selection of the parameters that can be changed. For the full set of parameters, see the :doc:`documentation <usage>`.

The problem can also be easily plotted by using the plot command, ::

    exp_sim.plot() #Plot the solution
   
If only the first solution component has to be plotted set

    exp_sim.plot(mask=[1,0])



Implicit ODE
=============

This example shows how to use Assimulo to solve an implicit ODE together with the solver IDA. 

(The full example can be found in the examples catalogue, test_Jacobian_IDA.py)

We first start by importing the necessary classes from Assimulo::

    import numpy as N
    from Assimulo.Implicit_ODE import IDA
    from Assimulo.Problem import Implicit_Problem

Here we imported an implicit solver, IDA, together with the base class for an implicit problem. 
Also we imported the numpy package for use of the numpy array.

Next, the residual have to be defined::

    def f(t,y,yd):
        
        res_0 = yd[0]-y[2]
        res_1 = yd[1]-y[3]
        res_2 = yd[2]+y[4]*y[0]
        res_3 = yd[3]+y[4]*y[1]+9.82
        res_4 = y[2]**2+y[3]**2-y[4]*(y[0]**2+y[1]**2)-y[1]*9.82

        return N.array([res_0,res_1,res_2,res_3,res_4])

In this example, we also provide the Jacobian. If the Jacobian is not provided it will be approximated (when using IDA or CVode). To define the Jacobian we type::

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

The Jacobian has to be defined as

.. math:: 

    J = \frac{dF}{dx} + c\cdot \frac{dF}{d\dot{x}}

where *c* is the inverse step-size. For more information about the Jacobian see SUNDIALS documentation for `IDA <http://computation.llnl.gov/casc/sundials/documentation/ida_guide/node3.html>`_. 

Now that we have defined the problem, we have to set-up the problem class for Assimulo. ::

    imp_mod = Implicit_Problem()
    imp_mod.f = f                          #Sets the residual function
    imp_mod.jac = jac                      #Sets the jacobian
    imp_mod.algvar = [1.0,1.0,1.0,1.0,0.0] #Set the algebraic components
    imp_mod.Problem_Name = 'Test IDA'      #Sets the name of the problem
    
Here we created an implicit problem, *imp_mod* which we provided the residual, the Jacobian, the algebraic components and also specified the name of the problem. When the problem specific methods have been passed to the problem we create our solver::

    y0 = [1.0,0.0,0.0,0.0,5] #Initial conditions
    yd0 = [0.0,0.0,0.0,-9.82,0.0] #Initial conditions
    imp_sim = IDA(imp_mod,y0,yd0) #Create a IDA solver

Here we created an IDA solver where we also provided the initial conditions. The initial conditions could just as well been set directly in the problem. ::

    imp_mod.y0 = y0
    imp_mod.yd0 = yd0

We then set the solver attributes, ::

    imp_sim.atol = 1e-6         #Absoulte tolerance. Default 1e-6
    imp_sim.rtol = 1e-6         #Relative tolerance. Default 1e-6
    imp_sim.suppress_alg = True #Suppress the algebraic variables on the error test in case of higher index DAEs.

This is just a selection of the parameters that can be changed. For the full set of parameters, see the :doc:`documentation <usage>`.

*IDA* also comes with two methods for calculating a consistent set of initial conditions, *IDA_YA_YDP_INIT* and *IDA_Y_INIT*. For more information about the methods, see SUNDIALS documentation for `IDA <http://computation.llnl.gov/casc/sundials/documentation/ida_guide/node3.html>`_. Here we calculate initial conditions by use of the *IDA_YA_YDP_INIT* option which calculates the differential parts of *yd0* and the algebraic parts of *y0* given the differential parts of *y0* and an initial guess for the algebraic parts of *y0*. ::

    imp_sim.make_consistency('IDA_YA_YDP_INIT')
    
Whats left is just to simulate the problem ::

    imp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points

and the problem can also be easily plotted by using the plot command, ::

    imp_sim.plot() #Plot the solution
    
