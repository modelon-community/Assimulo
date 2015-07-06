
Simulation in three steps
============================

Assimulo requires three steps for performing a simulation

   - definition of the problem
   
   - selection of a solver
   
   - executing the simulation
   
Problem
--------

A problem is the formulation of the differential equation to be solved
together with its initial values and special properties, e.g. events.

The problem class collects the complete mathematical description without 
referring to the solution process or solution control parameters.

The formulation of the problem might differ, it might be
an implicit problem (often also called a differential-algebraic system), 
an explicit problem, an overdetermined problem etc.

For each of these problem types special problem classes are provided 
which account for the different mathematical formulations.

More on problem classes see :doc:`problem`.


Solver
-------

Assimulo collects a variety of numerical ODE solvers. They are strongly 
related to the problem type. Some of them can only solve implicit problems.
The use of others might be restricted to explicit problems.

More on solver classes see :doc:`explicit`  and :doc:`implicit`.


.. note::

   Do not confound the notion implicit problem with the notion implicit solver.
   There are implicit solvers which are implemented for explicit methods.

A solver together with a problem is used to initiate a solver instance.

Execution:
-----------

The solver instance has a method simulate which has to be invoked to 
perform simulation, i.e. to numerically solve the differential equation.  
