
RungeKutta4
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from Assimulo.Explicit_ODE import RungeKutta4
    from Assimulo.Problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <Assimulo.Problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`maxsteps <Assimulo.Explicit_ODE.RungeKutta4.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finnish the simulation.
- :class:`problemname <Assimulo.Explicit_ODE.RungeKutta4.problemname>` Defines the name of the problem.
- :class:`verbosity <Assimulo.Explicit_ODE.RungeKutta4.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`RungeKutta4.simulate(tfinal, ncp) <Assimulo.Explicit_ODE.RungeKutta4.simulate>` 

Plot the solution.

    :class:`RungeKutta4.plot() <Assimulo.Explicit_ODE.RungeKutta4.plot>`

.. note::

    Only IDA and CVode supports discontinuous systems.