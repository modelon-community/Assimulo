
RungeKutta34
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from Assimulo.Explicit_ODE import RungeKutta34
    from Assimulo.Problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <Assimulo.Problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`initstep <Assimulo.Explicit_ODE.RungeKutta34.initstep>` This determines the initial step-size to be used in the integration.
- :class:`maxsteps <Assimulo.Explicit_ODE.RungeKutta34.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finnish the simulation.
- :class:`problemname <Assimulo.Explicit_ODE.RungeKutta34.problemname>` Defines the name of the problem.
- :class:`verbosity <Assimulo.Explicit_ODE.RungeKutta34.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`RungeKutta34.simulate(tfinal, ncp) <Assimulo.Explicit_ODE.RungeKutta34.simulate>` 

Plot the solution.

    :class:`RungeKutta34.plot() <Assimulo.Explicit_ODE.RungeKutta34.plot>`

.. note::

    Only IDA and CVode supports discontinuous systems.