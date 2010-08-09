
Explicit_Euler
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from Assimulo.Explicit_ODE import Explicit_Euler
    from Assimulo.Problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <Assimulo.Problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`maxsteps <Assimulo.Explicit_ODE.Explicit_Euler.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finnish the simulation.
- :class:`problemname <Assimulo.Explicit_ODE.Explicit_Euler.problemname>` Defines the name of the problem.
- :class:`verbosity <Assimulo.Explicit_ODE.Explicit_Euler.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`Explicit_Euler.simulate(tfinal, ncp) <Assimulo.Explicit_ODE.Explicit_Euler.simulate>` 

Plot the solution.

    :class:`Explicit_Euler.plot() <Assimulo.Explicit_ODE.Explicit_Euler.plot>`

.. note::

    Only IDA and CVode supports discontinuous systems.