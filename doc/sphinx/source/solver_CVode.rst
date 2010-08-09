
CVode
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from Assimulo.Explicit_ODE import CVode
    from Assimulo.Problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <Assimulo.Problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`maxord <Assimulo.Explicit_ODE.CVode.maxord>` This determines the maximal order that is be used by the solver.
- :class:`usejac <Assimulo.Explicit_ODE.CVode.usejac>` This sets the option to use the user defined jacobian.
- :class:`initstep <Assimulo.Explicit_ODE.CVode.initstep>` This determines the initial step-size to be used in the integration.
- :class:`iter <Assimulo.Explicit_ODE.CVode.iter>` This determines the iteration method that is be used by the solver.
- :class:`discr <Assimulo.Explicit_ODE.CVode.discr>` This determines the discretization method.
- :class:`is_disc <Assimulo.Explicit_ODE.CVode.is_disc>` Method to test if we are at an event.
- :class:`disc_info <Assimulo.Explicit_ODE.CVode.disc_info>` Attribute that returns information about an occured event.
- :class:`stats <Assimulo.Explicit_ODE.CVode.stats>` Attribute that returns the run-time statistics from the Integrator.
- :class:`maxh <Assimulo.Explicit_ODE.CVode.maxh>` Defines the maximal step-size that is to be used by the solver.
- :class:`rtol <Assimulo.Explicit_ODE.CVode.rtol>` Defines the relative tolerance that is to be used by the solver.
- :class:`atol <Assimulo.Explicit_ODE.CVode.atol>` Defines the absolute tolerance(s) that is to be used by the solver.
- :class:`maxsteps <Assimulo.Explicit_ODE.CVode.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finnish the simulation.
- :class:`problemname <Assimulo.Explicit_ODE.CVode.problemname>` Defines the name of the problem.
- :class:`verbosity <Assimulo.Explicit_ODE.CVode.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`CVode.simulate(tfinal, ncp) <Assimulo.Explicit_ODE.CVode.simulate>` 

Plot the solution.

    :class:`CVode.plot() <Assimulo.Explicit_ODE.CVode.plot>`

Information.

- :class:`CVode.print_statistics() <Assimulo.Explicit_ODE.CVode.print_statistics>` Prints the run-time statistics for the problem.
- :class:`CVode.print_event_info() <Assimulo.Explicit_ODE.CVode.print_event_info>` Prints the event information.
