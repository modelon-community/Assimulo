
IDA
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from Assimulo.Implicit_ODE import IDA
    from Assimulo.Problem import Implicit_Problem

Define the problem. 

    :class:`Implicit_Problem <Assimulo.Problem.Implicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`maxord <Assimulo.Implicit_ODE.IDA.maxord>` This determines the maximal order that is be used by the solver.
- :class:`lsoff <Assimulo.Implicit_ODE.IDA.lsoff>` Boolean value to turn OFF Sundials LineSearch when calculating initial conditions.
- :class:`tout1 <Assimulo.Implicit_ODE.IDA.tout1>` Sets the value used in the internal Sundials function for determine initial conditions.
- :class:`initstep <Assimulo.Implicit_ODE.IDA.initstep>` This determines the initial step-size to be used in the integration.
- :class:`suppress_alg <Assimulo.Implicit_ODE.IDA.suppress_alg>` A boolean flag which indicates that the error-tests are suppressed on algebraic variables.
- :class:`is_disc <Assimulo.Implicit_ODE.IDA.is_disc>` Method to test if we are at an event.
- :class:`usejac <Assimulo.Implicit_ODE.IDA.usejac>` This sets the option to use the user defined jacobian.
- :class:`algvar <Assimulo.Implicit_ODE.IDA.algvar>` A vector for defining which variables are differential and which are algebraic.
- :class:`disc_info <Assimulo.Implicit_ODE.IDA.disc_info>` Attribute that returns information about an occured event.
- :class:`stats <Assimulo.Implicit_ODE.IDA.stats>` Attribute that returns the run-time statistics from the Integrator.
- :class:`maxh <Assimulo.Implicit_ODE.IDA.maxh>` Defines the maximal step-size that is to be used by the solver.
- :class:`rtol <Assimulo.Implicit_ODE.IDA.rtol>` Defines the relative tolerance that is to be used by the solver.
- :class:`atol <Assimulo.Implicit_ODE.IDA.atol>` Defines the absolute tolerance(s) that is to be used by the solver.
- :class:`maxsteps <Assimulo.Implicit_ODE.IDA.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finnish the simulation.
- :class:`problemname <Assimulo.Implicit_ODE.IDA.problemname>` Defines the name of the problem.
- :class:`verbosity <Assimulo.Implicit_ODE.IDA.verbosity>` Defines the verbosity of the integrator.

Methods:

	:class:`make_consistency(method) <Assimulo.Implicit_ODE.IDA.make_consistency>` Directs IDA to try to calculate consistant initial conditions.

Simulate the problem.

    :class:`IDA.simulate(tfinal, ncp) <Assimulo.Implicit_ODE.IDA.simulate>` 

Plot the solution.

    :class:`IDA.plot() <Assimulo.Implicit_ODE.IDA.plot>`

Information.

- :class:`IDA.print_statistics() <Assimulo.Implicit_ODE.IDA.print_statistics>` Prints the run-time statistics for the problem.
- :class:`IDA.print_event_info() <Assimulo.Implicit_ODE.IDA.print_event_info>` Prints the event information.
