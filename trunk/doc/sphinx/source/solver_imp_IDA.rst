
IDA
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from assimulo.implicit_ode import IDA
    from assimulo.problem import Implicit_Problem

Define the problem. 

    :class:`Implicit_Problem <assimulo.problem.Implicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`lsoff <assimulo.implicit_ode.IDA.lsoff>` Boolean value to turn OFF Sundials LineSearch when calculating initial conditions.
- :class:`initstep <assimulo.implicit_ode.IDA.initstep>` This determines the initial step-size to be used in the integration.
- :class:`algvar <assimulo.implicit_ode.IDA.algvar>` A vector for defining which variables are differential and which are algebraic.
- :class:`usejac <assimulo.implicit_ode.IDA.usejac>` This sets the option to use the user defined jacobian.
- :class:`maxord <assimulo.implicit_ode.IDA.maxord>` This determines the maximal order that is be used by the solver.
- :class:`tout1 <assimulo.implicit_ode.IDA.tout1>` Sets the value used in the internal Sundials function for determine initial conditions.
- :class:`suppress_alg <assimulo.implicit_ode.IDA.suppress_alg>` A boolean flag which indicates that the error-tests are suppressed on algebraic variables.
- :class:`is_disc <assimulo.implicit_ode.IDA.is_disc>` Method to test if we are at an event.
- :class:`disc_info <assimulo.implicit_ode.IDA.disc_info>` Attribute that returns information about an occured event.
- :class:`stats <assimulo.implicit_ode.IDA.stats>` Attribute that returns the run-time statistics from the Integrator.
- :class:`maxh <assimulo.implicit_ode.IDA.maxh>` Defines the maximal step-size that is to be used by the solver.
- :class:`rtol <assimulo.implicit_ode.IDA.rtol>` Defines the relative tolerance that is to be used by the solver.
- :class:`atol <assimulo.implicit_ode.IDA.atol>` Defines the absolute tolerance(s) that is to be used by the solver.
- :class:`store_cont <assimulo.implicit_ode.IDA.store_cont>` Defines how the method handle_result is called.
- :class:`maxsteps <assimulo.implicit_ode.IDA.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finish the simulation.
- :class:`problem_name <assimulo.implicit_ode.IDA.problem_name>` Defines the name of the problem.
- :class:`verbosity <assimulo.implicit_ode.IDA.verbosity>` Defines the verbosity of the integrator.

Sensitivity parameters:

- :class:`pbar <assimulo.implicit_ode.IDA.pbar>` Specifies the order of magnitude for the parameters.
- :class:`DQrhomax <assimulo.implicit_ode.IDA.DQrhomax>` Specifies the selection parameters used in deciding switching between a simultaneous or separate approximation of the two terms in the sensitivity residual.
- :class:`DQtype <assimulo.implicit_ode.IDA.DQtype>` Specifies the difference quotient type in the sensitivity calculations and can be either 'IDA_CENTERED' or 'IDA_FORWARD'.
- :class:`usesens <assimulo.implicit_ode.IDA.usesens>` Specifies if the sensitivity calculations should be used or turned off.
- :class:`suppress_sens <assimulo.implicit_ode.IDA.suppress_sens>` Specifies whether sensitivity variables are included in the error test or not.
- :class:`sensmethod <assimulo.implicit_ode.IDA.sensmethod>` Specifies the sensitivity solution method.
- :class:`maxsensiter <assimulo.implicit_ode.IDA.maxsensiter>` Specifies the maximum number of nonlinear solver iterations for sensitivity variables per step.


Methods:

- :class:`make_consistent <assimulo.implicit_ode.IDA.make_consistent>`
- :class:`interpolate <assimulo.implicit_ode.IDA.interpolate>`
- :class:`interpolate_sensitivity <assimulo.implicit_ode.IDA.interpolate_sensitivity>`

Simulate the problem.

    :class:`IDA.simulate(tfinal, ncp) <assimulo.implicit_ode.IDA.simulate>` 

Plot the solution.

    :class:`IDA.plot() <assimulo.implicit_ode.IDA.plot>`

Information.

- :class:`IDA.print_statistics() <assimulo.implicit_ode.IDA.print_statistics>` Prints the run-time statistics for the problem.
- :class:`IDA.print_event_info() <assimulo.implicit_ode.IDA.print_event_info>` Prints the event information.
- :class:`IDA.echo_options() <assimulo.implicit_ode.IDA.echo_options>` Echo the current solver options.
