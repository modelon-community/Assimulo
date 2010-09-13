
CVode
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from assimulo.explicit_ode import CVode
    from assimulo.problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <assimulo.problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`maxord <assimulo.explicit_ode.CVode.maxord>` This determines the maximal order that is be used by the solver.
- :class:`usejac <assimulo.explicit_ode.CVode.usejac>` This sets the option to use the user defined jacobian.
- :class:`initstep <assimulo.explicit_ode.CVode.initstep>` This determines the initial step-size to be used in the integration.
- :class:`iter <assimulo.explicit_ode.CVode.iter>` This determines the iteration method that is be used by the solver.
- :class:`discr <assimulo.explicit_ode.CVode.discr>` This determines the discretization method.
- :class:`is_disc <assimulo.explicit_ode.CVode.is_disc>` Method to test if we are at an event.
- :class:`disc_info <assimulo.explicit_ode.CVode.disc_info>` Attribute that returns information about an occured event.
- :class:`stats <assimulo.explicit_ode.CVode.stats>` Attribute that returns the run-time statistics from the Integrator.
- :class:`maxh <assimulo.explicit_ode.CVode.maxh>` Defines the maximal step-size that is to be used by the solver.
- :class:`rtol <assimulo.explicit_ode.CVode.rtol>` Defines the relative tolerance that is to be used by the solver.
- :class:`atol <assimulo.explicit_ode.CVode.atol>` Defines the absolute tolerance(s) that is to be used by the solver.
- :class:`store_cont <assimulo.explicit_ode.CVode.store_cont>` Defines how the method handle_result is called.
- :class:`maxsteps <assimulo.explicit_ode.CVode.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finish the simulation.
- :class:`problem_name <assimulo.explicit_ode.CVode.problem_name>` Defines the name of the problem.
- :class:`verbosity <assimulo.explicit_ode.CVode.verbosity>` Defines the verbosity of the integrator.

Sensitivity parameters:

- :class:`pbar <assimulo.explicit_ode.CVode.pbar>` Specifies the order of magnitude for the parameters.
- :class:`dqrhomax <assimulo.explicit_ode.CVode.dqrhomax>` Specifies the selection parameters used in deciding switching between a simultaneous or separate approximation of the two terms in the sensitivity residual.
- :class:`dqtype <assimulo.explicit_ode.CVode.dqtype>` Specifies the difference quotient type in the sensitivity calculations and can be either 'IDA_CENTERED' or 'IDA_FORWARD'.
- :class:`usesens <assimulo.explicit_ode.CVode.usesens>` Specifies if the sensitivity calculations should be used or turned off.
- :class:`suppress_sens <assimulo.explicit_ode.CVode.suppress_sens>` Specifies whether sensitivity variables are included in the error test or not.
- :class:`sensmethod <assimulo.explicit_ode.CVode.sensmethod>` Specifies the sensitivity solution method.
- :class:`maxsensiter <assimulo.explicit_ode.CVode.maxsensiter>` Specifies the maximum number of nonlinear solver iterations for sensitivity variables per step.


Methods:

- :class:`interpolate <assimulo.explicit_ode.CVode.interpolate>`
- :class:`interpolate_sensitivity <assimulo.explicit_ode.CVode.interpolate_sensitivity>`

Simulate the problem.

    :class:`CVode.simulate(tfinal, ncp) <assimulo.explicit_ode.CVode.simulate>` 

Plot the solution.

    :class:`CVode.plot() <assimulo.explicit_ode.CVode.plot>`

Information.

- :class:`CVode.print_statistics() <assimulo.explicit_ode.CVode.print_statistics>` Prints the run-time statistics for the problem.
- :class:`CVode.print_event_info() <assimulo.explicit_ode.CVode.print_event_info>` Prints the event information.
- :class:`CVode.echo_options() <assimulo.explicit_ode.CVode.echo_options>` Echo the current solver options.
