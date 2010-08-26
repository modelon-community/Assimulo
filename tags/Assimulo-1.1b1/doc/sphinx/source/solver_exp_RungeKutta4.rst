
RungeKutta4
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from assimulo.explicit_ode import RungeKutta4
    from assimulo.problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <assimulo.problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`store_cont <assimulo.explicit_ode.RungeKutta4.store_cont>` Defines how the method handle_result is called.
- :class:`maxsteps <assimulo.explicit_ode.RungeKutta4.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finish the simulation.
- :class:`problem_name <assimulo.explicit_ode.RungeKutta4.problem_name>` Defines the name of the problem.
- :class:`verbosity <assimulo.explicit_ode.RungeKutta4.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`RungeKutta4.simulate(tfinal, ncp) <assimulo.explicit_ode.RungeKutta4.simulate>` 

Plot the solution.

    :class:`RungeKutta4.plot() <assimulo.explicit_ode.RungeKutta4.plot>`

Information.

- :class:`RungeKutta4.print_statistics() <assimulo.explicit_ode.RungeKutta4.print_statistics>` Prints the run-time statistics for the problem.

.. note::

    Only IDA and CVode supports discontinuous systems.