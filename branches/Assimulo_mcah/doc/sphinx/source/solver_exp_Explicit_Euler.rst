
Explicit_Euler
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from assimulo.explicit_ode import Explicit_Euler
    from assimulo.problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <assimulo.problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`store_cont <assimulo.explicit_ode.Explicit_Euler.store_cont>` Defines how the method handle_result is called.
- :class:`maxsteps <assimulo.explicit_ode.Explicit_Euler.maxsteps>` Determines the maximum number of steps the solver is allowed to take to finish the simulation.
- :class:`problem_name <assimulo.explicit_ode.Explicit_Euler.problem_name>` Defines the name of the problem.
- :class:`verbosity <assimulo.explicit_ode.Explicit_Euler.verbosity>` Defines the verbosity of the integrator.

Simulate the problem.

    :class:`Explicit_Euler.simulate(tfinal, ncp) <assimulo.explicit_ode.Explicit_Euler.simulate>` 

Plot the solution.

    :class:`Explicit_Euler.plot() <assimulo.explicit_ode.Explicit_Euler.plot>`

Information.

- :class:`Explicit_Euler.print_statistics() <assimulo.explicit_ode.Explicit_Euler.print_statistics>` Prints the run-time statistics for the problem.

.. note::

    Only IDA and CVode supports discontinuous systems.