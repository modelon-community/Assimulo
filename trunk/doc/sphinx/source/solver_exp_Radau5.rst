
Radau5
=================================

Usage
--------------

Import the solver together with correct problem:: 

    from assimulo.explicit_ode import Radau5
    from assimulo.problem import Explicit_Problem

Define the problem. 

    :class:`Explicit_Problem <assimulo.problem.Explicit_Problem>`

Create and modify the solver parameters.

Parameters:

- :class:`usejac <assimulo.explicit_ode.Radau5.usejac>` This sets the option to use the user defined jacobian.
- :class:`quot2 <assimulo.explicit_ode.Radau5.quot2>` If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
- :class:`initstep <assimulo.explicit_ode.Radau5.initstep>` This determines the initial step-size to be used in the integration.
- :class:`quot1 <assimulo.explicit_ode.Radau5.quot1>` If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
- :class:`fac2 <assimulo.explicit_ode.Radau5.fac2>` Parameters for step-size selection.
- :class:`fac1 <assimulo.explicit_ode.Radau5.fac1>` Parameters for step-size selection.
- :class:`safe <assimulo.explicit_ode.Radau5.safe>` The safety factor in the step-size prediction.
- :class:`thet <assimulo.explicit_ode.Radau5.thet>` Value for determine if the Jacobian is to be recomputed or not.
- :class:`fnewt <assimulo.explicit_ode.Radau5.fnewt>` Stopping criterion for Newton's method, usually chosen <1.
- :class:`newt <assimulo.explicit_ode.Radau5.newt>` Maximal number of Newton iterations.
- :class:`maxh <assimulo.explicit_ode.Radau5.maxh>` Defines the maximal step-size that is to be used by the solver.

Methods:

- :class:`interpolate <assimulo.explicit_ode.Radau5.interpolate>`

Simulate the problem.

    :class:`Radau5.simulate(tfinal, ncp) <assimulo.explicit_ode.Radau5.simulate>` 

Plot the solution.

    :class:`Radau5.plot() <assimulo.explicit_ode.Radau5.plot>`

Information.

- :class:`Radau5.print_statistics() <assimulo.explicit_ode.Radau5.print_statistics>` Prints the run-time statistics for the problem.
- :class:`Radau5.echo_options() <assimulo.explicit_ode.Radau5.echo_options>` Echo the current solver options.

.. note::

    Only IDA and CVode supports discontinuous systems.
