About
-------------------

Assimulo is a Cython/Python based simulation package that allows for simulation of both ordinary differential equations of the form $f(t,y)$ (explicit problems) and differential algebraic equations of the form $f(t, y, y_d)$ (implicit problems).
Assimulo currently supports Explicit Euler, adaptive Runge-Kutta of order 4 and Runge-Kutta of order 4.
Assimulo also wraps the popular [SUNDIALS](https://computation.llnl.gov/casc/sundials/main.html) solvers CVode (for explicit problems) and IDA (for implicit problems).
[Hairer's](http://www.unige.ch/~hairer/software.html) codes Radau5, Rodas and Dopri5 is also available.

CVode and IDA supports discontinuous systems of the form $f(t, y, s_w)$ and $f(t, y, y_d, s_w)$ with event(root)-functions defined as $g(t, y, s_w)$ and $g(t, y, y_d, s_w)$ where $s_w$ should be fixed during the integration.

The package comes with a Problem specifications class and subclasses corresponding to different types of problems. Choices are `Implicit_Problem` and `Explicit_Problem`.
To define and solve a problem, first import the solver of choice and the appropriate `Implicit`/`Explicit` class. Furthermore, define your function $f$ and the initial conditions to pass to the problem class constructor. Then, create your solver, set the attributes
(method, absolute/relative tolerance etc.) and use the simulate method to simulate.

For more information about Assimulo, documentation tutorial etc, see
the docs from the older [JModelica](http://www.jmodelica.org/assimulo) website.

Sundials Compliance
-------------------
Current Assimulo development aims for compliancy with Sundials [v2.7.0](https://github.com/LLNL/sundials/releases/tag/v2.7.0).

Some optional features in `Assimulo>=3.4.1`` are built on the modified Sundials version that can be accessed [here](https://github.com/modelon-community/sundials/releases/tag/v2.7.0-1).

Installation
-------------------
See the INSTALL file or the installation page on http://www.jmodelica.org/assimulo.

Contributing
-------------------
For information about contributing, see https://github.com/modelon/contributing.

Contact
-------------------
For contact, either via email or the easiest is immediately by filing an issue or pull request here on GitHub.
Email: `christian.winther` at `modelon.com`.
