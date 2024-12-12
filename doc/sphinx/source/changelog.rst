
==========
Changelog
==========
--- Assimulo-3.6.0---
    * Added get_sundials_version function (import from assimulo.solvers.sundials).
    * Fixed bug that Radau5ODE would not count state events in statistics.
    * Removed tests from the Assimulo installation.
    * Changed testing framework from `nose` to `pytest`.

--- Assimulo-3.5.2 ---
    * Allow to build without distutils for Python>=3.12 support

--- Assimulo-3.5.1 ---
    * Fixed build with sundials 7.x
    * Added "std=legacy" as default fortran compile flag. 
      This addresses compilation issues commonly experienced with more modern Fortran standards.

--- Assimulo-3.5.0 ---
    * Changed "numpy.float" to equivalent "numpy.float64" due to DeprecationWarnings in numpy >= 1.20.
    * Improved examples with sparse jacobians by omitting the zeros in the jacobians.
    * Upgraded to Cython >= 3.
    * Removed deprecated build option "sundials_with_superlu", which is checked automatically.
    * Added support for Sundials 7.x

--- Assimulo-3.4.3 ---
    * Improved compliance with newer scipy version by instead using corresponding numpy calls when more suitable.

--- Assimulo-3.4.2 ---
    * Updated an import statement in one of the examples to resolve compliance issues with scipy 1.10.1.

--- Assimulo-3.4.1 ---
    * Restored functionality of CVode with 'external_event_detection' = True, which broke with 3.4.0.
    * CVode now allows rtol = 0.
    * Added support for CVode (without sensitivities) with non-negative relative tolerance vector.
      Requires a SUNDIALS installation supporting this, see README.

--- Assimulo-3.4.0 ---
    * Removed Radau5ODE Fortran implementation. 
      This change should be seamless since the C implementation is numerically equivalent.
    * Radau5 C with SuperLU now requires compilation with '-D__RADAU5_WITH_SUPERLU'
    * Refactored Radau5 C implementation; moved to thirdparty/radau5/.
    * Replaced event localization for explicit problems with an equivalent C implementation.
    * Added missing counts on event indicator evaluations when setting problem data in solvers.
    * Fixed Radau5 bug where simulation would not stop when time limit has been exceeded.
    * Radau5(ODE) errors due to Exceptions in callbacks will now raise these Exceptions, rather than Radau5Error.
    * KeyboardInterrupt in Radau5ODE solout callback (e.g., event_indicator) will properly terminate simulation.
    * Added warning if a solver does not support setting time limits.
    * Added so that solver statistics is stored in case of early abort for CVode and Radau5.
    * Minor refactoring of Radau5 C implementation:
      Now uses actual return value to return error codes. idid is return of solout callback.
    * Separate fortran from C flags in setup.py
    * Changed "numpy.bool" to equivalent "bool" due to DeprecationWarnings for numpy >= 1.20.

--- Assimulo-3.3 ---
    * Sundials 6.x port
    * Radau5ODE and Radau5DAE now correctly terminate from unrecoverable errors in right-hand side and Jacobian calls.
    * Added support for sparse LU decompositons (via SuperLU) in Radau5 solver. 
      This is only available with the C version of the Radau5 solver.
    * Changed attribute name Radau5ODE.solver to Radau5ODE.implementation.
      Removed support for Radau5DAE.implementation = 'c'
    * Removed unused functionality from Radau5 C implementation: DAEs, 2nd order eq, banded matrix structures.

--- Assimulo-3.2.9 ---
   * Bugfix for C version of Radau5 solver

--- Assimulo-3.2.8 ---
    * Sundials 5.x port
    * Improved asserts in tests
    * Added C version of the Radau5 solver

--- Assimulo-3.2.7 ---
    * Resolved deprecation warnings visible when using numpy 1.20, related to deprecation of the alias numpy.float.
    * Resolved deprecation warnings visible when creating an ndarray from ragged nested sequences.

--- Assimulo-3.2.6 --
    * Fixed some tests that were not compatible with Python 3.

--- Assimulo-3.2.5 --
    * Added so that additional Fortran compile flags can be provided
      during setup.
    * Fixed minor issue that caused exceptions during setup when
      Sundials was not found.
      
--- Assimulo-3.2.4 --
    * Fixed issue with LSODAR (when being extended as a Cython problem)

--- Assimulo-3.2.3 --
    * Improved testing and made sure that event localization works
      even when a user reuses vectors in the problem class.

--- Assimulo-3.2.2 --
    * Fixed a regression with Radau5ODE from 3.2.1.

--- Assimulo-3.2.1 --
    * Improved detection of BLAS
    * Minor performance improvements

--- Assimulo-3.2 ---
    * Added support for Cython user extensions
    
--- Assimulo-3.1 ---
    * Added support for building Assimulo with Intel and MSVC on Windows (ticket: 434)
    * Added support for Sundials 4.1 (ticket:431)
    * Fixed so that *.pyf files are included in the source distro (ticket:430)
    * Moved ImportError output to stderr instead of stdout (ticket:429)
    * Disabled test dasp3_basic because problem appears to numerically unstable.

--- Assimulo-3.0 ---
    * Changed so that setuptools is used (support creating wheels) 
      (ticket:426)
    * Fixed so that sparse return type can be used from the jacobian
      method (ticket:423)
    * Delayed import of matplotlib
    * Fixed memory leaks in CVode and IDA (ticket:424)
    * Removed version check for numpy (it was only a problem for old 
      numpy version < 1.6.1 which we do no longer guard against) ticket:409)
    * Added license and changelog to the install folder (ticket:410)
    * Deprecated the setup option "sundials-with-superlu". SuperLU support
      using Sundials is not automatically checked (ticket:414)
    * Added support for Sundials 3.1 (ticket:418)
    * Renamed the option stablimit to stablimit (ticket:417)

--- Assimulo-2.9 ---
    * Added option to specify to use the 2-norm in CVode (ticket:401)
    * Added option to set max nonlinear iteration in CVode (ticket:400)
    * Renamed hmax to maxh in LSODAR for consistency (ticket:399)
    * Fixed version checking in setup for numpy (ticket:394)
    * Fixed bug with ncp list and backward integration (ticket:393)
    * Added method to retrieve current order in IDA (ticket:395)

--- Assimulo-2.8 ---
    * Added support for Sundials 2.6 (ticket:382)
    * Added support for sparse Jacobians (together with Sundials) 
      (ticket:383)
    * Added warning about chattering (ticket:387)
    * Added run-time status update (ticket:181)
    * Added option to set max-conv failures (ticket:386)
    * Removed warningen about comparison to None (ticket:381)
    * Update requirement for Cython to 0.18 (ticket:384)

--- Assimulo-2.7 ---
    * Added Python 3 support (ticket:296)
    * Fixed crash with atol as a matrix (ticket:351)
    * Added option for stability detection (ticket:355)
    * Fixed problem with event tolerance (ticket:367)

--- Assimulo-2.6 ---
    * Added version as an attribute (ticket:264)
    * Added more information is Sunials was found or not during install
      (ticket:265)
    * Fixed problem with storing event points (ticket:297)
    * Fixed wrong number of F-Evals in statistics in Radau/Rodas 
      (ticket:331)
    * Fixed problem with event detection in Euler (ticket:332)
    * Improved performance when using LSODAR (ticket:328)

--- Assimulo-2.5 ---
    * Added support for retrieving the last step in CVode (ticket:298)
    * Added support for retrieving the actual step in CVode (ticket:298)
    * Updated the documentation on the examples (ticket:316, ticket:315)
    * Added the name to the problem constructor (ticket:321)
    * Added option for timing a step (ticket:325)
    * Added an option to specify an upper bound on the integration time
      (ticket:289)
    * Added an option for a user specified J*v in the IDA case 
      (ticket:284)
    * Various bug fixes.

--- Assimulo-2.4 ---
    * Added support for simulating backward in time (ticket:267)
    * Added support event detection for Radau, Dopri, Explicit/Implicit
      Euler, Rodas, CVode (ticket:272)
    * Added the solver ODASSL.
    * Added the solver DASP3 (ticket:257)
    * Added basic Implicit Euler method (ticket:249)
    * Various bug fixes.

--- Assimulo-2.3 ---
    * Changed license to LGPL from GPL (ticket:261)
    * Fixed re_init problem with scalars (ticket:248)
    * Added a timer for measuring elapsed time of a step (ticket:260)
    * Added options to CVode to get order, weights and errors 
      (ticket:258, ticket:259)
    * Fixed problem with wrong dimensions when getting sensitivities in 
      CVode (ticket:255)
    * Added parameters when using Jac*Vec in CVode (ticket:250)
    * Added automatically creation of res function for explicit problems 
      (ticket:195)
    * Removed catching of exceptions in Explicit Euler (ticket:252)

--- Assimulo-2.2 ---
    * Added the solver LSODAR from ODEPACK (ticket:219)
    * Added number of state events to the statistics (ticket:224)
    * Fixed bug when storing result points together with events 
      (ticket:222)
    * Bug fixes.
    
--- Assimulo-2.1.1 ---
    * Fixed problem with binary distribution on Windows (ticket:213)

--- Assimulo-2.1 ---
    * Added support for passing in parameters when using Jacobians.
      (ticket:210)
    * Added warning when the solver does not support state events.
      (ticket:209)
    * Added RODAS by Hairer (ticket:207)
    * Added RADAU5 by Hairer (ticket:205)
    * Added DOPRI5 by Hairer (ticket:206)
    * Renamed the Python version of Radau with the prefix underscore.
      Radau5ODE -> _Radau5ODE, Radau5DAE -> _Radau5DAE

--- Assimulo-2.0 ---
    * Minor bug fixes in the setup script (ticket:191).
    * Fixed bug in type checking of switches (ticket:201). 

--- Assimulo-2.0b1 ---
    * Replaced setup_source.py and setup_binary.py with a single setup.py.
    * Base code migrated to Cython. 
    * Results are now returned from the simulate method.
    * Options and statistics are now stored in dictionaries.
    * Results are stored in variables appended with _sol. For example:
      y -> y_sol. (Also note that the result is now returned from simulate)
    * The current time and states (state derivative) have changed name from
      t_cur, y_cur, yd_cur -> t, y, yd.
    * Method in IDA make_consistency have been renamed to make_consistent.
    * Added a method get_support which returns a dictionary with 
      information about what the current solver supports.
    * Change name of the function in Explicit_Problem, f -> rhs
    * Change name of the function in Implicit_Problem, f -> res
    * Multiple name changes. (To be specified)
    * Improved the documentation
    * Speed improvements in the Sundials wrapper.
    * Fixed a couple of memory leaks in the Sundials wrapper.
    * Added support for specifying a list of output points.

--- Assimulo-1.4b3 ---
    * Fixed bug with t0 != 0 when using time events and step events
      (ticket:173)
    * Added support for specifying initial conditions for sensitivity
      variables (ticket:105)
    * Allowed pbar to specified in the problem (ticket:172)

--- Assimulo-1.4b2 ---
    * Fixed statistics for SPGMR (ticket:162).
    * Fixed bug when using fixed point iteration and jacobian related 
      calls (ticket:152).
    * Added options to terminate a simulation from handle_event via an 
      exception (ticket:163).
    * Fixed problem with atol and integers (ticket:161).

--- Assimulo-1.4b1 ---
    * Added option to use SPGMR in CVode (ticket:140).
    * Added new attributes in CVode, maxkrylov, pretype, linearsolver (ticket:140).
    * Added option to use a new method in Explicit_Problem, jacv (Jacobian*Vector)
      (ticket:144).
    * Fixed a bug with the completed simulation flag (ticket:133).
    * Fixed a bug when y0 is provided to Radau5 in the problem class (ticket:134).
    * Added an exception when the number of equations are zero (ticket:136).
    * Fixed a bug in the calling sequence of an event (ticket:138).
    * Added option to specify test attributes on tests (ticket:154).
    * Fixed various documentation inconsistencies.
    * Added a Kinsol wrapper (ticket:99)
    * Added a regularization technique (ticket:135).
    * Added SuperLU as a linear solver in Kinsol (ticket:153).
    * Fixed various bug related to Kinsol.
    * Updated the setup script to allow for specifying paths to SuperLU 
      and Blas (ticket:148).

--- Assimulo-1.3b1 ---
    * Improved the tolerance handling in RungeKutta34.
    * Improved information output from all the solvers.
    * Implemented basic support for calculating sensitivities using 
      IDAS.
    * Fixed a bug with the discretization method reseting the maximum
      order in CVode.
    * Minor bug fix in implicit Radau interpolate.
    * Changed the default value of pbar in CVodes and IDAs to the 
      absolute values of the parameters.

--- Assimulo-1.2b1 ---
    * Implemented basic support for calculating sensitivities using 
      CVodes.
    * Changed from using CVode to CVodes.
    * Added 'echo' methods used for viewing the current solver settings.
    * Fixed a bug with the reset method not resetting the statistics.
    * Fixed a bug which was exposed when overwriting the switches.
    * Added a custom error method in CVode and IDA.
    * Fixed a segmentation fault discovered on Mac when IDAS was used.
    * Renamed the test modules to lower-case.
    * Renamed the setup script to setup_from_binary (used when a
      pre-compiled binary is distributed)

--- Assimulo-1.1b1 ---
    * Fixed a bug with re-init resulting in resetting the options.
    * Moved the result handling to the problem class.
    * Renamed the event function to state_events.
    * Improved the information displayed after a simulation (mainly for 
      IDA and CVode).
    * Added support for step events (completed_step).
    * Added support for time events.
    * Implemented basic support for calculating sensitivities using 
      IDAS.
    * Renamed the modules to correspond to Python standard (all 
      lowercase). Classes starts with a capital letter.
    * Implemented Radau5 for both explicit and implicit problems.
    * Wrapped an interpolate method from Sundials (IDAGetDky, CVodeGetDky)
    * Changed from using IDA to IDAS
    * Changed assimulo.problem.Problem_Name to problem_name.
    * Changed assimulo.ODE.problemname to problem_name.
    * Fixed a bug when printing information when used FixedPoint.
    * Changed algvar to be more type independent.
    * Added **kwargs to the plotting functionality.

--- Assimulo-1.0b2 ---
    * Added an option to mask which variables that is to be plotted.
    * Added a .simulate function for use when simulating instead of
      __call__. Although __call__ can still be used.
    * Added a plotting functionality for plotting the step-size used
      together with the order used when the simulation have been
      run with one-step mode in either CVode or IDA.
    * Added so that when using IDA or CVode in one-step mode, the 
      current order and the last order are stored.
    * Added option to specify initial step-size in CVode.
    * Added support to switch between using the user defined Jacobian
      in CVode or not.
    * Added support to switch between using the user defined Jacobian
      in IDA or not.
    * Added support for user-defined Jacobians when using CVode.
    * Added support for user-defined Jacobians when using IDA.

--- Assimulo-1.0b1 ---
    * The rough first version.
