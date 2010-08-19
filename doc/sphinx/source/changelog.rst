
==========
Changelog
==========

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
    * Added support to switch between using the user defined jacobian
      in CVode or not.
    * Added support to switch between using the user defined jacobian
      in CVode or not.
    * Added support for user-defined Jacobians when using CVode.
    * Added support for user-defined Jacobians when using IDA.

--- Assimulo-1.0b1 ---
    * The rough first version.

