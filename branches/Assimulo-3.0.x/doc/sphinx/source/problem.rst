 
Problem Class
=============================

.. note::
   
   Here we describe the standard problem classes. 
   From these users might derive problem classes for their particular needs.
   This is mainly the case when the right hand side function (rhs / res) 
   has to be provided with additional information passed as instance attributes.
   For an example see :download:`cvode_with_disc.py <../../../examples/cvode_with_disc.py>`.

.. autoclass:: assimulo.problem.cProblem
   :members:
 

Explicit Problem
-----------------------
 
.. autoclass:: assimulo.problem.Explicit_Problem
    :members: 


  
Implicit Problem
-----------------------
 
.. autoclass:: assimulo.problem.Implicit_Problem
    :members:

Mechanical Problem 
-----------------------

.. autoclass:: assimulo.special_systems.cMechanical_System
    :members:
