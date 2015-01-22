

###############
Tutorial
###############

Introduction
===============

This tutorial is intended to give a short introduction to the Assimulo package for solving both explicit and implicit ordinary differential equations. The tutorial focuses on the solvers IDA and CVode which are part of the SUNDIALS package written in C. In Assimulo these solvers have been lifted to Python to provide an easy interface and an easy platform for experimenting.

The tutorial serves also a short introduction for students taking the course `FMNN05 <http://www.maths.lth.se/na/courses/FMNN05/>`_ at `Lund University <http://www.lu.se/>`_ , Lund, Sweden. 
   


Chapters
==========

.. toctree::
   :maxdepth: 2
   
   tutorial_philo
   tutorial_exp
   tutorial_imp
   tutorial_disc
   tutorial_sundials

Additional information
========================

The Assimulo package comes with a number of examples, showing how to use the different solvers on different types of problem. These examples can be found :doc:`here <examples>`.


.. note::

    If there are any questions about a method or class, it is recommended to have a look at the docstrings.  For example, when using IPython, the docstrings are displayed  by entering the name of the method/class followed by 
    a question mark(?). Example, ::
    
        CVode.atol?

