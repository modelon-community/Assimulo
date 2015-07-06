.. highlight:: rest

=============
Installation
=============

Dependencies:
    
- `Python-2.6 / 2.7 / 3.3 / 3.4 <http://www.python.org/>`_
- `Numpy <http://www.scipy.org/Download/>`_ (>1.6.1 for the solver DASP3)
- `Scipy <http://www.scipy.org/Download/>`_
- `Pylab <http://matplotlib.sourceforge.net/>`_

Additional dependencies for compiling from source:

- `Python-2.6 / 2.7 / 3.3 / 3.4 <http://www.python.org/>`_ with headers (python-dev package for Ubuntu)
- `Sundials-2.5/2.6 <http://computation.llnl.gov/casc/sundials/main.html>`_ (for 64bits machines, install Sundials using -fPIC)
- `Cython 0.18 <http://www.cython.org/>`_
- C compiler
- Fortran compiler
- BLAS (only needed for the solver GLIMDA)
- LAPACK (only needed for the solver GLIMDA)


Assimulo is found on the :doc:`download` page.


Installation flags
====================

When installing Assimulo from source there are a number of available flags that can be specified in order to point to dependencies and which should be provided after the install command::

    python setup.py install ...
    
The available flags are::

    - "--sundials-home=..." - Point to an Sundials installation
    - "--blas-home=..." - Point to an BLAS installation
    - "--lapack-home=..." - Point to an LAPACK installation

Example::

    python setup.py install --sundials-home=/home/chria/Sundials --blas-home=/home/chria/Blas

Ubuntu
==========

Once all the dependencies are satisfied an installation is done by::

    python setup.py install 
    
After a successful installation, the package will be located in Pythons dist-packages folder. Note, in case of 64bit systems
the section troubleshooting see :ref:`instTrouble` should be consulted before installation.

.. note::

    If Sundials has been installed on a different location then the default, use the sundials flag::
    
        --sundials-home=/path/to/sundials

.. note::

    To test Assimulo, go into the tests folder and type::
    
        nosetests
        
    Which requires python-nose.

Windows
==========

For installing on Windows it is recommended to download the binary installers from the :doc:`download` page which includes all the solvers available. The below instructions are for installing Assimulo on Windows from source.

.. note::

    Assimulo is also dependent on the Windows redistributable package for Windows.

Installing Sundials on Windows can be a bit tricky but here is a link for installing Sundials using cmake together with Visual Studio, http://sundials.wikidot.com/installation-cmake-vs . However I would recommend using Mingw instead of Visual Studio, here is link for installing Mingw on Windows and setting up Cython to use it, http://docs.cython.org/src/tutorial/appendix.html . If you would like to use Mingw instead of Visual Studio, just follow the above guide for installing Sundials until the step where Visual Studio is used. Instead of following those instructions, browse to Sundials build catalogue using the command prompt and type::

    make
    make install

Once Sundials and the rest of the packages are installed just install Assimulo by browsing to the folder in the command prompt and type::

    python setup.py install --sundials-home=/path/to/sundials
    
After a successful installation, the package will be located in pythons dist-packages folder.

.. note::

    To test Assimulo, go into the tests folder and type::
    
        nosetests
        
    Which requires python-nose.


.. _instTrouble:

Troubleshooting
================

Ubuntu 64bits
---------------
There have been some problems installing Assimulo on Ubuntu 64bits machines when Sundials has been installed with the default options. The problem generates the following error printout::

    /usr/bin/ld: /home/chria/sundialscode/lib/libsundials_cvodes.a(cvodes.o): relocation R_X86_64_32
    against `.rodata.str1.1' can not be used when making a shared object; recompile with -fPIC
    > /home/chria/sundialscode/lib/libsundials_cvodes.a: could not read symbols: Bad value
    > collect2: ld returned 1 exit status
    > error: command 'gcc' failed with exit status 1
    
To solve this problem, Sundials has to be installed with the flag 

    CFLAGS="-fPIC"
    
Consult the Sundials INSTALL_NOTES Sec. B.3 to see 
how this compiler flag has to be specified.
