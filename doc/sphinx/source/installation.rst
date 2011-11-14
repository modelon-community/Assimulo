.. highlight:: rest

=============
Installation
=============

Dependencies:
    
- `Python-2.6 / 2.7 <http://www.python.org/>`_ with headers (python-dev package for Ubuntu)
- `Cython <http://www.cython.org/>`_
- `Numpy <http://www.scipy.org/Download/>`_
- `Scipy <http://www.scipy.org/Download/>`_
- `Pylab <http://matplotlib.sourceforge.net/>`_
- `Sundials-2.4 <http://computation.llnl.gov/casc/sundials/main.html>`_


Assimulo is found on the :doc:`download` page.

Ubuntu
==========

Once all the dependencies are satisfied an installation is done by::

    python setup.py install 
    
After a successful installation, the package will be located in Pythons dist-packages folder. For troubleshooting see :ref:`instTrouble`.

.. note::

    If Sundials has been installed on a different location then the default, the argument::
    
        --sundials-home=/path/to/sundials
        
    should point to where Sundials is installed. Example, ::
    
        python setup.py install -sundials-home=/home/chria/Sundials

.. note::

    To test Assimulo, go into the tests folder and type::
    
        nosetests
        
    Which requires python-nose.

Windows
==========

Installing Sundials on Windows can be a bit tricky but here is a link for installing Sundials using cmake together with Visual Studio, http://sundials.wikidot.com/installation-cmake-vs . However I would recommend using Mingw instead of Visual Studio, here is link for installing Mingw on Windows and setting up Cython to use it, http://docs.cython.org/src/tutorial/appendix.html . If you would like to use Mingw instead of Visual Studio, just follow the above guide for installing Sundials until the step where Visual Studio is used. Instead of following those instructions, browse to Sundials build catalogue using the command prompt and type::

    make
    make install

Once Sundials and the rest of the packages are installed just install Assimulo by browsing to the folder in the command prompt and type::

    python setup.py install --sundials-home=/path/to/sundials
    
After a successful installation, the package will be located in pythons dist-packages folder.

.. note::

    The argument::
    
        --sundials-home 
        
    should point to where sundials is installed.

.. note::

    To test Assimulo, go into the tests folder and type::
    
        nosetests
        
    Which requires python-nose.


.. _instTrouble:

Troubleshooting
================

Ubuntu 64bits
---------------
There have been some problems installing Assimulo on Ubuntu 64bits machines when Sundials have been installed with the default options. The problem generates the following error printout::

    /usr/bin/ld: /home/chria/sundialscode/lib/libsundials_cvodes.a(cvodes.o): relocation R_X86_64_32
    against `.rodata.str1.1' can not be used when making a shared object; recompile with -fPIC
    > /home/chria/sundialscode/lib/libsundials_cvodes.a: could not read symbols: Bad value
    > collect2: ld returned 1 exit status
    > error: command 'gcc' failed with exit status 1
    
To solve this problem, Sundials have to be installed with the CFLAGS="-fPIC".
