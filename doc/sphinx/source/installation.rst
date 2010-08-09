.. highlight:: rest

=============
Installation
=============

Dependencies:
    
- `Python-2.6 (with headers) <http://www.python.org/>`_
- `Cython <http://www.cython.org/>`_
- `Numpy <http://numpy.scipy.org/>`_
- Pylab
- `Sundials-2.4 <http://computation.llnl.gov/casc/sundials/main.html>`_


Assimulo is found on the :doc:`download` page.

Ubuntu
==========

Once all the dependencies are satisfied an installation is done by::

    python setup_source.py install --sundials-home=/path/to/sundials
    
After a successful installation, the package will be located in Pythons dist-packages folder.

.. note::

    The argument::
    
        --sundials-home 
        
    should point to where Sundials is installed. If Sundials has been installed using the default path, this argument is not needed as the default path is /usr/local/. If it have been installed anywhere else, this path should point to its location.

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

    python setup_source.py install --sundials-home=/path/to/sundials
    
After a successful installation, the package will be located in pythons dist-packages folder.

.. note::

    The argument::
    
        --sundials-home 
        
    should point to where sundials is installed.

.. note::

    To test Assimulo, go into the tests folder and type::
    
        nosetests
        
    Which requires python-nose.

