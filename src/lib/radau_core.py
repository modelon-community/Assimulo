#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import pylab as P

class Radau_Exception(Exception):
    pass
    

class Radau_Common(object):
    """
    The common attributes for the Radau solvers.
    """
    def _get_h(self):
        """
        Sets the stepsize.
        """
        return self.__h
    
    def _set_h(self, h):
        """
        Sets the stepsize.
        """
        self.__h = h
        
    h = property(fget=_get_h,fset=_set_h)
    
    def print_statistics(self):
        """
        Prints the run-time statistics for the problem.
        """
        print 'Final Run Statistics: %s \n' % self.problem_name
        print 'Number of Steps                          = %s'%(self._nsteps)
        print 'Number of Function Evaluations           = %s'%(self._nfcn)
        print 'Number of Jacobian Evaluations           = %s'%(self._njac)
        print 'Number of F-Eval During Jac-Eval         = %s'%(self._njacfcn)
        print 'Number of Error Test Failures            = %s'%(self._errfail)
        print 'Number of Nonlinear Iterations           = %s'%(self._nniter)
        print 'Number of Nonlinear Convergence Failures = %s'%(self._nniterfail)
        print 'Number of LU decompositions              = %s'%(self._nlu)
    
    def plot_stepsize(self):
        """
        Plots the step-size.
        """
        P.semilogy(N.diff(self.t),drawstyle='steps-post')
        P.title(self.problem_name)
        P.ylabel('Step length')
        P.xlabel('Number of steps')
        P.show()
    
    def _set_newt(self, newt):
        """
        Maximal number of Newton iterations.
        
            Parameters::
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        try:
            newt = int(newt)
        except (ValueError, TypeError):
            raise Radau_Exception('The newt must be an integer or float.')
            
        self.__newt = newt
		
    def _get_newt(self):
        """
        Maximal number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        return self.__newt
		
    newt = property(_get_newt,_set_newt)
    
    def _set_fnewt(self, fnewt):
        """
        Stopping criterion for Newton's method, usually chosen <1.
        Smaller values of fnewt make the code slower, but safer.
        
            Parameters::
            
                fnewt
                        - Default min(0.03,rtol**0.5)
                        
                        - Should be a float.
                        
                            Example:
                                fnewt = 0.05
        """
        try:
            fnewt = float(fnewt)
        except (ValueError, TypeError):
            raise Radau_Exception('The fnewt must be an integer or float.')
            
        self.__fnewt = fnewt
        
    def _get_fnewt(self):
        """
        Stopping criterion for Newton's method, usually chosen <1.
        Smaller values of fnewt make the code slower, but safer.
        
            Parameters::
                fnewt
                        - Default min(0.03,rtol**0.5)
                        
                        - Should be a float.
                        
                            Example:
                                fnewt = 0.05
        """
        return self.__fnewt
        
    fnewt = property(_get_fnewt,_set_fnewt)
    
    def _set_safe(self, safe):
        """
        The safety factor in the step-size prediction.
        
            Parameters::
            
                safe
                        - Default '0.9'
                        
                        - Should be a float.
                        
                            Example:
                                safe = 0.8
        """
        try:
            safe = float(safe)
        except (ValueError, TypeError):
            raise Radau_Exception('The safe must be an integer or float.')
            
        self.__safe = safe
        
    def _get_safe(self):
        """
        The safety factor in the step-size prediction.
        
            Parameters::
            
                safe
                        - Default '0.9'
                        
                        - Should be float.
                        
                            Example:
                                safe = 0.8
        """
        return self.__safe
        
    safe = property(_get_safe, _set_safe)
    
    def _set_thet(self, thet):
        """
        Value for determine if the Jacobian is to be recomputed or not.
        Increasing thet makes the code compute new Jacobians more seldom.
        Negative thet forces the code to compute the Jacobian after every accepted step.
        
            Parameters::
            
                thet
                        - Default '0.003'
                        
                        - Should be float.
                        
                            Example:
                                thet = 0.01
        """
        try:
            thet = float(thet)
        except (ValueError, TypeError):
            raise Radau_Exception('The thet must be an integer or float.')
            
        self.__thet = thet
        
    def _get_thet(self):
        """
        Value for determine if the Jacobian is to be recomputed or not.
        Increasing thet makes the code compute new Jacobians more seldom.
        Negative thet forces the code to compute the Jacobian after every accepted step.
        
            Parameters::
            
                thet
                        - Default '0.003'
                        
                        - Should be float.
                        
                            Example:
                                thet = 0.01
        """
        return self.__thet
        
    thet = property(_get_thet, _set_thet)
    
    def _set_max_h(self,max_h):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default final time - current time.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        if not isinstance(max_h,float):
            raise Radau_Exception('Maximal stepsize must be a (scalar) float.')
        if max_h < 0.0:
            raise Radau_Exception('Maximal stepsize must be a positive (scalar) float.')
        self.__max_h=max_h
    
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default final time - current time.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.__max_h
        
    maxh=property(_get_max_h,_set_max_h)    
    
    def _set_initial_step(self, initstep):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Radau_Exception('The initial step must be an integer or float.')
        
        self.__initstep = initstep
    
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.__initstep
        
    initstep = property(_get_initial_step,_set_initial_step)
    
    
    def _set_quot1(self, quot1):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
            
                quot1
                        - Default 1.0
                        
                        - Should be a float.
                        
                            Example:
                                quot1 = 0.9
        """
        try:
            quot1 = float(quot1)
        except (ValueError, TypeError):
            raise Radau_Exception('The quot1 must be an integer or float.')
        
        self.__quot1 = quot1
    
    def _get_quot1(self):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
            
                quot1
                        - Default 1.0
                        
                        - Should be a float.
                        
                            Example:
                                quot1 = 0.9
        """
        return self.__quot1
        
    quot1 = property(_get_quot1, _set_quot1)
    
    def _set_quot2(self, quot2):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
        
                quot2
                        - Default 1.2
                        
                        - Should be a float.
                        
                            Example:
                                quot2 = 1.3
        """
        try:
            quot2 = float(quot2)
        except (ValueError, TypeError):
            raise Radau_Exception('The quot2 must be an integer or float.')
            
        self.__quot2 = quot2
    
    def _get_quot2(self):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
            
                quot2
                        - Default 1.2
                        
                        - Should be a float.
                        
                            Example:
                                quot2 = 1.2
        """
        return self.__quot2
        
    quot2 = property(_get_quot2, _set_quot2)
    
    
    def _set_fac1(self, fac1):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac1
                        - Default 0.2
                        
                        - Should be a float.
                        
                            Example:
                                fac1 = 0.1
        """
        try:
            fac1 = float(fac1)
        except (ValueError, TypeError):
            raise Radau_Exception('The fac1 must be an integer or float.')
            
        self.__fac1 = fac1
    def _get_fac1(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac1
                        - Default 0.2
                        
                        - Should be a float.
                        
                            Example:
                                fac1 = 0.1
        """
        return self.__fac1
    fac1 = property(_get_fac1, _set_fac1)
    
    def _set_fac2(self, fac2):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac2
                        - Default 8.0
                        
                        - Should be a float.
                        
                            Example:
                                fac2 = 10.0
        """
        try:
            fac2 = float(fac2)
        except (ValueError, TypeError):
            raise Radau_Exception('The fac2 must be an integer or float.')
            
        self.__fac2 = fac2
    def _get_fac2(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac2
                        - Default 8.0
                        
                        - Should be a float.
                        
                            Example:
                                fac2 = 10.0
        """
        return self.__fac2
    fac2 = property(_get_fac2, _set_fac2)
    
    def _set_usejac(self, jac):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        self.__usejac = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.__usejac
    
    usejac = property(_get_usejac,_set_usejac)
