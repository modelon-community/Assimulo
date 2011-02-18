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

from lib import sundials_core
import numpy as N
import pylab as P

class Sundials_Exception(Exception):
    pass
    

class Sundials:
    
    def __init__(self, y0, integrator):
        """
        Defines and sets up the integrator.
        """
        try:
            if isinstance(y0, int) or isinstance(y0, float):
                y0 = [y0]
            dim = len([N.array(y0, dtype=float)][0])
        except ValueError:
            dim = 0
            
        if integrator == 'IDA':
            self.Integrator = sundials_core.IDA_wrap(dim) #Creates a IDA integrator
        if integrator == 'CVode':
            self.Integrator = sundials_core.CVode_wrap(dim) #Creates a CVode integrator
    
    def _set_atol(self,atol):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        
        See SUNDIALS IDA documentation 4.5.2 for more details.
        """
        try:
            atol_arr = N.array(atol, dtype=float)
            if (atol_arr <= 0.0).any():
                raise Sundials_Exception('Absolute tolerance must be a positive float or a float vector.')
        except (ValueError,TypeError):
            raise Sundials_Exception('Absolute tolerance must be a positive float or a float vector.')
        if atol_arr.size == 1:
            self.Integrator.atol = atol_arr*N.ones(self.Integrator.dim)
            self.__atol = float(atol)
        elif atol_arr.size == self.Integrator.dim:
            self.Integrator.atol = atol_arr
            self.__atol = [float(x) for x in atol]
        else:
            raise Sundials_Exception('Absolute tolerance must be a float vector of same dimension as the problem or a scalar.')
    
    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        
        See SUNDIALS IDA documentation 4.5.2 for more details.
        """
        return self.__atol
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
                                
        """
        try:
            rtol = float(rtol)
        except (ValueError, TypeError):
            raise Sundials_Exception('Relative tolerance must be a (scalar) float.')
        if rtol <= 0.0:
            raise Sundials_Exception('Relative tolerance must be a positive (scalar) float.')
        self.Integrator.rtol=rtol
        self.__rtol=rtol
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
                                
        """
        return self.__rtol
        
    rtol=property(_get_rtol,_set_rtol)
    
    def _set_max_h(self,max_h):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0', which indicates that the maximal
                          step-size is infinity.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        if not isinstance(max_h,float):
            raise Sundials_Exception('Maximal stepsize must be a (scalar) float.')
        if max_h < 0.0:
            raise Sundials_Exception('Maximal stepsize must be a positive (scalar) float.')
        self.__max_h=max_h
        self.Integrator.maxh = self.__max_h
    
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0', which indicates that the maximal
                          step-size is infinity.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.__max_h
        
    maxh=property(_get_max_h,_set_max_h)    
    
    def plot_stepsize_order(self):
        """
        Plots the step-size used through-out the integration together with
        the order used. Only able when using one-step mode, 
        
        (number of communication points 'ncp' = 0).
        """
        if not self.Integrator.detailed_info:
            raise Sundials_Exception('There is no information. Either the problem' \
            ' was not simulated using the one-step mode or the problem have not ' \
            'been simulated.')
        
        P.subplot(211)
        P.semilogy(N.diff(self.t),drawstyle='steps-post')
        P.title(self.problemname)
        P.ylabel('Step length')
        P.xlabel('Number of steps')
        P.subplot(212)
        P.plot(self.Integrator.detailed_info['qlast'])
        #P.plot(self.Integrator.detailed_info['qcurrent'])
        P.ylabel('Order')
        P.xlabel('Steps')
        #P.legend(['Last order','Current order'], loc=4)
        P.show()
        
    
    @property
    def stats(self):
        """
        Returns the run-time statistics.
        
            Returns::
            
                A list of statistics.
                
            Example::
            
                stats = IDA/CVode.stats
                
                stats[0] # Number of Steps                         
                stats[1] # Number of Function Evaluations         
                stats[2] # Number of Jacobian Evaluations         
                stats[3] # Number of F-Eval During Jac-Eval        
                stats[4] # Number of Root Evaluations              
                stats[5] # Number of Error Test Failures           
                stats[6] # Number of Nonlinear Iterations          
                stats[7] # Number of Nonlinear Convergence Failures
        """
        return self.Integrator.solver_stats
    
    @property
    def disc_info(self):
        """
        Attribute that returns information about an occured event.
        """
        return [self.Integrator.event_time, self.Integrator.event_info]
        
    def _set_dqrhomax(self, dqrhomax):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        try:
            dqrhomax = float(dqrhomax)
            if dqrhomax < 0.0:
                raise Sundials_Exception('DQrhomax must be a positive float.')
        except (TypeError, ValueError):
            raise Sundials_Exception('DQrhomax must be convertable to a float.')
            
        self.Integrator.DQrhomax = dqrhomax
        
    def _get_dqrhomax(self):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        return self.Integrator.DQrhomax
    
    dqrhomax = property(_get_dqrhomax, _set_dqrhomax)
    
    def _set_usesens(self, usesens):
        """
        Specifies if the sensitivity calculations should be used or turned off.
        
            Parameters::
            
                usesens
                            - A boolean type.
                            - Default True
                            
            Returns::
            
                The current value of usesens (boolean)
        
        See SUNDIALS documentation '(IDA/CVode)SensToggleOff' 
        """
        self.Integrator.sensToggleOff = not bool(usesens)
        
    def _get_usesens(self):
        """
        Specifies if the sensitivity calculations should be used or turned off.
        
            Parameters::
            
                usesens
                            - A boolean type.
                            - Default True
                            
            Returns::
            
                The current value of usesens (boolean)
        
        See SUNDIALS documentation '(IDA/CVode)SensToggleOff' 
        """
        return not self.Integrator.sensToggleOff
        
    usesens = property(_get_usesens, _set_usesens)
    
    def _set_suppress_sens(self, suppress_sens):
        """
        Specifies whether sensitivity variables are included in the error test
        or not. True means that the variables are suppressed and not included 
        in the error test.
        
            Parameters::
            
                suppress_sens
                        - A boolean
                        - Default True
                        
            Returns::
                
                The current value of suppress_sens (boolean)
                
        See SUNDIALS documentation '(IDA/CVode)SetSensErrCon'.
        
        NOTE:: 
        
            That this method does the opposite of (IDA/CVode)SetSensERRCon to have 
            the same meaning as suppress_alg.
        """
        self.Integrator.errconS = not bool(suppress_sens)
        
    def _get_suppress_sens(self):
        """
        Specifies whether sensitivity variables are included in the error test
        or not. True means that the variables are suppressed and not included 
        in the error test.
        
            Parameters::
            
                suppress_sens
                        - A boolean
                        - Default True
                        
            Returns::
                
                The current value of suppress_sens (boolean)
                
        See SUNDIALS documentation '(IDA/CVode)SetSensErrCon'.
        
        NOTE:: 
        
            That this method does the opposite of (IDA/CVode)SetSensERRCon to have 
            the same meaning as suppress_alg.
        """
        return not self.Integrator.errconS
    
    suppress_sens = property(_get_suppress_sens, _set_suppress_sens)
    
    def _set_max_nonlin(self, maxsensiter):
        """
        Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables per step. (>0)
        
            Parameters::
            
                maxsensiter 
                            - An integer
                            - Default 3
                            
            Returns::
            
                The current value of maxsensiter.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensMaxNonlinIters'.
        """
        try:
            maxsensiter = int(maxsensiter)
            if maxsensiter < 1:
                raise Sundials_Exception('maxsensiter must be greater than zero.')
        except (TypeError, ValueError):
            raise Sundials_Exception('maxsensiter must be convertable to an integer.')
        
        self.Integrator.maxcorS = maxsensiter
        
    def _get_max_nonlin(self):
        """
        Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables per step. (>0)
        
            Parameters::
            
                maxsensiter 
                            - An integer
                            - Default 3
                            
            Returns::
            
                The current value of maxsensiter.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensMaxNonlinIters'.
        """
        return self.Integrator.maxcorS
    
    maxsensiter = property(_get_max_nonlin, _set_max_nonlin)

    def _set_pbar(self, pbar):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default absolute values of the parameters.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS documentation '(IDA/CVode)SetSensParams'
        """
        if len(pbar) != self.problem_data['dimSens']:
            raise Sundials_Exception('pbar must be of equal length as the parameters.')
        
        self.Integrator.pbar = pbar
    
    def _get_pbar(self):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default absolute values of the parameters.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS documentation '(IDA/CVode)SetSensParams'
        """
        return self.Integrator.pbar
    
    pbar = property(_get_pbar, _set_pbar)


    def _set_sensitivity_method(self, ism):
        """
        Specifies the sensitivity solution method. Can be either,
        'SIMULTANEOUS' or 'STAGGERED'.
            
        Parameters::
            
            ism
                    - A string of either 'SIMULTANEOUS' or 'STAGGERED'
                    - Default 'STAGGERED'
                        
        Returns::
            
            The current value of sensmethod (string)
        
        See SUNDIALS documentation '(IDA/CVode)SensInit'
        """
        if not isinstance(ism, str):
            raise Sundials_Exception('sensmethod must be string.')
        
        if ism.upper() == 'SIMULTANEOUS':
            self.Integrator.ism = 1
        elif ism.upper() == 'STAGGERED':
            self.Integrator.ism = 2
        else:
            raise Sundials_Exception('sensmethod must be either "SIMULTANEOUS" or "STAGGERED".')
        
    def _get_sensitivity_method(self):
        """
        Specifies the sensitivity solution method. Can be either,
        'SIMULTANEOUS' or 'STAGGERED'.
            
        Parameters::
            
            ism
                    - A string of either 'SIMULTANEOUS' or 'STAGGERED'
                    - Default 'STAGGERED'
                        
        Returns::
            
            The current value of sensmethod (string)
        
        See SUNDIALS documentation '(IDA/CVode)SensInit'
        """
        if self.Integrator.ism == 1:
            return 'SIMULTANEOUS'
        elif self.Integrator.ism == 2:
            return 'STAGGERED'
        else:
            raise Sundials_Exception('Unknown value of sensitivity solution method.')
    
    sensmethod = property(_get_sensitivity_method, _set_sensitivity_method)
    
    def _set_dqtype(self, dqtype):
        """
        Specifies the difference quotient type in the sensitivity calculations.
        Can be either, 'CENTERED' or 'FORWARD'.
        
        Parameters::
            
            DQtype 
                    - A string of either 'CENTERED' or 'FORWARD'
                    - Default 'CENTERED'
                        
        Returns::
            
            The current value of DQtype.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        if not isinstance(dqtype, str):
            raise Sundials_Exception('DQtype must be string.')
        
        if dqtype.upper() == 'CENTERED':
            self.Integrator.DQtype = 1
        elif dqtype.upper() == 'FORWARD':
            self.Integrator.DQtype = 2
        else:
            raise Sundials_Exception('DQtype must be either "CENTERED" or "FORWARD".')
            
    def _get_dqtype(self):
        """
        Specifies the difference quotient type in the sensitivity calculations.
        Can be either, 'CENTERED' or 'FORWARD'.
        
        Parameters::
            
            DQtype 
                    - A string of either 'CENTERED' or 'FORWARD'
                    - Default 'CENTERED'
                        
        Returns::
            
            The current value of DQtype.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        if self.Integrator.DQtype == 1:
            return 'CENTERED'
        elif self.Integrator.DQtype == 2:
            return 'FORWARD'
        else:
            raise Sundials_Exception('Unknown value of DQtype.')
    
    dqtype = property(_get_dqtype, _set_dqtype)

    def interpolate_sensitivity(self, t, k, p=-1):
        """        
        Calls Sundials internal function (IDA/CVodes)GetSensDky1 that computes the interpolated 
        values of the k-th derivative of the sensitivity p for any value of t in the 
        last internal step taken by IDAS or CVodes.
        
            Parameters::
            
                t
                    - Must be within tn − hu ≤ t ≤ tn  where tn denotes the current
                      internal time reached, and hu is the last internal step size used successfully.
                      
                    - Must be a float.
                      
                k
                    - Must be non-negative and samller than the last internal order used.
                    
                    - Must be an integer.
                    
                p
                    - An integer to determine of which parameter the solution is to be computed.
                    
                    - Default value -1 indicates that all is to be calculated.
                    
                    - Must be an integer.
                    
            Returns::
            
                A numpy array of the calculated sensitivites.
        """
        try:
            t = float(t)
        except (TypeError, ValueError):
            raise Sundials_Exception('t must be convertable to a float.')
        try:
            k = int(k)
        except (TypeError, ValueError):
            raise Sundials_Exception('k must be convertable to an integer.')
        try:
            p = int(p)
        except (TypeError, ValueError):
            raise Sundials_Exception('p must be convertable to an integer.')
            
        return self.Integrator.interpolate_sensitivity(t,k,p)

    def get_sensitivity_statistics(self):
        """
        Returns the sensitivity statistics.
        
            Returns::
            
                A list of statistics.
                
            Example::
            
                stats = IDA.get_sensitivity_statistics()
                
                stats[0] #Number of Sensitivity Calculations
                stats[1] #Number of F-Evals Due to Finite Approximation
                stats[2] #Number of Local Error Test Failures
                stats[3] #Number of Linear Setups            
                stats[4] #Number of Nonlinear iterations     
                stats[5] #Number of Nonlinear Convergance Failures
        """
        return self.Integrator.solver_sens_stats
