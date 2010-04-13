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
        """Sets the absolute tolerance(s)."""
        if not isinstance(atol,float):
            if not isinstance(atol, list) and not isinstance(atol, N.ndarray):
                raise Sundials_Exception('Absolute tolerance must be a float or a float vector.')
            if len(atol) != self.Integrator.dim:
                raise Sundials_Exception('Absolute tolerance must be a float vector of same dimension as the problem or a scalar.')
            for tol in atol:
                if not isinstance(tol,float):
                    raise Sundials_Exception('Absolute tolerance must be a float vector or a float scalar.')
                if tol <= 0.0:
                    raise Sundials_Exception('Absolute tolerance must be a positive (scalar) float or a positive (vector) float.')
            self.Integrator.abstol_ar=N.array(atol)
            self.__atol=atol
        else:
            if atol <= 0.0:
                raise Sundials_Exception('Absolute tolerance must be a positive (scalar) float or a positive (vector) float.')
            self.Integrator.abstol_ar=atol*N.ones(self.Integrator.dim)
            self.__atol=atol
    def _get_atol(self):
        """Returns the absolute tolerance(s)."""
        return self.__atol
    atol=property(_get_atol,_set_atol,doc='Absolute tolerance\n can be a float vector or a float scalar.')
    
    def _set_rtol(self,rtol):
        """Sets the relative tolerance."""
        if not isinstance(rtol,float):
            raise Sundials_Exception('Relative tolerance must be a (scalar) float.')
        if rtol <= 0.0:
            raise Sundials_Exception('Relative tolerance must be a positive (scalar) float.')
        self.Integrator.reltol=rtol
        self.__rtol=rtol
    def _get_rtol(self):
        """Returns the relative tolerance."""
        return self.__rtol
    rtol=property(_get_rtol,_set_rtol,doc='Absolute tolerance\n can only be a float scalar.')
    
    def _set_max_h(self,max_h):
        """Sets the maximal stepsize with the default value of infinity."""
        if not isinstance(max_h,float):
            raise Sundials_Exception('Maximal stepsize must be a (scalar) float.')
        if max_h < 0.0:
            raise Sundials_Exception('Maximal stepsize must be a positive (scalar) float.')
        self.Integrator.max_h=max_h
        self.__max_h=max_h
    def _get_max_h(self):
        """Returns the maximal stepsize."""
        return self.__max_h
    maxh=property(_get_max_h,_set_max_h,doc='Maximal stepsize')    
    
    def plot_stepsize_order(self):
        """
        Plots the step-size used throughout the integration together with
        the order used. Only able when using one-step mode.
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
        P.plot(self.Integrator.detailed_info['qcurrent'])
        P.ylabel('Order')
        P.xlabel('Steps')
        P.legend(['Last order','Current order'], loc=4)
        P.show()
        
    
    @property
    def stats(self):
        """Attribute that returns the run-time statistics from the Integrator."""
        return self.Integrator.stats
    
    @property
    def disc_info(self):
        """Attribute that returns information about an occured event."""
        return [self.Integrator.event_time, self.Integrator.event_info]
