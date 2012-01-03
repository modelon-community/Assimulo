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

import numpy as N

from assimulo.ode import *
from assimulo.explicit_ode import Explicit_ODE

from assimulo.exception import *

from assimulo.lib import dopri5

class Dopri5(Explicit_ODE):
    """
    Explicit Runge-Kutta method of order (4)5 with step-size control
    and continuous output. Based on the method by Dormand and Prince.
    
    Based on the FORTRAN code DOPRI5 by E.Hairer and G.Wanner, which can 
    be found here: http://www.unige.ch/~hairer/software.html
    
    Details about the implementation (FORTRAN) can be found in the book,::
    
        Solving Ordinary Differential Equations I,
        Nonstiff Problems
        
        Authors: E. Hairer, S. P. Norsett and G. Wanner
        Springer-Verlag, ISBN: 3-540-56670-8
    
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["safe"]     = 0.9 #Safety factor
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 10.0 #Parameters for step-size selection (upper bound)
        self.options["beta"]     = 0.04
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["inith"]    = 0.0
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["maxsteps"] = 100000
        
        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["errfail"]     = 0 #Number of step rejections
        self.statistics["nstepstotal"] = 0 #Number of total computed steps (may NOT be equal to nsteps+nerrfail)
        
        #Internal
        self._leny = len(self.y) #Dimension of the problem
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
        self._tlist = []
        self._ylist = []
        
    def _solout(self, nrsol, told, t, y, cont, lrc, irtrn):
        """
        This method is called after every successful step taken by Radau5
        """
        if self._opts["output_list"] == None:
            self._tlist.append(t)
            self._ylist.append(y.copy())
        else:
            output_list = self._opts["output_list"]
            output_index = self._opts["output_index"]
            try:
                while output_list[output_index] <= t:
                    self._tlist.append(output_list[output_index])
                    
                    yval = N.empty(self._leny)
                    for i in range(self._leny):
                        yval[i] = dopri5.contd5(i+1,output_list[output_index], cont, lrc)
                        
                    self._ylist.append(yval)

                    output_index = output_index+1
            except IndexError:
                pass
            self._opts["output_index"] = output_index
    
    def integrate(self, t, y, tf, opts):
        ITOL  = 1 #Both atol and rtol are vectors
        IOUT  = 2 #Dense out in solout
        WORK  = N.array([0.0]*(8*self.problem_info["dim"]+5*self.problem_info["dim"]+21))
        IWORK = N.array([0]*(self.problem_info["dim"]+21))
        
        #Setting work options
        WORK[1] = self.safe
        WORK[2] = self.fac1
        WORK[3] = self.fac2
        WORK[4] = self.beta
        WORK[5] = self.maxh
        WORK[6] = self.inith
        
        #Setting iwork options
        IWORK[0] = self.maxsteps
        IWORK[4] = self.problem_info["dim"] 
        
        #Store the opts
        self._opts = opts
        
        t, y, iwork, flag = dopri5.dopri5(self.problem.rhs, t, y.copy(), tf, self.rtol*N.ones(self.problem_info["dim"]), self.atol, ITOL, self._solout, IOUT, WORK, IWORK)
        
        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("Dopri5 failed with flag %d"%flag)
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[18]
        self.statistics["nfcn"]        += iwork[16]
        self.statistics["nstepstotal"] += iwork[17]
        self.statistics["errfail"]     += iwork[19]
        
        return flag, self._tlist, self._ylist
        
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["errfail"]),       verbose)
        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : Dopri5 ',          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
        
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise Dopri5_Exception("atol must be of length one or same as the dimension of the problem.")

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
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        try:
            self.options["rtol"] = float(rtol)
        except (ValueError, TypeError):
            raise Dopri5_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise Dopri5_Exception('Relative tolerance must be a positive (scalar) float.')
    
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
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)
    
    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 10000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise Dopri5_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def _set_fac1(self, fac1):
        try:
            self.options["fac1"] = float(fac1)
        except (ValueError, TypeError):
            raise Dopri5_Exception('The fac1 must be an integer or float.')
            
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
        return self.options["fac1"]
        
    fac1 = property(_get_fac1, _set_fac1)
    
    def _set_fac2(self, fac2):
        try:
            self.options["fac2"] = float(fac2)
        except (ValueError, TypeError):
            raise Dopri5_Exception('The fac2 must be an integer or float.')
            
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
        return self.options["fac2"]
        
    fac2 = property(_get_fac2, _set_fac2)
    
    def _set_safe(self, safe):
        try:
            self.options["safe"] = float(safe)
        except (ValueError, TypeError):
            raise Dopri5_Exception('The safe must be an integer or float.')

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
        return self.options["safe"]
        
    safe = property(_get_safe, _set_safe)
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Dopri5_Exception('The initial step must be an integer or float.')
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                inith    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    inith = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise Dopri5_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise Dopri5_Exception('Maximal stepsize must be a positiv (scalar) float.')
        
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
        return self.options["maxh"]
        
    maxh=property(_get_max_h,_set_max_h)
    
    def _set_beta(self, beta):
        self.options["beta"] = beta
        
    def _get_beta(self):
        """
        Option for stabilized step-size control.
        
            Parameters::
            
                beta
                        - Default 0.04
                        
                        - Should be a float.
        """
        return self.options["beta"]
    
    beta = property(_get_beta, _set_beta)

class RungeKutta34(Explicit_ODE):
    """
    Adaptive Runge-Kutta of order four.
    
    Obs. Step rejection not implemented.
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.                       
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["atol"] = 1.0e-6
        self.options["rtol"] = 1.0e-6
        self.options["inith"] = 0.01
        self.options["maxsteps"] = 10000
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self.Y4 = N.array([0.0]*len(self.y0))
        self.Z3 = N.array([0.0]*len(self.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
        
        #Internal values
        # - Statistic values
        self.statistics["nsteps"] = 0 #Number of steps
        self.statistics["nfcn"] = 0 #Number of function evaluations
    
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
    
    def _set_initial_step(self, initstep):
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
        self.options["inith"] = initstep
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                inith    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    inith = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_atol(self,atol):

        try:
            atol_arr = N.array(atol, dtype=float)
            if (atol_arr <= 0.0).any():
                raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        except (ValueError,TypeError):
            raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        if atol_arr.size == 1:
            self.options["atol"] = float(atol)
        elif atol_arr.size == len(self.y):
            self.options["atol"] = [float(x) for x in atol]
        else:
            raise Explicit_ODE_Exception('Absolute tolerance must be a float vector of same dimension as the problem or a scalar.')

    def _get_atol(self):
        """
        Sets the absolute tolerance to be used in the integration.
        
            Parameters::
            
                atol    
                            - Default 1.0e-6.
                            
                            - Should be float or an array/list of len(y)
                            
                                Example:
                                    atol=1.e5
                                    atol=[1.e-5,1.e-4]
        """
        return self.options["atol"]
    
    atol = property(_get_atol,_set_atol)
    
    def _set_rtol(self, rtol):
        try:
            rtol = float(rtol)
        except (TypeError,ValueError):
            raise Explicit_ODE_Exception('Relative tolerance must be a float.')
        if rtol <= 0.0:
            raise Explicit_ODE_Exception('Relative tolerance must be a positive (scalar) float.')
        self.options["rtol"] = rtol
            
    def _get_rtol(self):
        """
        The relative tolerance to be used in the integration.
        
            Parameters::
            
                rtol    
                            - Default 1.0e-6
                            
                            - Should be a float.
        """
        return self.options["rtol"]
    
    rtol = property(_get_rtol, _set_rtol)
    
    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 10000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise Explicit_ODE_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def step(self, t, y, tf, opts):
        initialize = opts["initialize"]
        
        if initialize:
            self.solver_iterator = self._iter(t,y,tf)

        return self.solver_iterator.next()
    
    def integrate(self, t, y, tf, opts):
        """
        Integrates (t,y) values until t > tf
        """
        [flags, tlist, ylist] = zip(*list(self._iter(t, y, tf)))
        
        return flags[-1], tlist, ylist
    
    def _iter(self,t,y,tf):
        maxsteps = self.options["maxsteps"]
        h = self.options["inith"]
        h = min(h, N.abs(tf-t))
        
        for i in range(maxsteps):
            if t+h < tf:
                t, y, error = self._step(t, y, h)
                self.statistics["nsteps"] += 1
                yield ID_PY_OK, t,y
                h=self.adjust_stepsize(h,error)
                h=min(h, N.abs(tf-t))
            else:
                break
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
            
        t, y, error = self._step(t, y, h)
        self.statistics["nsteps"] += 1
        yield ID_PY_COMPLETE, t, y

    def adjust_stepsize(self, h, error):
        """
        Adjusts the stepsize.
        """
        fac=min((1.0/error)**(1.0/4.0),2.)
        h *= fac
        
        return h
        
    def _step(self, t, y, h):
        """
        This calculates the next step in the integration.
        """
        self.statistics["nfcn"] += 5
        
        scaling = N.array(abs(y)*self.rtol + self.atol) # to normalize the error 
        f = self.f
        
        f(self.Y1, t, y)
        f(self.Y2, t + h/2., y + h*self.Y1/2.)
        f(self.Y3, t + h/2., y + h*self.Y2/2.)
        f(self.Z3, t + h, y - h*self.Y1 + 2.0*h*self.Y2)
        f(self.Y4, t + h, y + h*self.Y3)
        
        error = N.linalg.norm(h/6*(2*self.Y2 + self.Z3 - 2.0*self.Y3 - self.Y4)/scaling) #normalized 
        
        return t+h, y + h/6.0*(self.Y1 + 2.0*self.Y2 + 2.0*self.Y3 + self.Y4), error
    
    def print_statistics(self, verbose):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,                  verbose)
        self.log_message(' Number of Steps                : %s '%(self.statistics["nsteps"]),             verbose)
        self.log_message(' Number of Function Evaluations : %s '%(self.statistics["nfcn"]),               verbose)
        self.log_message('\nSolver options:\n',                                              verbose)
        self.log_message(' Solver             : RungeKutta4',                                verbose)
        self.log_message(' Solver type        : Adaptive',                                   verbose)
        self.log_message(' Relative tolerance : ' + str(self.options["rtol"]),        verbose)
        self.log_message(' Absolute tolerance : ' + str(self.options["atol"]) + '\n', verbose)
    
    
class RungeKutta4(Explicit_ODE):
    """
    This solver solves an explicit ordinary differential equation using 
    a Runge-Kutta method of order 4.
    
    We want to approximate the solution to the ordinary differential 
    equation of the form,
    
    .. math::

        \dot{y} = f(t,y), \quad y(t_0) = y_0 .
        
    Using a Runge-Kutta method of order 4, the approximation is defined as 
    follow,
    
    .. math::
    
        y_{n+1} = y_n + \\frac{1}{6}(k_1+2k_2+2k_3+k_4)
        
    where,
    
    .. math::
    
        k_1 = hf(t_n,y_n)
        
        k_2 = hf(t_n+\\frac{1}{2}h,y_n+\\frac{1}{2}k_1)
        
        k_3 = hf(t_n+\\frac{1}{2}h,y_n+\\frac{1}{2}k_2)
        
        k_4 = hf(t_n+h,y_n+k_3)
        
    with :math:`h` being the step-size and :math:`y_n` the previous 
    solution to the equation.
    """
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self.Y4 = N.array([0.0]*len(self.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
        
    def step(self, t, y, tf, opts):
        initialize = opts["initialize"]
        
        if initialize:
            self.solver_iterator = self._iter(t,y,tf)

        return self.solver_iterator.next()
    
    def integrate(self, t, y, tf, opts):
        """
        Integrates (t,y) values until t > tf
        """
        [flags, tlist, ylist] = zip(*list(self._iter(t, y, tf)))
        
        return flags[-1], tlist, ylist
    
    def _set_h(self,h):
        try:
            self.options["h"] = float(h)
        except:
            raise AssimuloException("Step-size must be a (scalar) float.")
    
    def _get_h(self):
        """
        Defines the step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0.01'.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["h"]
        
    h=property(_get_h,_set_h)
    
    def _iter(self,t,y,tf):
        h = self.options["h"]
        h = min(h, N.abs(tf-t))
        
        while t+h < tf:
            t, y = self._step(t, y, h)
            yield ID_PY_OK, t,y
            h=min(h, N.abs(tf-t))
        else:
            t, y = self._step(t, y, h)
            yield ID_PY_COMPLETE, t, y
    
    def _step(self, t, y, h):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        
        f(self.Y1, t, y)
        f(self.Y2, t + h/2., y + h*self.Y1/2.)
        f(self.Y3, t + h/2., y + h*self.Y2/2.)
        f(self.Y4, t + h, y + h*self.Y3)
        
        return t+h, y + h/6.*(self.Y1 + 2.*self.Y2 + 2.*self.Y3 + self.Y4)
        
    def print_statistics(self, verbose):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        self.log_message(' Step-length          : %s '%(self.options["h"]), verbose)
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : RungeKutta4',                       verbose)
        self.log_message(' Solver type       : Fixed step\n',                      verbose)
