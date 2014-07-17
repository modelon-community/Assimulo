#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N

from assimulo.exception import *
from assimulo.ode import *

from assimulo.implicit_ode import Implicit_ODE

try:
    from assimulo.lib.glimda import glimda
except ImportError:
    print("Could not find GLIMDA")

class GLIMDA(Implicit_ODE):
    """
    GLIMDA is a solver for nonlinear index-2 DAEs f(q'(t,x),x,t)=0.
    
    Details about the implementation (FORTRAN) can be found in the PhD
    dissertation::
    
        General Linear Methods for Integrated Circuit Design
        
        Author: Steffen Voigtmann
    """
    
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
        """
        Implicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["inith"]    = 0.01
        self.options["newt"]     = 5 #Maximum number of newton iterations
        self.options["maxord"]   = 3 #Maximum order used
        self.options["minord"]   = 1 #Minimum order used
        self.options["order"]    = 0 #Variable order (0) or fixed order (1-3)
        self.options["maxsteps"] = 100000 #Maximum number of steps
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["minh"]     = N.finfo(N.double).eps #Minimum step-size      
        self.options["maxretry"] = 15 #Maximum number of consecutive retries  
        
        #Solver support 
        self.supports["report_continuously"] = True 
        self.supports["interpolated_output"] = False
        self.supports["state_events"] = False 
        
        self._leny = len(self.y) #Dimension of the problem
        self._type = '(implicit)'
        
    def initialize(self):
        #Reset statistics
        self.statistics.reset()
            
        self._tlist = []
        self._ylist = []
        self._ydlist = []
        
    def _get_print_level(self):
        """
        Helper method to determine the print level of GLIMDA related
        to the verbosity setting.
        """
        val = int(self.verbosity/10)
        
        
        #val (1-5) 1=SCREAM
        #          2=LOUD
        #          3=NORMAL
        #          4=WHISPER
        #          5=QUIET
        
        if val == 3 or val == 4:
            return 1
        elif val == 5:
            return 0
        elif val == 2:
            return 2
        else:
            return 2
            
    def _solout(self, fsol, tph, h, p, y,yd):
        """
        tph = t+h
        h = stepsize
        p = order
        y = y
        """
        if self._opts["report_continuously"]:
            self.report_solution(tph, y, yd, self._opts)
        else:   
            self._tlist.append(tph)
            self._ylist.append(y.copy())
            self._ydlist.append(yd.copy())
    
    def integrate(self, t, y, yd, tf, opts):
        #F = f(y,x,t) #f(q'(x,t),x,t)=0
        #Q = q(x,t)
        
        ITOL  = 1 #Both atol and rtol are vectors
        INUMA = 1 #Evaluates A = df/dy numerically
        INUMD = 1 #Evaluates D = dq/dx numerically
        INUMB = 1 #Evaluates B = df/dx numerically
        IODE = 0 #Problem is a DAE (1 == ODE)
        IADCONST = 1 #The leading term (A,D) is constant 
        
        #Options vectors
        IOPT = N.array([0]*9) #Integer options
        ROPT = N.array([0.0]*11) #Real options
        
        #Setting work options
        IOPT[0] = self._get_print_level() #Print level (default 2)
        IOPT[1] = self.newt #Newton iterations
        IOPT[2] = self.maxord #Maximum order used
        IOPT[3] = self.minord #Minimum order used
        IOPT[4] = 1 #The particular order 3 method to use
        IOPT[5] = self.order #Use variable order
        IOPT[6] = 5 #Number of steps to keep before changing order
        IOPT[7] = self.maxsteps #Maximum number of steps
        IOPT[8] = self.maxretry #Maximum number of convergence failures

        #Setting iwork options
        ROPT[0] = 0.1 #Controls the error in the Newton iterations
        ROPT[1] = 0.1 #Controls the error in the Newton updates abs(dx_err) < tolX
        ROPT[2] = 0.1 #Controls whether a step is accepted
        ROPT[3] = 0.1 #Controls how much the stepsize is adjusted after conv fail 0.1*h
        ROPT[4] = 2.0 #Controls how much the stepsize can be increased
        ROPT[5] = 0.5 #Controls how much the stepsize can be decreased
        ROPT[6] = tf-t if tf-t < self.maxh else self.maxh #Controls how much the stepsize can be used
        ROPT[7] = self.minh #Minimal stepsize
        ROPT[8] = 1e-3 #Theshold for the convergence rate
        ROPT[9] = 0.8 #Threshould for changing the convergence rate when changing the order downwards
        ROPT[10]= 0.0 #Minimum condition number
        
        #Dummy methods
        dfdy_dummy = lambda t:x #df/dy
        dfdx_dummy = lambda t:x #df/dx
        dqdx_dummy = lambda t:x #dq/dx
        qeval_dummy = lambda x,t:x #q(x,t)
        res_dummy = lambda yd,y,t:self.problem.res(t,y,yd) #Needed to correct the order of the arguments

        #Store the opts
        self._opts = opts
        
        yret, ydret, istats, flag =  glimda(res_dummy, qeval_dummy, dfdy_dummy, 
                            dfdx_dummy,dqdx_dummy,t,tf,y.copy(),yd.copy(),self.inith,
                                self.atol, self.rtol*N.ones(self.problem_info["dim"]),ITOL,INUMA,INUMD,INUMB,
                                IODE, IADCONST, IOPT, ROPT, self._solout)
        #print self._tlist, self._ylist
        #Checking return
        if flag == 0:
            flag = ID_PY_COMPLETE
        else:
            err_mess = {-1:"Could not get method of order p",
                        -2:"Too many steps, increase maxsteps.",
                        -3:"Too many consecutive failures, increase ...",
                        -4:"Stepsize too small, decrease minh",
                        -7:"Could not compute the tolerances atol/rtol",
                        -8:"Could not evaluate partial derivatives.",
                        -9:"Could not evaluate q(x,t)."}
                        
            raise Exception("GLIMDA failed with flag %d."%flag+" "+err_mess[flag])
        
        #Retrieving statistics
        self.statistics["nsteps"]      += istats[0]
        self.statistics["nfcns"]        += istats[3]
        self.statistics["njacs"]        += istats[4]
        self.statistics["nnfails"]   += istats[2]
        self.statistics["nerrfails"]     += istats[1]
        self.statistics["nlus"]         += istats[5]
        
        return flag, self._tlist, self._ylist, self._ydlist
    
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        Implicit_ODE.print_statistics(self, verbose) #Calls the base class
        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : GLIMDA ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)

    def _set_newt(self, newt):
        try:
            self.options["newt"] = int(newt)
        except (ValueError, TypeError):
            raise GLIMDA_Exception('The newt must be an integer or float.')
        if self.options["newt"] < 0:
            raise GLIMDA_Exception("Maximum number of Newton iterations must be positive.")
        
    def _get_newt(self):
        """
        Maximum number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        return self.options["newt"]
        
    newt = property(_get_newt,_set_newt)
    
    def _set_maxord(self, maxord):
        try:
            self.options["maxord"] = int(maxord)
        except (ValueError, TypeError):
            raise GLIMDA_Exception('The maxord must be an integer or float.')
        if self.options["maxord"] < 1 or self.options["maxord"] > 3:
            raise GLIMDA_Exception('The maximum order must be between 1 and 3.')
        
    def _get_maxord(self):
        """
        Maximum order to be used (1-3).
        
            Parameters::
            
                maxord
                        - Default '3'.
                        
                        - Should be an integer.
                        
                            Example:
                                maxord = 2
        """
        return self.options["maxord"]
        
    maxord = property(_get_maxord,_set_maxord)
    
    def _set_minord(self, minord):
        try:
            self.options["minord"] = int(minord)
        except (ValueError, TypeError):
            raise GLIMDA_Exception('The minord must be an integer or float.')
        if self.options["minord"] < 1 or self.options["minord"] > 3:
            raise GLIMDA_Exception('The minimum order must be between 1 and 3.')
        
    def _get_minord(self):
        """
        Minimum order to be used (1-3).
        
            Parameters::
            
                minord
                        - Default '1'.
                        
                        - Should be an integer.
                        
                            Example:
                                maxord = 2
        """
        return self.options["minord"]
        
    minord = property(_get_minord,_set_minord)
    
    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 100000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise GLIMDA_Exception("Maximum number of steps must be a positive integer.")
        if max_steps < 0:
            raise GLIMDA_Exception("Maximum number of steps must be positive.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def _set_min_h(self,min_h):
        try:
            self.options["minh"] = float(min_h)
        except (ValueError,TypeError):
            raise GLIMDA_Exception('Minimum stepsize must be a (scalar) float.')
        if self.options["minh"] < 0:
            raise GLIMDA_Exception('Minimum stepsize must be a positiv (scalar) float.')
        
    def _get_min_h(self):
        """
        Defines the minimum step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default machine precision.
                          
                        - Should be a float.
                        
                            Example:
                                minh = 0.01
                                
        """
        return self.options["minh"]

    minh=property(_get_min_h,_set_min_h) 

    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise GLIMDA_Exception("atol must be of length one or same as the dimension of the problem.")

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
            raise GLIMDA_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise GLIMDA_Exception('Relative tolerance must be a positive (scalar) float.')
    
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
    
    def _set_order(self, order):
        try:
            self.options["order"] = int(order)
        except (ValueError, TypeError):
            raise GLIMDA_Exception('The order must be an integer or float.')
        if self.options["order"] < 0 or self.options["order"] > 3:
            raise GLIMDA_Exception('The order must be between 0 and 3.')
        
    def _get_order(self):
        """
        Determines if GLIMDA should use a variable order method (0) or 
        a fixed order method (1-3).
        
            Parameters::
            
                order
                        - Default '0'.
                        
                        - Should be an integer.
                        
                            Example:
                                order = 2
        """
        return self.options["order"]
        
    order = property(_get_order,_set_order)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise GLIMDA_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise GLIMDA_Exception('Maximal stepsize must be a positiv (scalar) float.')
        
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

    def _set_maxretry(self, maxretry):
        try:
            self.options["maxretry"] = int(maxretry)
        except (ValueError,TypeError):
            raise GLIMDA_Exception('Maximum number of consecutive retries must be an Integer.')
        if self.options["maxretry"] < 0:
            raise GLIMDA_Exception('Maximum number of consecutive retries must be an positive Integer.')

    def _get_maxretry(self):
        """
        Defines the maximum number of consecutive number of retries
        after a convergence failure.
        
            Parameters::
            
                maxretry
                        - Default '15'
                        
                        - Should be a Integer.
                        
                            Example:
                                maxretry = 5
                                
        """
        return self.options["maxretry"]
        
    maxretry = property(_get_maxretry, _set_maxretry)
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise GLIMDA_Exception('The initial step must be an integer or float.')
        
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
