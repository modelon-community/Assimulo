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
import scipy as S
import scipy.linalg as LIN
import copy

from assimulo.exception import *
from assimulo.ode import *

from assimulo.explicit_ode import Explicit_ODE
from assimulo.implicit_ode import Implicit_ODE

from assimulo.lib import radar5


class Radar_Exception(Exception):
    pass

class Radar5ODE(Explicit_ODE):
    """
    Radar5     
    """
    
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Delay_Explicit_Problem' class.
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["inith"]    = 0.01
        self.options["newt"]     = 7 #Maximum number of newton iterations
        self.options["thet"]     = 1.e-3 #Boundary for re-calculation of jac
#        self.options["fnewt"]    = 0.0 #Stopping critera for Newtons Method
        self.options["fnewt"]    = 0.03 #Stopping critera for Newtons Method
        self.options["quot1"]    = 1.0 #Parameters for changing step-size (lower bound)
        self.options["quot2"]    = 1.2 #Parameters for changing step-size (upper bound)
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 8.0 #Parameters for step-size selection (upper bound)
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000
        
        self.options["alpha"] = 0.0 # Parameter to tune the error control of dense output (smaller = stronger control)
        self.options["tckbp"] = 5.0 # Parameter for controlling the search for breaking points
        self.options["ieflag"] = 0 # Switch between different modes of error control
        self.options["mxst"] = 100 # The maximum number of stored dense output points
        self.options["usejaclag"]   = True if self.problem_info["jaclag_fcn"] else False
        
        SQ6 = N.sqrt(6.0)
        C1 = (4.0-SQ6)/10.0 
        C2 = (4.0+SQ6)/10.0 
        self.C1M1 = C1-1.0 
        self.C2M1 = C2-1.0 
        
        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["njac"]        = 0 #Number of Jacobian evaluations
        self.statistics["njacfcn"]     = 0 #Number of function evaluations when evaluating the jacobian
        self.statistics["errfail"]     = 0 #Number of step rejections
        self.statistics["nlu"]         = 0 #Number of LU decompositions
        self.statistics["nstepstotal"] = 0 #Number of total computed steps (may NOT be equal to nsteps+nerrfail)
        
        #Internal values
        self._leny = len(self.y) #Dimension of the problem
        self._type = '(explicit)'
        self._yDelayTemp = []
        self._ntimelags = len(self.problem.lagcompmap)
        for i in range(self._ntimelags):
            self._yDelayTemp.append(range(len(self.problem.lagcompmap[i])))
        flat_lagcompmap = []
        for comp in self.problem.lagcompmap:
            flat_lagcompmap.extend(comp)
        self._nrdens = len(N.unique(flat_lagcompmap))
        self._ipast = N.unique(flat_lagcompmap).tolist()+[0]
        self._grid = N.array([])
        
#        if hasattr(problem, 'pbar'):
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
        self._tlist = []
        self._ylist = []
        
    def _solout(self,nr, told, t, hold, y, cont,irtrn):
        """
        This method is called after every successful step taken by Radar5
        """
        #print "SOLOUT:", told, t, hold, y, cont
#        print cont
#        print told, t, told + hold
        if self._opts["output_list"] is None:
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
#                        yval[i] = radar5.contr5(i+1,self.problem_info["dim"],output_list[output_index],t,hold)
                        yval[i] = radar5.contr5(i+1,self.problem_info["dim"],output_list[output_index],cont,t,hold)
#                        yval[i] = radar5.contr5(i+1,output_list[output_index], cont)
                        
                    self._ylist.append(yval)

                    output_index = output_index+1
            except IndexError:
                pass
            self._opts["output_index"] = output_index
        
        return irtrn

    #def coutput(self,t):
        #Nx = self.problem_info["dim"]
        #y = N.zeros(Nx)
        
        #theta, pos = radar5.lagr5(10, t, None, self.arglag, self.past,  self.problem.phi,  self.problem.ipast)
        #for i in range(1,Nx+1):
            #y[i-1] = radar5.ylagr5(i, theta, pos, self.problem.phi,  self.past,  self.problem.ipast)
        #return y

    def coutput(self, t, i = -1):
        """
            Return the continous output solution at time t.
            
            t: time
            i: solution component (default -1 gives the whole vector)
        """
        Nx = self.problem_info["dim"]
        y = N.zeros(Nx)
        
        
        # t belongs to the interval (tk[ik], tk[ik+1])
        ik = N.searchsorted(self.tk, t) - 1
        
        I = self.idif*ik
        
        H = self.past[I+self.idif-1]
        theta = (t - (self.past[I] + H))/H
        
        
        if i == -1:
            # The line below this comment is what's effectively happening,
            # but it is unfortunately extremely slow compared to the
            # vectorized version below that doesn't use the cpoly function:
            #return N.array([self.cpoly(i, I, theta) for i in range(self.problem_info["dim"])])
            nrds = self._nrdens
            I = I + 1
            I2 = I + self.problem_info["dim"]
            return self.past[I:I2] + theta*(self.past[nrds+I:nrds+I2] + (theta-self.C2M1)*(self.past[2*nrds+I:2*nrds+I2] + (theta-self.C1M1)*(self.past[3*nrds+I:3*nrds+I2])))  
        elif i >= 0:
            return self.cpoly(i, I, theta)
        else:
            raise ValueError('i has to be either -1 or a positive integer <= the problem dimension')
            
    def cpoly(self, i, I, theta):
        """
            Evaluate the I:th dense output polynomial for component i at theta.
        """
        nrds = self._nrdens
        I = I + i + 1
        return self.past[I] + theta*(self.past[nrds+I] + (theta-self.C2M1)*(self.past[2*nrds+I] + (theta-self.C1M1)*(self.past[3*nrds+I])))  

    def arglag(self, i, t, y, past, ipast):
        return self.problem.time_lags(t,y)[i-1]
        

    def compute_ydelay(self, t, y, past, ipast):
        ydelay = self._yDelayTemp
        for i in range(1, self._ntimelags+1):
            theta, pos = radar5.lagr5(i, t, y, self.arglag, past,  self.problem.phi,  ipast)

            for j, val in enumerate(self.problem.lagcompmap[i-1]):
                ydelay[i-1][j] = radar5.ylagr5(val+1, theta, pos, self.problem.phi,  past,  ipast)

        return ydelay

    def F(self, t, y, past, ipast):
        #print 'F:', t, y, past, ipast
 
        # First find the correct place in the past vector for each time-lag
        # then evaluate all required solution components at that point
        ydelay = self.compute_ydelay(t,y, past,  ipast)

        # Now we can compute the right-hand-side
        return self.problem.rhs(t, y, ydelay)
    
    def Fjac(self, t, y, past, ipast):
        # First find the correct place in the past vector for each time-lag
        # then evaluate all required solution components at that point
        ydelay = self.compute_ydelay(t,y,  None, None,  past,  ipast)
        
        # Now we can compute the right-hand-side
        return self.problem.jac(t, y, ydelay)


    def integrate(self, t, y, tf, opts):
        ITOL  = 1 #Both atol and rtol are vectors
        IJAC  = 0 if self.usejac else 0 #Switch for the jacobian, 0==NO JACOBIAN
        MLJAC = self.problem_info["dim"] #The jacobian is full
        MUJAC = self.problem_info["dim"] #See MLJAC
        IMAS  = 0 #The mass matrix is the identity
        MLMAS = self.problem_info["dim"] #The mass matrix is full
        MUMAS = self.problem_info["dim"] #See MLMAS
        IOUT  = 1 #solout is called after every step
        WORK  = N.array([0.0]*30) #Work (double) vector
        IWORK = N.array([0]*30) #Work (integer) vector
        
        #Setting work options
        WORK[0] = N.finfo(N.double).eps        # Rounding unit
        WORK[1] = self.safe     
        WORK[2] = self.thet
        WORK[3] = self.fnewt
        WORK[4] = self.quot1
        WORK[5] = self.quot2
        WORK[6] = self.maxh
        WORK[7] = self.fac1
        WORK[8] = self.fac2
        WORK[9] = self.alpha
        WORK[10] = self.tckbp
        
        #Setting iwork options
        IWORK[1] = self.maxsteps
        IWORK[2] = self.newt
        IWORK[7] = 1
        IWORK[10] = self.ieflag
        IWORK[11] = self.mxst
        IWORK[12] = len(self.grid)
        IWORK[13] = 1
        IWORK[14] = self._nrdens
        
        self.idif = 4*self._nrdens + 2
        lrpast = self.mxst*self.idif
        past = N.zeros(lrpast)
        
        #past = N.zeros(self.mxst*(4*self.problem.nrdens+2))
#        print WORK
#        print IWORK
        
        #Dummy methods
        mas_dummy = lambda t:x
        jac_dummy = (lambda t:x) if not self.usejac else self.Fjac
        jaclag_dummy = (lambda t:x) if not self.usejaclag else self.problem.jaclag
        nlags = 0 if not self.usejaclag else self.problem.nlags
        njacl = 0 if not self.usejaclag else self.problem.njacl
        
        #Store the opts
        self._opts = opts
        #print "INIT", t,y,tf,self.inith, self.problem.ipast
        #print "GRID", self.problem.grid, self.problem.ngrid
        #t, y, h, iwork, flag, past = radar5.assimulo_radar5(self.F,            \
        a = radar5.assimulo_radar5(self.F,            \
                                       self.problem.phi,        \
                                       self.arglag,             \
                                       t,                       \
                                       y.copy(),                \
                                       tf,                      \
                                       self.inith,              \
                                       self.rtol*N.ones(self.problem_info["dim"]), \
                                       self.atol,               \
                                       ITOL,                    \
                                       jac_dummy,               \
                                       IJAC,                    \
                                       MLJAC,                   \
                                       MUJAC,                   \
                                       jaclag_dummy,            \
                                       nlags,                   \
                                       njacl,                   \
                                       IMAS,                    \
                                       self._solout,            \
                                       IOUT,                    \
                                       WORK,                    \
                                       IWORK,                   \
                                       self.grid.tolist()+[0.0],       \
                                       self._ipast,      \
                                       mas_dummy,               \
                                       MLMAS,                   \
                                       MUMAS,
                                       past),
#                                       lrpast)#,                   \
                                       #past,                    \
                                       #IWORK[14]+1,             \
##                                       IWORK[14],               \
                                       #self.problem.ngrid,      \
                                       #len(past)                \
                                       #)
        t, y, h, iwork, flag, past = a[0]
        #print a[0]
        #print len(a[0])
        #self.past = copy.deepcopy(past)
        self.past = past
        self.tk = N.trim_zeros(self.past[::self.idif], 'b')
        self.hk = N.trim_zeros(self.past[self.idif-1:-1:self.idif], 'b')
        
        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("Radar5 failed with flag %d"%flag)
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[16]
        self.statistics["nfcn"]        += iwork[13]
        self.statistics["njac"]        += iwork[14]
        self.statistics["nstepstotal"] += iwork[15]
        self.statistics["errfail"]     += iwork[17]
        self.statistics["nlu"]         += iwork[18]
        
        return flag, self._tlist, self._ylist
    
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of function evaluations           : '+str(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Jacobian evaluations           : '+ str(self.statistics["njac"]),    verbose)
        self.log_message(' Number of error test failures            : '+ str(self.statistics["errfail"]),       verbose)
        self.log_message(' Number of LU decompositions              : '+ str(self.statistics["nlu"]),       verbose)
        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : Radar5 ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
        
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
    
    
    def plot_stepsize(self):
        """
        Plots the step-size.
        """
        P.semilogy(N.diff(self.t),drawstyle='steps-post')
        P.title(self.problem.name)
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
            self.options["newt"] = int(newt)
        except (ValueError, TypeError):
            raise Radar_Exception('The newt must be an integer or float.')
        
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
        return self.options["newt"]
        
    newt = property(_get_newt,_set_newt)
    
    def _set_fnewt(self, fnewt):
        try:
            self.options["fnewt"] = float(fnewt)
        except (ValueError, TypeError):
            raise Radar_Exception('The fnewt must be an integer or float.')
        
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
        return self.options["fnewt"]
        
    fnewt = property(_get_fnewt,_set_fnewt)
    
    def _set_safe(self, safe):
        try:
            self.options["safe"] = float(safe)
        except (ValueError, TypeError):
            raise Radar_Exception('The safe must be an integer or float.')

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
    
    def _set_thet(self, thet):
        try:
            self.options["thet"] = float(thet)
        except (ValueError, TypeError):
            raise Radar_Exception('The thet must be an integer or float.')
        
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
        return self.options["thet"]
        
    thet = property(_get_thet, _set_thet)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise Radar_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise Radar_Exception('Maximal stepsize must be a positive (scalar) float.')
        
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
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Radar_Exception('The initial step must be an integer or float.')
        
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
    
    
    def _set_quot1(self, quot1):
        try:
            self.options["quot1"] = float(quot1)
        except (ValueError, TypeError):
            raise Radar_Exception('The quot1 must be an integer or float.')
    
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
        return self.options["quot1"]
        
    quot1 = property(_get_quot1, _set_quot1)
    
    def _set_quot2(self, quot2):
        try:
            self.options["quot2"] = float(quot2)
        except (ValueError, TypeError):
            raise Radar_Exception('The quot2 must be an integer or float.')
    
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
        return self.options["quot2"]
        
    quot2 = property(_get_quot2, _set_quot2)
    
    def _set_fac1(self, fac1):
        try:
            self.options["fac1"] = float(fac1)
        except (ValueError, TypeError):
            raise Radar_Exception('The fac1 must be an integer or float.')
            
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
            raise Radar_Exception('The fac2 must be an integer or float.')
        
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
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
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
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_usejaclag(self, jaclag):
        self.options["usejaclag"] = bool(jaclag)
    def _get_usejaclag(self, jaclag):
        return self.options["usejaclag"]
    usejaclag = property(_get_usejaclag,_set_usejaclag)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise Radar_Exception("atol must be of length one or same as the dimension of the problem.")

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
            raise Radar_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise Radar_Exception('Relative tolerance must be a positive (scalar) float.')
    
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
    
    def _get_grid(self):
        """
        TODO
        """
        return self._grid
    def _set_grid(self, grid):
        self._grid
    grid=property(_get_grid,_set_grid)
    
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
            raise Radar_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def _get_alpha(self):
        """
        Used to tune the error control of dense output. A small value means
        stricter control and should be used for problems with almost
        discontinuous solutions while alpha = 1.0 can be used for problems
        with fairly smooth solutions.
        
        
            Parameters::
            
                alpha
                            - Default 0
                            
                            - Should be in [0.0, 1.0]
        """
        return self.options["alpha"]
    
    def _set_alpha(self, alpha):
        try:
            self.options["alpha"] = float(alpha)
        except (TypeError, ValueError):
            raise Radar_Exception("alpha must be a float.")
        if (self.options["alpha"] < 0.0) or (self.options["alpha"] > 1.0):
            raise Radar_Exception("alpha must be in the interval [0,1].")
    
    alpha = property(_get_alpha, _set_alpha)
    
    
    def _get_tckbp(self):
        """
        Used to control the search for breaking points. If the error increases
        by a factor larger than tckbp in one step the routine that searches for
        breaking points is activated.
        
        
            Parameters::
            
                tckbp
                            - Default 5.0
                            
                            - Should be a positive float
        """
        return self.options["tckbp"]
    
    def _set_tckbp(self, tckbp):
        try:
            self.options["tckbp"] = float(tckbp)
        except (TypeError, ValueError):
            raise Radar_Exception("tckbp must be a float.")
        if (self.options["tckbp"] < 0.0):
            raise Radar_Exception("tckbp must be a positive float.")
    
    tckbp = property(_get_tckbp, _set_tckbp)
    

    def _get_ieflag(self):
        """
        Switch between different modes of error control.
        -1: pure control of the dense output (makes use of a quadratic and a
            linear polynomial interpolating the stage values)
         0: mixed control of dense output and discrete output
         1: simpler mixed control
         2: pure control of the discrete output (is provided by the routine ESTRAD)        
        
            Parameters::
            
                ieflag
                            - Default 0.0
                            
                            - Should be either -1, 0, 1 or 2
        """
        return self.options["ieflag"]
    
    def _set_ieflag(self, ieflag):
        try:
            self.options["ieflag"] = int(ieflag)
        except (TypeError, ValueError):
            raise Radar_Exception("ieflag must be an integer.")
        if (self.options["ieflag"] < -1) or (self.options["ieflag"] > 2):
            raise Radar_Exception("ieflag must be either -1, 0, 1 or 2.")
    
    ieflag = property(_get_ieflag, _set_ieflag)
    
    
    def _get_mxst(self):
        """
        The maximum number of steps stored in the dense output array.
        The dimension of this array will be nrdens*(4*self.nrdens+2).
        
            Parameters::
            
                mxst
                            - Default 100
                            
                            - Should be a positive integer
        """
        return self.options["mxst"]
    
    def _set_mxst(self, mxst):
        try:
            self.options["mxst"] = int(mxst)
        except (TypeError, ValueError):
            raise Radar_Exception("mxst must be a positive integer.")
        if (self.options["mxst"] < 1):
            raise Radar_Exception("mxst must be a positive integer.")

        mxst = property(_get_mxst, _set_mxst)
        
    def _set_usejaclag(self, jaclag):
        self.options["usejaclag"] = bool(jaclag)
    
    def _get_usejaclag(self):
        """
        This sets the option to use the user defined lagged jacobian. If a
        user provided lagged jacobian is implemented into the problem the
        default setting is to use that lagged jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejaclag  
                        - True - use user defined lagged jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejaclag = False
        """
        return self.options["usejaclag"]
    
    usejaclag = property(_get_usejaclag,_set_usejaclag)
