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
cimport numpy as N

from assimulo.support import set_type_shape_array

include "constants.pxi" #Includes the constants (textual include)

    
cdef class cProblem:
    
    name = "---"
    
    def __init__(self,  y0 = None, double t0 = 0.0, p0 = None, sw0 = None, name = None):
        
        if not y0  is None:
            self.y0  = set_type_shape_array(y0)
        if not p0  is None:
            self.p0 = set_type_shape_array(p0) 
        if not sw0 is None:
            self.sw0 = set_type_shape_array(sw0, bool)
        if name:
            self.name = name
        self.t0  = t0
        
    cdef public int _sensitivity_result
    
    cpdef initialize(self, solver):
        """
        Method for specializing initiation.
        """
        solver.log_message("No initialization defined for the problem.", LOUD)
    
    cpdef reset(self):
        """
        Resets a problem to its default values.
        """
        pass
        
    cpdef handle_event(self, object solver, event_info):
        """
        Defines how to handle a discontinuity. This functions gets called when
        a discontinuity has been found in the supplied event functions. The solver
        is the solver attribute while the event_info is a list of length 2 where
        the first element is a list containing information about state events and
        the second element is a Boolean for indicating if there has been a time
        event. If there has not been a state event the first element is an empty
        list. The state event list contains a set of integers of values (-1,0,1),
        the values indicates which state event has triggered (determined from 
        state_event(...) ) and the value indicates to where the state event is 'headed'.
        """
        solver.log_message("No event handling defined.", NORMAL)
       
    
    cpdef finalize(self,object solver):
        """
        Method for specifying the finalization options when the simulation have
        finished.
        """
        solver.log_message("No finalization defined for the problem.", LOUD)

cdef class cImplicit_Problem(cProblem):
    
    def __init__(self, object res=None, y0=None, yd0=None,double t0=0.0, 
                                                          p0=None, sw0=None, name = None):
        cProblem.__init__(self, y0, t0, p0, sw0, name)
        if res is not None:
            self.res = res
        if yd0 is not None:
            self.yd0 = set_type_shape_array(yd0)
        
    
    def handle_result(self, solver, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        """
        Method for specifying how the result is handled. By default the
        data is stored in three vectors: solver.(t/y/yd).
        """
        cdef int i = 0
        
        solver.t_sol.extend([t])
        solver.y_sol.extend([y])
        solver.yd_sol.extend([yd])
        
        #Store sensitivity result (variable _sensitivity_result are set from the solver by the solver)
        if self._sensitivity_result == 1:
            for i in range(solver.problem_info["dimSens"]):
                solver.p_sol[i] += [solver.interpolate_sensitivity(t, i=i)] 
        
    cpdef res_internal(self, N.ndarray[double, ndim=1] res, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        try:
            res[:] = self.res(t,y,yd)
        except:
            return ID_FAIL
        return ID_OK
        
cdef class cOverdetermined_Problem(cProblem):
    
    def __init__(self, object res=None, y0=None, yd0=None,double t0=0.0, 
                                                          p0=None, sw0=None, name=None):
        cProblem.__init__(self, y0, t0, p0, sw0, name)
        if res is not None:
            self.res = res
        if yd0 is not None:
            self.yd0 = set_type_shape_array(yd0)
        
    
    def handle_result(self, solver, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in three vectors: solver.(t/y/yd).
        """
        cdef int i = 0
        
        solver.t_sol.extend([t])
        solver.y_sol.extend([y])
        solver.yd_sol.extend([yd])
        
    cpdef res_internal(self, N.ndarray[double, ndim=1] res, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        try:
            res[:] = self.res(t,y,yd)
        except:
            return ID_FAIL
        return ID_OK
    
cdef class cExplicit_Problem(cProblem):
    
    def __init__(self, object rhs=None, y0=None,double t0=0.0, p0=None, sw0=None, name = None):
        
        cProblem.__init__(self, y0, t0, p0, sw0, name)        
        if rhs is not None:
            self.rhs = rhs
    def handle_result(self, solver, double t, N.ndarray[double, ndim=1] y):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in two vectors: solver.(t/y).
        """
        cdef int i = 0
        
        solver.t_sol.extend([t])
        solver.y_sol.extend([y])
        
        #Store sensitivity result (variable _sensitivity_result are set from the solver by the solver)
        if self._sensitivity_result == 1:
            for i in range(solver.problem_info["dimSens"]):
                solver.p_sol[i] += [solver.interpolate_sensitivity(t, i=i)] 
                
    cpdef int rhs_internal(self, N.ndarray[double, ndim=1] yd, double t, N.ndarray[double, ndim=1] y):
        try:
            yd[:] = self.rhs(t,y)
        except:
            return ID_FAIL
        return ID_OK
        
    cpdef N.ndarray res(self, t, y, yd, sw=None):
        if sw == None:
            return yd-self.rhs(t,y)
        else:
            return yd-self.rhs(t, y, sw)
        
cdef class cDelay_Explicit_Problem(cExplicit_Problem):
    def __init__(self, object rhs=None, y0=None, phi = None, arglag = None, lagcompmap = None, jaclag = None, nlags = None, njacl = None, double t0=0.0, p0=None, sw0=None):
        cExplicit_Problem.__init__(self, rhs, y0, t0, p0, sw0)
        
        if phi is not None:
            self.phi = phi
        if arglag is not None:
            self.arglag = arglag
        if jaclag is not None:
            self.jaclag = jaclag
            
        self.lagcompmap = lagcompmap

        self.nlags = nlags # Number of delay arguments to differentiate for
        self.njacl = njacl # Number of possible delay arguments that fcn can be differentiated with respect to

cdef class cSingPerturbed_Problem(cExplicit_Problem):
    
    def __init__(self, rhs1=None, rhs2=None, yy0=None, zz0=None, double t0=0.0, eps=None, name = None):
        if rhs1 is not None:
            self.rhs1 = rhs1
        if rhs2 is not None:
            self.rhs2 = rhs2
        if eps is not None:
            self.eps = eps
        if yy0 is not None:
            self.yy0 = set_type_shape_array(yy0)
        if zz0 is not None:
            self.zz0 = set_type_shape_array(zz0)
        # preparing things so that the problem also describes a 
        # classical explicit problem without exposing the structure
        # of a singularly perturbed problem
        if yy0 is not None and zz0 is not None:
            y0 = N.hstack((self.yy0,self.zz0))
        elif yy0 is not None:
            y0 = self.yy0
        elif zz0 is not None:
            y0 = self.zz0
        self.n = len(self.yy0) if yy0 is not None else 0
        self.m = len(self.zz0) if zz0 is not None else 0
        cExplicit_Problem.__init__(self, y0=y0, t0=t0, name=name)   
        
    def rhs(self,t,y):
        yy=y[:self.n]
        zz=y[self.n:]
        yydot=self.rhs1(t,yy,zz)
        zzdot=self.rhs2(t,yy,zz)
        # componentwise division by eps as it is the diagonal of 
        # a diagonal matrix
        if self.eps != None: 
            zzdot /= self.eps 
        return N.hstack((yydot,zzdot))
            
class Delay_Explicit_Problem(cDelay_Explicit_Problem):
    pass

class Implicit_Problem(cImplicit_Problem):
    """
        Problem for our implicit integrators (DAEs). A problem
        consists of the residual function and some initial conditions.
        
        Parameters ::
          
            res   
                Function that calculates the residual. Depending on
                the problem and the support of the solver, this function can
                have the following input parameters.
                
                    res(t,y,yd)   - Normal DAE
                    res(t,y,yd,sw)   - An DAE with different modes, sw is a list of
                                       switches (boolean list) which should be held
                                       constant during the integration and only be
                                       changed when an event have occured. Used together
                                       with event functions.
                    res(t,y,yd,p)   - An DAE with parameters for which sensitivities
                                       should be calculated.
                    res(t,y,yd,sw,p) - An DAE with both parameters and switches.
                    
                    Returns:
                        A numpy array of size len(y).
            y0
                Defines the starting values of y0.
            yd0
                Defines the starting values of yd0.
            t0
                Defines the starting time.
            p0 (Depending on if the solver supports sensitivity calculations)
                Parameters for which sensitivites are to be calculated
            sw0 (Depending on if the solver supports state events)
                Defines the starting values of the switches. 
                Should be a list of Booleans.
                
        Parameters (optionally contained in class) ::
        
            algvar
                Defines the differential and algebraic components of the problem.
                Should be a list of integers. For more information, see the
                property algvar in IDA.
        
        Signature of default or user provided methods. Their use is solver dependent. ::
        
            def state_events(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def time_events(self, t, y, yd, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, c, t, y, yd, sw)
                Defines the Jacobian, which should be of the form
                J = dF/dx + c*dF/dx'.
                
                Returns:
                    A numpy array of size len(y)*len(y).
                    
            def handle_result(self, solver, t, y, yd)
                Method for specifying how the result is  handled. 
                By default the data is stored in three vectors, solver.(t_sol/y_sol/yd_sol). 
                If the problem to be solved also involve sensitivities these results are
                stored in p_sol
                
            def handle_event(self, object solver, event_info):
                Defines how to handle a discontinuity. This functions gets called when
                a discontinuity has been found in the supplied event functions. The solver
                is the solver attribute. The event_info is a list of length 2 where
                the first element is a list containing information about state events and
                the second element is a Boolean  indicating if there has been a time
                event. If there has not been a state event the first element is an empty
                list. The state event list contains a set of integers of values (-1,0,1),
                the values indicate which state event has triggered (determined from 
                state_event(...) ) and the value indicates to where the state event is 'headed'.
    """
    pass
class Overdetermined_Problem(cOverdetermined_Problem):
    """
        Problem for integrators for overdetermined DAES (ODAEs). A problem
        consists of the residual function with more components than state variables 
        and some initial conditions.
        
        Parameters ::
          
            res   
                Function that calculates the residual. Depending on
                the problem and the support of the solver, this function can
                have the following input parameters.
                
                    res(t,y,yd)   - Normal ODAE
                    res(t,y,yd,sw)   - An ODAE with different modes, sw is a list of
                                       switches (boolean list) which should be held
                                       constant during the integration and only be
                                       changed when an event have occured. Used together
                                       with event functions.
                    
                    Returns:
                        A numpy array of size neq > len(y).
            y0
                Defines the starting values of y0.
            yd0
                Defines the starting values of yd0.
            t0
                Defines the starting time.
            sw0 (Depending on if the solver supports state events)
                Defines the starting values of the switches. 
                Should be a list of Booleans.
                
        Parameters (optionally contained in class) ::
        
            algvar
                Defines the differential and algebraic components of the problem.
                Should be a list of integers. For more information, see the
                property algvar in IDA.
        
        Available (optional) options (depending on the solver support)::
        
            def state_events(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def time_events(self, t, y, yd, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, c, t, y, yd, sw)
                Defines the Jacobian, which should be of the form
                J = dF/dx + c*dF/dx'.
                
                Returns:
                    A numpy array of size neq*len(y).
                    
            def handle_result(self, solver, t, y, yd)
                Method for specifying how the result is  handled. 
                By default the data is stored in three vectors, solver.(t_sol/y_sol/yd_sol). 
                If the problem to be solved also involve sensitivities these results are
                stored in p_sol
                
            def handle_event(self, object solver, event_info):
                Defines how to handle a discontinuity. This functions gets called when
                a discontinuity has been found in the supplied event functions. The solver
                is the solver attribute. The event_info is a list of length 2 where
                the first element is a list containing information about state events and
                the second element is a Boolean indicating if there has been a time
                event. If there have not been a state event the first element is an empty
                list. The state event list contains a set of integers of values (-1,0,1),
                the values indicate which state event has triggered (determined from 
                state_event(...) ) and the value indicates to where the state event is 'headed'.
    """
    pass
class Explicit_Problem(cExplicit_Problem):
    """
        Problem for our explicit integrators (ODEs). A problem
        consists of the right-hand-side and some initial conditions.
 
        Parameters::
            
            rhs 
                Function that calculates the right-hand-side. Depending on
                the problem and the support of the solver, this function has
                the following input parameters:
                
                    rhs(t,y)      - Normal ODE
                    rhs(t,y,sw)   - An ODE with different modes, sw is a list of
                                    switches (Boolean list) which should be held
                                    constant during the integration and only be
                                    changed when an event have occured. Used together
                                    with event functions.
                    rhs(t,y,p)  - An ODE with parameters for which sensitivities
                                    should be calculated.
                    rhs(t,y,sw,p) - An ODE with both parameters and switches.
                    
                    Returns:
                        A numpy array of size len(y).
            
            y0
                Defines the starting values 
            t0
                Defines the starting time
            p0 (Depending on if the solver supports sensitivity calculations)
                Parameters for which sensitivites are to be calculated
            sw0 (Depending on if the solver supports state events)
                Defines the starting values of the switches. 
                Should be a list of Booleans.
        
        Signature of default or user provided methods. Their use is solver dependent. ::
        
            def state_events(self ,t ,y, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                    
            def time_events(self, t, y, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, t, y, sw=None)
                Defines the jacobian. J=df/dx.
                
                Returns:
                    A numpy matrix of size len(y)*len(y).
                    
            def jacv(self, t, y, fy, v)
                Defines a Jacobian Vector product. df/dx*v.
                
                Returns:
                    A numpy vector of size len(y).
            
            def handle_result(self, solver, t, y)
                Method for specifying how the result is handled. 
                By default the data is stored in two vectors, solver.(t_sol/y_sol). If
                the problem to be solved also involve sensitivities these results are
                stored in p_sol
                
            def handle_event(self, object solver, event_info):
                Defines how to handle a discontinuity. This functions is called when
                a discontinuity has been found in the supplied event functions. The solver
                is the solver attribute. The event_info is a list of length 2 where
                the first element is a list containing information about state events and
                the second element is a Boolean  indicating if there has been a time
                event. If there has not been a state event the first element is an empty
                list. The state event list contains a set of integers of values (-1,0,1),
                the values indicate which state event has triggered (determined from 
                state_event(...) ) and the value indicates to where the state event is 'headed'.
    """
    pass
class SingPerturbed_Problem(cSingPerturbed_Problem):
    """
        Problem for singularly perturbed problems of the form
        The equation is in the form
        .. math::
           :nowrap:

           \begin{eqnarray}
           \dot{y} & = & \mathrm{rhs}_1(t,y,z) \\
           \vareps \dot{z} & = & \mathrm{rhs}_(t,y,z)
           \end{eqnarray}
 
        Parameters::
            
            rhs1 
                Function that calculates the 'slow' right-hand-side: f (in formula above)
                
                    rhs1(t,y,z)    - 'slow' ODE
                    Returns:
                        A numpy array of size len(y).
            rhs2 
                Function that calculates the 'fast' right-hand-side: g (in formula above)
                
                    rhs2(t,y,z)    - 'fast' ODE
                    Returns:
                        A numpy array of size len(z).
            eps diagonal of a len(z) x len(z) matrix with small numbers
                    A numpy array of size len(z)
            
            yy0
                Defines the starting values of y0
            zz0
                Defines the starting values of z0
            t0
                Defines the starting time
    """
    pass


cdef class cAlgebraic_Problem:
    
    name = "---"
    
    def __init__(self, object res, y0, y0_min = None, y0_max=None, 
                    y0_nominal=None, object jac=None, object jacv=None, 
                    object prec_solve=None, object prec_setup=None, name = None):
        
        if res != None:
            self.res = res
        if jac != None:
            self.jac = jac
        if jacv != None:
            self.jacv = jacv
        if prec_solve != None:
            self.prec_solve = prec_solve
        if prec_setup != None:
            self.prec_setup = prec_setup
        
        if not y0 is None:
            self.y0 = set_type_shape_array(y0)
        
        if not y0_min is None:
            self.y0_min = set_type_shape_array(y0_min)
        else:
            self.y0_min = set_type_shape_array([MIN_VALUE]*len(self.y0))
            
        if not y0_max is None:
            self.y0_max = set_type_shape_array(y0_max)
        else:
            self.y0_max = set_type_shape_array([MAX_VALUE]*len(self.y0))
            
        if not y0_nominal is None:
            self.y0_nominal = set_type_shape_array(y0_nominal)
        else:
            self.y0_nominal = set_type_shape_array([1.0]*len(self.y0))
            
        if name:
            self.name = name
    
    cpdef initialize(self, solver):
        """
        Method for specializing initiation.
        """
        solver.log_message("No initialization defined for the problem.", LOUD)
    
    cpdef finalize(self,object solver):
        """
        Method for specifying the finalization options when the simulation have
        finished.
        """
        solver.log_message("No finalization defined for the problem.", LOUD)

class Algebraic_Problem(cAlgebraic_Problem):
    """
        Problem class for solving systems of algebraic equations.
        
        .. math::
           :nowrap:
           
           0 = F(y)
           
    """
    pass
