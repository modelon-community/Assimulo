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

from ODE import *
from Problem import Explicit_Problem
from Sundials import Sundials, Sundials_Exception

class Explicit_ODE_Exception(Exception):
    """ An integrator exception. """
    pass

class Explicit_ODE(ODE):
    """
    Baseclass for our explicit ODE integrators.
    """
    
    def __init__(self, problem, y0=None, t0=None):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                
                t0          
                            - Default 'None'. The initial time. If 'None'. the
                              initial time are retrieved from the problem.t0.
                              If set it override problem.t0. If NOT set and NOT
                              defined in problem.t0, t0 is set to 0.0.
                            
                            - Should be a float.
                            
                                Example:
                                    t0 = 1.0
                                    
        """
        ODE.__init__(self) #Sets general attributes
        
        if problem == None:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem')
        
        if isinstance(problem, Explicit_Problem):
            self._problem = problem
            self.f = problem.f
            self.problemname = problem.Problem_Name
        else:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem.')
        
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Explicit_ODE_Exception('y0 must be specified. Either in the problem or in the initialization')

        try:
            if isinstance(y0, int) or isinstance(y0, float):
                y0 = [y0]
            self._problem.y0 = y0[:]
            self.y = [N.array(y0, dtype=float)]
        except ValueError:
            raise Explicit_ODE_Exception('Initial values must be a scalar/list/array of type int or float.')

        if t0 == None:
            if hasattr(problem, 't0'):
                t0 = problem.t0
            else:
                t0 = 0.0

        try:
            if isinstance(t0, list):
                self.t = [float(t0[0])]
            else:
                self.t = [float(t0)]
            self._problem.t0 = self.t[0]
        except ValueError:
            raise Explicit_ODE_Exception('Initial time must be an integer or float.')
            
    def reset(self):
        """
        
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0 and problem.y0.
        
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0)
    
    def integrate(self, t, y, tf, nt):
        pass 

    def __call__(self, tfinal, ncp=0):
        """
        Calls the integrator to perform the simulation over the given time-interval.
        
            Parameters::
            
                tfinal  
                        - Final time for the simulation
                
                        - Should be a float or integer greater than the initial time.
                        
                ncp     
                        - Default '0'. Number of communication points where the 
                          solution is returned. If '0', the integrator will return 
                          at its internal steps.
                          
                        - Should be an integer.
                          
                    Example:
                    
                        __call__(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
                 
        Returns the computed solution which is also saved in 'solver'.t
                                                             'solver'.y
                    
        """

        try:
            tfinal = float(tfinal)
        except ValueError:
            raise Explicit_ODE_Exception('Final time must be an integer or float.')
            
        if self.t[-1] > tfinal:
            raise Explicit_ODE_Exception('Final time must be greater than start time.')
        
        if not isinstance(ncp, int):
            raise Explicit_ODE_Exception('Number of communication points must be an integer')
        if ncp < 0:
            ncp = 0
            if self.verbosity > self.QUIET:
                print 'Number of communication points must be a positive integer, setting' \
                      ' nt = 0.'
        ncp_ori = ncp
        
        while self.t[-1] < tfinal:
            
            solution = list(self.integrate(self.t[-1], self.y[-1], tfinal,ncp))
        
            self.t.extend(q[0] for q in solution)
            self.y.extend(q[1] for q in solution)
            
            if self.is_disc: #Is discontinious?
                [tevent,event_info]=self.disc_info
                
                if ncp > 0:
                    ncp = ncp_ori-len(self.y)+1
                    if ncp < 0:
                        ncp = 0
                #Log the information
                self._log_event_info.append([self.t[-1], event_info])
                
                if self.verbosity > self.NORMAL:
                    print 'A discontinuity occured at t = %e.'%tevent
                if self.verbosity >= self.LOUD:
                    print 'Current switches: ', self.switches
                    print 'Event info: ', event_info
                    
                if self.verbosity >= self.SCREAM:
                    self.print_statistics() #Prints statistics
                    
                if self.verbosity > self.NORMAL:
                    print 'Calling problem specified event handling...'
                
                self._problem.handle_event(self, event_info) #self corresponds to the solver
                #self.event_iteration(event_info) #Event Iteration
            
            if self.verbosity >= self.NORMAL:
                self.print_statistics()
        
        return [self.t, self.y]
    
    def re_init(self,t0, y0):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  
                    - The initial time.
                y0  
                    - The initial values for the states
                
        See information in the __init__ method.
        """
        if len(self.y[-1]) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        Explicit_ODE.__init__(self, self._problem,y0,t0)
    
    def plot(self, mask=None):
        """
        Plot the computed solution.
        
            Parameters::
            
                mask    
                        - Default 'None'. Used to determine which variables that is to be plotted.
                          Used as a list of integers, ones represents the variable that is to be
                          plotted and zeros that is not. 
                        
                        - Should be a list of integers.
                        
                            Example:
                                mask = [1,0] , plots the first variable.
                        
        """
        if not mask:
            P.plot(self.t, self.y)
        else:
            if not isinstance(mask, list):
                raise Explicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.y[-1]):
                raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.y)[:,i])
        
        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problemname)
        P.show()
            
            
    
class Explicit_Euler(Explicit_ODE):
    """
    Explicit Euler.
    """
    def integrate(self, t, y, tf, nt):
        """
        Integrates (t,y) values until t > tf
        """
        if nt <= 0.0:
            raise Explicit_ODE_Exception('Explicit Euler is a fixed step-size method. Provide' \
                                         ' the number of communication points.')
        
        self.h = N.array((tf-self.t[-1])/nt)

        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            yield t,y
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        h = self.h
        return t + h, y + h*f(t, y) 
    
    
class RungeKutta34(Explicit_ODE):
    """
    Adaptive Runge-Kutta of order four.
    Obs. Step rejection not implemented.
    """
    def __init__(self, problem, y0=None, t0=None):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                
                t0          
                            - Default 'None'. The initial time. If 'None'. the
                              initial time are retrieved from the problem.t0.
                              If set it override problem.t0. If NOT set and NOT
                              defined in problem.t0, t0 is set to 0.0.
                            
                            - Should be a float.
                            
                                Example:
                                    t0 = 1.0
                                    
        """
        Explicit_ODE.__init__(self, problem, y0, t0) #Calls the base class
        
        #Default values
        self.initstep = 0.01
        
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
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
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
        
    
    def integrate(self, t, y, tf, nt):
        """
        Integrates (t,y) values until t > tf
        """
        self.h = self.initstep
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            yield t,y
            self.adjust_stepsize()
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def adjust_stepsize(self):
        """
        Adjusts the stepsize.
        """
        if self.error < self.atol*self.atol:
            self.error = self.atol*self.atol
        self.h *= (self.atol/self.error)**(1.0/4.0)
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        h = self.h
        Y1 = f(t, y)
        Y2 = f(t + h/2., y + h*Y1/2.)
        Y3 = f(t + h/2, y + h*Y2/2)
        Z3 = f(t + h, y - h*Y1 + 2*h*Y2)
        Y4 = f(t + h, y + h*Y3)
        self.error = N.linalg.norm(h/6*(2*Y2 + Z3 - 2*Y3 - Y4))
        return t+h, y + h/6*(Y1 + 2*Y2 + 2*Y3 + Y4)
     
    
    
class RungeKutta4(Explicit_ODE):
    """
    Runge-Kutta of order 4.
    """
    def integrate(self, t, y, tf, nt):
        """
        Integrates (t,y) values until t > tf
        """
        if nt <= 0.0:
            raise Explicit_ODE_Exception('RungeKutta4 is a fixed step-size method. Provide' \
                                         ' the number of communication points.')
        
        self.h = N.array((tf-self.t[-1])/nt)

        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            yield t,y
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        h = self.h
        Y1 = f(t, y)
        Y2 = f(t + h/2., y + h*Y1/2.)
        Y3 = f(t + h/2., y + h*Y2/2.)
        Y4 = f(t + h, y + h*Y3)
        return t+h, y + h/6.*(Y1 + 2.*Y2 + 2.*Y3 + Y4)
    
    
class CVode(Explicit_ODE, Sundials):
    """
    Sundials CVode.
    """
    def __init__(self, problem, y0=None, t0=None, switches0=None):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                    
                t0          
                            - Default 'None'. The initial time. If 'None'. the
                              initial time are retrieved from the problem.t0.
                              If set it override problem.t0. If NOT set and NOT
                              defined in problem.t0, t0 is set to 0.0.
                            
                            - Should be a float.
                            
                                Example:
                                    t0 = 1.0
                                    
                switches0   
                            - Default 'None'. Only used for hybrid (discontinuous)
                              systems. If 'None', the switches are retrieved from
                              the problem.switches0. If set, they override the
                              problem.switches0.
                            
                            - Should be a list of booleans.
                            
                                Example:
                                    switches0 = [True, False]
            
        """
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Explicit_ODE_Exception('y0 must not be None.')
            
        Sundials.__init__(self, y0, 'CVode') #Creates an integrator
        Explicit_ODE.__init__(self, problem, y0, t0) #Calls the base class
        
        #Default values
        self.discr = 'Adams' #Setting default discretization to Adams
        self.iter = 'FixedPoint' #Setting default iteration to FixedPoint
        self.maxord = 12 #Setting default maxord to maximum
        self.initstep = 0.0 #Setting the initial step to be estimated
        
        if hasattr(self._problem, 'switches0') and switches0 == None:
            switches0 = self._problem.switches0
        
        if isinstance(switches0, list):
            for x in switches0:
                if not isinstance(x, bool):
                    raise Explicit_ODE_Exception('Switches must be a list of booleans.')
        elif switches0 is not None:
            raise Explicit_ODE_Exception('Switches must be a list of booleans.')
        
        self._problem.switches0 = switches0
        self.switches = switches0
        
        #Determine if we have an event function and sets the integration data
        if hasattr(problem, 'event_fcn'):
            self.event_fcn = self._problem.event_fcn #problem.event_fcn
            self.Integrator.num_event_fcn=len(self.event_fcn(self._problem.t0,self._problem.y0,self._problem.switches0))
            self._ROOT = [self.event_fcn, self._problem.switches0]
        else:
            self.Integrator.num_event_fcn=0
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            if self.switches == None:
                trial = self._problem.jac(self.t[-1],self.y[-1])
            else:
                trial = self._problem.jac(self.t[-1],self.y[-1], self.switches)
            if trial.shape != (len(self.y[-1]),len(self.y[-1])):
                raise Explicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            
            self.jac = self._problem.jac    
            self.Integrator.jacobian = True
            self.usejac = True
            self._RHS = [self.f, self._problem.jac]
        else:
            self.Integrator.jacobian = False
            self.usejac = False
            self._RHS = [self.f]
        
        
        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RHS, self._ROOT]
        else:
            self.problem_spec = [self._RHS]
        
    
    def integrate(self,t,y,tfinal,nt):
        """
        Simulates the problem up until tfinal.
        """
        self.Integrator.cvinit(t,self.problem_spec,y,self.maxord,self.maxsteps,self.initstep)
        return self.Integrator.run(t,tfinal,nt)
    
    def _set_discr_method(self,discr='Adams'):
        """
        This determines the discretization method.
        
            Parameters::
            
                discr   
                        - Default 'Adams', which indicates the use
                          of the Adams method. Can also be set to
                          'BDF' which indicates the use of the BDF
                          method.
                
                    Example:
                        discr = 'BDF'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        if discr=='BDF':
            self.Integrator.discr=2
            self.maxord = 5
        elif discr=='Adams':
            self.Integrator.discr=1
            self.maxord = 12
        else:
            raise Sundials_Exception('Discretization method must be either Adams or BDF')
            
    def _get_discr_method(self):
        """
        This determines the discretization method.
        
            Parameters::
            
                discr   
                        - Default 'Adams', which indicates the use
                          of the Adams method. Can also be set to
                          'BDF' which indicates the use of the BDF
                          method.
                
                    Example:
                        discr = 'BDF'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        discr='Adams'
        if self.Integrator.discr==2:
            discr='BDF'
        return discr

    discr= property(_get_discr_method,_set_discr_method)
    
    def _set_initial_step(self, initstep):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
        self.__initstep = initstep
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.__initstep
        
    initstep = property(_get_initial_step,_set_initial_step)
    
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
        
        if not bool(jac):
            self.Integrator.jacobian = False
            self._RHS = [self.f]
        else:
            self.Integrator.jacobian = True
            if not hasattr(self, 'jac'):
                raise Explicit_ODE_Exception('No jacobian defined.')
            self._RHS = [self.f, self.jac]
            
        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RHS, self._ROOT]
        else:
            self.problem_spec = [self._RHS]
    
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
    
    def _set_iter_method(self,iter='FixedPoint'):
        """
        This determines the iteration method that is be used by the
        solver.
        
            Parameters::
            
                iter    
                        - Default 'FixedPoint', which indicates the
                          use of a fixedpoint iteration method. Can
                          also be set to 'Newton' which indicates
                          the use of a Newton method.
                          
                            Example:
                                iter = 'Newton'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        if iter=='Newton':
            self.Integrator.iter=2
        elif iter=='FixedPoint':
            self.Integrator.iter=1
        else:
            raise Sundials_Exception('Iteration method must be either FixedPoint or Newton')
    
    def _get_iter_method(self):
        """
        This determines the iteration method that is be used by the
        solver.
        
            Parameters::
            
                iter    
                        - Default 'FixedPoint', which indicates the
                          use of a fixedpoint iteration method. Can
                          also be set to 'Newton' which indicates
                          the use of a Newton method.
                          
                            Example:
                                iter = 'Newton'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        iter='FixedPoint'
        if self.Integrator.iter==2:
            iter='Newton'
        return iter
        
    iter = property(_get_iter_method,_set_iter_method)        
    
    def _set_max_ord(self,maxord):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '12', which is the maximum for the
                          Adams method, which is also default. For the
                          BDF method the maximum order is 5. 'maxord'
                          can be set in an interval from 1 to the
                          maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        if not isinstance(maxord,int):
            raise Sundials_Exception('The maximal order must be an integer.')
        if self.Integrator.discr == 1: #Adams
            if maxord > 12:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is greater than that of the Adams. ' \
                           'Setting the maximal order to Adams maximum.'
                self.__maxord=12
            elif maxord < 1:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is lower than the minimum of Adams. ' \
                          'Setting the maximal order to Adams minimum.'
                self.__maxord=1
            else:
                self.__maxord=maxord
        if self.Integrator.discr == 2: #BDF
            if maxord > 5:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is greater than that of the BDF. ' \
                           'Setting the maximal order to BDF maximum.'
                self.__maxord=5
            elif maxord < 1:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is lower than the minimum of BDF. ' \
                          'Setting the maximal order to BDF minimum.'
                self.__maxord=1
            else:
                self.__maxord=maxord
    
    def _get_max_ord(self):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '12', which is the maximum for the
                          Adams method, which is also default. For the
                          BDF method the maximum order is 5. 'maxord'
                          can be set in an interval from 1 to the
                          maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.__maxord

    maxord=property(_get_max_ord,_set_max_ord)
    
    def print_statistics(self):
        """
        Prints the run-time statistics for the problem.
        """
        print 'Final Run Statistics: %s \n' % self.problemname
        
        statistics = self.stats
        keys = statistics.keys()
        keys.sort()
        
        for x in keys:
            print '%s = %s'%(x, statistics[x])
    
    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return self.t[-1]==self.Integrator.event_time
