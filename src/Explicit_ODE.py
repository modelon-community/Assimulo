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
        Defines the problem and sets the initial values.
        
            f - The 'right-hand-side' function
            y0 - The initial starting values
            t0 - The initial starting time
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
                y0 = problem.y0[:]
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
        Resets the problem if defined in the problem.
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0)
        
    #def event_fcn_adjust(self, t, y, sw):
    #    """
    #    This function adjusts the event functions according to Martin Otter et al defined
    #    in (...)
    #    """
    #    r = N.array(self._problem.event_fcn(t,y,sw))
    #    rp = N.zeros(len(r))
    #    self.eps_adjust = N.zeros(len(r))
    #    
    #    for i in range(len(sw)):
    #        if sw[i]:
    #            self.eps_adjust[i]=self.eps
    #        else:
    #            self.eps_adjust[i]=-self.eps
    #    
    #    rp = r + self.eps_adjust
    #    
    #    return rp
    
    def integrate(self, t, y, tf, nt):
        pass 

    def __call__(self, tfinal, ncp=0):
        """
        Calls the integrator to perform the simulation over the given time-interval.
        
            tfinal - Final time for the simulation
            ncp - Number of communication points where the solution is returned.
                 If nt=0, the integrator will return it's internal steps.
                 
        Returns the computed solution.
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
                
                if self.verbosity >= self.NORMAL:
                    print 'A discontinuity occured at t = %e.'%tevent
                if self.verbosity >= self.LOUD:
                    print 'Current switches: ', self.switches
                    print 'Event info: ', event_info
                    
                if self.verbosity >= self.SCREAM:
                    self.print_statistics() #Prints statistics
                    
                if self.verbosity >= self.NORMAL:
                    print 'Calling problem specified event handling...'
                
                self._problem.handle_event(self, event_info) #self corresponds to the solver
                #self.event_iteration(event_info) #Event Iteration
        
        return [self.t, self.y]
        
    #def event_iteration(self, event_info):
    #    """
    #    Handles the event iteration.
    #    """
    #    nbr_iteration = 0
    #        
    #    while self.max_eIter > nbr_iteration: #Event Iteration
    #        self._problem.event_switch(self, event_info) #Turns the switches
    #        
    #        b_mode = self.event_fcn_adjust(self.t[-1], self.y[-1], self.switches)
    #        b_mode -= self.eps_adjust #Adjust for the appended epsilon
    #        self._problem.init_mode(self) #Pass in the solver to the problem specified init_mode
    #        
    #        a_mode = self.event_fcn_adjust(self.t[-1], self.y[-1], self.switches)
    #        a_mode -= self.eps_adjust #Adjust for the appended epsilon
    #        
    #        event_info = self.check_eIter(b_mode, a_mode)
    #        
    #        if self.verbosity >= self.SCREAM:
    #            print 'Event iteration?: ', event_info
    #            
    #        if not True in event_info: #Breaks the iteration loop
    #            break
    #        
    #        nbr_iteration += 1
    
    def re_init(self,t0, y0):
        """
        Re initiates the solver
        """
        if len(self.y[-1]) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        Explicit_ODE.__init__(self, self._problem,y0,t0)
    
    def plot(self):
        """
        Plot the computed solution.
        """
        P.plot(self.t, self.y)
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
    def integrate(self, t, y, tf, nt):
        """
        Integrates (t,y) values until t > tf
        """
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
        Defines the problem and sets the initial values.
        
            f - The 'right-hand-side' function
            y0 - Starting values for the none differentiated variables 
            t0 - Starting time
            event_fcn - The event function (To keep track of changes)
            switches0 - Sets the starting mode
        """
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Explicit_ODE_Exception('y0 must not be None.')
            
        Sundials.__init__(self, y0, 'CVode') #Creates an integrator
        Explicit_ODE.__init__(self, problem, y0, t0) #Calls the base class
        
        self.discr = 'Adams' #Setting default discretization to Adams
        self.iter = 'FixedPoint' #Setting default iteration to FixedPoint
        self.maxord = 12 #Setting default maxord to maximum
        
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
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            if self.switches == None:
                trial = self._problem.jac(self.t[-1],self.y[-1])
            else:
                trial = self._problem.jac(self.t[-1],self.y[-1], self.switches)
            if trial.shape != (len(self.y[-1]),len(self.y[-1])):
                raise Explicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
                
            self.Integrator.jacobian = True
            self.jac = self._problem.jac
            self._RHS = [self.f, self._problem.jac]
        else:
            self.Integrator.jacobian = False
            self._RHS = [self.f]
        
        #Determine if we have an event function and sets the integration data
        if hasattr(problem, 'event_fcn'):
            self.event_fcn = self._problem.event_fcn #problem.event_fcn
            self.Integrator.num_event_fcn=len(self.event_fcn(self._problem.t0,self._problem.y0,self._problem.switches0))
            self._ROOT = [self.event_fcn, self._problem.switches0]
            self.problem_spec=[self._RHS, self._ROOT]
        else:
            self.Integrator.num_event_fcn=0
            self.problem_spec=[self._RHS]
        
        
    
    def integrate(self,t,y,tfinal,nt):
        """
        Simulates the problem up until tfinal.
        """
        self.Integrator.cvinit(t,self.problem_spec,y,self.maxord,self.maxsteps)
        return self.Integrator.run(t,tfinal,nt)
    
    def _set_discr_method(self,discr='Adams'):
        """
        This sets the discretization method which can either be set
        to Adams (default) or BDF.
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
        """Returns the discretization method."""
        discr='Adams'
        if self.Integrator.discr==2:
            discr='BDF'
        return discr
    discrdocstring = 'Can be set to "BDF" or "Adams" (default)'
    discr= property(_get_discr_method,_set_discr_method,doc=discrdocstring)
    
    def _set_iter_method(self,iter='FixedPoint'):
        """
        This sets the iteration method which can either be set to
        FixedPoint (default) or Newton.
        """
        if iter=='Newton':
            self.Integrator.iter=2
        elif iter=='FixedPoint':
            self.Integrator.iter=1
        else:
            raise Sundials_Exception('Iteration method must be either FixedPoint or Newton')
    
    def _get_iter_method(self):
        """Returns the iteration method."""
        iter='FixedPoint'
        if self.Integrator.iter==2:
            iter='Newton'
        return iter
    iterdocstring = 'Can be set to "Newton" or "FixedPoint" (default)'
    iter = property(_get_iter_method,_set_iter_method,doc=iterdocstring)        
    
    def _set_max_ord(self,maxord):
        """
        Sets the maximal order of the method:
        defaults = maximal values:
        Adams:  maxord=12
        BDF  :  maxord= 5
        
        An input value greater than the default will result in the default value.
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
        Returns the used maximal order.
        """
        return self.__maxord
    maxorddocstring='Maxord: Maximal Order\n Adams  0 < maxord < 13\n BDF 0 < maxord < 6' \
                    '\n Defaults to None, which corresponds to the maximal values above.'
    maxord=property(_get_max_ord,_set_max_ord,doc=maxorddocstring)
    
    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return self.t[-1]==self.Integrator.event_time
