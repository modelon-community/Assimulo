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
from Problem import Implicit_Problem
from Sundials import Sundials, Sundials_Exception

class Implicit_ODE_Exception(Exception):
    """An integrator exception"""
    pass 

class Implicit_ODE(ODE):
    """
    Baseclass for our implicit ODE integrators.
    """
    
    def __init__(self, problem, y0=None, yd0=None, t0=None):
        """
        Defines the problem and sets the initial values.
        
            res - The residual function
            y0 - Starting values for the none differentiated variables 
            yd0 - Starting values for the differentiated variables
            t0 - Starting time
        """
        ODE.__init__(self) #Sets general attributes
        
        if problem == None:
            raise Implicit_ODE_Exception('The problem needs to be a subclass of a Implicit_Problem')
        
        if isinstance(problem, Implicit_Problem):
            self._problem = problem
            self.res_fcn = problem.f
            self.problemname = problem.Problem_Name
        else:
            raise Implicit_ODE_Exception('The problem needs to be a subclass of a Implicit_Problem.')
        
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0[:]
            else:
                raise Implicit_ODE_Exception('y0 must be specified. Either in the problem or in the initialization')
        
        if yd0 == None:
            if hasattr(problem, 'yd0'):
                yd0 = problem.yd0[:]
            else:
                raise Implicit_ODE_Exception('yd0 must be specified. Either in the problem or in the initialization')
        
        try:
            if isinstance(y0, int) or isinstance(y0, float):
                y0 = [y0]
            if isinstance(yd0, int) or isinstance(yd0, float):
                yd0 = [yd0]
            self._problem.y0 = y0[:]
            self._problem.yd0 = yd0[:]
            self.y = [N.array(y0, dtype=float)]
            self.yd = [N.array(yd0, dtype=float)]
        except ValueError:
            raise Implicit_ODE_Exception('Initial values must be a scalar/list/array of type int or float.')
        
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
            raise Implicit_ODE_Exception('Initial time must be an integer or float.')

    
    def integrate(self, t, y, yd, tf,nt):
        pass
        
    def reset(self):
        """
        Resets the problem if defined in the problem.
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0, self._problem.yd0)
        
    def re_init(self,t0, y0, yd0):
        """
        Re initiates the solver
        """
        if len(self.y[-1]) != len(y0) or len(self.yd[-1]) != len(yd0):
            raise Explicit_ODE_Exception('y0/yd0 must be of the same length as the original problem.')
        
        Implicit_ODE.__init__(self, self._problem,y0,yd0,t0)
            
    #def event_fcn_adjust(self, t, y, yd, sw):
    #    """
    #    This function adjusts the event functions according to Martin Otter et al defined
    #    in (...)
    #    """
    #    r = N.array(self._problem.event_fcn(t,y,yd,sw))
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

    def __call__(self, tfinal, ncp=0):
        """
        Calls the integrator to perform the simulation over the given time-interval.
        
            tfinal - Final time for the simulation
            nt - Number of communication points where the solution is returned.
                 If nt=0, the integrator will return at it's internal steps.
                 
        Returns the computed solution.
        """
        try:
            tfinal = float(tfinal)
        except ValueError:
            raise Implicit_ODE_Exception('Final time must be an integer or float.')
            
        if self.t[-1] > tfinal:
            raise Implicit_ODE_Exception('Final time must be greater than start time.')
        if not isinstance(ncp, int):
            raise Implicit_ODE_Exception('Number of communication points must be an integer.')
        if ncp < 0:
            ncp = 0
            if self.verbosity > self.QUIET:
                print 'Number of communication points must be a positive integer, setting' \
                      ' nt = 0.'
        ncp_ori = ncp
        
        while self.t[-1] < tfinal:
            solution = list(self.integrate(self.t[-1], self.y[-1], self.yd[-1], tfinal,ncp))
        
            self.t.extend(q[0] for q in solution)
            self.y.extend(q[1] for q in solution)
            self.yd.extend(q[2] for q in solution)
            
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
                #self.event_iteration(event_info) #Handles the event iteration
                
        
        return [self.t, self.y, self.yd]
        
    #def event_iteration(self, event_info):
    #    """
    #    Handles the event iteration.
    #    """
    #    nbr_iteration = 0
    #        
    #    while self.max_eIter > nbr_iteration: #Event Iteration
    #        self._problem.event_switch(self, event_info) #Turns the switches
    #        
    #        b_mode = self.event_fcn_adjust(self.t[-1], self.y[-1], self.yd[-1], self.switches)
    #        b_mode -= self.eps_adjust #Adjust for the appended epsilon
    #        self._problem.init_mode(self) #Pass in the solver to the problem specified init_mode
    #        
    #        a_mode = self.event_fcn_adjust(self.t[-1], self.y[-1], self.yd[-1], self.switches)
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
        
    def plot(self):
        """
        Plot the computed solution.
        """
        P.figure(1)
        P.plot(self.t, self.y)
        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problemname)
        #P.show()
        P.figure(2)
        P.plot(self.t, self.yd)
        P.xlabel('time')
        P.ylabel('state derivatives')
        P.title(self.problemname)
        P.show()

        
        
class IDA(Implicit_ODE, Sundials):
    """
    Sundials IDA.
    """
    
    def __init__(self, problem, y0=None, yd0=None, t0=None, switches0=None):
        """
        Defines the problem and sets the initial values.
        
            problem - Defines the problem
            y0 - Starting values for the none differentiated variables 
            yd0 - Starting values for the differentiated variables
            t0 - Starting time
            switches0 - Sets the starting mode
        """
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Implicit_ODE_Exception('y0 must not be None.')
                
        Sundials.__init__(self, y0, 'IDA') #Creates an integrator
        Implicit_ODE.__init__(self, problem, y0, yd0, t0) #Calls the base class

        if hasattr(self._problem, 'algvar'): #Check if the algebraic components are defined in the Problem specifications
            self.algvar = self._problem.algvar
        else:
            self.algvar = [1.0]*len(self.y[0]) #No algebraic variables are set
            
        self.maxord = 5 #Maximal order is set to max
        
        if hasattr(self._problem, 'switches0') and switches0 == None:
            switches0 = self._problem.switches0
        
        if isinstance(switches0, list):
            for x in switches0:
                if not isinstance(x, bool):
                    raise Implicit_ODE_Exception('Switches must be a list of booleans.')
        elif switches0 is not None:
            raise Implicit_ODE_Exception('Switches must be a list of booleans.')
        
        self.switches = switches0
        self._problem.switches0 = switches0
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            #try:
            if self.switches == None:
                trial = self._problem.jac(1,self.t[-1],self.y[-1],self.yd[-1])
            else:
                trial = self._problem.jac(1,self.t[-1],self.y[-1],self.yd[-1], self.switches)
            if trial.shape != (len(self.y[-1]),len(self.y[-1])):
                raise Implicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            #except:
            #    raise Implicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            self.Integrator.jacobian = True
            self.jac = self._problem.jac
            self._RES = [self.res_fcn, self._problem.jac]
        else:
            self.Integrator.jacobian = False
            self._RES = [self.res_fcn]
        
        #Determine if we have an event function and sets the integration data
        if hasattr(problem, 'event_fcn'):
            self.event_fcn = self._problem.event_fcn
            self.Integrator.num_event_fcn=len(self.event_fcn(self._problem.t0,self._problem.y0,self._problem.yd0,self._problem.switches0))
            self._ROOT = [self.event_fcn, self.switches]
            self.problem_spec=[self._RES,self._ROOT]
        else:
            self.Integrator.num_event_fcn=0
            self.problem_spec=[self._RES]
            
        #Default values
        self.tout1 = 0.001
        self.lsoff = False #Use LineSearch
        self.suppress_alg = False #Don't suppres algebraic variables
        
    
    def _set_calcIC_tout1(self, tout1):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        """
        try:
            tout1 = float(tout1)
        except (ValueError,TypeError):
            raise Implicit_ODE_Exception('tout1 must be an integer or float.')
        self.__tout1 = tout1
        
    def _get_calcIC_tout1(self):
        """
        Gets the internal directional value used in make_consistency.
        """
        return self.__tout1
    tout1docstring='Value used for determine the direction of the ' \
                    'integration in make_consistency'
    tout1=property(_get_calcIC_tout1,_set_calcIC_tout1,doc=tout1docstring)
    
    
    def _set_lsoff(self, lsoff):
        """
        Boolean value to turn of Sundials LineSearch when calculating
        initial conditions.
        
            lsoff - False/True (False = LineSearch on)
        """
        if not isinstance(lsoff, bool):
            raise Implicit_ODE_Exception('lsoff must be a boolean.')
        self.__lsoff = lsoff
        
    def _get_lsoff(self):
        """
        Gets the boolean value that determines if LineSearch is active.
        """
        return self.__lsoff
    
    lsoffdocstring='Value to turn of LineSearch'
    lsoff = property(_get_lsoff, _set_lsoff, doc=lsoffdocstring)
    
    def make_consistency(self, method):
        """
        Directs IDA to try to calculate consistant initial conditions.
        
        The options are:
            method = IDA_YA_YDP_INIT - This tries to calculate the algebraic
            componets y and the differential componets of yd given the
            differential componets of y. The algebraic componets of y must 
            have been specified with algvar.
            
            method = IDA_Y_INIT - This tries to calculate all componets of y
            given yd.
            
        See SUNDIALS IDA documentation 4.5.4 for more details.
        """
        self.Integrator.idinit(self.t[-1], self.problem_spec, self.y[-1], self.yd[-1], self.maxord, self.maxsteps)
        
        if method == 'IDA_YA_YDP_INIT':
            [flag, y, yd] = self.Integrator.calc_IC(method,self.tout1,self.lsoff)
        elif method == 'IDA_Y_INIT':
            [flag, y, yd] = self.Integrator.calc_IC(method,self.tout1,self.lsoff)
        else:
            raise Sundials_Exception('Make consistency must be one of the following:'\
                                            ' IDA_YA_YDP_INIT or IDA_Y_INIT')
        
        if flag < 0:
            raise Sundials_Exception('Calculation of initial conditions failed. IDA returned flag %d'%flag)
        else:
            self.y[-1] = y
            self.yd[-1] = yd
        
        return [self.y[-1], self.yd[-1]]
    
    def integrate(self,t,y,yd,tfinal,nt=0):
        """
        Simulates the problem up until tfinal.
        """
        self.Integrator.idinit(t,self.problem_spec,y,yd,self.maxord, self.maxsteps)
        return self.Integrator.run(t,tfinal,nt)
    
    def _set_max_ord(self,maxord):
        """
        Sets the maximal order of the method:
        
            defaults : maximal values:
            BDF  :  maxord= 5
        
        An input value greater than the maximal order will result in the maximum value.
        """
        if not isinstance(maxord,int):
            raise Sundials_Exception('The maximal order must be an integer.')
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
        """Returns the maximal order."""
        return self.__maxord
    maxorddocstring='Maxord: Maximal Order BDF 0 < maxord < 6\n Defaults to None, ' \
                     'which corresponds to the maximal values above.'
    maxord=property(_get_max_ord,_set_max_ord,doc=maxorddocstring)
    
    def _set_suppress_alg(self,suppress_alg):
        """
        Sets the boolean flag which indicates that error-tests are 
        suppressed on algebraic variables.
        """
        if not isinstance(suppress_alg, bool):
            raise Sundials_Exception('Option suppress algebraic variables' \
                                              ' must be defined by a boolean.')
        self.Integrator.suppress_alg=suppress_alg
    def _get_suppress_alg(self):
        """
        Returns the boolean flag which indicates that error-tests are 
        suppressed on algebraic variables.
        """
        return self.Integrator.suppress_alg    
    supprdocstring='True indicates that the error test on algebraic variables is suppressed.' \
              '\n Define the algebraic variables by the attribute algvar.'
    suppress_alg=property(_get_suppress_alg,_set_suppress_alg,doc=supprdocstring)
    
    def _set_algvar(self,algvar):
        """
        Sets the ndarray which indicates which variables are differential (1.0)
        and which are algebraic (0.0).
        """
        algvar = N.array(algvar)

        if isinstance(algvar,N.ndarray) and algvar.dtype=='float':
            if len(algvar) != self.Integrator.dim:
                raise Sundials_Exception('When setting the algebraic variables, the' \
                                    ' vector must be of the same size as the problem dimension.')
            for alg in algvar:
                if alg != 1.0 and alg != 0.0:
                    raise Sundials_Exception('The vector must consist of 1.0 or 0.0 .')
            self.Integrator.algvar=algvar
        else:
            raise Sundials_Exception('Argval vector must be of type float.')
    def _get_algvar(self):
        """
        Returns an ndarray which indicates which variables are differential (1.0)
        and which are algebraic (0.0).
        """
        return self.Integrator.algvar    
    algvardocstring='An ndarray which has an entry 1.0 for a differential variable, otherwise 0.0'
    algvar=property(_get_algvar,_set_algvar,doc=algvardocstring)
    
    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return self.t[-1]==self.Integrator.event_time
