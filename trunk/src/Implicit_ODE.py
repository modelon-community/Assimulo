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
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                    
                yd0         
                            - Default 'None'. The initial values for the state
                              derivatives. If 'None', the initial values are
                              retrieved from the problem.yd0. If set they
                              override problem.yd0.
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    yd0 = [0.0, 0.0]
                                
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
        
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0, problem.y0 and problem.yd0
        
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0, self._problem.yd0)
        
    def re_init(self,t0, y0, yd0):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  - The initial time.
                y0  - The initial values for the states
                yd0 - The initial values for the state derivatives.
                
        See information in the __init__ method.
        """
        if len(self.y[-1]) != len(y0) or len(self.yd[-1]) != len(yd0):
            raise Explicit_ODE_Exception('y0/yd0 must be of the same length as the original problem.')
        
        Implicit_ODE.__init__(self, self._problem,y0,yd0,t0)

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
                                                             'solver'.yd
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
        
        while N.abs(self.t[-1]-tfinal) > self._SAFETY*(N.abs(tfinal)+N.abs(self.t[-1]-tfinal)/(ncp+1.0)):
            
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
                #self.event_iteration(event_info) #Handles the event iteration
        
        if self.verbosity >= self.NORMAL:
            self.print_statistics()
        
        
        return [self.t, self.y, self.yd]
        
        
    def plot(self, mask=None, der=False, **kwargs):
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
                        
                der     
                        - Default 'False'. When 'True' plots the derivative variables also.
                        
                        - Should be a boolean.
                        
                            Example:
                                der = True
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        P.figure(1)
        if not mask:
            P.plot(self.t, self.y, **kwargs)
        else:
            if not isinstance(mask, list):
                raise Implicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.y[-1]):
                raise Implicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.y)[:,i], **kwargs)

        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problemname)

        
        if der and not mask:
            P.figure(2)
            P.plot(self.t, self.yd, **kwargs)
            P.xlabel('time')
            P.ylabel('state derivatives')
            P.title(self.problemname)
        elif mask and der:
            P.figure(2)
            if not isinstance(mask, list):
                raise Implicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.yd[-1]):
                raise Implicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.yd)[:,i], **kwargs)
                    
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
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                    
                yd0         
                            - Default 'None'. The initial values for the state
                              derivatives. If 'None', the initial values are
                              retrieved from the problem.yd0. If set they
                              override problem.yd0.
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    yd0 = [0.0, 0.0]
                                
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
                raise Implicit_ODE_Exception('y0 must not be None.')
                
        Sundials.__init__(self, y0, 'IDA') #Creates an integrator
        Implicit_ODE.__init__(self, problem, y0, yd0, t0) #Calls the base class

        if hasattr(self._problem, 'algvar'): #Check if the algebraic components are defined in the Problem specifications
            self.algvar = self._problem.algvar
        else:
            self.algvar = [1.0]*len(self.y[0]) #No algebraic variables are set
            
        
        
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
        
        #Determine if we have an event function and sets the integration data
        if hasattr(problem, 'event_fcn'):
            self.event_fcn = self._problem.event_fcn
            self.Integrator.num_event_fcn=len(self.event_fcn(self._problem.t0,self._problem.y0,self._problem.yd0,self._problem.switches0))
            self._ROOT = [self.event_fcn, self.switches]
        else:
            self.Integrator.num_event_fcn=0
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            if self.switches == None:
                trial = self._problem.jac(1,self.t[-1],self.y[-1],self.yd[-1])
            else:
                trial = self._problem.jac(1,self.t[-1],self.y[-1],self.yd[-1], self.switches)
            if trial.shape != (len(self.y[-1]),len(self.y[-1])):
                raise Implicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            
            self.jac = self._problem.jac
            self.Integrator.jacobian = True
            self.usejac = True
            self._RES = [self.res_fcn, self._problem.jac]
        else:
            self.Integrator.jacobian = False
            self.usejac = False
            self._RES = [self.res_fcn]
        
        
        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RES, self._ROOT]
        else:
            self.problem_spec = [self._RES]
        
        
        #Default values
        self.tout1 = 0.001
        self.lsoff = False #Use LineSearch
        self.suppress_alg = False #Don't suppres algebraic variables
        self.initstep = 0.0 #Setting the initial step to be estimated
        self.maxord = 5 #Maximal order is set to max
        self.maxh = 0.0 #Setting the maximum absolute step length to infinity
        
    
    def _set_calcIC_tout1(self, tout1):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        
            Parameters::
            
                tout1       
                            - Default '0.001'.
                            
                            - Should be a float.
                            
                                Example:
                                    tout1 = 0.01
        """
        try:
            tout1 = float(tout1)
        except (ValueError,TypeError):
            raise Implicit_ODE_Exception('tout1 must be an integer or float.')
        self.__tout1 = tout1
        
    def _get_calcIC_tout1(self):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        
            Parameters::
            
                tout1       
                            - Default '0.001'.
                            
                            - Should be a float.
                            
                                Example:
                                    tout1 = 0.01
        """
        return self.__tout1

    tout1=property(_get_calcIC_tout1,_set_calcIC_tout1)
    
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
            self._RES = [self.res_fcn]
        else:
            self.Integrator.jacobian = True
            if not hasattr(self, 'jac'):
                raise Implicit_ODE_Exception('No jacobian defined.')
            self._RES = [self.res_fcn, self.jac]
            
        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RES, self._ROOT]
        else:
            self.problem_spec = [self._RES]
    
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
    
    
    def _set_lsoff(self, lsoff):
        """
        Boolean value to turn OFF Sundials LineSearch when calculating
        initial conditions.
            
            Parameters::
            
                lsoff   
                        - Default 'False'. False indicates the use of
                          linesearch.

                        - Should be a boolean.
                        
                            Example:
                                lsoff = True
        """
        if not isinstance(lsoff, bool):
            raise Implicit_ODE_Exception('lsoff must be a boolean.')
        self.__lsoff = lsoff
        
    def _get_lsoff(self):
        """
        Boolean value to turn OFF Sundials LineSearch when calculating
        initial conditions.
            
            Parameters::
            
                lsoff   
                        - Default 'False'. False indicates the use of
                          linesearch.

                        - Should be a boolean.
                        
                            Example:
                                lsoff = True
        """
        return self.__lsoff
    
    lsoff = property(_get_lsoff, _set_lsoff)
    
    def make_consistency(self, method):
        """
        Directs IDA to try to calculate consistant initial conditions.
            
            Parameters::
            
                method  
                        - 'IDA_YA_YDP_INIT'
                                - This tries to calculate the
                                  algebraic components of y and the differential
                                  components of yd given the differential components
                                  of y. The algebraic components of y must have been
                                  specified with the property 'algvar'. The property
                                  'tout1' is also used in the calculations which should
                                  represent the the next output point.
                                
                        - 'IDA_Y_INIT' 
                                - This tries to calculate all components
                                  of y given yd.
            
        See SUNDIALS IDA documentation 4.5.4 for more details.
        """
        self.Integrator.idinit(self.t[-1], self.problem_spec, self.y[-1], self.yd[-1], self.maxord, self.maxsteps,self.initstep,self.maxh)
        
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
        self.Integrator.idinit(t,self.problem_spec,y,yd,self.maxord, self.maxsteps, self.initstep, self.maxh)
        return self.Integrator.run(t,tfinal,nt)
    
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
            raise Implicit_ODE_Exception('The initial step must be an integer or float.')
        
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
    
    def _set_max_ord(self,maxord):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '5', which is the maximum.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
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
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '5', which is the maximum.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.__maxord

    maxord=property(_get_max_ord,_set_max_ord)
    
    def _set_suppress_alg(self,suppress_alg):
        """
        A boolean flag which indicates that the error-tests are 
        suppressed on algebraic variables. The algebraic variables
        are defined by setting the property 'algvar'.
        
            Parameters::
            
                suppress_alg    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 4.5.7 'IDASetSuppressAlg' 
        for more details.
        """
        if not isinstance(suppress_alg, bool):
            raise Sundials_Exception('Option suppress algebraic variables' \
                                              ' must be defined by a boolean.')
        self.Integrator.suppress_alg=suppress_alg
    def _get_suppress_alg(self):
        """
        A boolean flag which indicates that the error-tests are 
        suppressed on algebraic variables. The algebraic variables
        are defined by setting the property 'algvar'.
        
            Parameters::
            
                suppress_alg    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 4.5.7 'IDASetSuppressAlg' 
        for more details.
        """
        return self.Integrator.suppress_alg    

    suppress_alg=property(_get_suppress_alg,_set_suppress_alg)
    
    def _set_algvar(self,algvar):
        """
        A vector for defining which variables are differential and
        which are algebraic.
        
            Parameters::
            
                algvar  
                        - The value True(1.0) indicates an differential
                          variable and the value False(0.0) indicates an
                          algebraic variable.
                          
                        - Should be a list or a numpy vector (ndarray)
                        
                            Example:
                                algvar = [1.0, 0.0, 1.0]
                                algvar = [True, False, True]
                                algvar = [1,0,1]
                                
        """
        try:
            algvar = N.array(algvar,dtype='float')
        except (ValueError, TypeError):
            raise Sundials_Exception('algvar needs to be a list or a numpy vector.')
            
        try:
            size = len(algvar)
        except TypeError:
            raise Sundials_Exception('algvar needs to be a list or a numpy vector.')
        
        if len(algvar) != self.Integrator.dim:
            raise Sundials_Exception('When setting the algebraic variables, the' \
                                ' vector must be of the same size as the problem dimension.')
        
        for alg in algvar:
            if alg != 1.0 and alg != 0.0:
                raise Sundials_Exception('The vector must consist of 1.0 or 0.0 .')
        self.Integrator.algvar=algvar

    def _get_algvar(self):
        """
        A vector for defining which variables are differential and
        which are algebraic.
        
            Parameters::
            
                algvar  
                        - The value 1.0 indicates an differential
                          variable and the value 0.0 indicates an
                          algebraic variable.
                          
                        - Should be a list or a numpy vector (ndarray)
                        
                            Example:
                                algvar = [1.0, 0.0, 1.0]
                                
        """
        return self.Integrator.algvar    

    algvar=property(_get_algvar,_set_algvar)
    
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
