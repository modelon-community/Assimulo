from lib import sundials_core
import numpy as N

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
        if max_h <= 0.0:
            raise Sundials_Exception('Maximal stepsize must be a positive (scalar) float.')
        self.Integrator.max_h=max_h
        self.__max_h=max_h
    def _get_max_h(self):
        """Returns the maximal stepsize."""
        return self.__max_h
    maxh=property(_get_max_h,_set_max_h,doc='Maximal stepsize')    

    def print_statistics(self):
        """Prints the run-time statistics for the problem."""
        print 'Final Run Statistics: %s \n' % self.problemname
        
        statistics = self.stats
        keys = statistics.keys()
        keys.sort()
        
        for x in keys:
            print '%s = %s'%(x, statistics[x])
    
    @property
    def stats(self):
        """Attribute that returns the run-time statistics from the Integrator."""
        return self.Integrator.stats
    
    @property
    def disc_info(self):
        """Attribute that returns information about an occured event."""
        return [self.Integrator.event_time, self.Integrator.event_info]
