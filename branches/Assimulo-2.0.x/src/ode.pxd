import numpy as N
cimport numpy as N

cdef class ODE:
    cdef public dict options, solver_options, problem_info
    cdef public dict supports
    
    cdef public list event_data
    
    cdef public object problem
    
    cdef public double t_cur
    cdef public N.ndarray y_cur,yd_cur
    
    cdef public list t,y,yd
        
    cpdef log_message(self, message, int level)
    cpdef log_event(self, double time, object event_info, int level)
    cpdef simulate(self, double tfinal, int ncp=*, object ncp_list=*)
    cpdef get_options(self)
