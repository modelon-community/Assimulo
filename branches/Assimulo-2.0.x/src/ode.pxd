import numpy as N

cdef class ODE:
    cpdef logg_message(self, message, int level)
    cpdef logg_event(self, double time, object event_info, int level)
    cpdef simulate(self, double tfinal, int ncp=*, object ncp_list=*)
    cpdef get_options(self)
