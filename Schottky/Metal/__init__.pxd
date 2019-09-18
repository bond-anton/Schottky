cdef class Metal(object):
    cdef:
        str __label
        double __work_function

    cpdef double work_function_boltzmann_t(self, double temperature)
    cpdef double[:] work_function_boltzmann(self, double[:] temperature)
