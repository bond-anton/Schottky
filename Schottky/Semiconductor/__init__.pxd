cdef class Semiconductor(object):
    cdef:
        str __label
        dict __reference
        list __dopants

    cpdef double[:] band_gap(self, double[:] temperature, bint electron_volts=*)
    cpdef double[:] n_c(self, double[:] temperature)
    cpdef double[:] n_v(self, double[:] temperature)
    cpdef double[:] v_e(self, double[:] temperature)
    cpdef double[:] v_h(self, double[:] temperature)
    cpdef double[:] bulk_charge(self, double[:] temperature, double mu)
