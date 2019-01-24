cdef class Semiconductor(object):
    cdef:
        str __label
        dict __reference
        list __dopants

    cdef double __band_gap_t(self, double temperature)
    cpdef double band_gap_t(self, double temperature, bint electron_volts=*)
    cpdef double[:] band_gap(self, double[:] temperature, bint electron_volts=*)
    cpdef double n_c_t(self, double temperature)
    cpdef double[:] n_c(self, double[:] temperature)
    cpdef double n_v_t(self, double temperature)
    cpdef double[:] n_v(self, double[:] temperature)
    cpdef double v_e_t(self, double temperature)
    cpdef double[:] v_e(self, double[:] temperature)
    cpdef double v_h_t(self, double temperature)
    cpdef double[:] v_h(self, double[:] temperature)
    cpdef double bulk_charge(self, double mu, double temperature)
    cpdef double el_chem_pot(self, double temperature)
