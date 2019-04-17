from Schottky.Trap cimport Trap
from Schottky.Helpers.Cache cimport Cache


cdef class Semiconductor(object):
    cdef:
        str __label
        dict __reference
        list __dopants
        Cache __band_gap_cache
        Cache __trap_eq_occupation_cache
        Cache __bulk_charge_cache
        Cache __el_chem_pot_cache

    cdef double __band_gap_t(self, double temperature)
    cpdef double band_gap_t(self, double temperature, bint electron_volts=*)
    cpdef double[:] band_gap(self, double[:] temperature, bint electron_volts=*)
    cpdef double n_c_t(self, double temperature)
    cpdef double[:] n_c(self, double[:] temperature)
    cpdef double n_v_t(self, double temperature)
    cpdef double[:] n_v(self, double[:] temperature)
    cpdef double n_e_t(self, double mu, double temperature)
    cpdef double[:] n_e(self, double[:] mu, double[:] temperature)
    cpdef double n_h_t(self, double mu, double temperature)
    cpdef double[:] n_h(self, double[:] mu, double[:] temperature)
    cpdef double v_e_t(self, double temperature)
    cpdef double[:] v_e(self, double[:] temperature)
    cpdef double v_h_t(self, double temperature)
    cpdef double[:] v_h(self, double[:] temperature)
    cpdef double trap_eq_occupation(self, Trap trap, double mu, double temperature,
                                    double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double bulk_charge(self, double mu, double temperature, double z=*,
                             double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double el_chem_pot_t(self, double temperature,
                               double f_threshold=*, int max_iter=*, bint verbose=*)
    cpdef double[:] el_chem_pot(self, double[:] temperature,
                                double f_threshold=*, int max_iter=*, bint verbose=*)
