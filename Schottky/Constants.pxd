cdef class Constants(object):
    cdef:
        double __c
        double __avogadro
        double __m_e
        double __q
        double __k
        double __k_by_q
        double __A_R
        double __epsilon_0

    cpdef double thermal_voltage_t(self, double temperature)
    cpdef double[:] thermal_voltage(self, double[:] temperature)

    cpdef double joule_to_ev_point(self, double joule)
    cpdef double[:] joule_to_ev(self, double[:] joule)

    cpdef double ev_to_joule_point(self, double ev)
    cpdef double[:] ev_to_joule(self, double[:] ev)

    cpdef double ev_to_boltzmann_point(self, double ev, double temperature)
    cpdef double[:] ev_to_boltzmann(self, double[:] ev, double temperature)

    cpdef double boltzmann_to_ev_point(self, double boltzmann, double temperature)
    cpdef double[:] boltzmann_to_ev(self, double[:] boltzmann, double temperature)

    cpdef double joule_to_boltzmann_point(self, double joule, double temperature)
    cpdef double[:] joule_to_boltzmann(self, double[:] joule, double temperature)

    cpdef double boltzmann_to_joule_point(self, double boltzmann, double temperature)
    cpdef double[:] boltzmann_to_joule(self, double[:] boltzmann, double temperature)

cdef Constants constant
