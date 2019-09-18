from scipy.constants import c, N_A, k, e, m_e, epsilon_0

from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

cdef class Constants(object):

    def __init__(self):
        self.__c = c
        self.__avogadro = N_A
        self.__m_e = m_e
        self.__q = e
        self.__k = k
        self.__k_by_q = self.__k / self.__q
        self.__A_R = 1.20173e6  # A*m-2*K-2
        self.__epsilon_0 = epsilon_0

    @property
    def c(self):
        return self.__c

    @property
    def avogadro(self):
        return self.__avogadro

    @property
    def q(self):
        return self.__q

    @property
    def m_e(self):
        return self.__m_e

    @property
    def A_R(self):
        return self.__A_R

    @property
    def k(self):
        return self.__k

    @property
    def epsilon_0(self):
        return self.__epsilon_0

    cpdef double thermal_voltage_t(self, double temperature):
        return self.__k_by_q * temperature

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] thermal_voltage(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.thermal_voltage_t(temperature[i])
        return result

    cpdef double joule_to_ev_point(self, double joule):
        return joule / self.__q

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] joule_to_ev(self, double[:] joule):
        cdef:
            int n = joule.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.joule_to_ev_point(joule[i])
        return result

    cpdef double ev_to_joule_point(self, double ev):
        return ev * self.__q

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] ev_to_joule(self, double[:] ev):
        cdef:
            int n = ev.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.ev_to_joule_point(ev[i])
        return result

    cpdef double ev_to_boltzmann_point(self, double ev, double temperature):
        return ev / self.thermal_voltage_t(temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] ev_to_boltzmann(self, double[:] ev, double temperature):
        cdef:
            int n = ev.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.ev_to_boltzmann_point(ev[i], temperature)
        return result

    cpdef double boltzmann_to_ev_point(self, double boltzmann, double temperature):
        return boltzmann * self.thermal_voltage_t(temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] boltzmann_to_ev(self, double[:] boltzmann, double temperature):
        cdef:
            int n = boltzmann.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.boltzmann_to_ev_point(boltzmann[i], temperature)
        return result

    cpdef double joule_to_boltzmann_point(self, double joule, double temperature):
        return joule / (self.__k * temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] joule_to_boltzmann(self, double[:] joule, double temperature):
        cdef:
            int n = joule.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.joule_to_boltzmann_point(joule[i], temperature)
        return result

    cpdef double boltzmann_to_joule_point(self, double boltzmann, double temperature):
        return boltzmann * self.__k * temperature

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] boltzmann_to_joule(self, double[:] boltzmann, double temperature):
        cdef:
            int n = boltzmann.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.boltzmann_to_joule_point(boltzmann[i], temperature)
        return result

constant = Constants()
