from __future__ import division, print_function

from libc.math cimport sqrt, exp, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from Schottky.Constants cimport constant
from Schottky.Dopant cimport Dopant


cdef class Semiconductor(object):

    def __init__(self, str label, dict reference, list dopants=None):
        self.__label = label
        self.__reference = reference
        self.__dopants = []
        if dopants is not None:
            for dopant in dopants:
                if isinstance(dopant, Dopant):
                    self.__dopants.append(dopant)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def reference(self):
        return self.__reference

    @reference.setter
    def reference(self, dict reference):
        self.__reference = reference

    @property
    def dopants(self):
        return self.__dopants

    @dopants.setter
    def dopants(self, list dopants):
        self.__dopants = []
        for dopant in dopants:
            if isinstance(dopant, Dopant):
                self.__dopants.append(dopant)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] band_gap(self, double[:] temperature, bint electron_volts=False):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = (self.__reference['band_gap']['E_g0'] -
                         self.__reference['band_gap']['alpha'] * temperature[i] ** 2 /
                         (temperature[i] + self.__reference['band_gap']['beta']))
            if electron_volts:
                result[i] /= constant.q
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_c(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__reference['N_c0'] * temperature[i] ** (3 / 2)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_v(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__reference['N_c0'] * temperature[i] ** (3 / 2)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] v_e(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = sqrt(3 * constant.k * temperature[i] /
                             (self.__reference['carrier_mass']['m_e_coeff'] * constant.m_e))
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] v_h(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = sqrt(3 * constant.k * temperature[i] /
                             (self.__reference['carrier_mass']['m_h_coeff'] * constant.m_e))
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double bulk_charge(self, double temperature, double mu):
        cdef:
            double f, ff, f_threshold = 1.0e-8
            double band_gap, n_c, n_v, v_e, v_h, result
            array[double] z = array('d'), t = array('d')
        z[0] = 1.0e5
        t[0] = temperature
        band_gap = self.band_gap(t)[0]
        n_c = self.n_c(t)[0]
        n_v = self.n_v(t)[0]
        v_e = self.v_e(t)[0]
        v_h = self.v_h(t)[0]
        result = -n_c * exp(-mu / (constant.k * temperature))
        result += n_v * exp((mu - band_gap) / (constant.k * temperature))
        for dopant in self.__dopants:
            f = 2.0
            ff = 0.0
            while fabs(f - ff) > f_threshold:
                f = ff
                ff = dopant.f_eq(temperature,
                             v_e, n_c * exp(-mu / (constant.k * temperature)), n_c,
                             v_h, n_v * exp((mu - band_gap) / (constant.k * temperature)), n_v,
                             f)
            result += dopant.n_t(z)[0] * ((dopant.charge_state[1] - dopant.charge_state[0]) * ff
                                             + dopant.charge_state[0])
        return result

    cpdef double el_chem_pot(self, double temperature):
        cdef:
            double mu = 0.0
        return mu


    def __str__(self):
        return 'Semiconductor: ' + self.__label
