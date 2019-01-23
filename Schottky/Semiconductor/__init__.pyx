from __future__ import division, print_function

from libc.math cimport sqrt, exp
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

    cpdef double[:] bulk_charge(self, double[:] temperature, double mu):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            double ff
            double[:] band_gap = self.band_gap(temperature)
            double[:] n_c = self.n_c(temperature)
            double[:] n_v = self.n_v(temperature)
            double[:] v_e = self.v_e(temperature)
            double[:] v_h = self.v_h(temperature)
            array[double] result, z = array('d'), template = array('d')
        result = clone(template, n, zero=False)
        z[0] = 1.0e5
        for i in range(n):
            result[i] = -n_c[i] * exp(-mu / (constant.k * temperature[i]))
            result[i] += n_v[i] * exp((mu - band_gap[i]) / (constant.k * temperature[i]))
            for dopant in self.__dopants:
                ff = dopant.f_eq(temperature[i],
                                 v_e[i], n_c[i] * exp(-mu / (constant.k * temperature[i])), n_c[i],
                                 v_h[i], n_v[i] * exp((mu - band_gap[i]) / (constant.k * temperature[i])), n_v[i],
                                 dopant.f(z)[0])
                result[i] += dopant.n_t(z)[0] * ((dopant.charge_state[1] - dopant.charge_state[0]) * ff
                                                 + dopant.charge_state[0])
        return result

    def __str__(self):
        return 'Semiconductor: ' + self.__label
