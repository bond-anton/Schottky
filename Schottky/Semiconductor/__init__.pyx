from __future__ import division, print_function

from libc.math cimport sqrt, exp, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from Schottky.Constants cimport constant
from Schottky.Dopant cimport Dopant
from Schottky.Trap cimport Trap


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

    cdef double __band_gap_t(self, double temperature):
        return (self.__reference['band_gap']['E_g0'] -
                self.__reference['band_gap']['alpha'] * temperature ** 2 /
                (temperature + self.__reference['band_gap']['beta']))

    cpdef double band_gap_t(self, double temperature, bint electron_volts=False):
        if electron_volts:
            return self.__band_gap_t(temperature) / constant.q
        return self.__band_gap_t(temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] band_gap(self, double[:] temperature, bint electron_volts=False):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.band_gap_t(temperature[i], electron_volts)
        return result

    cpdef double n_c_t(self, double temperature):
        return self.__reference['N_c0'] * temperature ** (3 / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_c(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_c_t(temperature[i])
        return result

    cpdef double n_v_t(self, double temperature):
        return self.__reference['N_v0'] * temperature ** (3 / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_v(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_v_t(temperature[i])
        return result

    cpdef double n_e_t(self, double mu, double temperature):
        return self.n_c_t(temperature) * exp(-mu / (constant.k * temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e(self, double[:] mu, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_t(mu[i], temperature[i])
        return result

    cpdef double n_h_t(self, double mu, double temperature):
        return self.n_v_t(temperature) * exp((mu - self.band_gap_t(temperature)) / (constant.k * temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h(self, double[:] mu, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_t(mu[i], temperature[i])
        return result

    cpdef double v_e_t(self, double temperature):
        return sqrt(3 * constant.k * temperature /
                    (self.__reference['carrier_mass']['m_e_coeff'] * constant.m_e))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] v_e(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.v_e_t(temperature[i])
        return result

    cpdef double v_h_t(self, double temperature):
        return sqrt(3 * constant.k * temperature /
                    (self.__reference['carrier_mass']['m_h_coeff'] * constant.m_e))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] v_h(self, double[:] temperature):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.v_h_t(temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double trap_eq_occupation(self, Trap trap, double mu, double temperature, int max_iter=100):
        cdef:
            int i = 0
            double fa, fb, fm=0.0, ff_a, ff_b, ff_m, f_threshold = 1.0e-8, dm
            double band_gap, n_c, n_v, v_e, v_h, n_e, n_h
        band_gap = self.__band_gap_t(temperature)
        n_c = self.n_c_t(temperature)
        n_v = self.n_v_t(temperature)
        v_e = self.v_e_t(temperature)
        v_h = self.v_h_t(temperature)
        n_e = self.n_e_t(mu, temperature)
        n_h = self.n_h_t(mu, temperature)
        if trap.cb_bound:
            trap.energy_v = band_gap - trap.energy_c
        else:
            trap.energy_c = band_gap - trap.energy_v
        fa = 0.0
        fb = 1.0
        ff_a = fa - trap.f_eq(temperature, v_e, n_e, n_c, v_h, n_h, n_v, fa)
        ff_b = fb - trap.f_eq(temperature, v_e, n_e, n_c, v_h, n_h, n_v, fb)
        if ff_a == 0:
            fm = fa
        elif ff_b == 0:
            fm = fb
        else:
            dm = fb - fa
            for i in range(max_iter):
                dm *= 0.5
                fm = fa + dm
                ff_m = fm - trap.f_eq(temperature, v_e, n_e, n_c, v_h, n_h, n_v, fm)
                if ff_m * ff_a >= 0:
                    fa = fm
                if ff_m == 0 or fabs(dm) < f_threshold:
                    break
        # if i == max_iter - 1:
        #     print('    trap_eq_occupation iters=', i, 'T=', temperature, 'K')
        return fm

    @boundscheck(False)
    @wraparound(False)
    cpdef double bulk_charge(self, double mu, double temperature, double z=1.0e5, int max_iter=100):
        cdef:
            double fm, result
            array[double] _z = clone(array('d'), 1, zero=False)
        _z[0] = z
        n_e = self.n_e_t(mu, temperature)
        n_h = self.n_h_t(mu, temperature)
        result = n_h - n_e
        for dopant in self.__dopants:
            fm = self.trap_eq_occupation(dopant, mu, temperature, max_iter=max_iter)
            result += dopant.n_t(_z)[0] * ((dopant.charge_state[1] - dopant.charge_state[0]) * fm
                                          + dopant.charge_state[0])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double el_chem_pot_t(self, double temperature, int max_iter=100):
        cdef:
            int i = 0
            double dm, xm=0.0, fm, fa, fb
            double tol, xtol = constant.k * temperature / 1000, rtol = 4.5e-16
            double xa = 0.0, xb = self.__band_gap_t(temperature)
        # if temperature < 8.0:
        #     fa = self.el_chem_pot_t(10.0)
        #     fb = self.el_chem_pot_t(8.0)
        #     return temperature * (fa - fb) / 2.0 + fa - 10.0 * (fa - fb) / 2.0
        tol = xtol + rtol*(fabs(xa) + fabs(xb))
        fa = self.bulk_charge(xa, temperature)
        fb = self.bulk_charge(xb, temperature)
        if fa * fb > 0:
            return -1.0
        if fa == 0.0:
            return xa
        if fb == 0:
            return xb
        dm = xb - xa
        for i in range(max_iter):
            dm *= 0.5
            xm = xa + dm
            fm = self.bulk_charge(xm, temperature, max_iter=max_iter)
            if fm * fa >= 0:
                xa = xm
            if fm == 0 or fabs(dm) < tol:
                break
        # if i == max_iter - 1:
        #     print('el_chem_pot_t iters=', i, 'T=', temperature, 'K')
        return xm

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] el_chem_pot(self, double[:] temperature, int max_iter=100):
        cdef:
            Py_ssize_t n = len(temperature)
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.el_chem_pot_t(temperature[i], max_iter=max_iter)
        return result

    def __str__(self):
        return 'Semiconductor: ' + self.__label
