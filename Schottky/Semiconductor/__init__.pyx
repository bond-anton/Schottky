from libc.math cimport sqrt, exp, log, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDFunction1D.Standard cimport Zero

from Schottky.Constants cimport constant
from Schottky.Dopant cimport Dopant
from Schottky.Trap cimport Trap
from Schottky.Helpers.Cache cimport Cache, hash_list


cdef class Semiconductor(object):

    def __init__(self, str label, dict reference, list dopants=None):
        self.__label = label
        self.__reference = reference
        self.__dopants = []
        self.__band_gap_cache = Cache('Band Gap Cache')
        self.__trap_eq_occupation_cache = Cache('Trap Eq Occupation Cache')
        self.__bulk_charge_cache = Cache('Bulk Charge Cache')
        self.__el_chem_pot_cache = Cache('El-Chem Potential Cache')
        self.__dopants_n = Zero()
        if dopants is not None:
            for dopant in dopants:
                if isinstance(dopant, Dopant):
                    self.__dopants.append(dopant)
                    self.__dopants_n += dopant.concentration

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
        self.__dopants_n = Zero()
        for dopant in dopants:
            if isinstance(dopant, Dopant):
                self.__dopants.append(dopant)
                self.__dopants_n += dopant.concentration

    @boundscheck(False)
    @wraparound(False)
    cdef double __band_gap_t(self, double temperature):
        return (self.__reference['band_gap']['E_g0'] -
                self.__reference['band_gap']['alpha'] * temperature ** 2 /
                (temperature + self.__reference['band_gap']['beta']))

    @boundscheck(False)
    @wraparound(False)
    cpdef double band_gap_t(self, double temperature):
        cdef:
            long key = hash(temperature)
        try:
            return self.__band_gap_cache[key]
        except KeyError:
            self.__band_gap_cache[key] = self.__band_gap_t(temperature)
        return self.__band_gap_cache[key]

    @boundscheck(False)
    @wraparound(False)
    cpdef double band_gap_ev_t(self, double temperature):
        cdef:
            long key = hash(temperature)
        try:
            return constant.joule_to_ev_point(self.__band_gap_cache[key])
        except KeyError:
            self.__band_gap_cache[key] = self.__band_gap_t(temperature)
        return constant.joule_to_ev_point(self.__band_gap_cache[key])

    @boundscheck(False)
    @wraparound(False)
    cpdef double band_gap_boltzmann_t(self, double temperature):
        cdef:
            long key = hash(temperature)
        try:
            return constant.joule_to_boltzmann_point(self.__band_gap_cache[key], temperature)
        except KeyError:
            self.__band_gap_cache[key] = self.__band_gap_t(temperature)
        return constant.joule_to_boltzmann_point(self.__band_gap_cache[key], temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] band_gap(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.band_gap_t(temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] band_gap_ev(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.band_gap_ev_t(temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] band_gap_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.band_gap_boltzmann_t(temperature[i])
        return result

    cpdef double n_c_t(self, double temperature):
        return self.__reference['N_c0'] * temperature ** (3 / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_c(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
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
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_v_t(temperature[i])
        return result

    cpdef double pn_t(self, double temperature):
        return self.__reference['N_c0'] * self.__reference['N_v0']\
               * exp(-self.band_gap_boltzmann_t(temperature))\
               * temperature ** 3

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] pn(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.pn_t(temperature[i])
        return result

    cpdef double n_i_t(self, double temperature):
        return sqrt(self.__reference['N_c0'] * self.__reference['N_v0'])\
               * exp(-self.band_gap_boltzmann_t(temperature) / 2)\
               * temperature ** 1.5

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_i(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_i_t(temperature[i])
        return result

    cpdef double one_by_n_i_t(self, double temperature):
        return sqrt(self.__reference['N_v0'] * self.__reference['N_c0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2)\
               * temperature ** (-1.5)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] one_by_n_i(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.one_by_n_i_t(temperature[i])
        return result

    cpdef double e_i_t(self, double temperature):
        return 0.5 * (
                self.band_gap_t(temperature) + constant.k * temperature
                * log(self.__reference['N_c0'] / self.__reference['N_v0'])
        )

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] e_i(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.e_i_t(temperature[i])
        return result

    cpdef double e_i_ev_t(self, double temperature):
        return 0.5 * (
                self.band_gap_ev_t(temperature) + constant.thermal_voltage_t(temperature)
                * log(self.__reference['N_c0'] / self.__reference['N_v0'])
        )

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] e_i_ev(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.e_i_ev_t(temperature[i])
        return result

    cpdef double e_i_boltzmann_t(self, double temperature):
        return 0.5 * (self.band_gap_boltzmann_t(temperature) + log(self.__reference['N_c0'] / self.__reference['N_v0']))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] e_i_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.e_i_boltzmann_t(temperature[i])
        return result

    cpdef double n_c_n_i_t(self, double temperature):
        return sqrt(self.__reference['N_c0'] / self.__reference['N_v0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_c_n_i(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_c_n_i_t(temperature[i])
        return result

    cpdef double n_v_n_i_t(self, double temperature):
        return sqrt(self.__reference['N_v0'] / self.__reference['N_c0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_v_n_i(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_v_n_i_t(temperature[i])
        return result

    cpdef double n_e_t(self, double mu, double temperature):
        return self.n_c_t(temperature) * exp(-constant.joule_to_boltzmann_point(mu, temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_t(mu[i], temperature[i])
        return result

    cpdef double n_e_ev_t(self, double mu_ev, double temperature):
        return self.n_c_t(temperature) * exp(-constant.ev_to_boltzmann_point(mu_ev, temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_ev(self, double[:] mu_ev, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_ev_t(mu_ev[i], temperature[i])
        return result

    cpdef double n_e_boltzmann_t(self, double mu, double temperature):
        return self.n_c_t(temperature) * exp(-mu)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_boltzmann(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_boltzmann_t(mu[i], temperature[i])
        return result

    cpdef double n_e_n_i_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_c0'] / self.__reference['N_v0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2 - constant.joule_to_boltzmann_point(mu, temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_n_i(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_n_i_t(mu[i], temperature[i])
        return result

    cpdef double n_e_n_i_ev_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_c0'] / self.__reference['N_v0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2 - constant.ev_to_boltzmann_point(mu, temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_n_i_ev(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_n_i_ev_t(mu[i], temperature[i])
        return result

    cpdef double n_e_n_i_boltzmann_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_c0'] / self.__reference['N_v0'])\
               * exp(self.band_gap_boltzmann_t(temperature) / 2 - mu)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_n_i_boltzmann(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_e_n_i_boltzmann_t(mu[i], temperature[i])
        return result

    cpdef double n_e_to_mu_t(self, double n_e, double temperature):
        return constant.__k * temperature * log (self.n_c_t(temperature) / n_e)

    cpdef double n_e_to_mu_ev_t(self, double n_e, double temperature):
        return constant.__k * temperature * log (self.n_c_t(temperature) / n_e) / constant.__q

    cpdef double n_e_to_mu_boltzmann_t(self, double n_e, double temperature):
        return log (self.n_c_t(temperature) / n_e)

    cpdef double n_h_t(self, double mu, double temperature):
        return self.n_v_t(temperature) * exp(
            constant.joule_to_boltzmann_point(mu, temperature) - self.band_gap_boltzmann_t(temperature)
        )

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_t(mu[i], temperature[i])
        return result

    cpdef double n_h_ev_t(self, double mu_ev, double temperature):
        return self.n_v_t(temperature) * exp(
            constant.ev_to_boltzmann_point(mu_ev, temperature) - self.band_gap_boltzmann_t(temperature)
        )

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_ev(self, double[:] mu_ev, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_ev_t(mu_ev[i], temperature[i])
        return result

    cpdef double n_h_boltzmann_t(self, double mu, double temperature):
        return self.n_v_t(temperature) * exp(mu - self.band_gap_boltzmann_t(temperature))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_boltzmann(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_boltzmann_t(mu[i], temperature[i])
        return result

    cpdef double n_h_n_i_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_v0'] / self.__reference['N_c0'])\
               * exp(constant.joule_to_boltzmann_point(mu, temperature) - self.band_gap_boltzmann_t(temperature) / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_n_i(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_n_i_t(mu[i], temperature[i])
        return result

    cpdef double n_h_n_i_ev_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_v0'] / self.__reference['N_c0'])\
               * exp(constant.ev_to_boltzmann_point(mu, temperature) - self.band_gap_boltzmann_t(temperature) / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_n_i_ev(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_n_i_ev_t(mu[i], temperature[i])
        return result

    cpdef double n_h_n_i_boltzmann_t(self, double mu, double temperature):
        return sqrt(self.__reference['N_v0'] / self.__reference['N_c0'])\
               * exp(mu - self.band_gap_boltzmann_t(temperature) / 2)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_n_i_boltzmann(self, double[:] mu, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n_h_n_i_boltzmann_t(mu[i], temperature[i])
        return result

    cpdef double n_h_to_mu_t(self, double n_h, double temperature):
        return constant.__k * temperature * log (n_h / self.n_v_t(temperature)) + self.band_gap_t(temperature)

    cpdef double n_h_to_mu_ev_t(self, double n_h, double temperature):
        return constant.__k * temperature * log (n_h / self.n_v_t(temperature)) / constant.__q \
               + self.band_gap_ev_t(temperature)

    cpdef double n_h_to_mu_boltzmann_t(self, double n_h, double temperature):
        return log (n_h / self.n_v_t(temperature)) + self.band_gap_boltzmann_t(temperature)

    cpdef double v_e_t(self, double temperature):
        return sqrt(3 * constant.k * temperature /
                    (self.__reference['carrier_mass']['m_e_coeff'] * constant.m_e))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] v_e(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
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
    cpdef double trap_eq_occupation(self, Trap trap, double temperature, double mu_e, double delta_mu_h_e=0.0,
                                    double f_threshold=1.0e-23, int max_iter=100,
                                    bint verbose=False):
        cdef:
            int i = 0
            double fa, fb, fm=0.0, ff_a, ff_b, ff_m, dm
            double band_gap, n_c, n_v, n_e, n_h
            double mu_h = mu_e + delta_mu_h_e
            long key = hash_list([trap, mu_e, mu_h, temperature, f_threshold, max_iter])
        try:
            return self.__trap_eq_occupation_cache[key]
        except KeyError:
            band_gap = self.__band_gap_t(temperature)
            n_c = self.n_c_t(temperature)
            n_v = self.n_v_t(temperature)
            n_e = self.n_e_t(mu_e, temperature)
            n_h = self.n_h_t(mu_h, temperature)
            if trap.__cb_bound:
                trap.energy_v = band_gap - trap.energy_c
            else:
                trap.energy_c = band_gap - trap.energy_v
            fa = 0.0
            fb = 1.0
            ff_a = fa - trap.f_eq(temperature, n_e, n_c, n_h, n_v, fa, verbose=verbose)
            ff_b = fb - trap.f_eq(temperature, n_e, n_c, n_h, n_v, fb, verbose=verbose)
            if ff_a == 0:
                fm = fa
            elif ff_b == 0:
                fm = fb
            else:
                dm = fb - fa
                for i in range(max_iter):
                    dm *= 0.5
                    fm = fa + dm
                    ff_m = fm - trap.f_eq(temperature, n_e, n_c, n_h, n_v, fm, verbose=verbose)
                    if ff_m * ff_a >= 0:
                        fa = fm
                    if ff_m == 0 or fabs(dm) < f_threshold:
                        break
            if verbose and i == max_iter - 1:
                print('    trap_eq_occupation (%s)' % trap.__label, 'reached max iters=', i + 1, 'T=', temperature, 'K')
            self.__trap_eq_occupation_cache[key] = fm
        return self.__trap_eq_occupation_cache[key]

    @boundscheck(False)
    @wraparound(False)
    cpdef double bulk_charge(self, double mu, double temperature, double z=1.0e5,
                             double f_threshold=1.0e-23, int max_iter=100,
                             bint verbose=False):
        cdef:
            double fm, result
            long key = hash_list([mu, temperature, z, f_threshold, max_iter])
            Dopant dopant
        try:
            return self.__bulk_charge_cache[key]
        except KeyError:
            result = self.n_h_t(mu, temperature) - self.n_e_t(mu, temperature)
            for dopant in self.__dopants:
                fm = self.trap_eq_occupation(dopant, temperature, mu,
                                             f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
                result += dopant.__concentration.evaluate_point(z) * (
                        (dopant.charge_state[1] - dopant.charge_state[0]) * fm
                        + dopant.charge_state[0])
            self.__bulk_charge_cache[key] = result
        return self.__bulk_charge_cache[key]

    @boundscheck(False)
    @wraparound(False)
    cpdef double el_chem_pot_t(self, double temperature,
                               double f_threshold=1.0e-23, int max_iter=100,
                               bint verbose=False):
        cdef:
            int i = 0
            double dm, xm=0.0, fm, fa, fb
            double tol, xtol = constant.__k * temperature * 0.001, rtol = 4.5e-16
            double xa = 0.0, xb = self.__band_gap_t(temperature)
            long key = hash_list([temperature, f_threshold, max_iter])
        try:
            return self.__el_chem_pot_cache[key]
        except KeyError:
            tol = xtol + rtol*(fabs(xa) + fabs(xb))
            fa = self.bulk_charge(xa, temperature, f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
            fb = self.bulk_charge(xb, temperature, f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
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
                fm = self.bulk_charge(xm, temperature, f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
                if fm * fa >= 0:
                    xa = xm
                if fm == 0 or fabs(dm) < tol:
                    break
            if verbose and i == max_iter - 1:
                print('el_chem_pot_t reached max iters=', i + 1, 'T=', temperature, 'K')
            self.__el_chem_pot_cache[key] = xm
        return self.__el_chem_pot_cache[key]

    @boundscheck(False)
    @wraparound(False)
    cpdef double el_chem_pot_ev_t(self, double temperature,
                                  double f_threshold=1.0e-23, int max_iter=100,
                                  bint verbose=False):
        return constant.joule_to_ev_point(self.el_chem_pot_t(temperature,
                                                             f_threshold=f_threshold,
                                                             max_iter=max_iter,
                                                             verbose=verbose))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] el_chem_pot(self, double[:] temperature,
                                double f_threshold=1.0e-23, int max_iter=100,
                                bint verbose=False):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.el_chem_pot_t(temperature[i],
                                           f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] el_chem_pot_ev(self, double[:] temperature,
                                   double f_threshold=1.0e-23, int max_iter=100,
                                   bint verbose=False):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.el_chem_pot_ev_t(temperature[i],
                                              f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double work_function_t(self, double temperature,
                                 double f_threshold=1.0e-23, int max_iter=100,
                                 bint verbose=False):
        return self.el_chem_pot_t(temperature,
                                  f_threshold=f_threshold,
                                  max_iter=max_iter,
                                  verbose=verbose) + self.__reference['affinity']

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] work_function(self, double[:] temperature,
                                  double f_threshold=1.0e-23, int max_iter=100,
                                  bint verbose=False):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.work_function_t(temperature[i],
                                             f_threshold=f_threshold,
                                             max_iter=max_iter, verbose=verbose)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double work_function_ev_t(self, double temperature,
                                    double f_threshold=1.0e-23, int max_iter=100,
                                    bint verbose=False):
        return constant.joule_to_ev_point(self.el_chem_pot_t(temperature,
                                                             f_threshold=f_threshold,
                                                             max_iter=max_iter,
                                                             verbose=verbose) + self.__reference['affinity'])

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] work_function_ev(self, double[:] temperature,
                                     double f_threshold=1.0e-23, int max_iter=100,
                                     bint verbose=False):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.work_function_ev_t(temperature[i],
                                                f_threshold=f_threshold,
                                                max_iter=max_iter, verbose=verbose)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double work_function_boltzmann_t(self, double temperature,
                                           double f_threshold=1.0e-23, int max_iter=100,
                                           bint verbose=False):
        return constant.joule_to_boltzmann_point(self.el_chem_pot_t(temperature,
                                                                    f_threshold=f_threshold,
                                                                    max_iter=max_iter,
                                                                    verbose=verbose) + self.__reference['affinity'],
                                                 temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] work_function_boltzmann(self, double[:] temperature,
                                            double f_threshold=1.0e-23, int max_iter=100,
                                            bint verbose=False):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.work_function_boltzmann_t(temperature[i],
                                                       f_threshold=f_threshold,
                                                       max_iter=max_iter, verbose=verbose)
        return result

    cpdef double mobility_e_point_t(self, double dopants_n, double field, double pn, double temperature):
        cdef:
            double mu_l, mu_i, mu_ccs, mu, x, _pn, field_coeff
        _pn = pn * 1e-12
        mu_l = self.__reference['mobility_e']['mu_L0'] * (
                (temperature / 300.0) ** (-self.__reference['mobility_e']['alpha'])
        )
        mu_i = (self.__reference['mobility_e']['A'] * (temperature ** 1.5) / (dopants_n * 1e-6)) \
               / (log(1 + self.__reference['mobility_e']['B'] * (temperature ** 2) / (dopants_n * 1e-6))
                  - self.__reference['mobility_e']['B'] * (temperature ** 2)
                  / (self.__reference['mobility_e']['B'] * (temperature ** 2) + (dopants_n * 1e-6)))
        if pn > 0:
            mu_ccs = 2e17 * (temperature ** 1.5) / sqrt(_pn) / log(1 + 8.28e8 * (temperature ** 2) * (_pn ** (-1 / 3)))
            x = sqrt(6 * mu_l * (mu_i + mu_ccs) / (mu_i * mu_ccs))
        else:
            x = sqrt(6 * mu_l / mu_i)
        mu = mu_l * (1.025 / (1 + ((x / 1.68) ** 1.43)) - 0.025)
        field_coeff = (1 +
                       (mu * field * 1e-2 / self.__reference['mobility_e']['v_s'])
                       ** self.__reference['mobility_e']['beta']
                      ) ** (-1 / self.reference['mobility_e']['beta'])
        return mu * field_coeff * 1e-4

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] mobility_e_point(self, double dopants_n, double field, double pn, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.mobility_e_point_t(dopants_n, field, pn, temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] mobility_e_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature):
        cdef:
            int n = dopants_n.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.mobility_e_point_t(dopants_n[i], field[i], pn[i], temperature)
        return result

    cpdef double mobility_h_point_t(self, double dopants_n, double field, double pn, double temperature):
        cdef:
            double mu_l, mu_i, mu_ccs, mu, x, _pn, field_coeff
        _pn = pn * 1e-12
        mu_l = self.__reference['mobility_h']['mu_L0'] * (
                (temperature / 300.0) ** (-self.__reference['mobility_h']['alpha'])
        )
        mu_i = (self.__reference['mobility_h']['A'] * (temperature ** 1.5) / (dopants_n * 1e-6)) \
               / (log(1 + self.__reference['mobility_h']['B'] * (temperature ** 2) / (dopants_n * 1e-6))
                  - self.__reference['mobility_h']['B'] * (temperature ** 2)
                  / (self.__reference['mobility_h']['B'] * (temperature ** 2) + (dopants_n * 1e-6)))
        if pn > 0:
            mu_ccs = 2e17 * (temperature ** 1.5) / sqrt(_pn) / log(1 + 8.28e8 * (temperature ** 2) * (_pn ** (-1 / 3)))
            x = sqrt(6 * mu_l * (mu_i + mu_ccs) / (mu_i * mu_ccs))
        else:
            x = sqrt(6 * mu_l / mu_i)
        mu = mu_l * (1.025 / (1 + ((x / 1.68) ** 1.43)) - 0.025)
        field_coeff = (1 +
                       (mu * field * 1e-2 / self.__reference['mobility_h']['v_s'])
                       ** self.__reference['mobility_h']['beta']
                      ) ** (-1 / self.reference['mobility_h']['beta'])
        return mu * field_coeff * 1e-4

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] mobility_h_point(self, double dopants_n, double field, double pn, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.mobility_h_point_t(dopants_n, field, pn, temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] mobility_h_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature):
        cdef:
            int n = dopants_n.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.mobility_h_point_t(dopants_n[i], field[i], pn[i], temperature)
        return result

    cpdef double diffusivity_e_point_t(self, double dopants_n, double field, double pn, double temperature):
        return constant.thermal_voltage_t(temperature) * self.mobility_e_point_t(dopants_n, field, pn, temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] diffusivity_e_point(self, double dopants_n, double field, double pn, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.diffusivity_e_point_t(dopants_n, field, pn, temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] diffusivity_e_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature):
        cdef:
            int n = dopants_n.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.diffusivity_e_point_t(dopants_n[i], field[i], pn[i], temperature)
        return result

    cpdef double diffusivity_h_point_t(self, double dopants_n, double field, double pn, double temperature):
        return constant.thermal_voltage_t(temperature) * self.mobility_h_point_t(dopants_n, field, pn, temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] diffusivity_h_point(self, double dopants_n, double field, double pn, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.diffusivity_h_point_t(dopants_n, field, pn, temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] diffusivity_h_t(self, double[:] dopants_n, double[:] field, double[:] pn, double temperature):
        cdef:
            int n = dopants_n.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.diffusivity_h_point_t(dopants_n[i], field[i], pn[i], temperature)
        return result

    cpdef double debye_length_point_t(self, double dopants_n, double temperature):
        return sqrt(constant.__epsilon_0 * self.__reference['epsilon'] * constant.k * temperature / dopants_n)\
               / constant.__q

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] debye_length_point(self, double dopants_n, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.debye_length_point_t(dopants_n, temperature[i])
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] debye_length_t(self, double[:] dopants_n, double temperature):
        cdef:
            int n = dopants_n.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.debye_length_point_t(dopants_n[i], temperature)
        return result

    def __str__(self):
        return 'Semiconductor: ' + self.__label

    def cache_info(self):
        self.__band_gap_cache.info()
        self.__trap_eq_occupation_cache.info()
        self.__bulk_charge_cache.info()
        self.__el_chem_pot_cache.info()
