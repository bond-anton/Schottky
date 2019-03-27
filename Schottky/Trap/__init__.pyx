from libc.math cimport fabs, exp

from Schottky.Potential.TrapPotential cimport TrapPotential, NullPotential
from Schottky.Constants cimport constant


cdef class Trap(object):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label, bint conduction_band_bound,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0,
                 TrapPotential e_potential=None, TrapPotential h_potential=None):
        '''
        Constructor
        '''
        self.__label = label
        self.__cb_bound = conduction_band_bound
        self.__energy_c = energy_c
        self.__energy_v = energy_v
        self.__e_cs0 = e_cs0
        self.__h_cs0 = h_cs0
        self.__e_cs_activation = e_cs_activation
        self.__h_cs_activation = h_cs_activation
        self.__charge_state = {0: 0, 1: -1}
        self.__g = {0: 1, 1: 2}
        self.__capture_barrier = {0: 0.0, 1: 0.0}
        if e_potential is None:
            self.__e_potential = NullPotential('electron null potential', self)
        else:
            self.__e_potential = e_potential
        if h_potential is None:
            self.__h_potential = NullPotential('hole null potential', self)
        else:
            self.__h_potential = h_potential

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def cb_bound(self):
        return self.__cb_bound

    @cb_bound.setter
    def cb_bound(self, bint conduction_band_bound):
        self.__cb_bound = conduction_band_bound

    @property
    def vb_bound(self):
        return not self.__cb_bound

    @vb_bound.setter
    def vb_bound(self, bint valence_band_bound):
        self.__cb_bound = not valence_band_bound

    @property
    def energy_c(self):
        return self.__energy_c

    @energy_c.setter
    def energy_c(self, double energy_c):
        self.__energy_c = energy_c

    @property
    def energy_v(self):
        return self.__energy_v

    @energy_v.setter
    def energy_v(self, double energy_v):
        self.__energy_v = energy_v

    @property
    def energy_c_ev(self):
        return self.__energy_c / constant.__q

    @energy_c_ev.setter
    def energy_c_ev(self, double energy_c_ev):
        self.__energy_c = energy_c_ev * constant.__q

    @property
    def energy_v_ev(self):
        return self.__energy_v / constant.__q

    @energy_v_ev.setter
    def energy_v_ev(self, double energy_v_ev):
        self.__energy_v = energy_v_ev * constant.__q

    @property
    def e_cs0(self):
        return self.__e_cs0

    @e_cs0.setter
    def e_cs0(self, double e_cs0):
        self.__e_cs0 = e_cs0

    @property
    def h_cs0(self):
        return self.__h_cs0

    @h_cs0.setter
    def h_cs0(self, double h_cs0):
        self.__h_cs0 = h_cs0

    @property
    def e_cs_activation(self):
        return self.__e_cs_activation

    @e_cs_activation.setter
    def e_cs_activation(self, double e_cs_activation):
        self.__e_cs_activation = e_cs_activation

    @property
    def h_cs_activation(self):
        return self.__h_cs_activation

    @h_cs_activation.setter
    def h_cs_activation(self, double h_cs_activation):
        self.__h_cs_activation = h_cs_activation

    @property
    def e_cs_activation_ev(self):
        return self.__e_cs_activation / constant.__q

    @e_cs_activation_ev.setter
    def e_cs_activation_ev(self, double e_cs_activation_ev):
        self.__e_cs_activation = e_cs_activation_ev * constant.__q

    @property
    def h_cs_activation_ev(self):
        return self.__h_cs_activation / constant.__q

    @h_cs_activation_ev.setter
    def h_cs_activation_ev(self, double h_cs_activation_ev):
        self.__h_cs_activation = h_cs_activation_ev * constant.__q

    @property
    def capture_barrier(self):
        return self.__capture_barrier

    @capture_barrier.setter
    def capture_barrier(self, dict capture_barrier):
        self.__capture_barrier[0] = capture_barrier[0]
        self.__capture_barrier[1] = capture_barrier[1]

    @property
    def capture_barrier_ev(self):
        return {0: self.__capture_barrier[0] / constant.__q, 1: self.__capture_barrier[1] / constant.__q}

    @capture_barrier_ev.setter
    def capture_barrier_ev(self, dict capture_barrier_ev):
        self.__capture_barrier[0] = capture_barrier_ev[0] * constant.__q
        self.__capture_barrier[1] = capture_barrier_ev[1] * constant.__q

    @property
    def e_potential(self):
        return self.__e_potential

    @e_potential.setter
    def e_potential(self, TrapPotential e_potential):
        self.__e_potential = e_potential

    @property
    def h_potential(self):
        return self.__h_potential

    @h_potential.setter
    def h_potential(self, TrapPotential h_potential):
        self.__h_potential = h_potential

    cpdef double e_cs(self, double temperature):
        return self.__e_cs0 * exp(-self.__e_cs_activation / (constant.__k * temperature))

    cpdef double h_cs(self, double temperature):
        return self.__h_cs0 * exp(-self.__h_cs_activation / (constant.__k * temperature))

    cpdef double e_c(self, double temperature, double v_e):
        return self.e_cs(temperature) * v_e

    cpdef double h_c(self, double temperature, double v_h):
        return self.h_cs(temperature) * v_h

    cpdef double e_cr(self, double temperature, double v_e, double n_e, double f):
        return self.e_c(temperature, v_e) * n_e * exp(-self.__capture_barrier[0] * f / (constant.__k * temperature))

    cpdef double h_cr(self, double temperature, double v_h, double n_h, double f):
        return self.h_c(temperature, v_h) * n_h \
               * exp(-self.__capture_barrier[1] * (1 - f) / (constant.__k * temperature))

    cpdef double e_er(self, double temperature, double v_e, double n_c, double f):
        cdef:
            double cr, exp_f, g, enhancement
        e_c = self.e_c(temperature, v_e)
        g = self.__g[0] / self.__g[1]
        exp_f = n_c * exp(-self.__energy_c /(constant.__k * temperature))
        enhancement = self.__e_potential.emission_rate_enhancement(temperature, f)
        return e_c * g * exp_f * enhancement

    cpdef double h_er(self, double temperature, double v_h, double n_v, double f):
        cdef:
            double cr, exp_f, g, enhancement
        h_c = self.h_c(temperature, v_h)
        g = self.__g[1] / self.__g[0]
        exp_f = n_v * exp(-self.__energy_v /(constant.__k * temperature))
        enhancement = self.__h_potential.emission_rate_enhancement(temperature, f)
        return h_c * g * exp_f * enhancement

    cpdef double f_eq(self, double temperature,
                      double v_e, double n_e, double n_c,
                      double v_h, double n_h, double n_v,
                      double f,
                      bint verbose=False):
        cdef:
            double e_c, e_e, h_c, h_e
        e_c = self.e_cr(temperature, v_e, n_e, f)
        e_e = self.e_er(temperature, v_e, n_c, f)
        h_c = self.h_cr(temperature, v_h, n_h, f)
        h_e = self.h_er(temperature, v_h, n_v, f)
        if fabs(e_c + e_e) > fabs(h_c + h_e):
            return e_c / (e_c + e_e)
        else:
            if fabs(h_c + h_e) > 0.0:
                return h_e / (h_c + h_e)
            if verbose:
                print('%s:' % self.__label, 'f_eq problem @ T =', temperature, 'K')
            return 0.0

    cpdef double df_dt(self, double temperature,
                       double v_e, double n_e, double n_c,
                       double v_h, double n_h, double n_v,
                       double f):
        cdef:
            double gain, loss
        gain = (self.e_cr(temperature, v_e, n_e, f) + self.h_er(temperature, v_h, n_v, f)) * (1 - f)
        loss = (self.e_er(temperature, v_e, n_c, f) + self.h_cr(temperature, v_h, n_h, f)) * f
        return gain - loss

    @property
    def charge_state(self):
        return self.__charge_state

    @charge_state.setter
    def charge_state(self, charge_state):
        self.__charge_state[0] = charge_state[0]
        self.__charge_state[1] = charge_state[1]

    @property
    def g(self):
        return self.__g

    @g.setter
    def g(self, g):
        self.__g[0] = g[0]
        self.__g[1] = g[1]

    def __str__(self):
        description = 'Trap: %s\n' % self.label
        if self.cb_bound:
            description += 'Ec-Et: %2.2f eV (%2.2g J)' % (self.energy_c_ev, self.energy_c)
        else:
            description += 'Et-Ev: %2.2f eV (%2.2g J)' % (self.energy_v_ev, self.energy_v)
        return description
