from __future__ import division, print_function

from libc.math cimport exp

from Schottky.Constants cimport constant


cdef class Trap(object):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0):
        '''
        Constructor
        '''
        self.__label = label
        self.__energy_c = energy_c
        self.__energy_v = energy_v
        self.__e_cs0 = e_cs0
        self.__h_cs0 = h_cs0
        self.__e_cs_activation = e_cs_activation
        self.__h_cs_activation = h_cs_activation
        self.__f = 0.0
        self.__charge_state = {0: 0, 1: -1}
        self.__g = {0: 1, 1: 2}

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

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

    cpdef double e_cs(self, double temperature):
        return self.__e_cs0 * exp(-self.__e_cs_activation / (constant.__k * temperature))

    cpdef double h_cs(self, double temperature):
        return self.__h_cs0 * exp(-self.__h_cs_activation / (constant.__k * temperature))

    cpdef double e_c(self, double temperature, double v_e):
        return self.e_cs(temperature) * v_e

    cpdef double h_c(self, double temperature, double v_h):
        return self.h_cs(temperature) * v_h

    cpdef double e_cr(self, double temperature, double v_e, double n_e):
        return self.e_c(temperature, v_e) * n_e

    cpdef double h_cr(self, double temperature, double v_h, double n_h):
        return self.h_c(temperature, v_h) * n_h

    cpdef double e_er(self, double temperature, double v_e, double n_c):
        cdef:
            double cr, exp_f, g
        e_c = self.e_c(temperature, v_e)
        g = self.__g[0] / self.__g[1]
        exp_f = n_c * exp(-self.__energy_c /(constant.__k * temperature))
        return e_c * g * exp_f

    cpdef double h_er(self, double temperature, double v_h, double n_v):
        cdef:
            double cr, exp_f, g
        h_c = self.h_c(temperature, v_h)
        g = self.__g[1] / self.__g[0]
        exp_f = n_v * exp(-self.__energy_v /(constant.__k * temperature))
        return h_c * g * exp_f

    cpdef f_eq(self, double temperature,
               double v_e, double n_e, double n_c,
               double v_h, double n_h, double n_v):
        cdef:
            double e_c, e_e, h_c, h_e
        e_c = self.e_cr(temperature, v_e, n_e)
        e_e = self.e_er(temperature, v_e, n_c)
        h_c = self.h_cr(temperature, v_h, n_h)
        h_e = self.h_er(temperature, v_h, n_v)
        if e_c + e_e > h_c + h_e:
            return e_c / (e_c + e_e)
        else:
            return h_e / (h_c + h_e)

    cpdef df_dt(self, double temperature,
                double v_e, double n_e, double n_c,
                double v_h, double n_h, double n_v):
        cdef:
            double gain, loss
        gain = (self.e_cr(temperature, v_e, n_e) + self.h_er(temperature, v_h, n_v)) * (1 - self.__f)
        loss = (self.e_er(temperature, v_e, n_c) + self.h_cr(temperature, v_h, n_h)) * self.__f
        return gain - loss

    cdef double __coerce_f(self, double f):
        if f > 1.0:
            return 1.0
        elif f < 0.0:
            return 0.0
        else:
            return f

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

    @property
    def f(self):
        return self.__f

    @f.setter
    def f(self, double f):
        self.__f = self.__coerce_f(f)

    def __str__(self):
        return 'Trap: %s\nEc-Et: %2.2f eV (%2.2g J)\nEt-Ev: %2.2f eV (%2.2g J)' % (self.label,
                                                                                   self.energy_c_ev,
                                                                                   self.energy_c,
                                                                                   self.energy_v_ev,
                                                                                   self.energy_v)
