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
        return self.e_c(temperature, v_e) * n_e * (1 - self.__f)

    cpdef double h_cr(self, double temperature, double v_h, double n_h):
        return self.h_c(temperature, v_h) * n_h * self.__f

    @property
    def f(self):
        return self.__f

    def __str__(self):
        return 'Trap: %s\nEc-Et: %2.2f eV (%2.2g J)\nEt-Ev: %2.2f eV (%2.2g J)' % (self.label,
                                                                                   self.energy_c_ev,
                                                                                   self.energy_c,
                                                                                   self.energy_v_ev,
                                                                                   self.energy_v)
