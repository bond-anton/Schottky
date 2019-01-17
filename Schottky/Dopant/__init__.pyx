from __future__ import division, print_function

from libc.math cimport exp

from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from Schottky.Trap cimport Trap


cdef class Dopant(Trap):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label, TreeMesh1DUniform concentration, TreeMesh1DUniform f,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0):
        '''
        Constructor
        '''
        self.__concentration = concentration
        self.__f = f
        super(Dopant, self).__init__(energy_c, energy_v,
                                     e_cs0, h_cs0,
                                     e_cs_activation, h_cs_activation)

    @property
    def concentration(self):
        return self.__concentration

    @concentration.setter
    def concentration(self, TreeMesh1DUniform concentration):
        self.__concentration = concentration

    @property
    def occupation(self):
        return self.__f

    @occupation.setter
    def occupation(self, TreeMesh1DUniform f):
        self.__f = f

    cpdef double n_t(self, double z):
        return 1.0

    cpdef double f(self, double z):
        return 1.0

    def __str__(self):
        description = 'Dopant: %s\n' % self.label
        description += super(Dopant, self).__str__()
        return description
