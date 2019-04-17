import warnings

from BDSpace.Field.SphericallySymmetric cimport HyperbolicPotentialSphericalConservativeField
from Schottky.Constants cimport constant
from libc.math cimport M_PI


cdef class PointChargeCoulombPotential(HyperbolicPotentialSphericalConservativeField):

    def __init__(self, str name, double q, double r, double epsilon=1.0):
        cdef:
            double a
            str field_type = 'Electric Field'
        if epsilon <= 0:
            warnings.warn('Epsilon must be grater than zero. Setting epsilon to 1.0 (Vacuum)')
            self.__epsilon = 1.0
        else:
            self.__epsilon = epsilon
        self.__q = q
        a = q / (4 * M_PI * constant.epsilon_0 * epsilon)
        super(PointChargeCoulombPotential, self).__init__(name, field_type, a, r)

    @property
    def epsilon(self):
        return self.__epsilon

    @epsilon.setter
    def epsilon(self, double epsilon):
        if epsilon <= 0:
            warnings.warn('Epsilon must be grater than zero. Setting epsilon to 1.0 (Vacuum)')
            self.__epsilon = 1.0
        else:
            self.__epsilon = epsilon

    @property
    def charge(self):
        return self.__q

    @charge.setter
    def charge(self, double q):
        self.__q = q
        self.__a = q / (4 * M_PI * constant.epsilon_0 * self.__epsilon)
