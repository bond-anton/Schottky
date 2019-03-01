import warnings
import numpy as np

from libc.math cimport M_PI

from BDSpace.Coordinates.transforms cimport unit_vector, vector_norm
from BDSpace.Field.Field cimport ConstantVectorConservativeField, HyperbolicPotentialSphericalConservativeField
from BDSpace.Field.CurveField cimport CurveField, HyperbolicPotentialCurveConservativeField
from BDSpace.Field.SuperposedField cimport SuperposedField
from Schottky.Constants cimport constant


cdef class ExternalField(ConstantVectorConservativeField):

    def __init__(self, str name, double[:] direction, double magnitude):
        cdef:
            double[:] potential = unit_vector(direction)
            int i, s = direction.shape[0]
            str field_type = 'Electric Field'
        for i in range(s):
            potential[i] *= magnitude
        super(ExternalField, self).__init__(name, field_type, potential)

    @property
    def magnitude(self):
        return vector_norm(self.__potential)

    @magnitude.setter
    def magnitude(self, double magnitude):
        cdef:
            double[:] potential = unit_vector(self.__potential)
            int i, s = self.__potential.shape[0]
        for i in range(s):
            potential[i] *= magnitude
        self.__potential = potential

    @property
    def direction(self):
        return unit_vector(self.__potential)

    @direction.setter
    def direction(self, double[:] direction):
        cdef:
            double magnitude = vector_norm(self.__potential)
            double[:] potential = unit_vector(direction)
            int i, s = direction.shape[0]
        for i in range(s):
            potential[i] *= magnitude
        self.__potential = potential


cdef class PointChargeCoulombPotential(HyperbolicPotentialSphericalConservativeField):

    def __init__(self, str name, double q, double r, double epsilon=1.0):
        cdef:
            double a = q / (4 * M_PI * constant.epsilon_0 * epsilon)
            str field_type = 'Electric Field'
        super(PointChargeCoulombPotential, self).__init__(name, field_type, a, r)
        if epsilon <= 0:
            warnings.warn('Epsilon must be grater than zero. Setting epsilon to 1.0 (Vacuum)')
            self.__epsilon = 1.0
        else:
            self.__epsilon = epsilon
        self.__q = q

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
