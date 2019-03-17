import warnings

from BDSpace.Field.SphericallySymmetric cimport SphericallySymmetric, HyperbolicPotentialSphericalConservativeField
from BDSpace.Field.SuperposedField cimport SuperposedField
from Schottky.Potential.ExternalField cimport ExternalField
from Schottky.Constants cimport constant
from libc.math cimport M_PI


cdef class PointLikeInExternalField(SuperposedField):

    def __init__(self, str name, SphericallySymmetric point_like, ExternalField external_field,
                 double r_min, double r_max):
        self.__point_like = point_like
        self.__external_field = external_field
        self.__r_min = r_min
        self.__r_max = r_max
        super(PointLikeInExternalField, self).__init__(name, [point_like, external_field])


    # cdef scalar_potential_polar(self, r, theta, phi):
        




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

