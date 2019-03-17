from BDSpace.Field.SphericallySymmetric cimport SphericallySymmetric, HyperbolicPotentialSphericalConservativeField
from BDSpace.Field.SuperposedField cimport SuperposedField
from Schottky.Potential.ExternalField cimport ExternalField


cdef class PointLikeInExternalField(SuperposedField):
    cdef:
        SphericallySymmetric __point_like
        ExternalField __external_field
        double __r_min
        double __r_max


cdef class PointChargeCoulombPotential(HyperbolicPotentialSphericalConservativeField):

    cdef:
        double __epsilon
        double __q
