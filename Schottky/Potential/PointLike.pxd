from BDSpace.Field.SphericallySymmetric cimport HyperbolicPotentialSphericalConservativeField


cdef class PointChargeCoulombPotential(HyperbolicPotentialSphericalConservativeField):

    cdef:
        double __epsilon
        double __q
