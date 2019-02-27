from BDSpace.Field.Field cimport ConstantVectorConservativeField, HyperbolicPotentialSphericalConservativeField


cdef class ExternalField(ConstantVectorConservativeField):
    pass


cdef class PointChargeCoulombPotential(HyperbolicPotentialSphericalConservativeField):

    cdef:
        double __epsilon
        double __q
