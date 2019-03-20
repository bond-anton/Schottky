from cpython.array cimport array, clone
from BDSpace.Field.Field cimport Field, ConstantScalarConservativeField
from BDSpace.Field.SuperposedField cimport SuperposedField
from BDSpace.Field.SphericallySymmetric cimport SphericallySymmetric, HyperbolicPotentialSphericalConservativeField
from Schottky.Potential.ExternalField cimport ExternalField
from Schottky.Trap cimport Trap


cdef class TrapPotential(SuperposedField):

    def __init__(self, str name, Trap trap, Field trap_field, ExternalField external_field=None):
        cdef:
            array[double] direction
        self.__trap = trap
        self.__trap_field = trap_field
        if external_field is None:
            direction = clone(array('d'), 3, zero=False)
            direction[0] = 0.0
            direction[1] = 0.0
            direction[2] = 1.0
            self.__external_field = ExternalField(name='External Field',
                                                  direction=direction,
                                                  magnitude=0.0)
        else:
            self.__external_field = external_field
        super(TrapPotential, self).__init__(name, [self.__trap_field, self.__external_field])

    cpdef double emission_rate_enhancement(self):
        return 1.0


cdef class NullPotential(TrapPotential):

    def __init__(self, str name, Trap trap, ExternalField external_field=None):
        super(NullPotential, self).__init__(name, trap,
                                            ConstantScalarConservativeField(name='Null Potential',
                                                                            type='Electric Field',
                                                                            potential=0.0),
                                            external_field)

    cpdef double emission_rate_enhancement(self):
        return 1.0


cdef class PointLikeInExternalField(TrapPotential):

    def __init__(self, str name, Trap trap, SphericallySymmetric point_like, ExternalField external_field=None,
                 double r_min=0.0, double r_max=1.0e15):
        self.__r_min = r_min
        self.__r_max = r_max
        super(PointLikeInExternalField, self).__init__(name, trap, point_like, external_field)


    cpdef double emission_rate_enhancement(self):
        return 1.0
