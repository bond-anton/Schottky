from BDSpace.Field.Field cimport Field
from BDSpace.Field.SuperposedField cimport SuperposedField
from Schottky.Potential.ExternalField cimport ExternalField
from Schottky.Trap cimport Trap


cdef class TrapPotential(SuperposedField):
    cdef:
        Trap __trap
        Field __trap_field
        ExternalField __external_field
    cpdef double emission_rate_enhancement(self)


cdef class NullPotential(TrapPotential):
    pass


cdef class PointLikeInExternalField(TrapPotential):
    cdef:
        double __r_min
        double __r_max