from BDSpace.Field.Field cimport Field


cdef class ElectricField(Field):
    pass

cdef class ConstantPotentialField(ElectricField):
    cdef:
        double __potential
