from BDSpace.Field.Field cimport ConstantVectorConservativeField


cdef class ExternalField(ConstantVectorConservativeField):
    cdef:
        double[:] __direction
