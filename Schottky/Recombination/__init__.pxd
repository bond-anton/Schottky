from BDPoisson1D.Function cimport Function


cdef class ConstantRecombinationFunction(Function):
    cdef:
        double __rate

cdef class ZerorecombinationFunction(Function):
    pass

cdef class ExponentialGenerationFunction(Function):
    cdef:
        double __rate
        double __absorption_length
