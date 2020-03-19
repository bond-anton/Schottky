from BDPoisson1D.Function cimport Function


cdef class ConstantGenerationFunction(Function):
    cdef:
        double __rate

cdef class ZeroGenerationFunction(Function):
    pass
