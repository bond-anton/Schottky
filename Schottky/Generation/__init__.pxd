from BDFunction1D cimport Function
from BDFunction1D.Standard cimport Constant, Zero


cdef class ConstantGenerationFunction(Constant):
    pass

cdef class ZeroGenerationFunction(Zero):
    pass

cdef class ExponentialGenerationFunction(Function):
    cdef:
        double __rate
        double __absorption_length
