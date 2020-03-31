from libc.math cimport exp

from cython cimport boundscheck, wraparound

from BDFunction1D cimport Function
from BDFunction1D.Standard cimport Constant, Zero


cdef class ConstantGenerationFunction(Constant):

    def __init__(self, double rate):
        super(ConstantGenerationFunction, self).__init__(rate)

    @property
    def rate(self):
        return self.c

    @rate.setter
    def rate(self, double rate):
        self.c = rate


cdef class ZeroGenerationFunction(Zero):

    @property
    def rate(self):
        return 0.0


cdef class ExponentialGenerationFunction(Function):

    def __init__(self, double rate, double absorption_length):
        super(ExponentialGenerationFunction, self).__init__()
        self.__rate = rate
        self.__absorption_length = absorption_length

    @property
    def rate(self):
        return self.__rate

    @rate.setter
    def rate(self, double rate):
        self.__rate = rate

    @property
    def absorption_length(self):
        return self.__absorption_length

    @absorption_length.setter
    def absorption_length(self, double absorption_length):
        self.__absorption_length = absorption_length

    @boundscheck(False)
    @wraparound(False)
    cpdef double evaluate_point(self, double x):
        return self.__rate * exp(-x / self.__absorption_length)
