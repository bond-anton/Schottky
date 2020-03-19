from libc.math cimport exp

from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDPoisson1D.Function cimport Function, Functional


cdef class ConstantGenerationFunction(Function):

    def __init__(self, double rate):
        super(ConstantGenerationFunction, self).__init__()
        self.__rate = rate

    @property
    def rate(self):
        return self.__rate

    @rate.setter
    def rate(self, double rate):
        self.__rate = rate

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] evaluate(self, double[:] x):
        cdef:
            array[double] result
            int i, n = x.shape[0]
        result = clone(array('d'), n, zero=False)
        for i in range(n):
            result[i] = self.__rate


cdef class ZeroGenerationFunction(Function):

    @property
    def rate(self):
        return 0.0

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] evaluate(self, double[:] x):
        return clone(array('d'), x.shape[0], zero=True)


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
    cpdef double[:] evaluate(self, double[:] x):
        cdef:
            array[double] result
            int i, n = x.shape[0]
        result = clone(array('d'), n, zero=False)
        for i in range(n):
            result[i] = self.__rate * exp(-x[i] / self.__absorption_length)
