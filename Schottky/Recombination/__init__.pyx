import numpy as np

from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDPoisson1D.Function cimport Function, Functional


cdef class ConstantRecombinationFunction(Function):

    def __init__(self, double rate):
        super(ConstantRecombinationFunction, self).__init__()
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


cdef class ZeroRecombinationFunction(Function):

    @property
    def rate(self):
        return 0.0

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] evaluate(self, double[:] x):
        return clone(array('d'), x.shape[0], zero=True)
