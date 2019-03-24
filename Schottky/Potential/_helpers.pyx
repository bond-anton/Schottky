from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone
from cython.parallel import prange


@boundscheck(False)
@wraparound(False)
cdef double trapz_1d(double[:] y, double[:] x) nogil:
    cdef:
        int nx = x.shape[0], i
        double result = 0.0
    for i in prange(nx - 1):
        result += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2
    return result


@boundscheck(False)
@wraparound(False)
cdef double[:] linspace(double start, double stop, int num):
    cdef:
        int i
        double step = (stop - start) / (num - 1)
        array[double] result, template = array('d')
    result = clone(template, num, zero=False)
    with nogil:
        for i in prange(num - 1):
            result[i] = start + i * step
    result[num - 1] = stop
    return result
