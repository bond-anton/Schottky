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


@boundscheck(False)
@wraparound(False)
cdef double[:] gradient1d(double[:] y, double[:] x):
    cdef:
        int i, n = x.shape[0]
        double a, b, c, dx1, dx2
        array[double] result, template = array('d')
    result = clone(template, n, zero=False)
    dx1 = x[1] - x[0]
    dx2 = x[2] - x[1]
    a = -(2. * dx1 + dx2)/(dx1 * (dx1 + dx2))
    b = (dx1 + dx2) / (dx1 * dx2)
    c = - dx1 / (dx2 * (dx1 + dx2))
    result[0] = a * y[0] + b * y[1] + c * y[2]
    dx1 = x[n - 2] - x[n - 3]
    dx2 = x[n - 1] - x[n - 2]
    a = dx2 / (dx1 * (dx1 + dx2))
    b = - (dx2 + dx1) / (dx1 * dx2)
    c = (2.0 * dx2 + dx1) / (dx2 * (dx1 + dx2))
    result[n - 1] = a * y[n - 3] + b * y[n - 2] + c * y[n - 1]
    for i in range(1, n - 1):
        result[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    return result
