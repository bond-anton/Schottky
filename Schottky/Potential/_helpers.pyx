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


cdef long hash_list(list seq):
    cdef long result = 0x345678
    cdef long mult = 1000003
    cdef long h
    cdef long l = 0
    try:
        l = len(seq)
    except TypeError:
        # NOTE: This probably means very short non-len-able sequences
        # will not be spread as well as they should, but I'm not
        # sure what else to do.
        l = 100
    for element in seq:
        try:
            h = hash(element)
        except TypeError:
            h = hash_list(element)
        result ^= h
        result *= mult
        mult += 82520 + l + l
    result += 97531
    return result ^ hash(type(seq))
