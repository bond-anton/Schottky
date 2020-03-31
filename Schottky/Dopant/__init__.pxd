from BDFunction1D cimport Function
from Schottky.Trap cimport Trap

cdef class Dopant(Trap):
    cdef:
        Function __concentration
        str __color
        str __linestyle
        str __marker
