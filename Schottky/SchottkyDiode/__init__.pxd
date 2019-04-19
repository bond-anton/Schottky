from Schottky.Metal cimport Metal
from Schottky.Semiconductor cimport Semiconductor

cdef class SchottkyDiode(object):
    cdef:
        str __label
        Metal __metal
        Semiconductor __semiconductor
        double __area
        double __length
        double __serial_resistance
