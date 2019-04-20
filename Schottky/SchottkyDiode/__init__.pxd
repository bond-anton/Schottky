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

    cpdef double built_in_voltage_t(self, double temperature, bint electron_volts=*)
    cpdef double[:] built_in_voltage(self, double[:] temperature, bint electron_volts=*)
