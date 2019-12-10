from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

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

    cpdef double phi_b_n_t(self, double temperature)
    cpdef double[:] phi_b_n(self, double[:] temperature)
    cpdef double phi_b_n_ev_t(self, double temperature)
    cpdef double[:] phi_b_n_ev(self, double[:] temperature)
    cpdef double phi_b_n_boltzmann_t(self, double temperature)
    cpdef double[:] phi_b_n_boltzmann(self, double[:] temperature)
    cpdef double phi_b_p_t(self, double temperature)
    cpdef double[:] phi_b_p(self, double[:] temperature)
    cpdef double phi_b_p_ev_t(self, double temperature)
    cpdef double[:] phi_b_p_ev(self, double[:] temperature)
    cpdef double phi_b_p_boltzmann_t(self, double temperature)
    cpdef double[:] phi_b_p_boltzmann(self, double[:] temperature)
    cpdef double built_in_voltage_t(self, double temperature)
    cpdef double[:] built_in_voltage(self, double[:] temperature)
    cpdef double built_in_voltage_ev_t(self, double temperature)
    cpdef double[:] built_in_voltage_ev(self, double[:] temperature)
    cpdef double built_in_voltage_boltzmann_t(self, double temperature)
    cpdef double[:] built_in_voltage_boltzmann(self, double[:] temperature)
