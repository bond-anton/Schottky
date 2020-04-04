from BDFunction1D cimport Function, Functional

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

        double __temperature
        double __bias

        Function __ep
        Function __ef
        Function __qf_e
        Function __qf_h
        Function __grad_qf_e
        Function __grad_qf_h
        Function __generation
        Function __recombination
        Function __ne
        Function __nh
        Function __pn
        Function __mu_e
        Function __mu_h
        Function __dopants_charge

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
    cpdef double n0_t(self, double temperature)
    cpdef double[:] n0(self, double[:] temperature)
    cpdef double p0_t(self, double temperature)
    cpdef double[:] p0(self, double[:] temperature)
    cpdef double thermionic_emission_current_e(self)
    cpdef double thermionic_emission_current_h(self)
    cpdef stationary_grad_qf_e_solver(self)
    cpdef stationary_grad_qf_h_solver(self)
    cpdef Function poisson_eq_solver(self)


cdef class TestSCR(Function):
    pass


cdef class QFeNe(Functional):
    cdef:
        SchottkyDiode __diode


cdef class QFhNh(Functional):
    cdef:
        SchottkyDiode __diode


cdef class NeQFe(Functional):
    cdef:
        SchottkyDiode __diode


cdef class NhQFh(Functional):
    cdef:
        SchottkyDiode __diode


cdef class MobilityE(Function):
    cdef:
        SchottkyDiode __diode


cdef class MobilityH(Function):
    cdef:
        SchottkyDiode __diode


cdef class DopantsEqCharge(Function):
    cdef:
        SchottkyDiode __diode


cdef class GradQFeF(Functional):
    cdef:
        Function __g
        Function __r
        Function __mu
        Function __n
        double __temperature
        double __qkt


cdef class GradQFhF(Functional):
    cdef:
        Function __g
        Function __r
        Function __mu
        Function __n
        double __temperature
        double __qkt
