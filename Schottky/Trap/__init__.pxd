from Schottky.Potential.TrapPotential cimport TrapPotential


cdef class Trap(object):
    cdef:
        str __label
        bint __cb_bound
        double __energy_c
        double __energy_v
        double __e_cs0
        double __h_cs0
        double __e_cs_activation
        double __h_cs_activation
        dict __charge_state
        dict __g
        dict __capture_barrier
        TrapPotential __e_potential
        TrapPotential __h_potential

    cpdef double e_cs(self, double temperature)
    cpdef double h_cs(self, double temperature)
    cpdef double e_c(self, double temperature, double v_e)
    cpdef double h_c(self, double temperature, double v_h)
    cpdef double e_cr(self, double temperature, double v_e, double n_e, double f)
    cpdef double h_cr(self, double temperature, double v_h, double n_h, double f)

    cpdef double e_er(self, double temperature, double v_e, double n_c, double f)
    cpdef double h_er(self, double temperature, double v_h, double n_v, double f)

    cpdef double f_eq(self, double temperature,
                      double v_e, double n_e, double n_c,
                      double v_h, double n_h, double n_v,
                      double f)

    cpdef double df_dt(self, double temperature,
                       double v_e, double n_e, double n_c,
                       double v_h, double n_h, double n_v,
                       double f)
