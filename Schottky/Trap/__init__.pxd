cdef class Trap(object):
    cdef:
        str __label
        double __energy_c
        double __energy_v
        double __e_cs
        double __h_cs
        double __e_cs_activation
        double __h_cs_activation
        double __f

    cpdef double e_cs(self, double temperature)
    cpdef double h_cs(self, double temperature)
    cpdef double e_c(self, double temperature, double v_e)
    cpdef double h_c(self, double temperature, double v_h)
    cpdef double e_cr(self, double temperature, double v_e, double n_e)
    cpdef double h_cr(self, double temperature, double v_h, double n_h)
