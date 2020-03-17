from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from Schottky.SchottkyDiode cimport SchottkyDiode


cdef class DCMeasurement(object):
    cdef:

        SchottkyDiode __diode

        str __label

        double __temperature
        double __bias

        double __initial_step

        TreeMesh1DUniform __ep
        TreeMesh1DUniform __qf_e
        TreeMesh1DUniform __qf_h
        TreeMesh1DUniform __n_e
        TreeMesh1DUniform __n_h
        TreeMesh1DUniform __generation
        TreeMesh1DUniform __recombination

    cpdef prepare_psi0(self)
    cpdef double[:] qf_e_to_n_e(self, double[:] qf_e)
    cpdef double[:] qf_h_to_n_h(self, double[:] qf_h)
    cpdef double[:] n_e_to_qf_e(self, double[:] n_e)
    cpdef double[:] n_h_to_qf_h(self, double[:] n_h)
    cpdef double thermionic_emission_current_e(self)
    cpdef double thermionic_emission_current_h(self)
