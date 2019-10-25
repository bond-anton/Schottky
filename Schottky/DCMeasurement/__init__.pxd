from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from Schottky.SchottkyDiode cimport SchottkyDiode


cdef class DCMeasurement(object):
    cdef:

        SchottkyDiode __diode

        str __label

        double __initial_step

        TreeMesh1DUniform __ep
        TreeMesh1DUniform __qfe
        TreeMesh1DUniform __qfh
        TreeMesh1DUniform __ne
        TreeMesh1DUniform __nh
        TreeMesh1DUniform __generation
        TreeMesh1DUniform __recombination

    cpdef prepare_psi0(self, double bias, double temperature)
