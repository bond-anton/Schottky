from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from Schottky.Trap cimport Trap

cdef class Dopant(Trap):
    cdef:
        TreeMesh1DUniform __concentration
        TreeMesh1DUniform __f

    cpdef double[:] n_t(self, double[:] z)
    cpdef double[:] f(self, double[:] z)
