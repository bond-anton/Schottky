from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.TreeMesh1D cimport TreeMesh1D
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from Schottky.Trap cimport Trap

cdef class Dopant(Trap):
    cdef:
        TreeMesh1DUniform __concentration
        TreeMesh1DUniform __f

    cpdef double[:] n_t(self, double[:] z)
    cpdef double[:] f(self, double[:] z)
    cpdef coerce_mesh_occupation(self, Mesh1D mesh)
    cpdef __coerce_mesh_tree_occupation(self, TreeMesh1D mesh)
