from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.TreeMesh1D cimport TreeMesh1D
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from Schottky.Trap cimport Trap
from Schottky.Potential.TrapPotential cimport TrapPotential


cdef class Dopant(Trap):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label, bint conduction_band_bound,
                 TreeMesh1DUniform concentration, TreeMesh1DUniform f,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0,
                 TrapPotential e_potential=None, TrapPotential h_potential=None):
        '''
        Constructor
        '''
        self.__concentration = concentration
        self.__f = f
        self.__coerce_mesh_tree_occupation(self.__f)
        self.__color = 'k'
        self.__linestyle = '-'
        self.__marker = ''
        super(Dopant, self).__init__(label, conduction_band_bound,
                                     energy_c, energy_v,
                                     e_cs0, h_cs0,
                                     e_cs_activation, h_cs_activation,
                                     e_potential, h_potential)

    @property
    def concentration(self):
        return self.__concentration

    @boundscheck(False)
    @wraparound(False)
    cdef void __coerce_mesh_occupation(self, Mesh1D mesh):
        cdef:
            Py_ssize_t n = mesh.num
            int i
            array[double] new_solution, template = array('d')
        new_solution = clone(template, n, zero=False)
        for i in range(n):
            if mesh.solution[i] > 1.0:
                new_solution[i] = 1.0
            elif mesh.solution[i] < 0.0:
                new_solution[i] = 0.0
            else:
                new_solution[i] = mesh.solution[i]
        mesh.solution = new_solution

    @boundscheck(False)
    @wraparound(False)
    cdef void __coerce_mesh_tree_occupation(self, TreeMesh1D mesh_tree):
        cdef:
            Py_ssize_t n
            int level
        for level in mesh_tree.levels:
            for i in range(len(mesh_tree.tree[level])):
                self.__coerce_mesh_occupation(mesh_tree.tree[level][i])

    @concentration.setter
    def concentration(self, TreeMesh1DUniform concentration):
        self.__concentration = concentration

    @property
    def occupation(self):
        return self.__f

    @occupation.setter
    def occupation(self, TreeMesh1DUniform f):
        self.__f = f
        self.__coerce_mesh_tree_occupation(self.__f)

    cpdef double[:] n_t(self, double[:] z):
        return self.__concentration.interpolate_solution(z)

    cpdef double[:] f(self, double[:] z):
        return self.__f.interpolate_solution(z)

    @property
    def color(self):
        return self.__color

    @color.setter
    def color(self, str color):
        self.__color = color

    @property
    def linestyle(self):
        return self.__linestyle

    @linestyle.setter
    def linestyle(self, str linestyle):
        self.__linestyle = linestyle

    @property
    def marker(self):
        return self.__marker

    @marker.setter
    def marker(self, str marker):
        self.__marker = marker

    def __str__(self):
        s = 'Dopant: %s\nEc-Et: %2.2f eV (%2.2g J)\nEt-Ev: %2.2f eV (%2.2g J)' % (self.label,
                                                                                  self.energy_c_ev,
                                                                                  self.energy_c,
                                                                                  self.energy_v_ev,
                                                                                  self.energy_v)
        return s
