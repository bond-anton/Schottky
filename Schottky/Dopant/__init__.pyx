from __future__ import division, print_function

import numpy as np

from cython cimport boundscheck, wraparound
from libc.math cimport exp

from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.TreeMesh1D cimport TreeMesh1D
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform
from Schottky.Trap cimport Trap


cdef class Dopant(Trap):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label, TreeMesh1DUniform concentration, TreeMesh1DUniform f,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0):
        '''
        Constructor
        '''
        self.__concentration = concentration
        self.__f = f
        super(Dopant, self).__init__(label,
                                     energy_c, energy_v,
                                     e_cs0, h_cs0,
                                     e_cs_activation, h_cs_activation)

    @property
    def concentration(self):
        return self.__concentration

    @boundscheck(False)
    @wraparound(False)
    cpdef coerce_mesh_occupation(self, Mesh1D mesh):
        cdef:
            Py_ssize_t n = mesh.num
            int i
            double[:] new_solution = np.zeros(n, dtype=np.double)
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
    cpdef __coerce_mesh_tree_occupation(self, TreeMesh1D mesh_tree):
        cdef:
            Py_ssize_t n
            int level
        for level in mesh_tree.levels:
            print(level)

    @concentration.setter
    def concentration(self, TreeMesh1DUniform concentration):
        self.__concentration = concentration

    @property
    def occupation(self):
        return self.__f

    @occupation.setter
    def occupation(self, TreeMesh1DUniform f):
        self.__f = f

    cpdef double[:] n_t(self, double[:] z):
        return self.__concentration.interpolate_solution(z)

    cpdef double[:] f(self, double[:] z):
        return self.__f.interpolate_solution(z)

    def __str__(self):
        s = 'Dopant: %s\nEc-Et: %2.2f eV (%2.2g J)\nEt-Ev: %2.2f eV (%2.2g J)' % (self.label,
                                                                                  self.energy_c_ev,
                                                                                  self.energy_c,
                                                                                  self.energy_v_ev,
                                                                                  self.energy_v)
        return s
