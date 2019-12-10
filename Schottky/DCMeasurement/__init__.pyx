import numpy as np
from BDMesh.Mesh1DUniform cimport Mesh1DUniform
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from Schottky.SchottkyDiode cimport SchottkyDiode


cdef class DCMeasurement(object):

    def __init__(self, str label, SchottkyDiode diode, double temperature=0, double initial_step=5e-7):
        self.__label = label
        self.__diode = diode
        self.__temperature = temperature
        self.__initial_step = initial_step
        self.__ep = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                      aligned=True)
        self.__qfe = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                       aligned=True)
        self.__qfh = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                       aligned=True)
        self.__ne = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                      aligned=True)
        self.__nh = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                      aligned=True)
        self.__generation = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                              aligned=True)
        self.__recombination = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                                 aligned=True)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, double temperature):
        self.__temperature = temperature

    @property
    def diode(self):
        return self.__diode

    @property
    def electric_potential(self):
        return self.__ep

    @property
    def quasi_fermi_e(self):
        return self.__qfe

    @property
    def quasi_fermi_h(self):
        return self.__qfh

    @property
    def n_e(self):
        return self.__ne

    @property
    def n_h(self):
        return self.__nh

    @property
    def generation(self):
        return self.__generation

    @property
    def recombination(self):
        return self.__recombination

    cpdef prepare_psi0(self, double bias):
        v_bi = self.__diode.built_in_voltage_ev_t(self.__temperature)
        v_bc = v_bi - bias
        psi_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        psi_0.solution = v_bc - np.asarray(psi_0.physical_nodes) * v_bc / psi_0.physical_nodes[-1]
        print(np.asarray(psi_0.physical_nodes))
        print(np.asarray(psi_0.solution))
        self.__ep = TreeMesh1DUniform(psi_0, aligned=True)
        qfe_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        qfe_0.solution = np.asarray(qfe_0.solution) + self.__diode.semiconductor.el_chem_pot_ev_t(self.__temperature)
        self.__qfe = TreeMesh1DUniform(qfe_0, aligned=True)
        qfh_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        qfh_0.solution = np.asarray(qfh_0.solution) + self.__diode.semiconductor.el_chem_pot_ev_t(self.__temperature)
        self.__qfh = TreeMesh1DUniform(qfh_0, aligned=True)
