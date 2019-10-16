from BDMesh.Mesh1DUniform import Mesh1DUniform
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from Schottky.SchottkyDiode cimport SchottkyDiode


cdef class DCMeasurement(object):

    def __init__(self, str label, SchottkyDiode diode, double initial_step=5e-7 ):
        self.__label = label
        self.__diode = diode
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
