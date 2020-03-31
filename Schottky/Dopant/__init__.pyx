from BDFunction1D cimport Function

from Schottky.Trap cimport Trap
from Schottky.Potential.TrapPotential cimport TrapPotential


cdef class Dopant(Trap):
    '''
    Charge carrier trap class
    '''

    def __init__(self, str label, bint conduction_band_bound,
                 Function concentration,
                 double energy_c, double energy_v,
                 double e_cs0, double h_cs0,
                 double e_cs_activation=0.0, double h_cs_activation=0.0,
                 TrapPotential e_potential=None, TrapPotential h_potential=None):
        '''
        Constructor
        '''
        self.__concentration = concentration
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

    @concentration.setter
    def concentration(self, Function concentration):
        self.__concentration = concentration

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
