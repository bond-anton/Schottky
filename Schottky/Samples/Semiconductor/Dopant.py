from __future__ import division

import numpy as np
import sympy as sym

from Schottky.Helpers import interp_Fn
from Schottky.Notation import q
from Schottky.Samples.Semiconductor.Trap_old import Trap


class Dopant(Trap):
    """
    Describes dopant in semiconductor
    """

    equipment = 'Semiconductor Dopant Simulator'
    description = 'Simulates properties of a dopant in bulk semiconductor'

    def __init__(self, name, concentration, charge_states=None,
                 energy_distribution_function='Single Level', energy_spread=0.3 * q,
                 electron_capture_cross_section=1e-17, electron_capture_cross_section_activation_energy=0,
                 hole_capture_cross_section=1e-17, hole_capture_cross_section_activation_energy=0):
        super(Dopant, self).__init__(name, charge_states, energy_distribution_function, energy_spread,
                                     electron_capture_cross_section, electron_capture_cross_section_activation_energy,
                                     hole_capture_cross_section, hole_capture_cross_section_activation_energy)
        self.concentration = self.__prepare_1dfunc(concentration)
        self.F = self.__prepare_1dfunc(0.0)
        self.dF = self.__prepare_1dfunc(0.0)

        self.equilibrium_f_memo = {}
        self.d_equilibrium_f_d_fermi_energy_memo = {}

    def set_F_func(self, F):
        self.F = self.__prepare_1dfunc(F)

    def set_dF_func(self, dF):
        self.dF = self.__prepare_1dfunc(dF)

    def set_F_interp(self, Z, F):
        '''
        z and F must be 1D arrays of equal size
        '''
        self.F = interp_Fn(Z, F, interp_type='last')

    def set_dF_interp(self, Z, dF):
        '''
        z and dF must be 1D arrays of equal size
        '''
        self.dF = interp_Fn(Z, dF, interp_type='last')

    def __str__(self):
        z = sym.symbols('z')
        descr = 'Dopant: ' + self.name + ', z -> ' + str(self.concentration(z)) + '\n'
        for state in self.charge_states:
            try:
                descr += ' Charge state: ' + str(state[0]) + '@%2.3f' % (state[1] / q) + 'eV'
            except:
                descr += ' Charge state: ' + str(state[0]) + '@' + str(state[1])
                pass
        return descr

    def __prepare_1dfunc(self, F):
        if hasattr(F, '__call__'):
            pre_f = F
        else:
            pre_f = lambda z: np.float(F)

        def final_f(z):
            if isinstance(z, np.ndarray):
                return np.array(map(pre_f, z))
            else:
                return pre_f(z)

        return final_f
