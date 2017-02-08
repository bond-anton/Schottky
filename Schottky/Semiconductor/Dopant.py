from __future__ import division

import sympy as sym

from Schottky.Notation import q
from Schottky.Semiconductor.Trap import Trap
from Schottky.Helpers import interp_Fn
from Schottky.Notation import q, k
from Schottky.Helpers import to_numeric
import numpy as np
import os.path



class Dopant(Trap):
    """
    Describes dopant in semiconductor
    """

    def __init__(self, name, concentration, charge_states=None,
                 energy_distribution_function='Single Level', energy_spread=0.3 * q,
                 electron_capture_cross_section=1e-17, electron_capture_cross_section_activation_energy=0,
                 hole_capture_cross_section=1e-17, hole_capture_cross_section_activation_energy=0,
                 trap_potential=None, poole_frenkel_lookup=None):
        super(Dopant, self).__init__(name, charge_states, energy_distribution_function, energy_spread,
                                     electron_capture_cross_section, electron_capture_cross_section_activation_energy,
                                     hole_capture_cross_section, hole_capture_cross_section_activation_energy,
                                     trap_potential)
        self.concentration = self.__prepare_1dfunc(concentration)
        self.F = self.__prepare_1dfunc(0.0)
        self.dF = self.__prepare_1dfunc(0.0)

        self.equilibrium_f_memo = {}
        self.d_equilibrium_f_d_fermi_energy_memo = {}
        self.poole_frenkel_lookup = poole_frenkel_lookup

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

    '''
    def generate_pf_lookup(self, Nl_range, F_range, T):
        if self.poole_frenkel_lookup is not None:
            if os.path.isfile(fname):
                return False
            kT = to_numeric(k * schottky_diode.T / q)
            poole_frenkel = np.ones((len(Nl_range), len(F_range)), dtype=np.float)
            theta_points = 100
            theta = np.linspace(0, np.pi, num=theta_points, endpoint=True)
            for i, Nl in enumerate(Nl_range):
                barrier_lowering = np.zeros((theta_points, len(F_range)), dtype=np.float)
                for j, F in enumerate(F_range):
                    local_electric_field_r = abs(F)
                    local_electric_field_theta = 0 if F >= 0 else np.pi
                    local_electric_field_3d = (local_electric_field_r, local_electric_field_theta, 0.0)
                    self.trap_potential.get_potential_by_name('Charged Dislocation') \
                        .set_linear_charge_density(Nl)
                    self.trap_potential.get_potential_by_name('External Field') \
                        .external_field = local_electric_field_3d
                    barrier_lowering[:, j] = np.array(
                        [self.trap_potential.barrier_lowering(theta_i)[0] for theta_i in theta])
                poole_frenkel[i, :] = 0.5 * np.trapz(np.exp(abs(barrier_lowering) / kT), theta, axis=0)

    '''

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
