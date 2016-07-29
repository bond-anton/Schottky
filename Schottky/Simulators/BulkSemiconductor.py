from __future__ import division, print_function
# import timeit
import numbers
import numpy as np
from scipy.optimize import fsolve, root, bisect

from Schottky.Samples.Semiconductor import Semiconductor
from Schottky.Simulators import Simulator
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap
from Schottky import constants


class BulkSemiconductor(Simulator):

    def __init__(self, client, semiconductor, description=None):
        assert isinstance(semiconductor, Semiconductor), 'Valid Semiconductor Sample object expected'
        self.semiconductor = semiconductor
        samples = [self.semiconductor]
        name = 'Bulk Semiconductor Simulator'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}
        measurement_types = [{'name': 'Bulk Semiconductor energetics',
                              'description': 'Measurement of Semiconductor energetics',
                              'children': []}]
        measurements = [
            {'name': 'Band Gap versus temperature',
             'description': 'Measurement of Semiconductor band gap versus temperature',
             'type': 'Bulk Semiconductor energetics'},
            {'name': 'Emission rate',
             'description': 'measure charge carriers emission rate',
             'type': 'Capture and Emission Kinetics'},
        ]
        parts = []
        if semiconductor.dopants:
            parts += [ChargeCarrierTrap(client=client, trap=dopant.trap) for dopant in semiconductor.dopants]
        if not parts:
            parts = None
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=measurements)

    def band_gap(self, temperature):
        band_gap = self.semiconductor.band_gap_parameters['Eg_0']
        alpha = self.semiconductor.band_gap_parameters['alpha']
        beta = self.semiconductor.band_gap_parameters['beta']
        shift = alpha * temperature ** 2 / (temperature + beta)
        return band_gap - shift

    def bands_density_of_states(self, temperature):
        t_3_2 = temperature ** (3 / 2)
        conduction_band = self.semiconductor.bands_density_of_states['Nc'] * t_3_2
        valence_band = self.semiconductor.bands_density_of_states['Nv'] * t_3_2
        return {'Nc': conduction_band, 'Nv': valence_band}

    def carriers_thermal_velocity(self, temperature):
        m_e = self.semiconductor.effective_mass['electron']
        m_h = self.semiconductor.effective_mass['hole']
        v_e = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_e) * 100
        v_h = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_h) * 100
        return {'electron': v_e, 'hole': v_h}

    def mobility(self, temperature, field=None):
        if isinstance(temperature, (list, tuple, np.ndarray)):
            temperature = np.array(temperature)
        elif isinstance(temperature, numbers.Number):
            temperature = np.array([np.float(temperature)])
        zero_temperature = np.where(temperature == 0)
        nonzero_temperature = np.where(temperature > 0)

        if field is None:
            field = 0.0
        result = {
            'electron':
                {'lattice': 0,
                 'impurities': 0,
                 'carrier-carrier': 0,
                 'total': 0},
            'hole':
                {'lattice': 0,
                 'impurities': 0,
                 'carrier-carrier': 0,
                 'total': 0}
        }

        np.seterr(divide='warn', invalid='warn')

        # Lattice scattering
        t_300 = temperature / 300.0
        # electrons
        alpha = self.semiconductor.electron_mobility_parameters['alpha']
        mu_l0 = self.semiconductor.electron_mobility_parameters['mu_L0']
        result['electron']['lattice'] = np.ones_like(temperature) * mu_l0
        result['electron']['lattice'][nonzero_temperature] *= (t_300[nonzero_temperature] ** (-alpha))
        result['electron']['lattice'][zero_temperature] = 1e200
        # holes
        alpha = self.semiconductor.hole_mobility_parameters['alpha']
        mu_l0 = self.semiconductor.hole_mobility_parameters['mu_L0']
        result['hole']['lattice'] = np.ones_like(temperature) * mu_l0
        result['hole']['lattice'][nonzero_temperature] *= (t_300[nonzero_temperature] ** (-alpha))
        result['hole']['lattice'][zero_temperature] = 1e200

        # Impurities scattering
        t_3_2 = temperature ** (3 / 2)
        t_2 = temperature ** 2
        dopants_concentration = 0
        for dopant in self.semiconductor.dopants:
            dopants_concentration += dopant.concentration
        # electrons
        a = self.semiconductor.electron_mobility_parameters['A']
        b = self.semiconductor.electron_mobility_parameters['B']
        numerator = a * t_3_2 / dopants_concentration
        denominator = np.log(1 + b * t_2 / dopants_concentration) - b * t_2 / (b * t_2 + dopants_concentration)
        result['electron']['impurities'] = np.ones_like(temperature) * numerator
        result['electron']['impurities'][nonzero_temperature] /= denominator[nonzero_temperature]
        # holes
        a = self.semiconductor.hole_mobility_parameters['A']
        b = self.semiconductor.hole_mobility_parameters['B']
        numerator = a * t_3_2 / dopants_concentration
        denominator = np.log(1 + b * t_2 / dopants_concentration) - b * t_2 / (b * t_2 + dopants_concentration)
        result['hole']['impurities'] = np.ones_like(temperature) * numerator
        result['hole']['impurities'][nonzero_temperature] /= denominator[nonzero_temperature]

        # Carrier-Carrier scattering
        result['electron']['total'] = np.zeros_like(temperature)
        inv_mobility = 1 / (result['electron']['lattice'][nonzero_temperature] * 0.88)
        inv_mobility += 1 / (result['electron']['impurities'][nonzero_temperature] * 0.66)
        result['electron']['total'][nonzero_temperature] = 1 / inv_mobility

        result['hole']['total'] = np.zeros_like(temperature)
        inv_mobility = 1 / (result['hole']['lattice'][nonzero_temperature] * 0.88)
        inv_mobility += 1 / (result['hole']['impurities'][nonzero_temperature] * 0.66)
        result['hole']['total'][nonzero_temperature] = 1 / inv_mobility

        # Electric field effect
        v_s = self.semiconductor.electron_mobility_parameters['v_s']
        beta = self.semiconductor.electron_mobility_parameters['beta']
        field_coefficient_e = (1 + (result['electron']['total'] * field / v_s) ** beta) ** (-1 / beta)
        result['electron']['total'] *= field_coefficient_e

        v_s = self.semiconductor.hole_mobility_parameters['v_s']
        beta = self.semiconductor.hole_mobility_parameters['beta']
        field_coefficient_h = (1 + (result['hole']['total'] * field / v_s) ** beta) ** (-1 / beta)
        result['hole']['total'] *= field_coefficient_h

        np.seterr(divide='warn', invalid='warn')
        return result

    def electrochemical_potential(self, temperature):
        assert isinstance(temperature, numbers.Number), 'Temperature must be a number'
        if temperature < 8:
            e1 = self.electrochemical_potential(temperature=8)
            e2 = self.electrochemical_potential(temperature=10)
            a = (e1 - e2) / (8 - 10)
            b = e1 - a*8
            return a * temperature + b
        energy_scale = constants['k'] * np.float(temperature)
        band_gap = self.band_gap(np.float(temperature))
        bands_density_of_states = self.bands_density_of_states(np.float(temperature))

        def equation(mu):
            np.seterr(divide='warn', invalid='warn', over='raise')
            # bands charge
            if abs(temperature) > 2 * np.finfo(np.float).eps:
                exp_term_cb = np.exp(-mu / energy_scale)
                exp_term_vb = np.exp((mu - band_gap) / energy_scale)
                cb_charge = -bands_density_of_states['Nc'] * exp_term_cb
                vb_charge = bands_density_of_states['Nv'] * exp_term_vb
                total_charge = cb_charge + vb_charge
            else:
                total_charge = 0
            # dopants charge
            for dopant in self.semiconductor.dopants:
                for part in self.parts.values():
                    if isinstance(part, ChargeCarrierTrap):
                        if part.trap.name == dopant.trap.name:

                            delta_energy = part.energy_level(band_gap) - mu
                            if abs(delta_energy) < 2 * np.finfo(np.float).eps:
                                f = 1 / 2
                            elif energy_scale < 2 * np.finfo(np.float).eps:
                                f = 1 if delta_energy > 0 else 0
                            else:
                                try:
                                    f = 1 / (1 + np.exp(-delta_energy / energy_scale))
                                except FloatingPointError:
                                    f = 1 if delta_energy > 0 else 0
                            q1 = dopant.trap.charge_state['full']
                            q0 = dopant.trap.charge_state['empty']
                            n = dopant.concentration
                            total_charge += n * (q1 - q0) * f + n * q0
            # print(mu, part.energy_level(band_gap), energy_scale, n, f, total_charge)
            np.seterr(divide='warn', invalid='warn', over='warn')
            return total_charge


        return bisect(equation, a=0, b=band_gap)
