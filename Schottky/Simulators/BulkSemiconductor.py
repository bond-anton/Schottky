from __future__ import division, print_function
# import timeit
import numpy as np

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
        np.seterr(divide='warn', invalid='warn')
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
        t_300 = temperature / 300.0

        result['electron']['lattice'] = self.semiconductor.electron_mobility_parameters['mu_L0']
        result['electron']['lattice'] *= (t_300 ** (-self.semiconductor.electron_mobility_parameters['alpha']))
        result['hole']['lattice'] = self.semiconductor.hole_mobility_parameters['mu_L0']
        result['hole']['lattice'] *= (t_300 ** (-self.semiconductor.hole_mobility_parameters['alpha']))

        t_3_2 = temperature ** (3 / 2)
        t_2 = temperature ** 2
        dopants_concentration = 0
        for dopant in self.semiconductor.dopants:
            dopants_concentration += dopant.concentration
        a = self.semiconductor.electron_mobility_parameters['A']
        b = self.semiconductor.electron_mobility_parameters['B']
        numerator = a * t_3_2 / dopants_concentration
        denominator = np.log(1 + b * t_2 / dopants_concentration) - b * t_2 / (b * t_2 + dopants_concentration)
        result['electron']['impurities'] = numerator / denominator

        a = self.semiconductor.hole_mobility_parameters['A']
        b = self.semiconductor.hole_mobility_parameters['B']
        numerator = a * t_3_2 / dopants_concentration
        denominator = np.log(1 + b * t_2 / dopants_concentration) - b * t_2 / (b * t_2 + dopants_concentration)
        result['hole']['impurities'] = numerator / denominator

        band_gap = self.band_gap(temperature)
        energy_scale = constants['k'] * temperature
        bands_density_of_states = self.bands_density_of_states(temperature)
        exp_factor = np.exp(-band_gap / energy_scale)

        pn = bands_density_of_states['Nc'] * bands_density_of_states['Nv'] * exp_factor
        mu_ccs = (2e17 * t_3_2 / np.sqrt(pn)) / (np.log(1 + 8.28e8 * t_2 * (pn ** (-1 / 3))))
        cross_factor_e = 6 * result['electron']['lattice'] * (result['electron']['impurities'] + mu_ccs)
        cross_factor_e /= result['electron']['impurities'] * mu_ccs
        cross_factor_e = np.sqrt(cross_factor_e)
        cross_factor_h = 6 * result['hole']['lattice'] * (result['hole']['impurities'] + mu_ccs)
        cross_factor_h /= result['hole']['impurities'] * mu_ccs
        cross_factor_h = np.sqrt(cross_factor_h)
        result['electron']['carrier-carrier'] = mu_ccs
        result['hole']['carrier-carrier'] = mu_ccs

        result['electron']['total'] = result['electron']['lattice'] * 1
        result['electron']['total'] *= (1.025 / (1 + ((cross_factor_e / 1.68) ** 1.43)) - 0.025)

        result['hole']['total'] = result['hole']['lattice'] * 1
        result['hole']['total'] *= (1.025 / (1 + ((cross_factor_h / 1.68) ** 1.43)) - 0.025)

        v_s = self.semiconductor.electron_mobility_parameters['v_s']
        beta = self.semiconductor.electron_mobility_parameters['beta']
        field_coefficient_e = (1 + (result['electron']['total'] * field * 1e-2 / v_s) ** beta) ** (-1 / beta)
        result['electron']['total'] *= field_coefficient_e

        v_s = self.semiconductor.hole_mobility_parameters['v_s']
        beta = self.semiconductor.hole_mobility_parameters['beta']
        field_coefficient_h = (1 + (result['hole']['total'] * field * 1e-2 / v_s) ** beta) ** (-1 / beta)
        result['hole']['total'] *= field_coefficient_h

        np.seterr(divide='warn', invalid='warn')
        return result
