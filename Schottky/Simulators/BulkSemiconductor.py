from __future__ import division, print_function
# import timeit
import numbers
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

    def _neutrality_equation(self, mu, temperature):
        if isinstance(temperature, (list, tuple, np.ndarray)):
            temperature = np.array(temperature)
        elif isinstance(temperature, numbers.Number):
            temperature = np.array([np.float(temperature)])
        band_gap = self.band_gap(temperature)
        energy_scale = constants['k'] * temperature
        exp_term = np.zeros_like(temperature)
        nonzero_temperature = np.where(temperature > 0)
        exp_term_cb = np.exp(-mu / energy_scale[nonzero_temperature])
        exp_term_vb = np.exp((mu - band_gap) / energy_scale[nonzero_temperature])
        bands_density_of_states = self.bands_density_of_states(temperature)
        cb_charge = -bands_density_of_states['Nc'] * exp_term_cb
        vb_charge = bands_density_of_states['Nv'] * exp_term_vb
        for dopant in self.semiconductor.dopants:
            for part in self.parts.values():
                if isinstance(part, ChargeCarrierTrap):
                    if part.trap.name == dopant.trap.name:
                        dopant_energy = part.energy_level(band_gap)
            print(dopant.trap.charge_state)
            print(dopant_energy)


