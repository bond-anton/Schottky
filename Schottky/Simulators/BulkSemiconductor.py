from __future__ import division, print_function
import numpy as np
from scipy.optimize import bisect

from Schottky.Samples.Semiconductor import Semiconductor
from Schottky.Simulators import Simulator
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap
from Schottky import constants
from ._helpers import storage_manager, prepare_array


class BulkSemiconductor(Simulator):

    def __init__(self, client, semiconductor, description=None):
        assert isinstance(semiconductor, Semiconductor), 'Valid Semiconductor Sample object expected'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {
                'name': 'Simulation',
                'description': 'Simulation software',
                'subcategory': None
            }
        }

        measurement_types = [
            {
                'name': 'Bulk Semiconductor energetics',
                'description': 'Measurement of Semiconductor energetics',
                'children': []
            }
        ]

        self.measurement_details = {
            'band gap': {'name': 'Band Gap temperature dependence',
                         'description': 'Band Gap temperature dependence',
                         'type': 'Bulk Semiconductor energetics'},
            'xi': {'name': 'Electrochemical potential temperature dependence',
                   'description': 'Electrochemical potential temperature dependence',
                   'type': 'Bulk Semiconductor energetics'},
            'dos': {'name': 'Bands effective density of states temperature dependence',
                    'description': 'Semiconductor Bands effective density of states temperature dependence',
                    'type': 'Bulk Semiconductor energetics'},
            'carrier velocity': {'name': 'Charge carrier thermal velocity temperature dependence',
                                 'description': 'Charge carrier thermal velocity temperature dependence',
                                 'type': 'Bulk Semiconductor energetics'},
            'mobility': {'name': 'Charge carrier mobility temperature dependence',
                         'description': 'Charge carrier mobility temperature dependence',
                         'type': 'Bulk Semiconductor energetics'},
            'emission_rate': {'name': 'Emission rate',
                              'description': 'measure charge carriers emission rate',
                              'type': 'Capture and Emission Kinetics'},
        }

        self.measurement_specs = {
            'band gap': {'parameters': None,
                         'input data': None,
                         'variables': [{
                             'name': 'Temperature',
                             'description': 'Sample temperature',
                             'units': 'K'
                         }],
                         'result': [{
                             'name': 'Band Gap',
                             'description': 'Semiconductor Band Gap',
                             'units': 'eV'
                         }]},
            'xi': {'parameters': None,
                   'input data': None,
                   'variables': [{
                       'name': 'Temperature',
                       'description': 'Sample temperature',
                       'units': 'K'
                   }],
                   'result': [{
                       'name': 'Electrochemical potential',
                       'description': 'Semiconductor Electrochemical potential, measured from Ec',
                       'units': 'eV'
                   }]},
            'dos': {'parameters': None,
                    'input data': None,
                    'variables': [{
                        'name': 'Temperature',
                        'description': 'Sample temperature',
                        'units': 'K'
                    }],
                    'result': [
                        {
                            'name': 'DOS C.band',
                            'description': 'Effective Density of States in Conduction band',
                            'units': 'cm^-3'
                        },
                        {
                            'name': 'DOS V.band',
                            'description': 'Effective Density of States in Conduction band',
                            'units': 'cm^-3'
                        }
                    ]},
            'carrier velocity': {'parameters': None,
                                 'input data': None,
                                 'variables': [{
                                     'name': 'Temperature',
                                     'description': 'Sample temperature',
                                     'units': 'K'
                                 }],
                                 'result': [
                                     {
                                         'name': 'electron',
                                         'description': 'electrons average thermal velocity',
                                         'units': 'cm/s'
                                     },
                                     {
                                         'name': 'hole',
                                         'description': 'holes average thermal velocity',
                                         'units': 'cm/s'
                                     }
                                 ]},
            'mobility': {'parameters': [{'name': 'field',
                                         'type': 'numeric',
                                         'default value': 0.0,
                                         'description': 'Electric field',
                                         'units': 'V/cm'}],
                         'input data': None,
                         'variables': [{
                             'name': 'Temperature',
                             'description': 'Sample temperature',
                             'units': 'K'
                         }],
                         'result': [
                             {
                                 'name': 'electron',
                                 'description': 'electrons mobility',
                                 'units': 'cm^2/(V*s)'
                             },
                             {
                                 'name': 'hole',
                                 'description': 'holes mobility',
                                 'units': 'cm^2/(V*s)'
                             }
                         ]},
            'emission_rate': {'name': 'Emission rate',
                              'description': 'measure charge carriers emission rate',
                              'type': 'Capture and Emission Kinetics'},
        }

        self.semiconductor = semiconductor
        samples = [self.semiconductor]
        parts = []
        if semiconductor.dopants:
            parts += [ChargeCarrierTrap(client=client, trap=dopant.trap) for dopant in semiconductor.dopants]
        Simulator.__init__(
            self,
            client=client, name='Bulk Semiconductor Simulator', description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=list(self.measurement_details.values()))

    @storage_manager('band gap', use_storage=True)
    def band_gap(self, temperature=0.0):
        temperature = prepare_array(temperature)
        band_gap = self.semiconductor.band_gap_parameters['Eg_0']
        alpha = self.semiconductor.band_gap_parameters['alpha']
        beta = self.semiconductor.band_gap_parameters['beta']
        shift = alpha * temperature ** 2 / (temperature + beta)
        band_gap = band_gap - shift
        return band_gap

    @storage_manager('dos', use_storage=True)
    def effective_bands_density_of_states(self, temperature=0.0):
        temperature = prepare_array(temperature)
        t_3_2 = temperature ** (3 / 2)
        conduction_band = self.semiconductor.effective_bands_density_of_states['Nc'] * t_3_2
        valence_band = self.semiconductor.effective_bands_density_of_states['Nv'] * t_3_2
        return {'DOS C.band': conduction_band, 'DOS V.band': valence_band}

    @storage_manager('carrier velocity', use_storage=True)
    def carriers_thermal_velocity(self, temperature):
        temperature = prepare_array(temperature)
        m_e = self.semiconductor.effective_mass['electron']
        m_h = self.semiconductor.effective_mass['hole']
        v_e = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_e) * 100
        v_h = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_h) * 100
        return {'electron': v_e, 'hole': v_h}

    @storage_manager('mobility', use_storage=True)
    def mobility(self, temperature=0.0, field=None):
        temperature = prepare_array(temperature)
        result = {'electron': np.zeros_like(temperature), 'hole': np.zeros_like(temperature)}

        zero_temperature = np.where(temperature == 0)
        nonzero_temperature = np.where(temperature > 0)

        if field is None:
            field = 0.0
        np.seterr(divide='warn', invalid='warn')

        # Lattice scattering
        t_300 = temperature / 300.0
        # electrons
        alpha = self.semiconductor.electron_mobility_parameters['alpha']
        mu_l0 = self.semiconductor.electron_mobility_parameters['mu_L0']
        e_lattice = np.ones_like(temperature) * mu_l0
        e_lattice[nonzero_temperature] *= (t_300[nonzero_temperature] ** (-alpha))
        e_lattice[zero_temperature] = 1e200
        # holes
        alpha = self.semiconductor.hole_mobility_parameters['alpha']
        mu_l0 = self.semiconductor.hole_mobility_parameters['mu_L0']
        h_lattice = np.ones_like(temperature) * mu_l0
        h_lattice[nonzero_temperature] *= (t_300[nonzero_temperature] ** (-alpha))
        h_lattice[zero_temperature] = 1e200

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
        e_impurities = np.ones_like(temperature) * numerator
        e_impurities[nonzero_temperature] /= denominator[nonzero_temperature]
        # holes
        a = self.semiconductor.hole_mobility_parameters['A']
        b = self.semiconductor.hole_mobility_parameters['B']
        numerator = a * t_3_2 / dopants_concentration
        denominator = np.log(1 + b * t_2 / dopants_concentration) - b * t_2 / (b * t_2 + dopants_concentration)
        h_impurities = np.ones_like(temperature) * numerator
        h_impurities[nonzero_temperature] /= denominator[nonzero_temperature]

        # Carrier-Carrier scattering
        inv_mobility = 1 / (e_lattice[nonzero_temperature] * 0.88)
        inv_mobility += 1 / (e_impurities[nonzero_temperature] * 0.66)
        result['electron'][nonzero_temperature] = 1 / inv_mobility

        inv_mobility = 1 / (h_lattice[nonzero_temperature] * 0.88)
        inv_mobility += 1 / (h_impurities[nonzero_temperature] * 0.66)
        result['hole'][nonzero_temperature] = 1 / inv_mobility

        # Electric field effect
        v_s = self.semiconductor.electron_mobility_parameters['v_s']
        beta = self.semiconductor.electron_mobility_parameters['beta']
        field_coefficient_e = (1 + (result['electron'] * field / v_s) ** beta) ** (-1 / beta)
        result['electron'] *= field_coefficient_e

        v_s = self.semiconductor.hole_mobility_parameters['v_s']
        beta = self.semiconductor.hole_mobility_parameters['beta']
        field_coefficient_h = (1 + (result['hole'] * field / v_s) ** beta) ** (-1 / beta)
        result['hole'] *= field_coefficient_h
        np.seterr(divide='warn', invalid='warn')
        return result

    @storage_manager('xi', use_storage=True)
    def electrochemical_potential(self, temperature=0.0):
        temperature = prepare_array(temperature)
        energy_scale = constants['k'] * temperature
        band_gap = self.band_gap(temperature=temperature, use_storage=use_storage)
        bands_density_of_states = self.effective_bands_density_of_states(temperature=temperature, use_storage=use_storage)
        mu = np.zeros_like(temperature)
        for i in range(len(temperature)):
            if temperature[i] < 8:
                e1 = self.electrochemical_potential(temperature=8)
                e2 = self.electrochemical_potential(temperature=10)
                a = (e1 - e2) / (8 - 10)
                b = e1 - a * 8
                mu[i] = a * temperature[i] + b
                continue

            def equation(mu):
                np.seterr(divide='warn', invalid='warn', over='raise')
                # bands charge
                if abs(temperature[i]) > 2 * np.finfo(np.float).eps:
                    exp_term_cb = np.exp(-mu / energy_scale[i])
                    exp_term_vb = np.exp((mu - band_gap[i]) / energy_scale[i])
                    cb_charge = -bands_density_of_states['DOS C.band'][i] * exp_term_cb
                    vb_charge = bands_density_of_states['DOS V.band'][i] * exp_term_vb
                    total_charge = cb_charge + vb_charge
                else:
                    total_charge = 0
                # dopants charge
                for dopant in self.semiconductor.dopants:
                    for part in self.parts.values():
                        if isinstance(part, ChargeCarrierTrap):
                            if part.trap.name == dopant.trap.name:

                                delta_energy = part.energy_level(band_gap[i]) - mu
                                if abs(delta_energy) < 2 * np.finfo(np.float).eps:
                                    f = 1 / 2
                                elif energy_scale[i] < 2 * np.finfo(np.float).eps:
                                    f = 1 if delta_energy > 0 else 0
                                else:
                                    try:
                                        f = 1 / (1 + np.exp(-delta_energy / energy_scale[i]))
                                    except FloatingPointError:
                                        f = 1 if delta_energy > 0 else 0
                                q1 = dopant.trap.charge_state['full']
                                q0 = dopant.trap.charge_state['empty']
                                n = dopant.concentration
                                total_charge += n * (q1 - q0) * f + n * q0
                np.seterr(divide='warn', invalid='warn', over='warn')
                return total_charge

            mu[i] = bisect(equation, a=0, b=band_gap[i])

        return mu

    def get_type(self, temperature=0.0):
        temperature = prepare_array(temperature)
        xi = self.electrochemical_potential(temperature=temperature)
        bg2 = self.band_gap(temperature=temperature) / 2
        p = np.where(xi - bg2 > 0)
        n = np.where(xi - bg2 < 0)
        return p, n
