from __future__ import division, print_function
import timeit
import numbers
import numpy as np
from scipy.optimize import fsolve, root, bisect

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
            'xi': {'name': 'Electrochemical potential temperature dependence',
                   'description': 'Electrochemical potential temperature dependence',
                   'type': 'Bulk Semiconductor energetics'},
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
                            'units': 'm^-3'
                        },
                        {
                            'name': 'DOS V.band',
                            'description': 'Effective Density of States in Conduction band',
                            'units': 'm^-3'
                        }
                    ]},
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

    def carriers_thermal_velocity(self, temperature, use_storage=False):
        start_time = timeit.default_timer()
        if isinstance(temperature, (list, tuple, np.ndarray)):
            temperature = np.array(temperature)
        elif isinstance(temperature, numbers.Number):
            temperature = np.array([np.float(temperature)])
        else:
            raise TypeError('Unsupported temperature type')
        measurement_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        if use_storage:
            measurement_details = {
                'name': 'Charge carrier thermal velocity temperature dependence',
                'description': 'Measurement of Semiconductor Charge carrier thermal velocity temperature dependence',
                'type': 'Bulk Semiconductor energetics'}
            record = 'Starting Measurement "%s"' % (measurement_details['name'])
            self.client.log_manager.log_record(record=record, category='Information')
            parameters = None
            measurement = self.register_measurement(measurement_details=measurement_details, parameters=parameters,
                                                    input_data=None, force_new=False)
            if measurement.progress == 100:
                temperature_channels = self.client.measurement_manager.get_data_channels(measurement=measurement,
                                                                                         name='Temperature')
                for temperature_channel in temperature_channels:
                    stored = self.client.measurement_manager.get_data_points_array(temperature_channel)[:, 0]
                    if stored.size == temperature.size and (stored == temperature).all():
                        e_channel = self.load_create_data_channel(
                            channel_name='electron',
                            measurement=measurement,
                            description='Thermal velocity of electrons', unit_name='m/s')
                        v_e = self.client.measurement_manager.get_data_points_array(e_channel)[:, 0]
                        h_channel = self.load_create_data_channel(
                            channel_name='hole',
                            measurement=measurement,
                            description='Thermal velocity of holes', unit_name='m/s')
                        v_h = self.client.measurement_manager.get_data_points_array(h_channel)[:, 0]
                        db_time = timeit.default_timer() - start_time
                        record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                                 (measurement_details['name'], db_time, 0.0, db_time)
                        self.client.log_manager.log_record(record=record, category='Information')
                        return {'electron': v_e, 'hole': v_h}
            temperature_channel = self.load_create_data_channel(channel_name='Temperature', measurement=measurement,
                                                                description='Temperature', unit_name='K')
            e_channel = self.load_create_data_channel(
                channel_name='electron',
                measurement=measurement,
                description='Thermal velocity of electrons', unit_name='m/s')
            h_channel = self.load_create_data_channel(
                channel_name='hole',
                measurement=measurement,
                description='Thermal velocity of holes', unit_name='m/s')
        db_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        m_e = self.semiconductor.effective_mass['electron']
        m_h = self.semiconductor.effective_mass['hole']
        v_e = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_e) * 100
        v_h = np.sqrt(3 * constants['k'] * temperature * constants['q'] / m_h) * 100
        measurement_time += timeit.default_timer() - start_time
        if use_storage:
            start_time = timeit.default_timer()
            self.client.measurement_manager.create_data_points(channel=temperature_channel,
                                                               float_value=temperature)
            self.client.measurement_manager.create_data_points(channel=e_channel,
                                                               float_value=v_e)
            self.client.measurement_manager.create_data_points(channel=h_channel,
                                                               float_value=v_h)
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            db_time += timeit.default_timer() - start_time
            record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                     (measurement_details['name'], db_time + measurement_time,
                      measurement_time, db_time)
            self.client.log_manager.log_record(record=record, category='Information')
        return {'electron': v_e, 'hole': v_h}

    def mobility(self, temperature, field=None, use_storage=False):
        start_time = timeit.default_timer()
        if isinstance(temperature, (list, tuple, np.ndarray)):
            temperature = np.array(temperature)
        elif isinstance(temperature, numbers.Number):
            temperature = np.array([np.float(temperature)])
        else:
            raise TypeError('Unsupported temperature type')
        result = {
            'electron':
                {'lattice': 0,
                 'impurities': 0,
                 'total': 0},
            'hole':
                {'lattice': 0,
                 'impurities': 0,
                 'total': 0}
        }
        measurement_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        if use_storage:
            measurement_details = {
                'name': 'Charge carrier mobility temperature dependence',
                'description': 'Measurement of Semiconductor Charge carrier mobility temperature dependence',
                'type': 'Bulk Semiconductor energetics'}
            record = 'Starting Measurement "%s"' % (measurement_details['name'])
            self.client.log_manager.log_record(record=record, category='Information')
            parameters = None
            measurement = self.register_measurement(measurement_details=measurement_details, parameters=parameters,
                                                    input_data=None, force_new=False)
            if measurement.progress == 100:
                temperature_channels = self.client.measurement_manager.get_data_channels(measurement=measurement,
                                                                                         name='Temperature')
                for temperature_channel in temperature_channels:
                    stored = self.client.measurement_manager.get_data_points_array(temperature_channel)[:, 0]
                    if stored.size == temperature.size and (stored == temperature).all():
                        e_lattice = self.load_create_data_channel(
                            channel_name='electron lattice',
                            measurement=measurement,
                            description='Lattice component of electrons mobility', unit_name='m^2/(V*s)')
                        e_impurities = self.load_create_data_channel(
                            channel_name='electron impurities',
                            measurement=measurement,
                            description='Impurities component of electrons mobility', unit_name='m^2/(V*s)')
                        e_total = self.load_create_data_channel(
                            channel_name='electron total',
                            measurement=measurement,
                            description='Total electrons mobility', unit_name='m^2/(V*s)')
                        result['electron']['lattice'] = self.client.measurement_manager.get_data_points_array(
                            e_lattice)[:, 0]
                        result['electron']['impurities'] = self.client.measurement_manager.get_data_points_array(
                            e_impurities)[:, 0]
                        result['electron']['total'] = self.client.measurement_manager.get_data_points_array(
                            e_total)[:, 0]
                        h_lattice = self.load_create_data_channel(
                            channel_name='hole lattice',
                            measurement=measurement,
                            description='Lattice component of holes mobility', unit_name='m^2/(V*s)')
                        h_impurities = self.load_create_data_channel(
                            channel_name='hole impurities',
                            measurement=measurement,
                            description='Impurities component of holes mobility', unit_name='m^2/(V*s)')
                        h_total = self.load_create_data_channel(
                            channel_name='hole total',
                            measurement=measurement,
                            description='Total holes mobility', unit_name='m^2/(V*s)')
                        result['hole']['lattice'] = self.client.measurement_manager.get_data_points_array(
                            h_lattice)[:, 0]
                        result['hole']['impurities'] = self.client.measurement_manager.get_data_points_array(
                            h_impurities)[:, 0]
                        result['hole']['total'] = self.client.measurement_manager.get_data_points_array(
                            h_total)[:, 0]
                        db_time = timeit.default_timer() - start_time
                        record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                                 (measurement_details['name'], db_time, 0.0, db_time)
                        self.client.log_manager.log_record(record=record, category='Information')
                        return result
            temperature_channel = self.load_create_data_channel(channel_name='Temperature', measurement=measurement,
                                                                description='Temperature', unit_name='K')
            e_lattice = self.load_create_data_channel(
                channel_name='electron lattice',
                measurement=measurement,
                description='Lattice component of electrons mobility', unit_name='m^2/(V*s)')
            e_impurities = self.load_create_data_channel(
                channel_name='electron impurities',
                measurement=measurement,
                description='Impurities component of electrons mobility', unit_name='m^2/(V*s)')
            e_total = self.load_create_data_channel(
                channel_name='electron total',
                measurement=measurement,
                description='Total electrons mobility', unit_name='m^2/(V*s)')
            h_lattice = self.load_create_data_channel(
                channel_name='hole lattice',
                measurement=measurement,
                description='Lattice component of holes mobility', unit_name='m^2/(V*s)')
            h_impurities = self.load_create_data_channel(
                channel_name='hole impurities',
                measurement=measurement,
                description='Impurities component of holes mobility', unit_name='m^2/(V*s)')
            h_total = self.load_create_data_channel(
                channel_name='hole total',
                measurement=measurement,
                description='Total holes mobility', unit_name='m^2/(V*s)')
        db_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
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
        measurement_time += timeit.default_timer() - start_time
        if use_storage:
            start_time = timeit.default_timer()
            self.client.measurement_manager.create_data_points(channel=temperature_channel,
                                                               float_value=temperature)
            self.client.measurement_manager.create_data_points(channel=e_lattice,
                                                               float_value=result['electron']['lattice'])
            self.client.measurement_manager.create_data_points(channel=h_lattice,
                                                               float_value=result['hole']['lattice'])
            self.client.measurement_manager.create_data_points(channel=e_impurities,
                                                               float_value=result['electron']['impurities'])
            self.client.measurement_manager.create_data_points(channel=h_impurities,
                                                               float_value=result['hole']['impurities'])
            self.client.measurement_manager.create_data_points(channel=e_total,
                                                               float_value=result['electron']['total'])
            self.client.measurement_manager.create_data_points(channel=h_total,
                                                               float_value=result['hole']['total'])
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            db_time += timeit.default_timer() - start_time
            record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                     (measurement_details['name'], db_time + measurement_time,
                      measurement_time, db_time)
            self.client.log_manager.log_record(record=record, category='Information')
        return result

    def electrochemical_potential(self, temperature, use_storage=False):
        start_time = timeit.default_timer()
        if isinstance(temperature, (list, tuple, np.ndarray)):
            temperature = np.array(temperature)
        elif isinstance(temperature, numbers.Number):
            temperature = np.array([np.float(temperature)])
        else:
            raise TypeError('Unsupported temperature type')
        measurement_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        if use_storage:
            measurement_details = {
                'name': 'Electrochemical potential temperature dependence',
                'description': 'Measurement of Semiconductor Electrochemical potential temperature dependence',
                'type': 'Bulk Semiconductor energetics'}
            record = 'Starting Measurement "%s"' % (measurement_details['name'])
            self.client.log_manager.log_record(record=record, category='Information')
            parameters = None
            measurement = self.register_measurement(measurement_details=measurement_details, parameters=parameters,
                                                    input_data=None, force_new=False)
            if measurement.progress == 100:
                temperature_channels = self.client.measurement_manager.get_data_channels(measurement=measurement,
                                                                                         name='Temperature')
                for temperature_channel in temperature_channels:
                    stored = self.client.measurement_manager.get_data_points_array(temperature_channel)[:, 0]
                    if stored.size == temperature.size and (stored == temperature).all():
                        mu_channel = self.load_create_data_channel(channel_name='Electrochemical potential',
                                                                   measurement=measurement,
                                                                   description='Electrochemical potential',
                                                                   unit_name='eV')
                        mu = self.client.measurement_manager.get_data_points_array(mu_channel)[:, 0]
                        db_time = timeit.default_timer() - start_time
                        record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                                 (measurement_details['name'], db_time, 0.0, db_time)
                        self.client.log_manager.log_record(record=record, category='Information')
                        return mu
            temperature_channel = self.load_create_data_channel(channel_name='Temperature', measurement=measurement,
                                                                description='Temperature', unit_name='K')
            mu_channel = self.load_create_data_channel(channel_name='Electrochemical potential',
                                                       measurement=measurement,
                                                       description='Electrochemical potential',
                                                       unit_name='eV')
        db_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        energy_scale = constants['k'] * temperature
        band_gap = self.band_gap(temperature, use_storage=True)
        bands_density_of_states = self.effective_bands_density_of_states(temperature, use_storage=True)
        mu = np.zeros_like(temperature)
        for i in range(len(temperature)):
            if temperature[i] < 8:
                e1 = self.electrochemical_potential(temperature=8, use_storage=True)
                e2 = self.electrochemical_potential(temperature=10, use_storage=True)
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
                    cb_charge = -bands_density_of_states['Nc'][i] * exp_term_cb
                    vb_charge = bands_density_of_states['Nv'][i] * exp_term_vb
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
                # print(mu, part.energy_level(band_gap), energy_scale, n, f, total_charge)
                np.seterr(divide='warn', invalid='warn', over='warn')
                return total_charge
            mu[i] = bisect(equation, a=0, b=band_gap[i])
        measurement_time += timeit.default_timer() - start_time
        if use_storage:
            start_time = timeit.default_timer()
            self.client.measurement_manager.create_data_points(channel=temperature_channel,
                                                               float_value=temperature)
            self.client.measurement_manager.create_data_points(channel=mu_channel,
                                                               float_value=mu)
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            db_time += timeit.default_timer() - start_time
            record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                     (measurement_details['name'], db_time + measurement_time,
                      measurement_time, db_time)
            self.client.log_manager.log_record(record=record, category='Information')
        return mu
