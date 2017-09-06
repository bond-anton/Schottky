from __future__ import division, print_function
import timeit
import numpy as np

from Schottky import constants


def storage_manager(spec_key, **outer_kwargs):
    def measurement_timer(measurement_function):
        def wrapper(self, *args, **kwargs):
            start_time = timeit.default_timer()
            if 'use_storage' in kwargs:
                use_storage = kwargs.pop('use_storage')
            elif 'use_storage' in outer_kwargs:
                use_storage = outer_kwargs['use_storage']
            else:
                use_storage = False
            print('use_storage:', use_storage)
            if use_storage:
                record = 'Starting Measurement "%s"' % (self.measurement_details[spec_key]['name'])
                parameters = []
                if self.measurement_specs[spec_key]['parameters'] is not None:
                    for parameter_spec in self.measurement_specs[spec_key]['parameters']:
                        if parameter_spec['name'] in kwargs:
                            if kwargs[parameter_spec['name']] is None:
                                parameter_value = parameter_spec['default value']
                            else:
                                parameter_value = kwargs[parameter_spec['name']]
                            if parameter_spec['type'] == 'numeric':
                                parameter = self.client.parameter_manager.create_numeric_parameter(
                                    parameter_spec['name'], parameter_value,
                                    unit_name=parameter_spec['units'],
                                    description=parameter_spec['description'],
                                    parent=None, commit=False)
                                parameters.append(parameter)
                self.client.log_manager.log_record(record=record, category='Information')
                measurements = self.measurement_lookup(self.measurement_details[spec_key],
                                                       parameters=parameters,
                                                       input_data=self.measurement_specs[spec_key]['input data'],
                                                       progress_threshold=100)
                match_found = False
                for measurement in measurements:
                    for channel_spec in self.measurement_specs[spec_key]['variables']:
                        arg = prepare_array(kwargs[channel_spec['name'].lower()])
                        channels = self.client.measurement_manager.get_data_channels(measurement=measurement,
                                                                                     name=channel_spec['name'])
                        match_found = False
                        for channel in channels:
                            stored_size = self.client.measurement_manager.get_data_points_num(channel)
                            if stored_size != arg.size:
                                continue
                            stored = self.client.measurement_manager.get_data_points_array(channel)[:, 0]
                            if stored.size == arg.size and (stored == arg).all():
                                match_found = True
                                break
                        if not match_found:
                            break
                    if match_found:
                        break
                if match_found:
                    result = {}
                    for channel_spec in self.measurement_specs[spec_key]['result']:
                        channel = self.load_create_data_channel(channel_name=channel_spec['name'],
                                                                measurement=measurement,
                                                                description=channel_spec['description'],
                                                                unit_name=channel_spec['units'])
                        result[channel_spec['name']] = self.client.measurement_manager.get_data_points_array(
                            channel)[:, 0].astype(np.float64)
                    db_time = timeit.default_timer() - start_time
                    record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' %\
                             (self.measurement_details[spec_key]['name'], db_time, 0.0, db_time)
                    self.client.log_manager.log_record(record=record, category='Information')
                    if len(result) == 1:
                        return result[channel_spec['name']]
                    return result
            db_time = timeit.default_timer() - start_time
            start_time = timeit.default_timer()
            result = measurement_function(self, *args, **kwargs)
            measurement_time = timeit.default_timer() - start_time
            if use_storage:
                start_time = timeit.default_timer()
                for parameter in parameters:
                    # TODO: add commit_parameter method in parameter manager for better logging
                    self.client.session.add(parameter)
                    self.client.session.commit()
                measurement = self.measurement_new(self.measurement_details[spec_key],
                                                   parameters=parameters,
                                                   input_data=self.measurement_specs[spec_key]['input data'])
                for channel_spec in self.measurement_specs[spec_key]['variables']:
                    arg = prepare_array(kwargs[channel_spec['name'].lower()])
                    channel = self.load_create_data_channel(channel_name=channel_spec['name'], measurement=measurement,
                                                            description=channel_spec['description'],
                                                            unit_name=channel_spec['units'])
                    self.client.measurement_manager.create_data_points(channel=channel, float_value=arg)

                for channel_spec in self.measurement_specs[spec_key]['result']:
                    channel = self.load_create_data_channel(channel_name=channel_spec['name'],
                                                            measurement=measurement,
                                                            description=channel_spec['description'],
                                                            unit_name=channel_spec['units'])
                    if len(self.measurement_specs[spec_key]['result']) == 1:
                        self.client.measurement_manager.create_data_points(channel=channel,
                                                                           float_value=result)
                    else:
                        self.client.measurement_manager.create_data_points(channel=channel,
                                                                           float_value=result[channel_spec['name']])
                self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                            progress=100)
                db_time += timeit.default_timer() - start_time
                record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                         (self.measurement_details[spec_key]['name'], db_time + measurement_time, measurement_time, db_time)
                self.client.log_manager.log_record(record=record, category='Information')
            return result
        return wrapper
    return measurement_timer


def prepare_array(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return np.array(x, dtype=np.float64)
    else:
        return np.array([x], dtype=np.float64)


def fermi(energy, fermi_energy, temperature, g_ratio=1):
    """
    fermi distribution function
    :param energy: Energy level measured from Ec towards Ev
    :param fermi_energy: Fermi level energy measured from Ec towards Ev
    :param temperature: Temperature in K
    :param g_ratio: Degeneracy ratio
    :return: Fermi distrubution function
    """
    energy_scale = constants['k'] * temperature
    return 1 / (1 + g_ratio * np.exp((fermi_energy - energy) / energy_scale))


def d_fermi_d_delta_fermi_energy(energy, fermi_energy, temperature, g_ratio=1):
    """
    Derivative of fermi distribution function by delta_fermi_energy at delta_fermi_energy = 0
    :param energy: Energy level
    :param fermi_energy: Fermi level enrgy
    :param temperature: Temperature in K
    :param g_ratio: Degeneracy ratio
    :return: sympy symbolic expression
    """
    energy_scale = constants['k'] * temperature
    numerator = g_ratio * np.exp((energy - fermi_energy) / energy_scale) / energy_scale
    denominator = (1 + g_ratio * np.exp((energy - fermi_energy) / energy_scale)) ** 2
    return numerator / denominator
