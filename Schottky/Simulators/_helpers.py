from __future__ import division, print_function
import timeit
import numpy as np


def storage_manager(spec_key, **outer_kwargs):
    def measurement_timer(measurement_function):
        def wrapper(self, *args, **kwargs):
            start_time = timeit.default_timer()
            if 'use_storage' in outer_kwargs:
                use_storage = outer_kwargs['use_storage']
            else:
                use_storage = False
            if use_storage:
                record = 'Starting Measurement "%s"' % (self.measurement_details[spec_key]['name'])
                self.client.log_manager.log_record(record=record, category='Information')
                measurements = self.measurement_lookup(self.measurement_details[spec_key],
                                                       parameters=self.measurement_specs[spec_key]['parameters'],
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
                            # TODO: Add channel size retrieval function
                            stored = self.client.measurement_manager.get_data_points_array(channel)[:, 0]
                            #print(channel_spec['name'], arg, stored)
                            #print(stored.size == arg.size, stored == arg, (stored == arg).all())
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
                        result[channel_spec['name']] = self.client.measurement_manager.get_data_points_array(channel)[:, 0]
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
            start_time = timeit.default_timer()
            if use_storage:
                measurement = self.measurement_new(self.measurement_details[spec_key],
                                                   parameters=self.measurement_specs[spec_key]['parameters'],
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