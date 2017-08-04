from __future__ import division, print_function
import timeit
import numbers
import numpy as np

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators import Simulator
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky import constants


class SchottkyDiodeSimulator(Simulator):

    def __init__(self, client, diode, description=None):
        assert isinstance(diode, SchottkyDiode), 'Valid SchottkyDiode Sample object expected'
        self.diode = diode
        samples = [self.diode]
        name = 'Schottky Diode Simulator'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}
        measurement_types = [
            {'name': 'Diode energetics',
             'description': 'Measurement of diode energetics',
             'children': []},
            {'name': 'Current Measurement',
             'description': 'Measurement of electrical current',
             'children': []}]
        measurements = [
            {'name': 'Built-in potential temperature dependence',
             'description': 'Measurement of Diode Built-in potential temperature dependence',
             'type': 'Diode energetics'},
            {'name': 'IV',
             'description': 'measure I-V dependence',
             'type': 'Current Measurement'},
            {'name': 'Emission rate',
             'description': 'measure charge carriers emission rate',
             'type': 'Capture and Emission Kinetics'},
        ]
        parts = [BulkSemiconductor(client=client, semiconductor=diode.semiconductor)]
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=measurements)

    def v_bi(self, temperature, use_storage=False):
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
                'name': 'Built-in potential temperature dependence',
                'description': 'Measurement of Diode Built-in potential temperature dependence',
                'type': 'Diode energetics'}
            record = 'Starting Measurement "%s"' % (measurement_details['name'])
            self.client.log_manager.log_record(record=record, category='Information')
            measurements = self.client.measurement_manager.get_measurements(name=measurement_details['name'])
            for measurement in measurements:
                if measurement.progress != 100:
                    continue
                temperature_channels = self.client.measurement_manager.get_data_channels(measurement=measurement,
                                                                                         name='Temperature')
                for temperature_channel in temperature_channels:
                    stored = self.client.measurement_manager.get_data_points_array(temperature_channel)[:, 0]
                    if stored.size == temperature.size and (stored == temperature).all():
                        v_bi_channel = self.load_create_data_channel(channel_name='Built-in potential',
                                                                     measurement=measurement,
                                                                     description='Built-in potential of the diode',
                                                                     unit_name='eV')
                        bi_voltage = self.client.measurement_manager.get_data_points_array(v_bi_channel)[:, 0]
                        db_time = timeit.default_timer() - start_time
                        record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                                 (measurement_details['name'], db_time, 0.0, db_time)
                        self.client.log_manager.log_record(record=record, category='Information')
                        return bi_voltage
            parameters = None
            measurement = self.register_measurement(measurement_details=measurement_details,
                                                    parameters=parameters, input_data=None,
                                                    force_new=True)
            temperature_channel = self.load_create_data_channel(channel_name='Temperature', measurement=measurement,
                                                                description='Temperature', unit_name='K')
            v_bi_channel = self.load_create_data_channel(channel_name='Built-in potential',
                                                         measurement=measurement,
                                                         description='Built-in potential of the diode',
                                                         unit_name='eV')
        db_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        metal_wf = self.diode.metal.work_function
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                  use_storage=True)
        semiconductor_wf = xi + self.diode.semiconductor.affinity
        bi_voltage = metal_wf - semiconductor_wf
        measurement_time += timeit.default_timer() - start_time
        if use_storage:
            start_time = timeit.default_timer()
            self.client.measurement_manager.create_data_points(channel=temperature_channel,
                                                               float_value=temperature)
            self.client.measurement_manager.create_data_points(channel=v_bi_channel,
                                                               float_value=bi_voltage)
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            db_time += timeit.default_timer() - start_time
            record = 'Measurement "%s" complete in %.3f s (measurement: %.3f s, db: %.3f s)' % \
                     (measurement_details['name'], db_time + measurement_time,
                      measurement_time, db_time)
            self.client.log_manager.log_record(record=record, category='Information')
        return bi_voltage

    def thermionic_emission_current(self, bias, temperature):
        pass
