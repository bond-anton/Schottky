from __future__ import division, print_function
import timeit
import numpy as np

from Space.Coordinates import transforms as gt
from Space.Field import Field
from Schottky.Samples import Sample
from Schottky.Simulators import Simulator


class FieldSimulator(Simulator):

    def __init__(self, client, field, description=None):
        assert isinstance(field, Sample), 'Valid Field Sample object expected'
        assert isinstance(field, Field), 'Valid Field Sample object expected'
        self.field = field
        name = 'Field Simulator'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}
        measurement_types = [{'name': 'Space fields measurement',
                              'description': 'Measurement strength and potentials of fields in space',
                              'children': []}]
        measurements = [
            {'name': 'Field values on cartesian grid',
             'description': 'measure vector and scalar field on given cartesian grid',
             'type': 'Space fields measurement'},
            {'name': 'Field values on cylindrical grid',
             'description': 'measure vector and scalar field on given cylindrical grid',
             'type': 'Space fields measurement'},
            {'name': 'Field values on spherical grid',
             'description': 'measure vector and scalar field on given spherical grid',
             'type': 'Space fields measurement'},
            ]
        super(FieldSimulator, self).__init__(
            client=client, name=name, description=description,
            samples=[self.field], parts=None,
            category=category,
            measurement_types=measurement_types,
            measurements=measurements)

    def measure_field_cartesian_coordinates(self, x_range, y_range, z_range,
                                            length_unit='m', force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on cartesian grid',
            'description': 'measure vector and scalar field on given cartesian grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=None, input_data=None,
                                                force_new=force_recalculate)
        x_channel = self.client.measurement_manager.create_data_channel(
            name='X points', measurement=measurement,
            description='X coordinate points', unit_name=length_unit)
        y_channel = self.client.measurement_manager.create_data_channel(
            name='Y points', measurement=measurement,
            description='Y coordinate points', unit_name=length_unit)
        z_channel = self.client.measurement_manager.create_data_channel(
            name='Z points', measurement=measurement,
            description='Z coordinate points', unit_name=length_unit)
        scalar_field_channel = self.client.measurement_manager.create_data_channel(
            name='Scalar field', measurement=measurement,
            description='Scalar field values', unit_name=None)
        vector_field_channel = self.client.measurement_manager.create_data_channel(
            name='Vector field', measurement=measurement,
            description='Vector field values', unit_name=None)
        matches = np.array([False, False, False])
        stored = self.client.measurement_manager.get_data_points_array(x_channel)[:, 0]
        if stored.size == x_range.size and (stored == x_range).all():
            matches[0] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=x_channel)
            self.client.measurement_manager.create_data_points(channel=x_channel, float_value=x_range)
        stored = self.client.measurement_manager.get_data_points_array(y_channel)[:, 0]
        if stored.size == y_range.size and (stored == y_range).all():
            matches[1] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=y_channel)
            self.client.measurement_manager.create_data_points(channel=y_channel, float_value=y_range)
        stored = self.client.measurement_manager.get_data_points_array(z_channel)[:, 0]
        if stored.size == z_range.size and (stored == z_range).all():
            matches[2] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=z_channel)
            self.client.measurement_manager.create_data_points(channel=z_channel, float_value=z_range)
        x_grid, y_grid, z_grid = np.meshgrid(x_range, y_range, z_range)
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(scalar_field_channel)[:, 0]
            scalar_field = stored.reshape((len(y_range), len(x_range), len(z_range))).astype(np.float)
            stored = self.client.measurement_manager.get_data_points_array(vector_field_channel)[:, 0]
            vector_field = stored.reshape((len(y_range), len(x_range), len(z_range), 3)).astype(np.float)
        else:
            self.client.measurement_manager.delete_data_points(channel=scalar_field_channel)
            self.client.measurement_manager.delete_data_points(channel=vector_field_channel)
            xyz = np.vstack([x_grid.ravel(), y_grid.ravel(), z_grid.ravel()]).T
            scalar_field = self.field.scalar_field(xyz).reshape((len(y_range), len(x_range), len(z_range)))
            self.client.measurement_manager.create_data_points(channel=scalar_field_channel,
                                                               float_value=scalar_field.ravel())
            vector_field = self.field.vector_field(xyz).reshape((len(y_range), len(x_range), len(z_range), 3))
            self.client.measurement_manager.create_data_points(channel=vector_field_channel,
                                                               float_value=vector_field.ravel())
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return x_grid, y_grid, z_grid, scalar_field, vector_field

    def measure_field_cylindrical_coordinates(self, r_range, phi_range, z_range,
                                              length_unit='m', force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on cylindrical grid',
            'description': 'measure vector and scalar field on given cylindrical grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=None, input_data=None,
                                                force_new=force_recalculate)
        r_channel = self.client.measurement_manager.create_data_channel(
            name='R points', measurement=measurement,
            description='R coordinate points', unit_name=length_unit)
        phi_channel = self.client.measurement_manager.create_data_channel(
            name='Phi points', measurement=measurement,
            description='Phi coordinate points', unit_name='rad')
        z_channel = self.client.measurement_manager.create_data_channel(
            name='Z points', measurement=measurement,
            description='Z coordinate points', unit_name=length_unit)
        scalar_field_channel = self.client.measurement_manager.create_data_channel(
            name='Scalar field', measurement=measurement,
            description='Scalar field values', unit_name=None)
        vector_field_channel = self.client.measurement_manager.create_data_channel(
            name='Vector field', measurement=measurement,
            description='Vector field values', unit_name=None)
        matches = np.array([False, False, False])
        stored = self.client.measurement_manager.get_data_points_array(r_channel)[:, 0]
        if stored.size == r_range.size and (stored == r_range).all():
            matches[0] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=r_channel)
            self.client.measurement_manager.create_data_points(channel=r_channel, float_value=r_range)
        stored = self.client.measurement_manager.get_data_points_array(phi_channel)[:, 0]
        if stored.size == phi_range.size and (stored == phi_range).all():
            matches[1] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=phi_channel)
            self.client.measurement_manager.create_data_points(channel=phi_channel, float_value=phi_range)
        stored = self.client.measurement_manager.get_data_points_array(z_channel)[:, 0]
        if stored.size == z_range.size and (stored == z_range).all():
            matches[2] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=z_channel)
            self.client.measurement_manager.create_data_points(channel=z_channel, float_value=z_range)
        r_grid, phi_grid, z_grid = np.meshgrid(r_range, phi_range, z_range)
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(scalar_field_channel)[:, 0]
            scalar_field = stored.reshape((len(phi_range), len(r_range), len(z_range))).astype(np.float)
            stored = self.client.measurement_manager.get_data_points_array(vector_field_channel)[:, 0]
            vector_field = stored.reshape((len(phi_range), len(r_range), len(z_range), 3)).astype(np.float)
        else:
            self.client.measurement_manager.delete_data_points(channel=scalar_field_channel)
            self.client.measurement_manager.delete_data_points(channel=vector_field_channel)
            positions = np.vstack([r_grid.ravel(), phi_grid.ravel(), z_grid.ravel()]).T
            xyz = gt.cylindrical_to_cartesian(positions)
            scalar_field = self.field.scalar_field(xyz).reshape((len(phi_range), len(r_range), len(z_range)))
            self.client.measurement_manager.create_data_points(channel=scalar_field_channel,
                                                               float_value=scalar_field.ravel())
            vector_field = self.field.vector_field(xyz).reshape((len(phi_range), len(r_range), len(z_range), 3))
            self.client.measurement_manager.create_data_points(channel=vector_field_channel,
                                                               float_value=vector_field.ravel())
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return r_grid, phi_grid, z_grid, scalar_field, vector_field

    def measure_field_spherical_coordinates(self, r_range, theta_range, phi_range,
                                            length_unit='m', force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on spherical grid',
            'description': 'measure vector and scalar field on given spherical grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=None, input_data=None,
                                                force_new=force_recalculate)
        r_channel = self.client.measurement_manager.create_data_channel(
            name='R points', measurement=measurement,
            description='R coordinate points', unit_name=length_unit)
        theta_channel = self.client.measurement_manager.create_data_channel(
            name='Theta points', measurement=measurement,
            description='Theta coordinate points', unit_name='rad')
        phi_channel = self.client.measurement_manager.create_data_channel(
            name='Phi points', measurement=measurement,
            description='Phi coordinate points', unit_name='rad')
        scalar_field_channel = self.client.measurement_manager.create_data_channel(
            name='Scalar field', measurement=measurement,
            description='Scalar field values', unit_name=None)
        vector_field_channel = self.client.measurement_manager.create_data_channel(
            name='Vector field', measurement=measurement,
            description='Vector field values', unit_name=None)
        matches = np.array([False, False, False])
        stored = self.client.measurement_manager.get_data_points_array(r_channel)[:, 0]
        if stored.size == r_range.size and (stored == r_range).all():
            matches[0] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=r_channel)
            self.client.measurement_manager.create_data_points(channel=r_channel, float_value=r_range)
        stored = self.client.measurement_manager.get_data_points_array(theta_channel)[:, 0]
        if stored.size == theta_range.size and (stored == theta_range).all():
            matches[1] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=theta_channel)
            self.client.measurement_manager.create_data_points(channel=theta_channel, float_value=theta_range)
        stored = self.client.measurement_manager.get_data_points_array(phi_channel)[:, 0]
        if stored.size == phi_range.size and (stored == phi_range).all():
            matches[2] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=phi_channel)
            self.client.measurement_manager.create_data_points(channel=phi_channel, float_value=phi_range)
        r_grid, theta_grid, phi_grid = np.meshgrid(r_range, theta_range, phi_range)
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(scalar_field_channel)[:, 0]
            scalar_field = stored.reshape((len(theta_range), len(r_range), len(phi_range))).astype(np.float)
            stored = self.client.measurement_manager.get_data_points_array(vector_field_channel)[:, 0]
            vector_field = stored.reshape((len(theta_range), len(r_range), len(phi_range), 3)).astype(np.float)
        else:
            self.client.measurement_manager.delete_data_points(channel=scalar_field_channel)
            self.client.measurement_manager.delete_data_points(channel=vector_field_channel)
            positions = np.vstack([r_grid.ravel(), theta_grid.ravel(), phi_grid.ravel()]).T
            xyz = gt.spherical_to_cartesian(positions)
            scalar_field = self.field.scalar_field(xyz).reshape((len(theta_range), len(r_range), len(phi_range)))
            self.client.measurement_manager.create_data_points(channel=scalar_field_channel,
                                                               float_value=scalar_field.ravel())
            vector_field = self.field.vector_field(xyz).reshape((len(theta_range), len(r_range), len(phi_range), 3))
            self.client.measurement_manager.create_data_points(channel=vector_field_channel,
                                                               float_value=vector_field.ravel())
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return r_grid, theta_grid, phi_grid, scalar_field, vector_field