from __future__ import division, print_function
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

    def measure_field_cylindrical_coordinates(self, r_range, phi_range, z_range,
                                              length_unit='m', force_recalculate=False):
        measurement_details = {
            'name': 'Field values on cylindrical grid',
            'description': 'measure vector and scalar field on given cylindrical grid',
            'type': 'Space fields measurement'}
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
            self.client.measurement_manager.create_data_points(channel=r_channel, float_value=r_range)
        stored = self.client.measurement_manager.get_data_points_array(phi_channel)[:, 0]
        if stored.size == phi_range.size and (stored == phi_range).all():
            matches[1] = True
        else:
            self.client.measurement_manager.create_data_points(channel=phi_channel, float_value=phi_range)
        stored = self.client.measurement_manager.get_data_points_array(z_channel)[:, 0]
        if stored.size == z_range.size and (stored == z_range).all():
            matches[2] = True
        else:
            self.client.measurement_manager.create_data_points(channel=z_channel, float_value=z_range)
        r_grid, phi_grid, z_grid = np.meshgrid(r_range, phi_range, z_range)
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(scalar_field_channel)[:, 0]
            scalar_field = stored.reshape((len(phi_range), len(r_range), len(z_range))).astype(np.float)
            stored = self.client.measurement_manager.get_data_points_array(vector_field_channel)[:, 0]
            vector_field = stored.reshape((len(phi_range), len(r_range), len(z_range), 3)).astype(np.float)
        else:
            positions = np.vstack([r_grid.ravel(), phi_grid.ravel(), z_grid.ravel()]).T
            xyz = gt.cylindrical_to_cartesian(positions)
            scalar_field = self.field.scalar_field(xyz).reshape((len(phi_range), len(r_range), len(z_range)))
            self.client.measurement_manager.create_data_points(channel=scalar_field_channel,
                                                               float_value=scalar_field.ravel())
            vector_field = self.field.vector_field(xyz).reshape((len(phi_range), len(r_range), len(z_range), 3))
            self.client.measurement_manager.create_data_points(channel=vector_field_channel,
                                                               float_value=vector_field.ravel())
        return r_grid, phi_grid, z_grid, scalar_field, vector_field
