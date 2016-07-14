from __future__ import division, print_function
import timeit
import numpy as np

from Quaternions import Rotation
from Space.Coordinates import Cartesian, transforms as gt
from Space.Field import Field
from Schottky.Samples import Sample
from Schottky.Samples.Fields import SuperpositionField
from Schottky.Simulators import Simulator


class FieldSimulator(Simulator, Field):

    def __init__(self, client, field, description=None):
        assert isinstance(field, Sample), 'Valid Field Sample object expected'
        assert isinstance(field, Field), 'Valid Field Sample object expected'
        self.field = field
        samples = [self.field]
        if isinstance(self.field, SuperpositionField):
            for field_component in self.field.fields:
                samples.append(field_component)
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
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=None,
            category=category,
            measurement_types=measurement_types,
            measurements=measurements)
        Field.__init__(self, name=name, field_type=self.field.type)

    def scalar_field(self, xyz, length_unit='cm'):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Scalar Field values at given XYZ points',
            'description': 'measure scalar field on given points array',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        parameters = None
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=parameters, input_data=None,
                                                force_new=False)
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
        matches = np.array([False, False, False])
        stored = self.client.measurement_manager.get_data_points_array(x_channel)[:, 0]
        if stored.size == xyz[:, 0].size and (stored == xyz[:, 0]).all():
            matches[0] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=x_channel)
            self.client.measurement_manager.create_data_points(channel=x_channel, float_value=xyz[:, 0])
        stored = self.client.measurement_manager.get_data_points_array(y_channel)[:, 0]
        if stored.size == xyz[:, 1].size and (stored == xyz[:, 1]).all():
            matches[1] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=y_channel)
            self.client.measurement_manager.create_data_points(channel=y_channel, float_value=xyz[:, 1])
        stored = self.client.measurement_manager.get_data_points_array(z_channel)[:, 0]
        if stored.size == xyz[:, 2].size and (stored == xyz[:, 2]).all():
            matches[2] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=z_channel)
            self.client.measurement_manager.create_data_points(channel=z_channel, float_value=xyz[:, 2])
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(scalar_field_channel)[:, 0]
            scalar_field = stored.reshape(xyz[:, 0].shape).astype(np.float)
        else:
            self.client.measurement_manager.delete_data_points(channel=scalar_field_channel)
            scalar_field = self.field.scalar_field(xyz)
            self.client.measurement_manager.create_data_points(channel=scalar_field_channel,
                                                               float_value=scalar_field.ravel())
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            self.client.measurement_manager.finish_measurement(measurement=measurement)
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return scalar_field

    def vector_field(self, xyz, length_unit='cm'):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Vector Field values at given XYZ points',
            'description': 'measure vector field on given points array',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        parameters = None
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=parameters, input_data=None,
                                                force_new=False)
        x_channel = self.client.measurement_manager.create_data_channel(
            name='X points', measurement=measurement,
            description='X coordinate points', unit_name=length_unit)
        y_channel = self.client.measurement_manager.create_data_channel(
            name='Y points', measurement=measurement,
            description='Y coordinate points', unit_name=length_unit)
        z_channel = self.client.measurement_manager.create_data_channel(
            name='Z points', measurement=measurement,
            description='Z coordinate points', unit_name=length_unit)
        vector_field_channel = self.client.measurement_manager.create_data_channel(
            name='Vector field', measurement=measurement,
            description='Vector field values', unit_name=None)
        matches = np.array([False, False, False])
        stored = self.client.measurement_manager.get_data_points_array(x_channel)[:, 0]
        if stored.size == xyz[:, 0].size and (stored == xyz[:, 0]).all():
            matches[0] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=x_channel)
            self.client.measurement_manager.create_data_points(channel=x_channel, float_value=xyz[:, 0])
        stored = self.client.measurement_manager.get_data_points_array(y_channel)[:, 0]
        if stored.size == xyz[:, 1].size and (stored == xyz[:, 1]).all():
            matches[1] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=y_channel)
            self.client.measurement_manager.create_data_points(channel=y_channel, float_value=xyz[:, 1])
        stored = self.client.measurement_manager.get_data_points_array(z_channel)[:, 0]
        if stored.size == xyz[:, 2].size and (stored == xyz[:, 2]).all():
            matches[2] = True
        else:
            self.client.measurement_manager.delete_data_points(channel=z_channel)
            self.client.measurement_manager.create_data_points(channel=z_channel, float_value=xyz[:, 2])
        if matches.all():
            stored = self.client.measurement_manager.get_data_points_array(vector_field_channel)[:, 0]
            vector_field = stored.reshape(xyz.shape).astype(np.float)
        else:
            self.client.measurement_manager.delete_data_points(channel=vector_field_channel)
            vector_field = self.field.vector_field(xyz)
            self.client.measurement_manager.create_data_points(channel=vector_field_channel,
                                                               float_value=vector_field.ravel())
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            self.client.measurement_manager.finish_measurement(measurement=measurement)
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return vector_field

    def measure_field_cartesian_coordinates(self, x_range, y_range, z_range, length_unit='m',
                                            force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on cartesian grid',
            'description': 'measure vector and scalar field on given cartesian grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        parameters = None
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=parameters, input_data=None,
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
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            self.client.measurement_manager.finish_measurement(measurement=measurement)
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return x_grid, y_grid, z_grid, scalar_field, vector_field

    def measure_field_cylindrical_coordinates(self, r_range, phi_range, z_range, length_unit='m',
                                              frame_of_view=None,
                                              force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on cylindrical grid',
            'description': 'measure vector and scalar field on given cylindrical grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')

        if frame_of_view is None:
            frame_of_view = Cartesian()
        else:
            if not isinstance(frame_of_view, Cartesian):
                frame_of_view = Cartesian()
        fov_parameter = self.coordinate_system_to_parameter(frame_of_view, parameter_name='frame_of_view')
        parameters = [fov_parameter]
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=parameters, input_data=None,
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
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            self.client.measurement_manager.finish_measurement(measurement=measurement)
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return r_grid, phi_grid, z_grid, scalar_field, vector_field

    def measure_field_spherical_coordinates(self, r_range, theta_range, phi_range, length_unit='m',
                                            force_recalculate=False):
        start_time = timeit.default_timer()
        measurement_details = {
            'name': 'Field values on spherical grid',
            'description': 'measure vector and scalar field on given spherical grid',
            'type': 'Space fields measurement'}
        record = 'Starting Measurement "%s"' % (measurement_details['name'])
        self.client.log_manager.log_record(record=record, category='Information')
        parameters = None
        measurement = self.register_measurement(measurement_details=measurement_details,
                                                parameters=parameters, input_data=None,
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
            self.client.measurement_manager.update_measurement_progress(measurement=measurement,
                                                                        progress=100)
            self.client.measurement_manager.finish_measurement(measurement=measurement)
        elapsed = timeit.default_timer() - start_time
        record = 'Measurement "%s" complete in %3.3f s' % (measurement_details['name'], elapsed)
        self.client.log_manager.log_record(record=record, category='Information')
        return r_grid, theta_grid, phi_grid, scalar_field, vector_field

    def coordinate_system_to_parameter(self, coordinate_system, parameter_name='coordinate_system'):
        assert isinstance(coordinate_system, Cartesian)
        rotation = Rotation(euler_angles_convention=coordinate_system.euler_angles_convention['title'])
        rotation.rotation_matrix = coordinate_system.basis.T
        quaternion = rotation.quadruple
        origin = coordinate_system.origin
        parameter = self.client.parameter_manager.create_dict_parameter(name=parameter_name,
                                                                        commit=False)
        origin_parameter = self.client.parameter_manager.create_dict_parameter(
            name='Origin',
            description='Origin point of Field CS',
            parent=parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='x',
                                                               value=float(origin[0]),
                                                               description='X component of origin point',
                                                               parent=origin_parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='y',
                                                               value=float(origin[1]),
                                                               description='Y component of origin point',
                                                               parent=origin_parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='z',
                                                               value=float(origin[2]),
                                                               description='Z component of origin point',
                                                               parent=origin_parameter, commit=False)
        quaternion_parameter = self.client.parameter_manager.create_dict_parameter(
            name='Quaternion',
            description='Rot. quaternion of Field CS', parent=parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='identity',
                                                               value=float(quaternion[0]),
                                                               description='1 component of quaternion',
                                                               parent=quaternion_parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='i',
                                                               value=float(quaternion[1]),
                                                               description='i component of quaternion',
                                                               parent=quaternion_parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='j',
                                                               value=float(quaternion[2]),
                                                               description='j component of quaternion',
                                                               parent=quaternion_parameter, commit=False)
        self.client.parameter_manager.create_numeric_parameter(name='k',
                                                               value=float(quaternion[3]),
                                                               description='k component of quaternion',
                                                               parent=quaternion_parameter, commit=False)
        return parameter
