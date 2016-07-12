from __future__ import division, print_function
import numpy as np

from Quaternions import Rotation
from Space.Field import Field, SuperposedField

from Schottky.Samples import Sample


class SampleField(Sample, Field):

    def __init__(self, client, name, description=None, field_type=None, orientation=None):
        Sample.__init__(self, client=client, name=name, description=description)
        self.load_create_sample()
        self.field_type = None
        self._read_in_field_type(field_type=field_type)
        Field.__init__(self, name=name, field_type=self.field_type)
        self.orientation = {'origin': [0, 0, 0],
                            'rotation': [1, 0, 0, 0]}
        self._read_in_orientation(orientation)

    def _read_in_field_type(self, field_type):
        try:
            self.field_type = self.parameters['Field type'].string_value
        except KeyError:
            pass
        if field_type is not None:
            if self.field_type != field_type:
                self.set_field_type(field_type)

    def set_field_type(self, field_type):
        assert isinstance(field_type, str), 'Field type must be string'
        try:
            self.parameters['Field type'].string_value = field_type
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_string_parameter(
                name='Field type',
                value=field_type,
                description='Field type keyword')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.field_type = field_type

    def _read_in_orientation(self, orientation):
        force_set = False
        try:
            origin_components = self.parameters['Origin'].children
            for origin_component in origin_components:
                if origin_component.name == 'x':
                    self.orientation['origin'][0] = origin_component.float_value
                elif origin_component.name == 'y':
                    self.orientation['origin'][1] = origin_component.float_value
                elif origin_component.name == 'z':
                    self.orientation['origin'][2] = origin_component.float_value
        except KeyError:
            force_set = True
            pass
        try:
            quaternion_components = self.parameters['Quaternion'].children
            for quaternion_component in quaternion_components:
                if quaternion_component.name == 'identity':
                    self.orientation['rotation'][0] = quaternion_component.float_value
                elif quaternion_component.name == 'i':
                    self.orientation['rotation'][1] = quaternion_component.float_value
                elif quaternion_component.name == 'j':
                    self.orientation['rotation'][2] = quaternion_component.float_value
                elif quaternion_component.name == 'k':
                    self.orientation['rotation'][3] = quaternion_component.float_value
        except KeyError:
            force_set = True
            pass
        if orientation is not None or force_set:
            if self.orientation != orientation or force_set:
                if orientation is None:
                    orientation = self.orientation
                self.set_orientation(orientation)
        else:
            self._apply_orientation(self.orientation)

    def set_orientation(self, orientation):
        assert isinstance(orientation, dict), 'Orientation must be dictionary with origin and rotation quaternion'
        assert 'origin' in orientation.keys(), 'Orientation must be dictionary with origin and rotation quaternion'
        assert 'rotation' in orientation.keys(), 'Orientation must be dictionary with origin and rotation quaternion'
        assert isinstance(orientation['origin'], (list, tuple, np.ndarray)), 'Origin must be 3 component iterable'
        assert isinstance(orientation['rotation'], (list, tuple, np.ndarray)), 'Rotation must be 4 component iterable'
        assert len(orientation['origin']) == 3, 'Origin must be 3 component iterable'
        assert len(orientation['rotation']) == 4, 'Rotation must be 4 component iterable'
        quaternion = orientation['rotation']
        origin = orientation['origin']
        try:
            origin_components = self.parameters['Origin'].children
            for origin_component in origin_components:
                if origin_component.name == 'x':
                    origin_component.float_value = float(origin[0])
                elif origin_component.name == 'y':
                    origin_component.float_value = float(origin[1])
                elif origin_component.name == 'z':
                    origin_component.float_value = float(origin[2])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Origin',
                                                                            description='Origin point of Field CS')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='x',
                                                                   value=float(origin[0]),
                                                                   description='X component of origin point',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='y',
                                                                   value=float(origin[1]),
                                                                   description='Y component of origin point',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='z',
                                                                   value=float(origin[2]),
                                                                   description='Z component of origin point',
                                                                   parent=parameter)
            self.load_create_sample()
        try:
            quaternion_components = self.parameters['Quaternion'].children
            for quaternion_component in quaternion_components:
                if quaternion_component.name == 'identity':
                    quaternion_component.float_value = float(quaternion[0])
                elif quaternion_component.name == 'i':
                    quaternion_component.float_value = float(quaternion[1])
                elif quaternion_component.name == 'j':
                    quaternion_component.float_value = float(quaternion[2])
                elif quaternion_component.name == 'k':
                    quaternion_component.float_value = float(quaternion[3])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Quaternion',
                                                                            description='Rot. quaternion of Field CS')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='identity',
                                                                   value=float(quaternion[0]),
                                                                   description='1 component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='i',
                                                                   value=float(quaternion[1]),
                                                                   description='i component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='j',
                                                                   value=float(quaternion[2]),
                                                                   description='j component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='k',
                                                                   value=float(quaternion[3]),
                                                                   description='k component of quaternion',
                                                                   parent=parameter)
            self.load_create_sample()
        self.orientation = orientation
        self._apply_orientation(self.orientation)

    def _update_orientation(self):
        rotation = Rotation(euler_angles_convention=self.coordinate_system.euler_angles_convention['title'])
        rotation.rotation_matrix = self.coordinate_system.basis.T
        quaternion = rotation.quadruple
        origin = self.coordinate_system.origin
        orientation = {'origin': origin,
                       'rotation': quaternion}
        self.set_orientation(orientation)

    def _apply_orientation(self, orientation):
        rotation = Rotation(quadruple=orientation['rotation'],
                            euler_angles_convention=self.coordinate_system.euler_angles_convention['title'])
        self.coordinate_system.origin = orientation['origin']
        self.coordinate_system.basis = rotation.rotation_matrix.T

    def rotate(self, rotation, rot_center=None):
        self.coordinate_system.rotate(rotation=rotation, rot_center=rot_center)
        self._update_orientation()

    def rotate_axis_angle(self, axis, theta, rot_center=None):
        self.coordinate_system.rotate_axis_angle(axis=axis, theta=theta, rot_center=rot_center)
        self._update_orientation()

    def rotate_euler_angles(self, euler_angles, rot_center=None):
        self.coordinate_system.rotate_euler_angles(euler_angles=euler_angles, rot_center=rot_center)
        self._update_orientation()

    def set_origin(self, origin):
        self.coordinate_system.origin = origin
        self._update_orientation()


class SuperpositionField(Sample, SuperposedField):

    def __init__(self, client, name, fields=None, description=None, orientation=None):
        Sample.__init__(self, client=client, name=name, description=description, parameters=None)
        self.load_create_sample()
        self.fields = []
        self._read_in_fields(fields=fields)
        SuperposedField.__init__(self, name=name, fields=self.fields)
        self.orientation = {'origin': [0, 0, 0],
                            'rotation': [1, 0, 0, 0]}
        self._read_in_orientation(orientation)

    def _read_in_fields(self, fields):
        try:
            field_components = self.parameters['Fields'].children
            for field in field_components:
                field_module_name, field_class_name = field.string_value.split('::')
                field_name = field.name
                field_id = int(field.float_value)
                field_module = __import__(field_module_name, fromlist=[field_class_name])
                field_class = getattr(field_module, field_class_name)
                field_sample = field_class(client=self.client.session_manager, name=field_name)
                if field_sample.sample.id == field_id:
                    self.fields.append(field_sample)
                else:
                    print('Field IDs do not match')
        except KeyError:
            pass
        if fields is not None:
            if self.fields != fields:
                self.set_fields(fields)

    def set_fields(self, fields):
        assert isinstance(fields, (list, tuple)), 'Expected a list of fields'
        for field in fields:
            assert isinstance(field, Field), 'Expected a list of fields'
            assert isinstance(field, Sample), 'Expected a list of fields'
        try:
            parameter = self.parameters['Fields']
            for field in fields:
                matched = False
                for field_parameter in parameter.children:
                    if field.sample.id == int(field_parameter.float_value):
                        matched = True
                        field_parameter.name = field.name
                        field_parameter.string_value = field.__class__.__module__ + '::' + field.__class__.__name__
                        field_parameter.description = field.description
                        break
                if not matched:
                    self.client.parameter_manager.create_generic_parameter(
                        name=field.name,
                        float_value=float(field.sample.id),
                        string_value=field.__class__.__module__ + '::' + field.__class__.__name__,
                        description=field.description,
                        parent=parameter)
            field_sample_ids = [field.sample.id for field in fields]
            for field_parameter in parameter.children:
                if int(field_parameter.float_value) not in field_sample_ids:
                    self.client.parameter_manager.delete_parameter(field_parameter)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Fields',
                                                                            description='Fields dictionary')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            for field in fields:
                self.client.parameter_manager.create_generic_parameter(
                    name=field.name,
                    float_value=float(field.sample.id),
                    string_value=field.__class__.__module__ + '::' + field.__class__.__name__,
                    description=field.description,
                    parent=parameter)
        self.parameters = {}
        self.load_create_sample()
        self.fields = fields

    def _read_in_orientation(self, orientation):
        force_set = False
        try:
            origin_components = self.parameters['Origin'].children
            for origin_component in origin_components:
                if origin_component.name == 'x':
                    self.orientation['origin'][0] = origin_component.float_value
                elif origin_component.name == 'y':
                    self.orientation['origin'][1] = origin_component.float_value
                elif origin_component.name == 'z':
                    self.orientation['origin'][2] = origin_component.float_value
        except KeyError:
            force_set = True
            pass
        try:
            quaternion_components = self.parameters['Quaternion'].children
            for quaternion_component in quaternion_components:
                if quaternion_component.name == 'identity':
                    self.orientation['rotation'][0] = quaternion_component.float_value
                elif quaternion_component.name == 'i':
                    self.orientation['rotation'][1] = quaternion_component.float_value
                elif quaternion_component.name == 'j':
                    self.orientation['rotation'][2] = quaternion_component.float_value
                elif quaternion_component.name == 'k':
                    self.orientation['rotation'][3] = quaternion_component.float_value
        except KeyError:
            force_set = True
            pass
        if orientation is not None or force_set:
            if self.orientation != orientation or force_set:
                if orientation is None:
                    orientation = self.orientation
                self.set_orientation(orientation)
        else:
            self._apply_orientation(self.orientation)

    def set_orientation(self, orientation):
        assert isinstance(orientation, dict), 'Orientation must be dictionary with origin and rotation quaternion'
        assert 'origin' in orientation.keys(), 'Orientation must be dictionary with origin and rotation quaternion'
        assert 'rotation' in orientation.keys(), 'Orientation must be dictionary with origin and rotation quaternion'
        assert isinstance(orientation['origin'], (list, tuple, np.ndarray)), 'Origin must be 3 component iterable'
        assert isinstance(orientation['rotation'], (list, tuple, np.ndarray)), 'Rotation must be 4 component iterable'
        assert len(orientation['origin']) == 3, 'Origin must be 3 component iterable'
        assert len(orientation['rotation']) == 4, 'Rotation must be 4 component iterable'
        quaternion = orientation['rotation']
        origin = orientation['origin']
        try:
            origin_components = self.parameters['Origin'].children
            for origin_component in origin_components:
                if origin_component.name == 'x':
                    origin_component.float_value = float(origin[0])
                elif origin_component.name == 'y':
                    origin_component.float_value = float(origin[1])
                elif origin_component.name == 'z':
                    origin_component.float_value = float(origin[2])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Origin',
                                                                            description='Origin point of Field CS')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='x',
                                                                   value=float(origin[0]),
                                                                   description='X component of origin point',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='y',
                                                                   value=float(origin[1]),
                                                                   description='Y component of origin point',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='z',
                                                                   value=float(origin[2]),
                                                                   description='Z component of origin point',
                                                                   parent=parameter)
            self.load_create_sample()
        try:
            quaternion_components = self.parameters['Quaternion'].children
            for quaternion_component in quaternion_components:
                if quaternion_component.name == 'identity':
                    quaternion_component.float_value = float(quaternion[0])
                elif quaternion_component.name == 'i':
                    quaternion_component.float_value = float(quaternion[1])
                elif quaternion_component.name == 'j':
                    quaternion_component.float_value = float(quaternion[2])
                elif quaternion_component.name == 'k':
                    quaternion_component.float_value = float(quaternion[3])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_dict_parameter(name='Quaternion',
                                                                            description='Rot. quaternion of Field CS')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='identity',
                                                                   value=float(quaternion[0]),
                                                                   description='1 component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='i',
                                                                   value=float(quaternion[1]),
                                                                   description='i component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='j',
                                                                   value=float(quaternion[2]),
                                                                   description='j component of quaternion',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='k',
                                                                   value=float(quaternion[3]),
                                                                   description='k component of quaternion',
                                                                   parent=parameter)
            self.load_create_sample()
        self.orientation = orientation
        self._apply_orientation(self.orientation)

    def _update_orientation(self):
        rotation = Rotation(euler_angles_convention=self.coordinate_system.euler_angles_convention['title'])
        rotation.rotation_matrix = self.coordinate_system.basis.T
        quaternion = rotation.quadruple
        origin = self.coordinate_system.origin
        orientation = {'origin': origin,
                       'rotation': quaternion}
        self.set_orientation(orientation)

    def _apply_orientation(self, orientation):
        rotation = Rotation(quadruple=orientation['rotation'],
                            euler_angles_convention=self.coordinate_system.euler_angles_convention['title'])
        self.coordinate_system.origin = orientation['origin']
        self.coordinate_system.basis = rotation.rotation_matrix.T

    def rotate(self, rotation, rot_center=None):
        self.coordinate_system.rotate(rotation=rotation, rot_center=rot_center)
        self._update_orientation()

    def rotate_axis_angle(self, axis, theta, rot_center=None):
        self.coordinate_system.rotate_axis_angle(axis=axis, theta=theta, rot_center=rot_center)
        self._update_orientation()

    def rotate_euler_angles(self, euler_angles, rot_center=None):
        self.coordinate_system.rotate_euler_angles(euler_angles=euler_angles, rot_center=rot_center)
        self._update_orientation()

    def set_origin(self, origin):
        self.coordinate_system.origin = origin
        self._update_orientation()
