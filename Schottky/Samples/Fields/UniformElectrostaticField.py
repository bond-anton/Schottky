from __future__ import division, print_function
import numpy as np
import numbers

from Space.Coordinates import transforms as gt
from Space.Field import Field

from Schottky.Samples import Sample


class UniformElectrostaticField(Field, Sample):

    def __init__(self, client, strength=None, direction=None,
                 name='Uniform electrostatic field', description=None, field_type='electrostatic'):
        assert isinstance(strength, numbers.Number)
        Field.__init__(self, name=name, field_type=field_type)
        Sample.__init__(self, client=client, name=name, description=description)
        self.strength = None
        self.direction = [None, None, None]
        self._read_in_strength(strength=strength)
        self._read_in_direction(direction=direction)

    def _read_in_strength(self, strength):
        try:
            self.strength = self.parameters['Field strength'].float_value
        except KeyError:
            pass
        if self.strength != strength and strength is not None:
            self.set_strength(strength)

    def set_strength(self, strength):
        assert isinstance(strength, numbers.Number), 'Field strength spread must be a number'
        try:
            self.parameters['Field strength'].float_value = float(strength)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Field strength',
                                                                               value=float(strength),
                                                                               unit_name='V/cm',
                                                                               description='Strength of electric field')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.strength = float(strength)

    def _read_in_direction(self, direction):
        try:
            field_direction_components = self.parameters['Field direction'].children
            for field_direction in field_direction_components:
                if field_direction.name == 'x':
                    self.direction[0] = field_direction.float_value
                elif field_direction.name == 'y':
                    self.direction[1] = field_direction.float_value
                elif field_direction.name == 'z':
                    self.direction[2] = field_direction.float_value
                self.direction = gt.unit_vector(self.direction)
        except KeyError:
            pass
        if direction is not None:
            if (self.direction != gt.unit_vector(direction)).any():
                self.set_field_direction(direction)

    def set_field_direction(self, direction):
        assert isinstance(direction, (np.ndarray, list, tuple)), 'Direction must be iterable of size 3'
        assert len(direction) == 3, 'Direction must be iterable of size 3'
        direction = gt.unit_vector(direction)
        try:
            field_direction_components = self.parameters['Field direction'].children
            for field_direction in field_direction_components:
                if field_direction.name == 'x':
                    field_direction.float_value = float(direction[0])
                elif field_direction.name == 'y':
                    field_direction.float_value = float(direction[1])
                elif field_direction.name == 'z':
                    field_direction.float_value = float(direction[1])
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager._create_parameter(name='Field direction',
                                                                        description='Field direction vector',
                                                                        parameter_type='Dictionary')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='x',
                                                                   value=float(direction[0]),
                                                                   description='X component of field direction vector',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='y',
                                                                   value=float(direction[1]),
                                                                   description='Y component of field direction vector',
                                                                   parent=parameter)
            self.client.parameter_manager.create_numeric_parameter(name='z',
                                                                   value=float(direction[2]),
                                                                   description='Z component of field direction vector',
                                                                   parent=parameter)
            self.load_create_sample()
        self.direction = direction

    def scalar_field(self, xyz):
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            distance = np.dot(xyz, self.direction)
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                distance = np.dot(xyz, self.direction)
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        return -self.strength * distance

    def vector_field(self, xyz):
        xyz = np.array(xyz, dtype=np.float)
        if xyz.size == 3:
            field = np.array([1, 1, 1]) * self.direction
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                field = np.ones_like(xyz) * self.direction
            else:
                raise ValueError('N-points array shape must be (N, 3)')
        else:
            raise ValueError('at least 3 coordinates are needed for point')
        return self.strength * field
