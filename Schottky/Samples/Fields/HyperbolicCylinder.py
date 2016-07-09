from __future__ import division, print_function
import numpy as np
import numbers

from Space.Coordinates import transforms as gt
from Space.Field import Field

from Schottky.Samples import Sample


class HyperbolicCylinder(Field, Sample):

    def __init__(self, client, name, coefficient=None, radius=None, description=None):
        Field.__init__(self, name=name, field_type='electrostatic')
        Sample.__init__(self, client=client, name=name, description=description)
        self.radius = None
        self.coefficient = None
        self._read_in_radius(radius=radius)
        self._read_in_coefficient(coefficient=coefficient)

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Cylinder radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.set_radius(radius)

    def set_radius(self, radius):
        assert isinstance(radius, numbers.Number), 'Cylinder radius must be a number'
        try:
            self.parameters['Cylinder radius'].float_value = float(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Cylinder radius',
                                                                               value=float(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of the cylinder')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.radius = float(radius)

    def _read_in_coefficient(self, coefficient):
        try:
            self.coefficient = self.parameters['HyperbolicCylinder coefficient'].float_value
        except KeyError:
            pass
        if self.coefficient != coefficient and coefficient is not None:
            self.set_coefficient(coefficient)

    def set_coefficient(self, coefficient):
        assert isinstance(coefficient, numbers.Number), 'HyperbolicCylinder coefficient must be a number'
        try:
            self.parameters['HyperbolicCylinder coefficient'].float_value = float(coefficient)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='HyperbolicCylinder coefficient',
                                                                               value=float(coefficient),
                                                                               description='HyperbolicCylinder coeff.')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.coefficient = float(coefficient)

    def scalar_field(self, xyz):
        rpz = gt.cartesian_to_cylindrical(xyz)
        rpz[np.where(rpz[:, 0] < self.radius), 0] = self.radius
        return self.coefficient / rpz[:, 0]

    def vector_field(self, xyz):
        field = gt.cartesian_to_cylindrical(xyz)
        field[np.where(field[:, 0] < self.radius), 0] = self.radius
        field[:, 0] = abs(self.coefficient / (field[:, 0]) ** 2)
        if self.coefficient < 0:
            field[:, 1] += np.pi
        field[:, 1] = gt.reduce_angle(field[:, 1])
        field[:, 2] = 0
        return gt.cylindrical_to_cartesian(field)
