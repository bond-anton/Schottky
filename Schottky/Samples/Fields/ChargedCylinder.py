from __future__ import division, print_function
import time
import numpy as np
import numbers

from Space.Coordinates import transforms as gt
from Space.Field import Field

from Schottky import constants
from Schottky.Samples import Sample


class ChargedCylinder(Sample, Field):

    def __init__(self, client, name, charge_density=None, radius=None, epsilon=None, description=None):
        Sample.__init__(self, client=client, name=name, description=description)
        self.load_create_sample()
        self.radius = None
        self.charge_density = None
        self.epsilon = None
        self.client.session.commit()
        self._read_in_radius(radius=radius)
        self._read_in_charge_density(charge_density=charge_density)
        self._read_in_epsilon(epsilon=epsilon)
        Field.__init__(self, name=name, field_type='electrostatic')

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Charged cylinder radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.set_radius(radius)

    def set_radius(self, radius):
        assert isinstance(radius, numbers.Number), 'Charged cylinder radius must be a number'
        try:
            self.parameters['Charged cylinder radius'].float_value = float(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged cylinder radius',
                                                                               value=float(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of charged cylinder')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.radius = float(radius)

    def _read_in_charge_density(self, charge_density):
        try:
            self.charge_density = self.parameters['Charged cylinder charge density'].float_value
        except KeyError:
            pass
        if self.charge_density != charge_density and charge_density is not None:
            self.set_charge_density(charge_density)

    def set_charge_density(self, charge_density):
        assert isinstance(charge_density, numbers.Number), 'Charged cylinder charge density must be a number'
        try:
            self.parameters['Charged cylinder charge density'].float_value = float(charge_density)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged cylinder charge density',
                                                                               value=float(charge_density),
                                                                               unit_name='1/cm',
                                                                               description='Linear density of traps')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.charge_density = float(charge_density)

    def _read_in_epsilon(self, epsilon):
        try:
            self.epsilon = self.parameters['epsilon'].float_value
        except KeyError:
            pass
        if self.epsilon != epsilon and epsilon is not None:
            self.set_epsilon(epsilon)

    def set_epsilon(self, epsilon):
        assert isinstance(epsilon, numbers.Number), 'epsilon must be a number'
        try:
            self.parameters['epsilon'].float_value = float(epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='epsilon',
                                                                               value=float(epsilon),
                                                                               description='Permittivity of space')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.epsilon = float(epsilon)

    def scalar_field(self, xyz):
        coefficient = constants['q'] * self.charge_density / (2 * np.pi * constants['epsilon_0'] * self.epsilon)
        rpz = gt.cartesian_to_cylindrical(xyz)
        rpz[np.where(rpz[:, 0] < self.radius), 0] = self.radius
        return -coefficient * np.log(rpz[:, 0] / self.radius)

    def vector_field(self, xyz):
        coefficient = constants['q'] * self.charge_density / (2 * np.pi * constants['epsilon_0'] * self.epsilon)
        field = gt.cartesian_to_cylindrical(xyz)
        field[np.where(field[:, 0] < self.radius), 0] = self.radius
        field[:, 0] = abs(coefficient / field[:, 0])
        if coefficient < 0:
            field[:, 1] += np.pi
        field[:, 1] = gt.reduce_angle(field[:, 1])
        field[:, 2] = 0
        return gt.cylindrical_to_cartesian(field)
