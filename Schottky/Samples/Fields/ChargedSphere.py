from __future__ import division, print_function
import numpy as np
import numbers

from Space.Coordinates import transforms as gt
from Space.Field import Field

from Schottky import constants
from Schottky.Samples import Sample


class ChargedSphere(Sample, Field):

    def __init__(self, client, name, charge=None, radius=None, epsilon=None, description=None):
        Sample.__init__(self, client=client, name=name, description=description)
        self.load_create_sample()
        self.radius = None
        self.charge = None
        self.epsilon = None
        self._read_in_radius(radius=radius)
        self._read_in_charge(charge=charge)
        self._read_in_epsilon(epsilon=epsilon)
        Field.__init__(self, name=name, field_type='electrostatic')

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Charged sphere radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.set_radius(radius)

    def set_radius(self, radius):
        assert isinstance(radius, numbers.Number), 'Charged sphere radius must be a number'
        try:
            self.parameters['Charged sphere radius'].float_value = float(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged sphere radius',
                                                                               value=float(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of charged sphere')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.radius = float(radius)

    def _read_in_charge(self, charge):
        try:
            self.charge = self.parameters['Charged sphere charge'].float_value
        except KeyError:
            pass
        if self.charge != charge and charge is not None:
            self.set_charge(charge)

    def set_charge(self, charge):
        assert isinstance(charge, numbers.Number), 'Charged sphere charge must be a number'
        try:
            self.parameters['Charged sphere charge'].float_value = float(charge)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged sphere charge',
                                                                               value=float(charge),
                                                                               unit_name='C',
                                                                               description='Charge of charged sphere')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.load_create_sample()
        self.charge = float(charge)

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
        coefficient = self.charge / (4 * np.pi * constants['epsilon_0'] * self.epsilon)
        rtp = gt.cartesian_to_spherical(xyz)
        rtp[np.where(rtp[:, 0] < self.radius), 0] = self.radius
        return coefficient / rtp[:, 0]

    def vector_field(self, xyz):
        coefficient = self.charge / (4 * np.pi * constants['epsilon_0'] * self.epsilon)
        rtp = gt.cartesian_to_spherical(xyz)
        rtp[np.where(rtp[:, 0] < self.radius), 0] = self.radius
        r = rtp[:, 0] ** 2
        r = np.array([r, r, r]).T
        return coefficient * xyz / r
