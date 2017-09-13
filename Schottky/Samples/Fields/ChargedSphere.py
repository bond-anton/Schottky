from __future__ import division, print_function
import numpy as np

from BDSpace.Coordinates import transforms as gt

from Schottky import constants
from Schottky.Samples.Fields.General import SampleField


class ChargedSphere(SampleField):

    def __init__(self, client, name, charge=None, radius=None, epsilon=None, description=None, orientation=None):
        SampleField.__init__(self, client=client, name=name, description=description,
                             field_type='electrostatic', orientation=orientation)
        self.__radius = None
        self.__charge = None
        self.__epsilon = None
        self._read_in_radius(radius=radius)
        self._read_in_charge(charge=charge)
        self._read_in_epsilon(epsilon=epsilon)

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, radius):
        try:
            self.parameters['Charged sphere radius'].float_value = np.float64(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged sphere radius',
                                                                               value=np.float64(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of charged sphere')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__radius = np.float64(radius)

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Charged sphere radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.radius = radius

    @property
    def charge(self):
        return self.__charge

    @charge.setter
    def charge(self, charge):
        try:
            self.parameters['Charged sphere charge'].float_value = np.float64(charge)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged sphere charge',
                                                                               value=np.float64(charge),
                                                                               unit_name='C',
                                                                               description='Charge of charged sphere')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__charge = np.float64(charge)

    def _read_in_charge(self, charge):
        try:
            self.charge = self.parameters['Charged sphere charge'].float_value
        except KeyError:
            pass
        if self.charge != charge and charge is not None:
            self.charge = charge

    @property
    def epsilon(self):
        return self.__epsilon

    @epsilon.setter
    def epsilon(self, epsilon):
        try:
            self.parameters['epsilon'].float_value = np.float64(epsilon)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='epsilon',
                                                                               value=np.float64(epsilon),
                                                                               description='Permittivity of space')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__epsilon = np.float64(epsilon)

    def _read_in_epsilon(self, epsilon):
        try:
            self.epsilon = self.parameters['epsilon'].float_value
        except KeyError:
            pass
        if self.epsilon != epsilon and epsilon is not None:
            self.epsilon = epsilon

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
