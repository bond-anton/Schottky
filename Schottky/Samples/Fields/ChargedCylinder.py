from __future__ import division, print_function
import numpy as np

from BDSpace.Coordinates import transforms as gt

from Schottky import constants
from Schottky.Samples.Fields.General import SampleField


class ChargedCylinder(SampleField):

    def __init__(self, client, name, charge_density=None, radius=None, epsilon=None,
                 description=None, orientation=None):
        SampleField.__init__(self, client=client, name=name, description=description,
                             field_type='electrostatic', orientation=orientation)
        self.__radius = None
        self.__charge_density = None
        self.__epsilon = None
        self.client.session.commit()
        self._read_in_radius(radius=radius)
        self._read_in_charge_density(charge_density=charge_density)
        self._read_in_epsilon(epsilon=epsilon)

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, radius):
        try:
            self.parameters['Charged cylinder radius'].float_value = np.float64(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged cylinder radius',
                                                                               value=np.float64(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of charged cylinder')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__radius = np.float64(radius)

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Charged cylinder radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.radius = radius

    @property
    def charge_density(self):
        return self.__charge_density

    @charge_density.setter
    def charge_density(self, charge_density):
        try:
            self.parameters['Charged cylinder charge density'].float_value = np.float64(charge_density)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Charged cylinder charge density',
                                                                               value=np.float64(charge_density),
                                                                               unit_name='1/cm',
                                                                               description='Linear density of traps')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__charge_density = np.float64(charge_density)

    def _read_in_charge_density(self, charge_density):
        try:
            self.charge_density = self.parameters['Charged cylinder charge density'].float_value
        except KeyError:
            pass
        if self.charge_density != charge_density and charge_density is not None:
            self.charge_density = charge_density

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
