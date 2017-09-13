from __future__ import division, print_function
import numpy as np

from BDSpace.Coordinates import transforms as gt

from Schottky.Samples.Fields.General import SampleField


class HyperbolicCylinder(SampleField):

    def __init__(self, client, name, coefficient=None, radius=None, description=None, orientation=None):
        SampleField.__init__(self, client=client, name=name, description=description,
                             field_type='electrostatic', orientation=orientation)
        self.load_create_sample()
        self.__radius = None
        self.__coefficient = None
        self._read_in_radius(radius=radius)
        self._read_in_coefficient(coefficient=coefficient)

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, radius):
        try:
            self.parameters['Cylinder radius'].float_value = np.float64(radius)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='Cylinder radius',
                                                                               value=np.float64(radius),
                                                                               unit_name='cm',
                                                                               description='Radius of the cylinder')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__radius = np.float64(radius)

    def _read_in_radius(self, radius):
        try:
            self.radius = self.parameters['Cylinder radius'].float_value
        except KeyError:
            pass
        if self.radius != radius and radius is not None:
            self.radius = radius

    @property
    def coefficient(self):
        return self.__coefficient

    @coefficient.setter
    def coefficient(self, coefficient):
        try:
            self.parameters['HyperbolicCylinder coefficient'].float_value = np.float64(coefficient)
            self.save_sample_changes()
        except KeyError:
            parameter = self.client.parameter_manager.create_numeric_parameter(name='HyperbolicCylinder coefficient',
                                                                               value=np.float64(coefficient),
                                                                               description='HyperbolicCylinder coeff.')
            self.client.sample_manager.add_parameter_to_sample(sample=self.sample,
                                                               parameter=parameter)
            self.reload_parameters()
        self.__coefficient = np.float64(coefficient)

    def _read_in_coefficient(self, coefficient):
        try:
            self.coefficient = self.parameters['HyperbolicCylinder coefficient'].float_value
        except KeyError:
            pass
        if self.coefficient != coefficient and coefficient is not None:
            self.coefficient = coefficient

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
