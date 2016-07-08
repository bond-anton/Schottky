from __future__ import division, print_function
import numpy as np
import numbers

from Space.Coordinates import transforms as gt
from Space.Field import Field

from Schottky import constants


class UniformField(Field):
    def __init__(self, strength, direction, name='Uniform electric field', field_type='electrostatic'):
        assert isinstance(strength, numbers.Number)
        self.strength = np.float(strength)
        self.direction = gt.unit_vector(direction)
        super(UniformField, self).__init__(name, field_type)

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


class ChargedSphere(Field):

    def __init__(self, q, r, epsilon=1, name='Spherical Coulomb trap', field_type='electrostatic'):
        self.r = r
        self.q = q / (4 * np.pi * constants['epsilon_0'] * epsilon)
        super(ChargedSphere, self).__init__(name, field_type)

    def scalar_field(self, xyz):
        rtp = gt.cartesian_to_spherical(xyz)
        rtp[np.where(rtp[:, 0] < self.r), 0] = self.r
        return self.q / rtp[:, 0]

    def vector_field(self, xyz):
        rtp = gt.cartesian_to_spherical(xyz)
        rtp[np.where(rtp[:, 0] < self.r), 0] = self.r
        r = rtp[:, 0] ** 2
        r = np.array([r, r, r]).T
        return self.q * xyz / r


class ChargedCylinder(Field):

    def __init__(self, l, r, epsilon=1, name='Cylindrical Coulomb trap', field_type='electrostatic'):
        self.r = r  # radius of cylinder
        self.b = constants['q'] * l / (2 * np.pi * constants['epsilon_0'] * epsilon)
        super(ChargedCylinder, self).__init__(name, field_type)

    def scalar_field(self, xyz):
        rpz = gt.cartesian_to_cylindrical(xyz)
        rpz[np.where(rpz[:, 0] < self.r), 0] = self.r
        return self.b * np.log(rpz[:, 0] / self.r)

    def vector_field(self, xyz):
        field = gt.cartesian_to_cylindrical(xyz)
        field[np.where(field[:, 0] < self.r), 0] = self.r
        field[:, 0] = self.b / field[:, 0]
        field[:, 2] = 0
        return gt.cylindrical_to_cartesian(field)


class ReciprocalCylinder(Field):

    def __init__(self, a, r, name='1/r Cylindrical potential', field_type='electrostatic'):
        self.r = r  # radius of cylinder
        self.a = constants['q'] * a
        super(ReciprocalCylinder, self).__init__(name, field_type)

    def scalar_field(self, xyz):
        rpz = gt.cartesian_to_cylindrical(xyz)
        rpz[np.where(rpz[:, 0] < self.r), 0] = self.r
        return -self.a / rpz[:, 0]

    def vector_field(self, xyz):
        field = gt.cartesian_to_cylindrical(xyz)
        field[np.where(field[:, 0] < self.r), 0] = self.r
        field[:, 0] = self.a / (field[:, 0]) ** 2
        field[:, 2] = 0
        return -1*gt.cylindrical_to_cartesian(field)
