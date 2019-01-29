from __future__ import division, print_function
import numpy as np

from BDSpace.Field.Field cimport Field


cdef class ElectricField(Field):

    def __init__(self, str name):
        field_type = 'Electric Field'
        super(ElectricField, self).__init__(name, field_type)


cdef class ConstantPotentialField(ElectricField):

    def __init__(self, str name, double potential):
        self.__potential = potential
        super(ConstantPotentialField, self).__init__(name)

    def scalar_field(self, xyz):
        if xyz.size == 3:
            return self.__potential
        elif xyz.size > 3:
            if len(xyz.shape) == 2 and xyz.shape[1] == 3:
                return np.ones(xyz.shape[0]) * self.__potential
            else:
                raise ValueError('N-points array expected')
        else:
            raise ValueError('N-points array expected')
