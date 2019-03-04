from BDSpace.Coordinates.transforms cimport unit_vector, vector_norm
from BDSpace.Field.Field cimport ConstantVectorConservativeField, HyperbolicPotentialSphericalConservativeField


cdef class ExternalField(ConstantVectorConservativeField):

    def __init__(self, str name, double[:] direction, double magnitude):
        cdef:
            double[:] potential = unit_vector(direction)
            str field_type = 'Electric Field'
        self.__direction = unit_vector(direction)
        potential[0] = self.__direction[0] * magnitude
        potential[1] = self.__direction[1] * magnitude
        potential[2] = self.__direction[2] * magnitude
        super(ExternalField, self).__init__(name, field_type, potential)

    @property
    def magnitude(self):
        return vector_norm(self.__potential)

    @magnitude.setter
    def magnitude(self, double magnitude):
        self.__potential[0] = self.__direction[0] * magnitude
        self.__potential[1] = self.__direction[1] * magnitude
        self.__potential[2] = self.__direction[2] * magnitude

    @property
    def direction(self):
        return self.__direction

    @direction.setter
    def direction(self, double[:] direction):
        cdef:
            double magnitude = vector_norm(self.__potential)
        self.__direction = unit_vector(direction)
        self.__potential[0] = self.__direction[0] * magnitude
        self.__potential[1] = self.__direction[1] * magnitude
        self.__potential[2] = self.__direction[2] * magnitude
