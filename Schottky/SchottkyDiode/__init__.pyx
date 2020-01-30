from libc.math cimport M_PI, sqrt, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from Schottky.Metal cimport Metal
from Schottky.Semiconductor cimport Semiconductor
from Schottky.Constants cimport constant

cdef class SchottkyDiode(object):

    def __init__(self, str label, Metal metal, Semiconductor semiconductor,
                 double length=1.0e-5,
                 double area=1.0, double serial_resistance=0.0):
        self.__label = label
        self.__metal = metal
        self.__semiconductor = semiconductor
        self.__length = fabs(length)
        self.__area = fabs(area)
        self.__serial_resistance = fabs(serial_resistance)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def metal(self):
        return self.__metal

    @property
    def semiconductor(self):
        return  self.__semiconductor

    @property
    def length(self):
        return self.__length

    @length.setter
    def length(self, double length):
        self.__length = fabs(length)

    @property
    def area(self):
        return self.__area

    @area.setter
    def area(self, double area):
        self.__area = fabs(area)

    @property
    def contact_diameter(self):
        return sqrt(self.__area / M_PI) * 2

    @contact_diameter.setter
    def contact_diameter(self, double diameter):
        self.__area = M_PI * diameter * diameter / 4

    @property
    def contact_radius(self):
        return sqrt(self.__area / M_PI)

    @contact_radius.setter
    def contact_radius(self, double radius):
        self.__area = M_PI * radius * radius

    @property
    def serial_resistance(self):
        return self.__serial_resistance

    @serial_resistance.setter
    def serial_resistance(self, double serial_resistance):
        self.__serial_resistance = fabs(serial_resistance)

    cpdef double phi_b_n_t(self, double temperature):
        return self.__metal.__work_function - self.__semiconductor.__reference['affinity']

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_n(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_n_t(temperature[i])
        return result

    cpdef double phi_b_n_ev_t(self, double temperature):
        return constant.joule_to_ev_point(self.__metal.__work_function
                                          - self.__semiconductor.__reference['affinity'])

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_n_ev(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_n_ev_t(temperature[i])
        return result

    cpdef double phi_b_n_boltzmann_t(self, double temperature):
        return constant.joule_to_boltzmann_point(self.__metal.__work_function
                                                 - self.__semiconductor.__reference['affinity'],
                                                 temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_n_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_n_boltzmann_t(temperature[i])
        return result

    cpdef double phi_b_p_t(self, double temperature):
        return self.__semiconductor.__reference['affinity'] \
               + self.__semiconductor.band_gap_t(temperature) - self.__metal.__work_function

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_p(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_p_t(temperature[i])
        return result

    cpdef double phi_b_p_ev_t(self, double temperature):
        return constant.joule_to_ev_point(self.__semiconductor.__reference['affinity'] \
               + self.__semiconductor.band_gap_t(temperature) - self.__metal.__work_function)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_p_ev(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_p_ev_t(temperature[i])
        return result

    cpdef double phi_b_p_boltzmann_t(self, double temperature):
        return constant.joule_to_boltzmann_point(self.__semiconductor.__reference['affinity'] \
                                                 + self.__semiconductor.band_gap_t(temperature) \
                                                 - self.__metal.__work_function,
                                                 temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] phi_b_p_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.phi_b_p_boltzmann_t(temperature[i])
        return result

    cpdef double built_in_voltage_t(self, double temperature):
        return self.__metal.__work_function - self.__semiconductor.work_function_t(temperature,
                                                                                   f_threshold=1.0e-23,
                                                                                   max_iter=100, verbose=False)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] built_in_voltage(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.built_in_voltage_t(temperature[i])
        return result

    cpdef double built_in_voltage_ev_t(self, double temperature):
        return constant.joule_to_ev_point(
            self.__metal.__work_function - self.__semiconductor.work_function_t(temperature,
                                                                                f_threshold=1.0e-23,
                                                                                max_iter=100, verbose=False))

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] built_in_voltage_ev(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.built_in_voltage_ev_t(temperature[i])
        return result

    cpdef double built_in_voltage_boltzmann_t(self, double temperature):
        return constant.joule_to_boltzmann_point(
            self.__metal.__work_function - self.__semiconductor.work_function_t(temperature,
                                                                                f_threshold=1.0e-23,
                                                                                max_iter=100, verbose=False),
            temperature
        )

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] built_in_voltage_boltzmann(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.built_in_voltage_boltzmann_t(temperature[i])
        return result

    cpdef double n0_t(self, double temperature):
        return self.__semiconductor.n_e_t(self.phi_b_n_t(temperature), temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n0(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.n0_t(temperature[i])
        return result

    cpdef double p0_t(self, double temperature):
        return self.__semiconductor.n_h_t(self.phi_b_n_t(temperature), temperature)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] p0(self, double[:] temperature):
        cdef:
            int n = temperature.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.p0_t(temperature[i])
        return result
