from libc.math cimport M_PI, sqrt, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDFunction1D cimport Function, Functional
from BDFunction1D.Standard cimport LineThroughPoints
from BDFunction1D.Differentiation import NumericGradient

from Schottky.Metal cimport Metal
from Schottky.Semiconductor cimport Semiconductor
from Schottky.Constants cimport constant

from Schottky.Generation cimport ZeroGenerationFunction
from Schottky.Recombination cimport ZeroRecombinationFunction, ConstantRecombinationFunction


cdef class SchottkyDiode(object):

    def __init__(self, str label, Metal metal, Semiconductor semiconductor,
                 double length=1.0e-5,
                 double area=1.0, double serial_resistance=0.0,
                 double temperature=300, double bias=0.0,
                 Function electric_potential=None, Function quasi_fermi_e=None, Function quasi_fermi_h=None,
                 Function generation=None, Function recombination=None,):
        self.__label = label
        self.__metal = metal
        self.__semiconductor = semiconductor
        self.__length = fabs(length)
        self.__area = fabs(area)
        self.__serial_resistance = fabs(serial_resistance)
        self.__temperature = temperature
        self.__bias = bias

        if electric_potential is None:
            self.__ep = LineThroughPoints(0.0,
                                          self.built_in_voltage_ev_t(self.__temperature) - self.__bias,
                                          self.__length,
                                          0.0)
        else:
            self.__ep = electric_potential
        self.__ef = -NumericGradient(self.__ep)
        if quasi_fermi_e is None:
            self.__qf_e = LineThroughPoints(0.0,
                                            self.phi_b_n_ev_t(self.__temperature),
                                            self.__length,
                                            self.__semiconductor.el_chem_pot_ev_t(self.__temperature))
        else:
            self.__qf_e = quasi_fermi_e
        if quasi_fermi_h is None:
            self.__qf_h = LineThroughPoints(0.0,
                                            self.phi_b_n_ev_t(self.__temperature),
                                            self.__length,
                                            self.__semiconductor.el_chem_pot_ev_t(self.__temperature))
        else:
            self.__qf_h = quasi_fermi_h
        self.__ne = QFeNe(self, self.__qf_e)
        self.__nh = QFhNh(self, self.__qf_h)
        self.__pn = self.__ne * self.__nh
        self.__mu_e = MobilityE(self)
        self.__mu_h = MobilityH(self)
        if generation is None:
            self.__generation = ZeroGenerationFunction()
        else:
            self.__generation = generation
        if recombination is None:
            self.__recombination = ZeroRecombinationFunction()
        else:
            self.__recombination = recombination

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

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, double temperature):
        self.__temperature = temperature

    @property
    def bias(self):
        return self.__bias

    @bias.setter
    def bias(self, double bias):
        self.__bias = bias
        self.__ep += LineThroughPoints(0.0, -self.__bias, self.__length, 0.0)

    @property
    def v_bi(self):
        return self.built_in_voltage_ev_t(self.__temperature)

    @property
    def xi(self):
        return self.__semiconductor.el_chem_pot_ev_t(self.__temperature)

    @property
    def electric_potential(self):
        return self.__ep

    @property
    def quasi_fermi_e(self):
        return self.__qf_e

    @quasi_fermi_e.setter
    def quasi_fermi_e(self, Function quasi_fermi_e):
        self.__qf_e = quasi_fermi_e
        if isinstance(self.__ne, QFeNe):
            self.__ne.f = self.__qf_e
        else:
            self.__ne = QFeNe(self.__qf_e)
        self.__pn = self.__ne * self.__nh

    @property
    def quasi_fermi_h(self):
        return self.__qf_h

    @quasi_fermi_h.setter
    def quasi_fermi_h(self, Function quasi_fermi_h):
        self.__qf_h = quasi_fermi_h
        if isinstance(self.__nh, QFhNh):
            self.__nh.f = self.__qf_h
        else:
            self.__nh = QFhNh(self.__qf_h)
        self.__pn = self.__ne * self.__nh

    @property
    def n_e(self):
        return self.__ne

    @n_e.setter
    def n_e(self, Function n_e):
        self.__ne = n_e
        if isinstance(self.__qf_e, NeQFe):
            self.__qf_e.f = self.__ne
        else:
            self.__qf_e = NeQFe(self.__ne)
        self.__pn = self.__ne * self.__nh

    @property
    def n_h(self):
        return self.__nh

    @n_h.setter
    def n_h(self, Function n_h):
        self.__nh = n_h
        if isinstance(self.__qf_h, NhQFh):
            self.__qf_h.f = self.__nh
        else:
            self.__qf_h = NeQFe(self.__nh)
        self.__pn = self.__ne * self.__nh

    @property
    def generation(self):
        return self.__generation

    @generation.setter
    def generation(self, Function generation):
        self.__generation = generation

    @property
    def recombination(self):
        return self.__recombination

    @recombination.setter
    def recombination(self, Function recombination):
        self.__recombination = recombination

    @property
    def pn(self):
        return self.__pn

    @property
    def mu_e(self):
        return self.__mu_e

    @property
    def mu_h(self):
        return self.__mu_h

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

    @boundscheck(False)
    @wraparound(False)
    cpdef double thermionic_emission_current_e(self):
        cdef:
            double n0, nb, a, nc, t2
        n0 = self.__ne.evaluate_point(0.0)
        nb = self.n0_t(self.__temperature)
        nc = self.__semiconductor.n_c_t(self.__temperature)
        a = self.__semiconductor.__reference['thermionic_emission']['A_R_coeff_n'] * constant.__A_R
        t2 = self.__temperature * self.__temperature
        return a * t2 * (n0 - nb) / nc

    @boundscheck(False)
    @wraparound(False)
    cpdef double thermionic_emission_current_h(self):
        cdef:
            double p0, pb, a, nv, t2
        p0 = self.__nh.evaluate_point(0.0)
        pb = self.p0_t(self.__temperature)
        nv = self.__semiconductor.n_v_t(self.__temperature)
        a = self.__semiconductor.__reference['thermionic_emission']['A_R_coeff_p'] * constant.__A_R
        t2 = self.__temperature * self.__temperature
        return a * t2 * (pb - p0) / nv



cdef class QFeNe(Functional):

    def __init__(self, SchottkyDiode diode, Function f):
        self.__diode = diode
        super(QFeNe, self).__init__(f)

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.n_e_ev_t(self.__f.evaluate_point(x), self.__diode.__temperature)


cdef class QFhNh(Functional):

    def __init__(self, SchottkyDiode diode, Function f):
        self.__diode = diode
        super(QFhNh, self).__init__(f)

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.n_h_ev_t(self.__f.evaluate_point(x), self.__diode.__temperature)


cdef class NeQFe(Functional):

    def __init__(self, SchottkyDiode diode, Function f):
        self.__diode = diode
        super(NeQFe, self).__init__(f)

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.n_e_to_mu_ev_t(self.__f.evaluate_point(x), self.__diode.__temperature)


cdef class NhQFh(Functional):

    def __init__(self, SchottkyDiode diode, Function f):
        self.__diode = diode
        super(NhQFh, self).__init__(f)

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.n_h_to_mu_ev_t(self.__f.evaluate_point(x), self.__diode.__temperature)


cdef class MobilityE(Function):

    def __init__(self, SchottkyDiode diode):
        self.__diode = diode
        super(MobilityE, self).__init__()

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.mobility_e_point_t(
            self.__diode.__semiconductor.__dopants_n.evaluate_point(x),
            self.__diode.__ef.evaluate_point(x),
            self.__diode.__pn.evaluate_point(x),
            self.__diode.__temperature)


cdef class MobilityH(Function):

    def __init__(self, SchottkyDiode diode):
        self.__diode = diode
        super(MobilityH, self).__init__()

    cpdef double evaluate_point(self, double x):
        return self.__diode.__semiconductor.mobility_h_point_t(
            self.__diode.__semiconductor.__dopants_n.evaluate_point(x),
            self.__diode.__ef.evaluate_point(x),
            self.__diode.__pn.evaluate_point(x),
            self.__diode.__temperature)
