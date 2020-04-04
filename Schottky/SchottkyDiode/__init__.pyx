from libc.math cimport M_PI, sqrt, fabs
from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDMesh.Mesh1DUniform cimport Mesh1DUniform

from BDFunction1D cimport Function, Functional
from BDFunction1D.Functional cimport ScaledFunction
from BDFunction1D.Standard cimport Zero, Constant, LineThroughPoints
from BDFunction1D.Differentiation cimport NumericGradient

from BDPoisson1D.FirstOrderLinear cimport cauchy_first_order_solver_mesh
from BDPoisson1D.FirstOrderNonLinear cimport dirichlet_non_linear_first_order_solver_mesh
from BDPoisson1D.FirstOrderNonLinear cimport dirichlet_non_linear_first_order_solver_recurrent_mesh
from BDPoisson1D.DirichletLinear cimport dirichlet_poisson_solver_mesh, dirichlet_poisson_solver_amr

from Schottky.Metal cimport Metal
from Schottky.Semiconductor cimport Semiconductor
from Schottky.Dopant cimport Dopant
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
        self.__grad_qf_e = NumericGradient(self.__qf_e)
        if quasi_fermi_h is None:
            self.__qf_h = LineThroughPoints(0.0,
                                            self.phi_b_n_ev_t(self.__temperature),
                                            self.__length,
                                            self.__semiconductor.el_chem_pot_ev_t(self.__temperature))
        else:
            self.__qf_h = quasi_fermi_h
        self.__grad_qf_h = NumericGradient(self.__qf_h)
        self.__ne = QFeNe(self, self.__qf_e)
        self.__nh = QFhNh(self, self.__qf_h)
        self.__pn = self.__ne * self.__nh
        self.__mu_e = MobilityE(self)
        self.__mu_h = MobilityH(self)
        self.__dopants_charge = DopantsEqCharge(self)
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
        ep_0 = self.__ep.evaluate_point(0)
        if ep_0 != self.v_bi - self.__bias:
            delta = self.v_bi - self.__bias - ep_0
            self.__ep += LineThroughPoints(0.0, delta, self.__length, 0.0)
            self.__ef = -NumericGradient(self.__ep)

    @property
    def v_bi(self):
        return self.built_in_voltage_ev_t(self.__temperature)

    @property
    def xi(self):
        return self.__semiconductor.el_chem_pot_ev_t(self.__temperature)

    @property
    def dopants_charge(self):
        return self.__dopants_charge

    @property
    def electric_potential(self):
        return self.__ep

    @property
    def electric_field(self):
        return self.__ef

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

    cpdef stationary_grad_qf_e_solver(self):
        cdef:
            double qkt, qee, bc1, bc2
            Mesh1DUniform mesh
            Function p, sol
            Functional f, df_dy
        qkt = constant.__q / (constant.__k * self.__temperature)
        bc1 = 0*self.thermionic_emission_current_e() / (constant.__q *
                                                       self.__ne.evaluate_point(0) * self.__mu_e.evaluate_point(0))
        bc2 = -self.__ef.evaluate_point(self.length)*0
        print('BC:', bc1, bc2)
        p = -qkt * self.__ef
        f = GradQFeF(self.__grad_qf_e, self.__generation, self.__recombination,
                     self.__mu_e, self.__ne, self.__temperature)
        df_dy = ScaledFunction(self.__grad_qf_e, 2 * qkt)
        mesh = Mesh1DUniform(0.0, self.__length, bc1, bc2, 1.0e-8)
        # self.__grad_qf_e = dirichlet_non_linear_first_order_solver_mesh(mesh, self.__grad_qf_e, p, f, df_dy, w=0.7)
        self.__grad_qf_e = dirichlet_non_linear_first_order_solver_recurrent_mesh(mesh, self.__grad_qf_e, p, f, df_dy, w=0.7)
        mesh.boundary_condition_2 = self.xi
        self.quasi_fermi_e = cauchy_first_order_solver_mesh(mesh, self.__grad_qf_e, ic=0.0, idx=-1)

    cpdef stationary_grad_qf_h_solver(self):
        cdef:
            double qkt, bc1, bc2
            Mesh1DUniform mesh
            Function p, sol
            Functional f, df_dy
        qkt = constant.__q / (constant.__k * self.__temperature)
        bc1 = -self.thermionic_emission_current_h() / (constant.__q *
                                                       self.__nh.evaluate_point(0) * self.__mu_h.evaluate_point(0))
        bc2 = -self.__ef.evaluate_point(self.length)
        print('BC:', bc1, bc2)
        p = qkt * self.__ef
        f = GradQFhF(self.__grad_qf_h, self.__generation, self.__recombination,
                     self.__mu_h, self.__nh, self.__temperature)
        df_dy = ScaledFunction(self.__grad_qf_h, -2 * qkt)
        mesh = Mesh1DUniform(0.0, self.__length, bc1, bc2, 1.0e-8)
        # self.__grad_qf_h = dirichlet_non_linear_first_order_solver_mesh(mesh, self.__grad_qf_h, p, f, df_dy, w=0.7)
        mesh.__boundary_condition_2 = self.xi
        # self.quasi_fermi_h = cauchy_first_order_solver_mesh(mesh, self.__grad_qf_h, ic=0.0, idx=-1)
        self.__grad_qf_h = self.__grad_qf_e
        self.quasi_fermi_h = self.quasi_fermi_e

    cpdef Function poisson_eq_solver(self):
        cdef:
            double qee, bc1, bc2
        qee = constant.__q / (constant.__epsilon_0 * self.__semiconductor.__reference['epsilon'])
        # rho = qee * (self.__dopants_charge - self.__ne + self.__nh)
        rho = qee * (self.__dopants_charge - self.__ne)
        # rho = qee * TestSCR()
        # rho1 = TestSCR()
        # rho1 = self.__dopants_charge - self.__ne + self.__nh
        rho1 = self.__dopants_charge - self.__ne
        bc1 = self.v_bi - self.bias
        bc2 = 0.0
        mesh = Mesh1DUniform(0.0, self.__length, bc1, bc2, 1.0e-8)
        self.__ep = dirichlet_poisson_solver_mesh(mesh, rho)
        self.__ef = -NumericGradient(self.__ep)
        return rho1


cdef class TestSCR(Function):

    cpdef double evaluate_point(self, double x):
        if x < 1e-6:
            return 1.0e21
        return 0.0


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


cdef class DopantsEqCharge(Function):

    def __init__(self, SchottkyDiode diode):
        self.__diode = diode
        super(DopantsEqCharge, self).__init__()

    cpdef double evaluate_point(self, double x):
        cdef:
            double mu_e, mu_h, f_e, f_h, result
            Dopant dopant
        mu_e = self.__diode.__qf_e.evaluate_point(x)
        mu_h = self.__diode.__qf_h.evaluate_point(x)
        result = 0.0
        for dopant in self.__diode.__semiconductor.__dopants:
            f_e = self.__diode.__semiconductor.trap_eq_occupation(
                dopant, mu_e, self.__diode.__temperature,
                f_threshold=1.0e-23, max_iter=100, verbose=False)
            # f_h = self.__diode.__semiconductor.trap_eq_occupation(
            #     dopant, mu_h, self.__diode.__temperature,
            #     f_threshold=1.0e-23, max_iter=100, verbose=False)
            result += dopant.__concentration.evaluate_point(x) * (
                    (dopant.charge_state[1] - dopant.charge_state[0]) * f_e + dopant.charge_state[0])
        return result


cdef class GradQFeF(Functional):

    def __init__(self, Function y, Function g, Function r, Function mu, Function n, double temperature):
        super(GradQFeF, self).__init__(y)
        self.__g = g
        self.__r = r
        self.__mu = mu
        self.__n = n
        self.__temperature = temperature
        self.__qkt = constant.__q / (constant.__k * self.__temperature)

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, double temperature):
        self.__temperature = temperature
        self.__qkt = constant.__q / (constant.__k * self.__temperature)

    cpdef double evaluate_point(self, double x):
        cdef:
            double gr = 0.0
        if not isinstance(self.__g, Zero):
            if isinstance(self.__g, Constant):
                gr += self.__g.c
            else:
                gr += self.__g.evaluate_point(x)
        if not isinstance(self.__r, Zero):
            if isinstance(self.__r, Constant):
                gr -= self.__r.c
            else:
                gr -= self.__r.evaluate_point(x)
        if gr != 0.0:
            gr /= self.__mu.evaluate_point(x) * self.__n.evaluate_point(x)
        return self.__qkt * self.__f.evaluate_point(x) * self.__f.evaluate_point(x) + gr


cdef class GradQFhF(Functional):

    def __init__(self, Function y, Function g, Function r, Function mu, Function n, double temperature):
        super(GradQFhF, self).__init__(y)
        self.__g = g
        self.__r = r
        self.__mu = mu
        self.__n = n
        self.__temperature = temperature
        self.__qkt = constant.__q / (constant.__k * self.__temperature)

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, double temperature):
        self.__temperature = temperature
        self.__qkt = constant.__q / (constant.__k * self.__temperature)

    cpdef double evaluate_point(self, double x):
        cdef:
            double gr = 0.0
        if not isinstance(self.__g, Zero):
            if isinstance(self.__g, Constant):
                gr -= self.__g.c
            else:
                gr -= self.__g.evaluate_point(x)
        if not isinstance(self.__r, Zero):
            if isinstance(self.__r, Constant):
                gr += self.__r.c
            else:
                gr += self.__r.evaluate_point(x)
        if gr != 0.0:
            gr /= self.__mu.evaluate_point(x) * self.__n.evaluate_point(x)
        return -self.__qkt * self.__f.evaluate_point(x) * self.__f.evaluate_point(x) + gr