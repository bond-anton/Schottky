import numpy as np

from cython cimport boundscheck, wraparound
from cpython.array cimport array, clone

from BDMesh.Mesh1D cimport Mesh1D
from BDMesh.Mesh1DUniform cimport Mesh1DUniform
from BDMesh.TreeMesh1DUniform cimport TreeMesh1DUniform

from Schottky.SchottkyDiode cimport SchottkyDiode
from Schottky.Dopant cimport Dopant
from Schottky.Constants cimport constant


cdef class DCMeasurement(object):

    def __init__(self, str label, SchottkyDiode diode,
                 double temperature=300, double bias=300, double initial_step=5e-7):
        self.__label = label
        self.__diode = diode
        self.__temperature = temperature
        self.__bias = bias
        self.__initial_step = initial_step
        self.__ep = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                      aligned=True)
        self.__qf_e = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                        aligned=True)
        self.__qf_h = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                        aligned=True)
        self.__n_e = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                       aligned=True)
        self.__n_h = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                       aligned=True)
        self.__generation = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                              aligned=True)
        self.__recombination = TreeMesh1DUniform(Mesh1DUniform(0.0, self.__diode.__length, physical_step=initial_step),
                                                 aligned=True)
        self.prepare_psi0()

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, double temperature):
        self.__temperature = temperature
        self.prepare_psi0()

    @property
    def bias(self):
        return self.__bias

    @bias.setter
    def bias(self, double bias):
        self.__bias = bias
        self.prepare_psi0()

    @property
    def v_bi(self):
        return self.__diode.built_in_voltage_ev_t(self.__temperature)

    @property
    def mu(self):
        return self.__diode.semiconductor.el_chem_pot_ev_t(self.__temperature)

    @property
    def diode(self):
        return self.__diode

    @property
    def electric_potential(self):
        return self.__ep

    @property
    def quasi_fermi_e(self):
        return self.__qf_e

    @property
    def quasi_fermi_h(self):
        return self.__qf_h

    @property
    def n_e(self):
        return self.__n_e

    @property
    def n_h(self):
        return self.__n_h

    @property
    def generation(self):
        return self.__generation

    @property
    def recombination(self):
        return self.__recombination

    cpdef prepare_psi0(self):
        v_bi = self.__diode.built_in_voltage_ev_t(self.__temperature)
        v_bc = v_bi - self.__bias
        phi_b_n = self.diode.phi_b_n_ev_t(self.__temperature)
        xi = self.__diode.__semiconductor.el_chem_pot_ev_t(self.__temperature)
        phi_b_p = self.diode.phi_b_p_ev_t(self.__temperature)
        psi_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        psi_0.solution = v_bc - np.asarray(psi_0.physical_nodes) * v_bc / psi_0.physical_nodes[-1]
        self.__ep = TreeMesh1DUniform(psi_0, aligned=True)
        qf_e_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        qf_e_0.solution = phi_b_n + np.asarray(qf_e_0.physical_nodes) * (xi - phi_b_n) / qf_e_0.physical_nodes[-1]
        self.__qf_e = TreeMesh1DUniform(qf_e_0, aligned=True)
        qf_h_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        qf_h_0.solution = phi_b_n + np.asarray(qf_h_0.physical_nodes) * (xi - phi_b_n) / qf_h_0.physical_nodes[-1]
        self.__qf_h = TreeMesh1DUniform(qf_h_0, aligned=True)
        n_e_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        n_e_0.solution = self.qf_e_to_n_e(qf_e_0.solution)
        self.__n_e = TreeMesh1DUniform(n_e_0, aligned=True)
        n_h_0 = Mesh1DUniform(0.0, self.__diode.__length, physical_step=self.__initial_step)
        n_h_0.solution = self.qf_h_to_n_h(qf_h_0.solution)
        self.__n_h = TreeMesh1DUniform(n_h_0, aligned=True)

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] qf_e_to_n_e(self, double[:] qf_e):
        cdef:
            int n = qf_e.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__diode.__semiconductor.n_e_ev_t(qf_e[i], self.__temperature)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] qf_h_to_n_h(self, double[:] qf_h):
        cdef:
            int n = qf_h.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__diode.__semiconductor.n_h_ev_t(qf_h[i], self.__temperature)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_e_to_qf_e(self, double[:] n_e):
        cdef:
            int n = n_e.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__diode.__semiconductor.n_e_to_mu_ev_t(n_e[i], self.__temperature)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double[:] n_h_to_qf_h(self, double[:] n_h):
        cdef:
            int n = n_h.shape[0]
            int i
            array[double] result, template = array('d')
        result = clone(template, n, zero=False)
        for i in range(n):
            result[i] = self.__diode.__semiconductor.n_h_to_mu_ev_t(n_h[i], self.__temperature)
        return result

    @boundscheck(False)
    @wraparound(False)
    cpdef double thermionic_emission_current_e(self):
        cdef:
            double n0, nb, a, nc, t2
        n0 = self.__n_e.flatten().solution[0]
        nb = self.__diode.n0_t(self.__temperature)
        nc = self.__diode.__semiconductor.n_c_t(self.__temperature)
        a = self.__diode.__semiconductor.__reference['thermionic_emission']['A_R_coeff_n'] * constant.__A_R
        t2 = self.__temperature * self.__temperature
        return a * t2 * (n0 - nb) / nc

    @boundscheck(False)
    @wraparound(False)
    cpdef double thermionic_emission_current_h(self):
        cdef:
            double p0, pb, a, nv, t2
        p0 = self.__n_h.flatten().solution[0]
        pb = self.__diode.p0_t(self.__temperature)
        nv = self.__diode.__semiconductor.n_v_t(self.__temperature)
        a = self.__diode.__semiconductor.__reference['thermionic_emission']['A_R_coeff_p'] * constant.__A_R
        t2 = self.__temperature * self.__temperature
        return a * t2 * (pb - p0) / nv

    @boundscheck(False)
    @wraparound(False)
    cpdef double grad_qf_e_bc(self):
        cdef:
            double n0, p0, mu, j, nd = 0, ef, h
            Mesh1D ep = self.__ep.flatten()
            Dopant dopant
        n0 = self.__n_e.flatten().solution[0]
        p0 = self.__n_h.flatten().solution[0]
        h = ep.jacobian * (ep.local_nodes[1] - ep.local_nodes[0])
        ef = (ep.solution[0] - ep.solution[1]) / h
        for dopant in self.__diode.__semiconductor.__dopants:
            nd += dopant.__concentration.flatten().solution[0]
        j = self.thermionic_emission_current_e()
        mu = self.__diode.__semiconductor.mobility_e_point_t(nd, ef, p0*n0, self.__temperature)
        return -j / (constant.__q * n0 * mu)

    @boundscheck(False)
    @wraparound(False)
    cpdef double grad_qf_h_bc(self):
        cdef:
            double n0, p0, mu, j, nd = 0, ef, h
            Mesh1D ep = self.__ep.flatten()
            Dopant dopant
        n0 = self.__n_e.flatten().solution[0]
        p0 = self.__n_h.flatten().solution[0]
        h = ep.jacobian * (ep.local_nodes[1] - ep.local_nodes[0])
        ef = (ep.solution[0] - ep.solution[1]) / h
        for dopant in self.__diode.__semiconductor.__dopants:
            nd += dopant.__concentration.flatten().solution[0]
        j = self.thermionic_emission_current_h()
        mu = self.__diode.__semiconductor.mobility_h_point_t(nd, ef, p0*n0, self.__temperature)
        return -j / (constant.__q * p0 * mu)
