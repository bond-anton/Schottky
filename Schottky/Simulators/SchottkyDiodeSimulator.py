from __future__ import division, print_function
from functools import partial
import numpy as np
from scipy.special import lambertw
from scipy.interpolate import interp1d

from BDPoisson1D import dirichlet_non_linear_poisson_solver_amr

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators import Simulator
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky import constants
from ._helpers import storage_manager, prepare_array, fermi, d_fermi_d_delta_fermi_energy


class SchottkyDiodeSimulator(Simulator):

    def __init__(self, client, diode, description=None):
        assert isinstance(diode, SchottkyDiode), 'Valid SchottkyDiode Sample object expected'
        self.diode = diode
        samples = [self.diode]
        name = 'Schottky Diode Simulator'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}

        measurement_types = [
            {'name': 'Diode energetics',
             'description': 'Measurement of diode energetics',
             'children': []},
            {'name': 'Carriers concentration',
             'description': 'Measurement of carriers concentration in the Schottky diode',
             'children': []},
            {'name': 'Current Measurement',
             'description': 'Measurement of electrical current',
             'children': []}]

        self.measurement_details = {
            'v_bi': {
                'name': 'Built-in potential temperature dependence',
                'description': 'Measurement of Diode Built-in potential temperature dependence',
                'type': 'Diode energetics'},
            'free carriers concentration': {
                'name': 'Free carriers concentration profile',
                'description': 'Measurement of free carriers concentration profile',
                'type': 'Carriers concentration'},
            'iv': {
                'name': 'IV',
                'description': 'measure I-V dependence',
                'type': 'Current Measurement'},
        }

        parts = [BulkSemiconductor(client=client, semiconductor=diode.semiconductor)]
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=list(self.measurement_details.values()))

        self.measurement_specs = {
            'v_bi': {'parameters': None,
                     'input data': None,
                     'variables': [{
                         'name': 'Temperature',
                         'description': 'Sample temperature',
                         'units': 'K'
                     }],
                     'result': [{
                         'name': 'Built-in potential',
                         'description': 'Diode\'s Built-in potential',
                         'units': 'eV'
                         }]},
            'free carriers concentration': {'parameters': [{'name': 'temperature',
                                                            'type': 'numeric',
                                                            'default value': 0.0,
                                                            'description': 'Temperature',
                                                            'units': 'K'}],
                                            'input data': None,
                                            'variables': [{'name': 'z',
                                                           'description': 'Z coordinate',
                                                           'units': 'cm'}],
                                            'result': [
                                                {
                                                    'name': 'electrons',
                                                    'description': 'electrons concentration',
                                                    'units': 'cm^-3'
                                                },
                                                {
                                                    'name': 'holes',
                                                    'description': 'holes concentration',
                                                    'units': 'cm^-3'
                                                }
                                            ]},
            'dopants occupation level': {'parameters': [{'name': 'temperature',
                                                         'type': 'numeric',
                                                         'default value': 0.0,
                                                         'description': 'Temperature',
                                                         'units': 'K'}],
                                         'input data': None,
                                         'variables': [{'name': 'z',
                                                        'description': 'Z coordinate',
                                                        'units': 'cm'}],
                                         'result': [
                                             {
                                                 'name': dopant.trap.name,
                                                 'description': 'dopant occupation',
                                                 'units': ''
                                             } for dopant in self.diode.semiconductor.dopants]}
        }

    @storage_manager('v_bi', use_storage=True)
    def v_bi(self, temperature=0.0):
        temperature = prepare_array(temperature)
        metal_wf = self.diode.metal.work_function
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        semiconductor_wf = xi + self.diode.semiconductor.affinity
        bi_voltage = metal_wf - semiconductor_wf
        return bi_voltage

    @storage_manager('free carriers concentration', use_storage=True)
    def free_carriers_concentration(self, z=0.0, temperature=0.0):
        psi = lambda x: np.zeros_like(x)
        return self._free_carriers_concentration(z=z, psi=psi, temperature=temperature)

    def thermionic_emission_current(self, bias=0.0, area=None, temperature=0.0):
        """
        Thermionic emission current simulation.
        :param bias: 1D array or float value in Volts.
        :param area: Diode's area in cm^2
        :param temperature: Temperature in K
        :return: Diode's current as 1D array of the same shape as bias
        """
        if area is None:
            area = self.diode.area
        bias = prepare_array(bias)
        temperature = prepare_array(temperature)
        energy_scale = constants['k'] * temperature
        diode_type = self._diode_type(temperature=temperature)
        band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        metal_wf = self.diode.metal.work_function
        if diode_type == 'n':
            phi_bn = metal_wf - self.diode.semiconductor.affinity
            a_r = self.diode.semiconductor.thermo_emission_coefficient['electron'] * constants['A_R']
        elif diode_type == 'p':
            phi_bn = self.diode.semiconductor.affinity + band_gap - metal_wf
            print(self.diode.semiconductor.affinity, band_gap, metal_wf)
            a_r = self.diode.semiconductor.thermo_emission_coefficient['hole'] * constants['A_R']
        else:
            raise ValueError('Cannot determine conductivity type of diode')
        j_s = a_r * (temperature ** 2) * np.exp(-phi_bn / energy_scale)
        j = -j_s + energy_scale / (area * self.diode.serial_resistance) *\
            lambertw((area * self.diode.serial_resistance * j_s / energy_scale) *
                     np.exp((bias + area * j_s * self.diode.serial_resistance) / energy_scale))
        return np.real(j)

    def potential(self, bias=0.0, temperature=0.0, psi=None):
        band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature)
        diode_type = self._diode_type(temperature=temperature)
        v_bi = self.v_bi(temperature=temperature)[0]
        if diode_type == 'p':
            bias_coeff = 1
        elif diode_type == 'n':
            bias_coeff = -1
        else:
            raise ValueError('Can not determine diode type')
        v_diode = bias_coeff * bias
        if v_bi + v_diode <= 0.0:
            v_diode = -v_bi * 0.9
        rho = partial(self._semiconductor_charge, temperature=temperature,
                      band_gap=band_gap, xi=xi, dos=dos, diode_type=diode_type)
        d_rho_d_psi = partial(self._d_semiconductor_charge_d_psi, temperature=temperature,
                              band_gap=band_gap, xi=xi, dos=dos, diode_type=diode_type)
        if psi is None:
            psi = lambda x: (bias + v_diode) - x * (bias + v_diode) / self.diode.thickness
        converged = False
        while not converged:
            meshes = dirichlet_non_linear_poisson_solver_amr(0.0, self.diode.thickness, 1.0e-6,
                                                             psi, rho, d_rho_d_psi, v_bi + v_diode, 0.0,
                                                             max_iter=1000, residual_threshold=1.5e-3,
                                                             int_residual_threshold=1.5e-4,
                                                             max_level=5, mesh_refinement_threshold=1e-7)

            flat_grid = meshes.flatten()
            j = self.thermionic_emission_current(bias=bias, temperature=temperature)[0]
            v_serial = j * self.diode.area * self.diode.serial_resistance
            v_diode_new = bias_coeff * (bias - v_serial)
            print(v_diode, v_diode_new, v_serial, abs(v_diode - v_diode_new))
            if abs(v_diode - v_diode_new) < 1e-6:
                converged = True
            else:
                v_diode = v_diode_new
                psi = interp1d(flat_grid.physical_nodes, flat_grid.solution, bounds_error=False,
                               fill_value=(v_bi + v_diode, 0))
        return flat_grid

    def _diode_type(self, temperature=0.0):
        p, n = self.parts['Bulk Semiconductor Simulator'].get_type(temperature=temperature)
        p_type = n_type = False
        for region in p:
            if region.size > 0:
                p_type = True
        for region in n:
            if region.size > 0:
                n_type = True
        if p_type and not n_type:
            diode_type = 'p'
        elif n_type and not p_type:
            diode_type = 'n'
        else:
            diode_type = 'n/d'
        return diode_type

    def _semiconductor_charge(self, z=0.0, psi=None, temperature=0.0,
                              band_gap=None, xi=None, dos=None, diode_type=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature)
        total_charge = 0
        # bands charge
        carriers_concentration = self._free_carriers_concentration(z=z, psi=psi, temperature=temperature,
                                                                   band_gap=band_gap, xi=xi, dos=dos,
                                                                   diode_type=diode_type)
        total_charge += carriers_concentration['holes'] - carriers_concentration['electrons']
        # dopants charge
        dopants_occupation = self._dopants_equilibrium_occupation(z=z, psi=psi, temperature=temperature,
                                                                  band_gap=band_gap, xi=xi)
        for dopant in self.diode.semiconductor.dopants:
            n = dopant.concentration
            q1 = dopant.trap.charge_state['full']
            q0 = dopant.trap.charge_state['empty']
            total_charge += n * (q1 - q0) * dopants_occupation[dopant.trap.name] + n * q0
        return total_charge * constants['q'] / constants['epsilon_0'] / self.diode.semiconductor.epsilon

    def _free_carriers_concentration(self, z=0.0, psi=None, temperature=0.0,
                                     band_gap=None, xi=None, dos=None, diode_type=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        energy_scale = constants['k'] * temperature
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature)
        fermi_level = psi(z) + xi
        n = dos['DOS C.band'] * np.exp(-fermi_level / energy_scale)
        p = dos['DOS V.band'] * np.exp((fermi_level - band_gap) / energy_scale)
        if diode_type == 'n':
            p = np.zeros_like(z)
        elif diode_type == 'p':
            n = np.zeros_like(z)
        return {'electrons': n, 'holes': p}

    def _dopants_equilibrium_occupation(self, z=0.0, psi=None, temperature=0.0,
                                        band_gap=None, xi=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        fermi_level = psi(z) + xi
        result = {}
        for part_name in self.parts['Bulk Semiconductor Simulator'].parts:
            if not 'Charge Carrier Trap Simulator' in part_name:
                continue
            dopant = self.parts['Bulk Semiconductor Simulator'].parts[part_name]
            energy = dopant.energy_level(band_gap=band_gap)
            result[dopant.trap.name] = fermi(energy=energy, fermi_energy=fermi_level,
                                             temperature=temperature, g_ratio=1)
        return result

    def _d_free_carriers_concentration_d_psi(self, z=0.0, psi=None, temperature=0.0,
                                             band_gap=None, xi=None, dos=None, diode_type=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        energy_scale = constants['k'] * temperature
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature)
        fermi_level = psi(z) + xi
        dn_dpsi = -dos['DOS C.band'] * np.exp(-fermi_level / energy_scale) / energy_scale
        dp_dpsi = dos['DOS V.band'] * np.exp((fermi_level - band_gap) / energy_scale) / energy_scale
        if diode_type == 'n':
            dp_dpsi = np.zeros_like(z)
        elif diode_type == 'p':
            dn_dpsi = np.zeros_like(z)
        return {'electrons': dn_dpsi, 'holes': dp_dpsi}

    def _d_dopants_equilibrium_occupation_d_psi(self, z=0.0, psi=None, temperature=0.0,
                                                band_gap=None, xi=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        fermi_level = psi(z) + xi
        result = {}
        for part_name in self.parts['Bulk Semiconductor Simulator'].parts:
            if not 'Charge Carrier Trap Simulator' in part_name:
                continue
            dopant = self.parts['Bulk Semiconductor Simulator'].parts[part_name]
            energy = dopant.energy_level(band_gap=band_gap)
            result[dopant.trap.name] = d_fermi_d_delta_fermi_energy(energy=energy, fermi_energy=fermi_level,
                                                                    temperature=temperature, g_ratio=1)
        return result

    def _d_semiconductor_charge_d_psi(self, z=0.0, psi=None, temperature=0.0,
                                      band_gap=None, xi=None, dos=None, diode_type=None):
        if psi is None:
            psi = lambda x: np.zeros_like(x)
        z = prepare_array(z)
        temperature = prepare_array(temperature)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature)
        d_total_charge_d_psi = 0
        # bands charge
        d_carriers_concentration_d_psi = self._d_free_carriers_concentration_d_psi(z=z, psi=psi,
                                                                                   temperature=temperature,
                                                                                   band_gap=band_gap, xi=xi, dos=dos,
                                                                                   diode_type=diode_type)
        d_total_charge_d_psi += d_carriers_concentration_d_psi['holes'] - d_carriers_concentration_d_psi['electrons']
        # dopants charge
        d_dopants_occupation_d_psi = self._d_dopants_equilibrium_occupation_d_psi(z=z, psi=psi,
                                                                                  temperature=temperature,
                                                                                  band_gap=band_gap, xi=xi)
        for dopant in self.diode.semiconductor.dopants:
            n = dopant.concentration
            q1 = dopant.trap.charge_state['full']
            q0 = dopant.trap.charge_state['empty']
            d_total_charge_d_psi += n * (q1 - q0) * d_dopants_occupation_d_psi[dopant.trap.name] + n * q0
        return d_total_charge_d_psi * constants['q'] / constants['epsilon_0'] / self.diode.semiconductor.epsilon
