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
            {'name': 'Traps occupation',
             'description': 'Measurement of traps occupation in the Schottky diode',
             'children': []},
            {'name': 'Current measurement',
             'description': 'Measurement of electrical current',
             'children': []}]
        measurement_specs = {
            'v_bi': {'name': 'Built-in potential temperature dependence',
                     'description': 'Measurement of Diode Built-in potential temperature dependence',
                     'type': 'Diode energetics',
                     'parameters': None,
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
            'potential': {'name': 'Electron potential in the diode',
                          'description': 'Measurement of Diode potential depth-profile',
                          'type': 'Diode energetics',
                          'parameters':
                              [
                                  {'name': 'temperature',
                                           'type': 'numeric',
                                           'default value': 0.0,
                                           'description': 'Temperature',
                                           'units': 'K'},
                                  {'name': 'bias',
                                           'type': 'numeric',
                                           'default value': 0.0,
                                           'description': 'Voltage bias',
                                           'units': 'V'}
                              ],
                          'input data': None,
                          'variables': None,
                          'result': [
                              {
                                  'name': 'z coordinate',
                                  'description': 'Z coordinate',
                                  'units': 'cm'
                              },
                              {
                                  'name': 'potential',
                                  'description': 'potential',
                                  'units': 'V'
                              },
                              {
                                  'name': 'residual error',
                                  'description': 'Calculation residual error',
                                  'units': 'V/cm^2'
                              },
                              {
                                  'name': 'current density',
                                  'description': 'Current density through the diode',
                                  'units': 'A/cm^2'
                              },
                              {
                                  'name': 'diode voltage',
                                  'description': 'Voltage drop on a diode',
                                  'units': 'V'
                              },
                              {
                                  'name': 'diode voltage error',
                                  'description': 'Voltage drop on a diode error',
                                  'units': 'V'
                              }]},
            'free carriers concentration': {'name': 'Free carriers concentration profile',
                                            'description': 'Measurement of free carriers concentration profile',
                                            'type': 'Carriers concentration',
                                            'parameters':
                                                [
                                                    {'name': 'temperature',
                                                             'type': 'numeric',
                                                             'default value': 0.0,
                                                             'description': 'Temperature',
                                                             'units': 'K'},
                                                    {'name': 'bias',
                                                             'type': 'numeric',
                                                             'default value': 0.0,
                                                             'description': 'Voltage bias',
                                                             'units': 'V'}
                                                ],
                                            'input data': None,
                                            'variables': [{'name': 'z_coordinate',
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
                                                }]},
            'dopants equilibrium occupation': {'name': 'Dopants equilibrium occupation',
                                               'description': 'Dopants equilibrium occupation profile',
                                               'type': 'Traps occupation',
                                               'parameters':
                                                   [
                                                       {'name': 'temperature',
                                                                'type': 'numeric',
                                                                'default value': 0.0,
                                                                'description': 'Temperature',
                                                                'units': 'K'},
                                                       {'name': 'bias',
                                                                'type': 'numeric',
                                                                'default value': 0.0,
                                                                'description': 'Voltage bias',
                                                                'units': 'V'}
                                                   ],
                                               'input data': None,
                                               'variables': [{'name': 'z_coordinate',
                                                              'description': 'Z coordinate',
                                                              'units': 'cm'}],
                                               'result': [
                                                   {
                                                       'name': dopant.trap.name,
                                                       'description': 'dopant occupation',
                                                       'units': ''
                                                   } for dopant in self.diode.semiconductor.dopants]
                                               },
            'iv': {'name': 'IV',
                   'description': 'measure I-V dependence',
                   'type': 'Current measurement',
                   'parameters':
                       [
                           {'name': 'temperature',
                            'type': 'numeric',
                            'default value': 0.0,
                            'description': 'Temperature',
                            'units': 'K'}
                       ],
                   'input data': None,
                   'variables': [
                       {'name': 'bias',
                        'description': 'Voltage bias',
                        'units': 'V'}],
                   'result': [
                       {
                           'name': 'current density',
                           'description': 'Current density through the diode',
                           'units': 'A/cm^2'
                       }]}
        }

        assert isinstance(diode, SchottkyDiode), 'Valid SchottkyDiode Sample object expected'
        self.__diode = diode
        samples = [self.diode]
        parts = [BulkSemiconductor(client=client, semiconductor=diode.semiconductor)]

        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurement_specs=measurement_specs)

    @property
    def diode(self):
        return self.__diode

    @storage_manager('v_bi', use_storage=True)
    def v_bi(self, temperature=0.0):
        temperature = prepare_array(temperature)
        metal_wf = self.diode.metal.work_function
        if temperature.size > 1:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=True)
        else:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        semiconductor_wf = xi + self.diode.semiconductor.affinity
        bi_voltage = metal_wf - semiconductor_wf
        return bi_voltage

    @storage_manager('free carriers concentration', use_storage=True)
    def free_carriers_concentration(self, z_coordinate=0.0, bias=0.0, temperature=0.0):
        psi_grid = self.potential(bias=bias, temperature=temperature)
        psi = interp1d(psi_grid['z coordinate'], psi_grid['potential'], bounds_error=True)
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                  use_storage=False)
        band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                       use_storage=False)
        dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                           use_storage=False)
        diode_type = self._diode_type(temperature=temperature, use_storage=False)
        return self._free_carriers_concentration(z=z_coordinate, psi=psi, temperature=temperature,
                                                 band_gap=band_gap, xi=xi, dos=dos, diode_type=diode_type)

    @storage_manager('dopants equilibrium occupation', use_storage=True)
    def dopants_equilibrium_occupation(self, z_coordinate=0.0, bias=0.0, temperature=0.0):
        psi_grid = self.potential(bias=bias, temperature=temperature)
        psi = interp1d(psi_grid['z_coordinate'], psi_grid['potential'], bounds_error=True)
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                  use_storage=False)
        band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                       use_storage=False)
        return self._dopants_equilibrium_occupation(z=z_coordinate, psi=psi, temperature=temperature,
                                                    band_gap=band_gap, xi=xi)

    @storage_manager('iv', use_storage=True)
    def thermionic_emission_current(self, bias=0.0, temperature=0.0):
        """
        Thermionic emission current simulation.
        :param bias: 1D array or float value in Volts.
        :param temperature: Temperature in K
        :return: Diode's current as 1D array of the same shape as bias
        """
        bias = prepare_array(bias)
        j = []
        for v in bias:
            psi_grid = self.potential(bias=v, temperature=temperature)
            j.append(psi_grid['current density'][0])
        return np.array(j)

    @storage_manager('potential', use_storage=True)
    def potential(self, bias=0.0, temperature=0.0, psi=None):
        band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                       use_storage=False)
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                  use_storage=False)
        dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                           use_storage=False)
        diode_type = self._diode_type(temperature=temperature, use_storage=False)
        v_bi = self.v_bi(temperature=temperature, use_storage=False)[0]
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
            def psi(x): return (bias + v_diode) - x * (bias + v_diode) / self.diode.thickness
        converged = False
        while not converged:
            meshes = dirichlet_non_linear_poisson_solver_amr(0.0, self.diode.thickness, 1.0e-6,
                                                             psi, rho, d_rho_d_psi, v_bi + v_diode, 0.0,
                                                             max_iter=1000, residual_threshold=1.5e-3,
                                                             int_residual_threshold=1.5e-4,
                                                             max_level=5, mesh_refinement_threshold=1e-7)

            flat_grid = meshes.flatten()
            phi_bn = self._phi_bn(psi_nodes=flat_grid.solution, v_diode=v_diode, temperature=temperature,
                                  xi=xi, band_gap=band_gap, diode_type=diode_type)
            j = self._thermionic_emission_current(bias=bias, temperature=temperature,
                                                  diode_type=diode_type, phi_bn=phi_bn, area=self.diode.area)[0]
            v_serial = j * self.diode.area * self.diode.serial_resistance
            v_diode_new = bias_coeff * (bias - v_serial)
            print(v_diode, v_diode_new, v_serial, abs(v_diode - v_diode_new))
            if abs(v_diode - v_diode_new) < 1e-6:
                converged = True
            else:
                v_diode = v_diode_new
                psi = interp1d(flat_grid.physical_nodes, flat_grid.solution, bounds_error=False,
                               fill_value=(v_bi + v_diode, 0))
        return {'z coordinate': flat_grid.physical_nodes,
                'potential': flat_grid.solution,
                'residual error': flat_grid.residual,
                'current density': np.array([j]),
                'diode voltage': np.array([v_diode]),
                'diode voltage error': np.array([abs(v_diode - v_diode_new)])}

    def _phi_bn(self, psi_nodes=0.0, v_diode=0.0, temperature=0.0, xi=None, band_gap=None, diode_type=None):
        psi_nodes = prepare_array(psi_nodes)
        v_diode = prepare_array(v_diode)
        temperature = prepare_array(temperature)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature,
                                          use_storage=False)
        bias_sign = 1
        if diode_type == 'n':
            v_b = max(psi_nodes)
            bias_sign = -1
        elif diode_type == 'p':
            v_b = min(psi_nodes)
            xi = band_gap - xi
        else:
            raise ValueError('Diode type is undefined')
        return abs(v_b) + bias_sign * v_diode + xi

    def _thermionic_emission_current(self, bias=0.0, temperature=0.0, diode_type=None, phi_bn=None, area=None):
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
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature, use_storage=False)
        if phi_bn is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature, use_storage=False)
            metal_wf = self.diode.metal.work_function
        if diode_type == 'n':
            if phi_bn is None:
                phi_bn = metal_wf - self.diode.semiconductor.affinity
            a_r = self.diode.semiconductor.thermo_emission_coefficient['electron'] * constants['A_R']
        elif diode_type == 'p':
            if phi_bn is None:
                phi_bn = self.diode.semiconductor.affinity + band_gap - metal_wf
            a_r = self.diode.semiconductor.thermo_emission_coefficient['hole'] * constants['A_R']
        else:
            raise ValueError('Cannot determine conductivity type of a diode')
        j_s = a_r * (temperature ** 2) * np.exp(-phi_bn / energy_scale)
        j = -j_s + energy_scale / (area * self.diode.serial_resistance) *\
            lambertw((area * self.diode.serial_resistance * j_s / energy_scale) *
                     np.exp((bias + area * j_s * self.diode.serial_resistance) / energy_scale))
        return np.real(j)

    def _conduction_band(self, z=0.0, psi_nodes=None, v_diode=0.0, temperature=0.0, xi=None):
        z = prepare_array(z)
        v_diode = prepare_array(v_diode)
        temperature = prepare_array(temperature)
        if psi_nodes is None:
            psi_nodes = np.zeros_like(z)
        else:
            psi_nodes = prepare_array(psi_nodes)
        if xi is None:
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        return psi_nodes + v_diode + xi

    def _valence_band(self, z=0.0, psi_nodes=None, v_diode=0.0, temperature=0.0, band_gap=None, xi=None):
        temperature = prepare_array(temperature)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        return self._conduction_band(z=z, psi_nodes=psi_nodes, v_diode=v_diode,
                                     temperature=temperature, xi=xi) - band_gap

    def _diode_type(self, temperature=0.0, use_storage=True):
        p, n = self.parts['Bulk Semiconductor Simulator'].get_type(temperature=temperature, use_storage=use_storage)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                               use_storage=False)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature,
                                          use_storage=False)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                               use_storage=False)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature, use_storage=False)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                               use_storage=False)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature,
                                          use_storage=False)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
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
            xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature,
                                                                                      use_storage=False)
        if band_gap is None:
            band_gap = self.parts['Bulk Semiconductor Simulator'].band_gap(temperature=temperature,
                                                                           use_storage=False)
        if dos is None:
            dos = self.parts['Bulk Semiconductor Simulator'].effective_bands_density_of_states(temperature=temperature,
                                                                                               use_storage=False)
        if diode_type is None:
            diode_type = self._diode_type(temperature=temperature,
                                          use_storage=False)
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
