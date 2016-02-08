from __future__ import division
from functools import partial

import warnings
import sympy as sym
import numpy as np
import mpmath as mp
from scipy.optimize import root

from Schottky.Notation import k, q
from Schottky.Helpers import to_numeric, fermi, d_fermi_d_delta_fermi_energy, centered_linspace
from Potential.Potential_1D import ConstantPotential, LinearPotential, SuperposedPotential

energy_distribution_functions = {'Single Level': ['single level', 'monoenergetic level', 'single', 'monoenergetic'],
                                 'Gaussian Level': ['gaussian level', 'gaussian'],
                                 'Rectangular Level': ['rectangular level', 'rectangular', 'box']}
np.seterr(over='ignore', under='ignore')


class Trap(object):
    """
    Describes Carrier Trap
    """

    def __init__(self, name, charge_states=None, energy_distribution_function='Single Level', energy_spread=0.3 * q,
                 electron_capture_cross_section=1e-17, electron_capture_cross_section_activation_energy=0,
                 hole_capture_cross_section=1e-17, hole_capture_cross_section_activation_energy=0, trap_potential=None):
        if not charge_states:
            charge_states = [[+1, sym.symbols('E_d'), 1], [0, sym.symbols('E_d'), 1]]
        if not trap_potential:
            ef = LinearPotential('External Field', -1e6 * 0)
            np = ConstantPotential('No potential', 0)
            trap_potential = SuperposedPotential('Superposed', [np, ef])
        self.name = str(name)
        self.charge_states = self.__check_charge_states(charge_states)
        self.energy_distribution_function = energy_distribution_function.lower()
        self.energy_spread = energy_spread
        self.energy_distribution = self.__build_energy_distribution_function(self.energy_distribution_function,
                                                                             self.energy_spread)
        self.electron_capture_cross_section = electron_capture_cross_section
        self.electron_capture_cross_section_activation_energy = electron_capture_cross_section_activation_energy
        self.hole_capture_cross_section = hole_capture_cross_section
        self.hole_capture_cross_section_activation_energy = hole_capture_cross_section_activation_energy
        self.trap_potential = trap_potential

    @staticmethod
    def __build_energy_distribution_function(energy_distribution_function, energy_spread):
        conduction_band_energy, energy, trap_central_energy_level = sym.symbols('Ecb E Et')
        if energy_distribution_function in energy_distribution_functions['Single Level']:
            energy_distribution = sym.DiracDelta(energy - conduction_band_energy + trap_central_energy_level)
        elif energy_distribution_function in energy_distribution_functions['Gaussian Level']:
            energy_distribution = sym.exp(-((energy - conduction_band_energy + trap_central_energy_level) ** 2)
                                          / (2 * energy_spread ** 2))
            energy_distribution /= (sym.sqrt(sym.pi) * energy_spread)
            energy_distribution *= sym.sqrt(2) / 2
        elif energy_distribution_function in energy_distribution_functions['Rectangular Level']:
            energy_offset = conduction_band_energy - trap_central_energy_level
            energy_distribution = (sym.sign(energy - energy_offset + energy_spread / 2) + 1) / 2
            energy_distribution *= (sym.sign(energy_offset + energy_spread / 2 - energy) + 1) / 2
            energy_distribution /= energy_spread
        else:
            raise Exception('The distribution function supplied is not supported!')
        return energy_distribution

    @staticmethod
    def __check_charge_states(charge_states):
        """
        Charge states are described in the form of the list of two lists:
        empty (no electron) full (electron trapped)

        [[+1, sym.symbols('E_d'), 1], [0, sym.symbols('E_d'), 1]]

        for each state there is a list of [state_charge, energy distance from Ec, degeneracy]
        """
        # TODO: add charge states validation code
        return charge_states

    def energy_level(self, temperature, semiconductor, charge_state_idx=0, electron_volts=False):
        """
        Trap energy level measured from Ec (Conduction band) towards Ev (Valence band)
        """
        if isinstance(temperature, (tuple, list, np.ndarray)):
            energy_level_wrap = partial(self.energy_level, semiconductor=semiconductor,
                                        charge_state_idx=charge_state_idx, electron_volts=electron_volts)
            return np.array(map(energy_level_wrap, temperature))
        band_gap = semiconductor.band_gap(temperature, symbolic=False, electron_volts=False)
        try:
            level = self.charge_states[charge_state_idx][1].subs('Eg', band_gap)
        except AttributeError:
            level = self.charge_states[charge_state_idx][1]
        if electron_volts:
            level /= q
        return to_numeric(level)

    def equilibrium_f(self, temperature, semiconductor, fermi_level_from_conduction_band,
                      electron_volts=False, use_mpmath=False, debug=False):
        """
        Calculates trap equilibrium filling for given Fermi level position
        :param temperature: Temperature, K
        :param semiconductor: Semiconductor object
        :param fermi_level_from_conduction_band: distance from Conduction band to Fermi level
        :param electron_volts: if True assume all energy values to be in eV
        :param use_mpmath: if True integration is done using mpmath.quad function instead of numpy.trapz (default)
        :param debug: if True prints out some debug information
        :return: equilibrium f between 0 and 1
        """
        if electron_volts:
            energy_coefficient = to_numeric(q)
            energy_unit = 'eV'
        else:
            energy_coefficient = 1
            energy_unit = 'J'
        energy_scale = to_numeric(k * temperature) / energy_coefficient
        conduction_band = 0
        fermi_level = conduction_band - fermi_level_from_conduction_band
        trap_energy_level = self.energy_level(temperature, semiconductor, charge_state_idx=0,
                                              electron_volts=electron_volts)
        if debug:
            print 'Et = %2.2g ' % ((conduction_band - trap_energy_level)) + energy_unit
            print 'Ef =', fermi_level, energy_unit
        g_ratio = self.charge_states[0][2] / self.charge_states[1][2]

        if self.energy_distribution_function in energy_distribution_functions['Single Level']:
            fermi_level_grid, = np.meshgrid(fermi_level)
            exp_arg = (np.float(conduction_band - trap_energy_level) - fermi_level_grid) / energy_scale
            exp_term = np.exp(exp_arg)
            f = 1 / (1 + g_ratio * exp_term)
        else:
            energy = sym.symbols('E')
            energy_distribution = self.energy_distribution.subs([('Et', trap_energy_level * energy_coefficient),
                                                                 ('Ecb', conduction_band * energy_coefficient)])
            energy_distribution = energy_distribution.subs(q, to_numeric(q))
            if not use_mpmath:
                if debug:
                    print 'Numeric integration (numpy.trapz)'
                energy_range = centered_linspace(conduction_band - trap_energy_level,
                                                 10 * to_numeric(self.energy_spread) / energy_coefficient,
                                                 to_numeric(self.energy_spread) / energy_coefficient / 1000)
                energy_range_grid, fermi_level_grid = np.meshgrid(energy_range, fermi_level)
                fermi_function = 1 / (1 + g_ratio * np.exp((energy_range_grid - fermi_level_grid) / energy_scale))
                energy_distribution_function = sym.lambdify(energy, energy_distribution, 'numpy')
                integrand_array = energy_distribution_function(energy_range_grid * energy_coefficient) * fermi_function
                f = np.trapz(integrand_array, energy_range_grid * energy_coefficient, axis=1)
            else:
                if debug:
                    print 'Numeric integration (mpmath.quad)'
                fermi_level_grid, = np.meshgrid(fermi_level)
                f = np.zeros_like(fermi_level_grid)
                for i, fermi_level_i in enumerate(fermi_level_grid):
                    fermi_function = fermi(energy, fermi_level_i * energy_coefficient, temperature, g_ratio)
                    fermi_function = fermi_function.subs(k, to_numeric(k))
                    integrand = sym.lambdify(energy, energy_distribution * fermi_function)
                    f[i] = mp.quad(integrand,
                                   [to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                               - 10 * self.energy_spread),
                                    to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                               - 0.5 * self.energy_spread),
                                    to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                               + 0.5 * self.energy_spread),
                                    to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                               + 10 * self.energy_spread)])
        if debug:
            print 'F =', f
        if f.size == 1:
            f = f[0]
        return f

    def d_equilibrium_f_d_fermi_energy(self, temperature, semiconductor, fermi_level_from_conduction_band,
                                       electron_volts=False, use_mpmath=False, debug=False):
        """
        Calculates trap equilibrium filling derivative for given Fermi level position
        on small Fermi level change
        :param temperature: Temperature, K
        :param semiconductor: Semiconductor object
        :param fermi_level_from_conduction_band: distance from Conduction band to Fermi level
        :param electron_volts: if True assume all energy values to be in eV
        :param use_mpmath: if True integration is done using mpmath.quad function instead of numpy.trapz (default)
        :param debug: if True prints out some debug information
        :return: equilibrium f between 0 and 1
        """
        if electron_volts:
            energy_coefficient = to_numeric(q)
            energy_unit = 'eV'
        else:
            energy_coefficient = 1
            energy_unit = 'J'
        energy_scale = to_numeric(k * temperature) / energy_coefficient
        conduction_band = 0
        fermi_level = conduction_band - fermi_level_from_conduction_band
        trap_energy_level = self.energy_level(temperature, semiconductor, charge_state_idx=0,
                                              electron_volts=electron_volts)
        if debug:
            print 'Et = %2.2g ' % ((conduction_band - trap_energy_level)) + energy_unit
            print 'Ef =,', fermi_level, energy_unit
        g_ratio = self.charge_states[0][2] / self.charge_states[1][2]

        if self.energy_distribution_function in energy_distribution_functions['Single Level']:
            fermi_level_grid, = np.meshgrid(fermi_level)
            exp_arg = (np.float(conduction_band - trap_energy_level) - fermi_level_grid) / energy_scale
            exp_term = np.exp(exp_arg)
            exp_term = np.array(map(mp.exp, exp_arg))
            a = g_ratio / (energy_coefficient * energy_scale)
            #print exp_term #/ (1 + g_ratio * exp_term)
            b = exp_term / (1 + g_ratio * exp_term) ** 2
            d_f = a * b
            # d_f = g_ratio / (energy_coefficient * energy_scale) * exp_term / (1 + g_ratio * exp_term) ** 2
            d_f[np.where(d_f == np.nan)] = 0
        else:
            energy = sym.symbols('E')
            energy_distribution = self.energy_distribution.subs([('Et', trap_energy_level * energy_coefficient),
                                                                 ('Ecb', conduction_band * energy_coefficient)])
            energy_distribution = energy_distribution.subs(q, to_numeric(q))
            if not use_mpmath:
                if debug:
                    print 'Numeric integration (numpy.trapz)'
                energy_range = centered_linspace(conduction_band - trap_energy_level,
                                                 10 * to_numeric(self.energy_spread) / energy_coefficient,
                                                 to_numeric(self.energy_spread) / energy_coefficient / 1000)
                energy_range_grid, fermi_level_grid = np.meshgrid(energy_range, fermi_level)
                exp_term = np.exp((energy_range_grid - fermi_level_grid) / energy_scale)
                fore_factor = g_ratio / (energy_coefficient * energy_scale)
                d_fermi_function = fore_factor * exp_term / (1 + g_ratio * exp_term) ** 2
                d_fermi_function[np.where(d_fermi_function == np.nan)] = 0
                energy_distribution_function = sym.lambdify(energy, energy_distribution, 'numpy')
                integrand_array = energy_distribution_function(energy_range_grid * energy_coefficient) * d_fermi_function
                d_f = np.trapz(integrand_array, energy_range_grid * energy_coefficient, axis=1)
            else:
                if debug:
                    print 'Numeric integration (mpmath.quad)'
                fermi_level_grid, = np.meshgrid(fermi_level)
                d_f = np.zeros_like(fermi_level_grid)
                for i, fermi_level_i in enumerate(fermi_level_grid):
                    d_fermi_function = d_fermi_d_delta_fermi_energy(energy, fermi_level_i * energy_coefficient,
                                                                    temperature, g_ratio)
                    d_fermi_function = d_fermi_function.subs(k, to_numeric(k))
                    integrand = sym.lambdify(energy, energy_distribution * d_fermi_function)
                    d_f[i] = mp.quad(integrand,
                                     [to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                                 - 10 * self.energy_spread),
                                      to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                                 - 0.5 * self.energy_spread),
                                      to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                                 + 0.5 * self.energy_spread),
                                      to_numeric(energy_coefficient * (conduction_band - trap_energy_level)
                                                 + 10 * self.energy_spread)])
        if debug:
            print 'dF =', d_f
        if d_f.size == 1:
            d_f = d_f[0]
        return d_f

    def f_to_equilibrium_fermi_level(self, temperature, semiconductor, f, electron_volts=False,
                                     use_mpmath=False, parallel=False, debug=False):
        """
        Calculates equilibrium Fermi level position for given occupation of the Trap F
        :param temperature: Temperature, K
        :param semiconductor: Semiconductor object
        :param f: Trap occupation from 0.0 to 1.0
        :param electron_volts: if True assume all energy values to be in eV
        :param use_mpmath: if True integration is done using mpmath.quad function instead of numpy.trapz (default)
        :param debug: if True prints out some debug information
        :return: Fermi level as distance from Conduction band to Fermi level

        In calculation we use eV since solver has much better stability in smaller order numbers
        """
        energy_unit = 'eV'
        energy_coefficient = to_numeric(q)
        trap_energy_level = self.energy_level(temperature, semiconductor, charge_state_idx=0,
                                              electron_volts=True)
        f_grid, = np.meshgrid(f)
        #print f_grid
        def equation(fermi_level_from_conduction_band, f):
            test_f = self.equilibrium_f(temperature, semiconductor, fermi_level_from_conduction_band,
                                        electron_volts=True, use_mpmath=use_mpmath, debug=debug)
            residual = f - test_f
            if debug:
                print 'Fermi level type:', type(fermi_level_from_conduction_band)
                print 'Test F =', test_f
                print 'F residual =', residual
            return residual

        def solver(args):
            equation, lower_boundary, upper_boundary, initial_guess, f, use_mpmath = args
            if not use_mpmath:
                warnings.filterwarnings('ignore')
                solution = root(equation, initial_guess, args=f, method='hybr')
                solution = solution.x[0]
                warnings.resetwarnings()
            else:
                equation_wrap = partial(equation, f=f)
                try:
                    solution = mp.findroot(equation_wrap, (lower_boundary, upper_boundary),
                                           maxsteps=1000, solver='anderson', tol=5e-16)
                except ValueError as err:
                    print err
                    print 'Lowering tolerance to 5e-6'
                    solution = mp.findroot(equation_wrap, (lower_boundary, upper_boundary),
                                           maxsteps=1000, solver='anderson', tol=5e-6)
                solution = np.float(solution)
            return solution

        fermi_level_lower_boundary = abs(to_numeric(trap_energy_level - 2 * self.energy_spread / energy_coefficient))
        fermi_level_upper_boundary = abs(to_numeric(trap_energy_level + 2 * self.energy_spread / energy_coefficient))
        if debug:
            print 'Fermi level lower boundary = %2.2g ' % fermi_level_lower_boundary + energy_unit
            print 'Fermi level upper boundary = %2.2g ' % fermi_level_upper_boundary + energy_unit
        args = np.empty((len(f_grid), 6), dtype=object)
        args[:, 0] = equation
        args[:, 1] = fermi_level_lower_boundary
        args[:, 2] = fermi_level_upper_boundary
        args[:, 3] = (fermi_level_lower_boundary + fermi_level_upper_boundary) / 2
        args[:, 4] = f_grid
        args[:, 5] = use_mpmath
        if parallel:
            try:
                from pathos.pools import ProcessPool as Pool
                pool = Pool()
                solutions = np.array(pool.map(solver, args))
            except ImportError:
                print 'Parallel calculation needs pathos! Using standard map() instead.'
                solutions = np.array(map(solver, args))
        else:
            solutions = np.array(map(solver, args))
        if not electron_volts:
            solutions *= energy_coefficient
        return solutions

    def capture_cross_sections(self, temperature):
        """
        Calculates temperature dependent capture cross section
        :param temperature: Temperature in K
        :return: Temperature dependent capture cross sections
        """
        energy_scale = to_numeric(k * temperature)
        exp_term_e = np.exp(-self.electron_capture_cross_section_activation_energy / energy_scale)
        exp_term_h = np.exp(-self.hole_capture_cross_section_activation_energy / energy_scale)
        return self.electron_capture_cross_section * exp_term_e, self.hole_capture_cross_section * exp_term_h

    def capture_rate(self, temperature, semiconductor, f, n, p, debug=False):
        """
        Calculate carriers capture rate for bot electrons and holes
        :param temperature: Temperature, K
        :param semiconductor: Semiconductor object
        :param f: Trap occupation from 0.0 to 1.0
        :param n: concentration of electrons at trap coordinate (1/m^3)
        :param p: concentration of holes at trap coordinate (1/m^3)
        :param debug: if True prints out some debug information
        :return: capture_e, capture_h

        Function does not account for barrier-limited capture and other fenomena
        """
        v_e = np.float(semiconductor.v_T('e', temperature, symbolic=False))
        v_h = np.float(semiconductor.v_T('h', temperature, symbolic=False))
        if debug:
            print '<v_e> =', v_e, 'm/s'
            print '<v_h> =', v_h, 'm/s'
        sigma_n, sigma_p = self.capture_cross_sections(temperature)
        c_n = sigma_n * v_e * n
        c_p = sigma_p * v_h * p
        capture_rate_e = np.copy(c_n)
        capture_rate_e = capture_rate_e.reshape(capture_rate_e.size)
        capture_rate_e[np.where(capture_rate_e == 0)] = 1e-250
        capture_rate_h = np.copy(c_p)
        capture_rate_h = capture_rate_h.reshape(capture_rate_h.size)
        capture_rate_h[np.where(capture_rate_h == 0)] = 1e-250
        capture_time_constant_e = 1 / capture_rate_e
        capture_time_constant_h = 1 / capture_rate_h
        if capture_time_constant_e.size == 1:
            capture_time_constant_e = capture_time_constant_e[0]
        if capture_time_constant_h.size == 1:
            capture_time_constant_h = capture_time_constant_h[0]
        if debug:
            print 'c_n =', c_n
            print 'c_p =', c_p
        capture_e = c_n * (1 - f)
        capture_h = c_p * f
        if debug:
            print 'capt_e =', capture_e
            print 'capt_h =', capture_h
        return capture_e, capture_h, capture_time_constant_e, capture_time_constant_h

    def emission_rate(self, temperature, semiconductor, f, poole_frenkel_e=1.0, poole_frenkel_h=1.0,
                      barrier_lowering_e=None, barrier_lowering_h=None, use_mpmath=False, debug=False):
        """
        Calculate carriers emission rate for bot electrons and holes
        :param temperature: Temperature, K
        :param semiconductor: Semiconductor object
        :param f: Trap occupation from 0.0 to 1.0
        :param poole_frenkel_e: emission rate boost due to Poole-Frenkel effect for electron
        :param poole_frenkel_h: emission rate boost due to Poole-Frenkel effect for electron
        :param barrier_lowering_e: lowering of activation energy for electrons
        :param barrier_lowering_h: lowering of activation energy for holes
        :param use_mpmath: if True integration is done using mpmath.quad function instead of numpy.trapz (default)
        :param debug: if True prints out some debug information
        :return: emission_e, emission_h
        """
        if barrier_lowering_e is None:
            barrier_lowering_e = np.zeros_like(f, dtype=np.float)
        if barrier_lowering_h is None:
            barrier_lowering_h = np.zeros_like(f, dtype=np.float)
        energy_scale = to_numeric(k * temperature)
        conduction_band = 0
        band_gap = semiconductor.band_gap(temperature, symbolic=False, electron_volts=False)
        valence_band = conduction_band - band_gap
        trap_energy_level_e = np.float(self.energy_level(temperature, semiconductor,
                                                         charge_state_idx=1, electron_volts=False))
        trap_energy_level_h = np.float(self.energy_level(temperature, semiconductor,
                                                         charge_state_idx=0, electron_volts=False))
        trap_energy_level_e_positive = conduction_band - trap_energy_level_e
        trap_energy_level_h_positive = conduction_band - trap_energy_level_h

        g_ratio_e = self.charge_states[0][2] / self.charge_states[1][2]
        g_ratio_h = self.charge_states[1][2] / self.charge_states[0][2]

        v_e = np.float(semiconductor.v_T('e', temperature, symbolic=False))
        v_h = np.float(semiconductor.v_T('h', temperature, symbolic=False))
        if debug:
            print '<v_e> =', v_e, 'm/s'
            print '<v_h> =', v_h, 'm/s'
        sigma_n, sigma_p = self.capture_cross_sections(temperature)
        fore_factor_n = np.float(sigma_n * v_e * semiconductor.Nc(temperature, symbolic=False) * g_ratio_e)
        fore_factor_p = np.float(sigma_p * v_h * semiconductor.Nv(temperature, symbolic=False) * g_ratio_h)
        fore_factor_n *= poole_frenkel_e
        fore_factor_p *= poole_frenkel_h
        if debug:
            print 'factor_n =', fore_factor_n
            print 'factor_p =', fore_factor_p

        if self.energy_distribution_function in energy_distribution_functions['Single Level']:
            activation_energy_e = conduction_band - trap_energy_level_e_positive - barrier_lowering_e
            activation_energy_h = trap_energy_level_h_positive - valence_band - barrier_lowering_h
            emission_rate_e = fore_factor_n * np.exp(-activation_energy_e / energy_scale)
            emission_rate_h = fore_factor_p * np.exp(-activation_energy_h / energy_scale)
            emission_time_constant_e = 1 / emission_rate_e
            emission_time_constant_h = 1 / emission_rate_h
            emission_e = emission_rate_e * f
            emission_h = emission_rate_h * (1 - f)
        else:
            quasi_fermi_level = self.f_to_equilibrium_fermi_level(temperature, semiconductor, f,
                                                                  electron_volts=False,
                                                                  use_mpmath=use_mpmath, debug=False)
            quasi_fermi_level = conduction_band - quasi_fermi_level
            if debug:
                print 'Eqf =', quasi_fermi_level / to_numeric(q), 'eV'
            energy = sym.symbols('E')
            energy_distribution_e = self.energy_distribution.subs([('Et', trap_energy_level_e),
                                                                   ('Ecb', conduction_band)])
            energy_distribution_e = energy_distribution_e.subs(q, to_numeric(q))
            energy_distribution_h = self.energy_distribution.subs([('Et', trap_energy_level_h),
                                                                   ('Ecb', conduction_band)])
            energy_distribution_h = energy_distribution_h.subs(q, to_numeric(q))
            energy_range_e = centered_linspace(conduction_band - trap_energy_level_e,
                                               10 * to_numeric(self.energy_spread),
                                               to_numeric(self.energy_spread) / 1000)
            energy_range_h = centered_linspace(conduction_band - trap_energy_level_h,
                                               10 * to_numeric(self.energy_spread),
                                               to_numeric(self.energy_spread) / 1000)
            energy_distribution_function_e = sym.lambdify(energy, energy_distribution_e, 'numpy')
            energy_distribution_function_h = sym.lambdify(energy, energy_distribution_h, 'numpy')
            energy_range_grid_e, barrier_lowering_grid_e = np.meshgrid(energy_range_e, barrier_lowering_e)
            energy_range_grid_h, barrier_lowering_grid_h = np.meshgrid(energy_range_h, barrier_lowering_h)
            exp_term_e = np.exp(-(conduction_band - energy_range_grid_e + barrier_lowering_grid_e) / energy_scale)
            exp_term_h = np.exp(-(energy_range_grid_h - barrier_lowering_grid_h - valence_band) / energy_scale)
            emission_rate_e = energy_distribution_function_e(energy_range_e) * exp_term_e * fore_factor_n
            emission_rate_h = energy_distribution_function_h(energy_range_h) * exp_term_h * fore_factor_p
            emission_rate_e_max = np.max(emission_rate_e, axis=1)
            emission_rate_h_max = np.max(emission_rate_h, axis=1)
            emission_time_constant_e = 1 / emission_rate_e_max
            emission_time_constant_h = 1 / emission_rate_h_max
            if not use_mpmath:
                if debug:
                    print 'Numeric integration (numpy.trapz)'
                energy_range_grid_e, fermi_level_grid_e = np.meshgrid(energy_range_e, quasi_fermi_level)
                energy_range_grid_h, fermi_level_grid_h = np.meshgrid(energy_range_h, quasi_fermi_level)
                fermi_function_e = 1 / (1 + g_ratio_e * np.exp((energy_range_grid_e - fermi_level_grid_e) / energy_scale))
                fermi_function_h = 1 - 1 / (1 + g_ratio_h * np.exp((energy_range_grid_h - fermi_level_grid_h) / energy_scale))
                exp_term_e = np.exp(-(conduction_band - energy_range_grid_e + barrier_lowering_e) / energy_scale)
                exp_term_h = np.exp(-(energy_range_grid_h - barrier_lowering_h - valence_band) / energy_scale)
                emission_rate_e = energy_distribution_function_e(energy_range_grid_e) * exp_term_e
                emission_rate_h = energy_distribution_function_h(energy_range_grid_h) * exp_term_h
                integrand_array_e = emission_rate_e * fermi_function_e
                integrand_array_h = emission_rate_h * fermi_function_h
                emission_e = np.trapz(integrand_array_e, energy_range_grid_e, axis=1) * fore_factor_n
                emission_h = np.trapz(integrand_array_h, energy_range_grid_h, axis=1) * fore_factor_p
            else:
                if debug:
                    print 'Numeric integration (mpmath.quad)'
                exp_term_e = sym.exp(-(conduction_band - energy + barrier_lowering_e) / energy_scale)
                exp_term_h = sym.exp(-(energy - barrier_lowering_h - valence_band) / energy_scale)
                quasi_fermi_level_grid, = np.meshgrid(quasi_fermi_level)
                emission_e = np.zeros_like(quasi_fermi_level_grid)
                emission_h = np.zeros_like(quasi_fermi_level_grid)
                for i, quasi_fermi_level_i in enumerate(quasi_fermi_level_grid):
                    fermi_function_e = fermi(energy, quasi_fermi_level_i, temperature, g_ratio_e)
                    fermi_function_e = fermi_function_e.subs(k, to_numeric(k))
                    fermi_function_h = fermi(quasi_fermi_level_i, energy, temperature, g_ratio_h)
                    fermi_function_h = fermi_function_h.subs(k, to_numeric(k))
                    integrand_e = sym.lambdify(energy, energy_distribution_e * fermi_function_e * exp_term_e)
                    integrand_h = sym.lambdify(energy, energy_distribution_h * fermi_function_h * exp_term_h)
                    emission_integral_e = mp.quad(integrand_e,
                                                  [to_numeric(trap_energy_level_e_positive - 10 * self.energy_spread),
                                                   to_numeric(trap_energy_level_e_positive - 0.5 * self.energy_spread),
                                                   to_numeric(trap_energy_level_e_positive + 0.5 * self.energy_spread),
                                                   to_numeric(trap_energy_level_e_positive + 10 * self.energy_spread)])
                    emission_integral_h = mp.quad(integrand_h,
                                                  [to_numeric(trap_energy_level_h_positive - 10 * self.energy_spread),
                                                   to_numeric(trap_energy_level_h_positive - 0.5 * self.energy_spread),
                                                   to_numeric(trap_energy_level_h_positive + 0.5 * self.energy_spread),
                                                   to_numeric(trap_energy_level_h_positive + 10 * self.energy_spread)])
                    emission_e[i] = np.float(emission_integral_e * fore_factor_n)
                    emission_h[i] = np.float(emission_integral_h * fore_factor_p)
        if isinstance(emission_e, np.ndarray):
            if emission_e.size == 1:
                emission_e = emission_e[0]
                emission_h = emission_h[0]
                emission_time_constant_e = emission_time_constant_e[0]
                emission_time_constant_h = emission_time_constant_h[0]
        if debug: print 'emission_e =', emission_e
        if debug: print 'emission_h =', emission_h
        if debug: print 'emission_tau_e =', emission_time_constant_e
        if debug: print 'emission_tau_h =', emission_time_constant_e
        return emission_e, emission_h, emission_time_constant_e, emission_time_constant_h

    def df_dt(self, temperature, semiconductor, f, n, p,
              poole_frenkel_e=1.0, poole_frenkel_h=1.0,
              barrier_lowering_e=None, barrier_lowering_h=None,
              use_mpmath=False, debug=False):
        capture_e, capture_h, capture_tau_e, capture_tau_h = self.capture_rate(temperature, semiconductor, f,
                                                                               n, p, debug=debug)
        emission_e, emission_h, emission_tau_e, emission_tau_h = self.emission_rate(temperature, semiconductor, f,
                                                                                    poole_frenkel_e,
                                                                                    poole_frenkel_h,
                                                                                    barrier_lowering_e,
                                                                                    barrier_lowering_h,
                                                                                    use_mpmath=use_mpmath, debug=debug)
        capture_tau_e = np.copy(capture_tau_e)
        capture_tau_e = capture_tau_e.reshape(capture_tau_e.size)
        capture_tau_h = np.copy(capture_tau_h)
        capture_tau_h = capture_tau_h.reshape(capture_tau_h.size)
        emission_tau_e = np.copy(emission_tau_e)
        emission_tau_e = emission_tau_e.reshape(emission_tau_e.size)
        emission_tau_h = np.copy(emission_tau_h)
        emission_tau_h = emission_tau_h.reshape(emission_tau_h.size)
        min_tau = min(min(capture_tau_e), min(capture_tau_h), min(emission_tau_e), min(emission_tau_h))
        return capture_e - capture_h - emission_e + emission_h, min_tau

    def energy_distribution_diagram(self, ax, temperature, semiconductor,
                                    trap_concentration=1, trap_concentration_units='',
                                    fermi_level_from_conduction_band = None,
                                    electron_volts=True, fancy_labels=False):
        if electron_volts:
            energy_coefficient = to_numeric(q)
            energy_unit = 'eV'
        else:
            energy_coefficient = 1
            energy_unit = 'J'
        energy_scale = to_numeric(k * temperature) / energy_coefficient
        conduction_band = 0
        g_ratio = self.charge_states[0][2] / self.charge_states[1][2]
        band_gap = semiconductor.band_gap(temperature, symbolic=False, electron_volts=electron_volts)
        energy_range = np.linspace(-np.float(band_gap), 0, num=1001, endpoint=True)
        charge_state_idx = 0
        trap_energy_level = self.energy_level(temperature, semiconductor, charge_state_idx, electron_volts)
        if fermi_level_from_conduction_band is None:
            fermi_level_from_conduction_band = trap_energy_level
        fermi_level = conduction_band - fermi_level_from_conduction_band
        if self.energy_distribution_function in energy_distribution_functions['Single Level']:
            ax.plot(np.zeros_like(energy_range), energy_range, linewidth=2, color='black', linestyle='-')
            ax_max = ax.get_xlim()[1]
            if trap_concentration > ax_max:
                head_length = 0.03 * trap_concentration
            else:
                head_length = 0.03 * ax_max
            ax.arrow(0, conduction_band - trap_energy_level, trap_concentration - head_length, 0,
                     linewidth=2, head_width=head_length, head_length=head_length, fc='black', ec='black')
            f = self.equilibrium_f(temperature, semiconductor, fermi_level_from_conduction_band, electron_volts)
            ax.arrow(0, conduction_band - trap_energy_level, trap_concentration * f, 0,
                     linewidth=2, head_width=0, head_length=0, fc='red', ec='red')
        else:
            energy = sym.symbols('E')
            energy_distribution = self.energy_distribution.subs([('Et', trap_energy_level * energy_coefficient),
                                                                 ('Ecb', conduction_band * energy_coefficient)])
            energy_distribution = energy_distribution.subs(q, to_numeric(q))
            energy_distribution_function = sym.lambdify(energy, energy_distribution, 'numpy')
            fermi_function = 1 / (1 + g_ratio * np.exp((energy_range - fermi_level) / energy_scale))
            energy_states_distribution = energy_distribution_function(energy_range * energy_coefficient)
            energy_states_distribution *= trap_concentration * energy_coefficient
            trapped_carriers_distribution = energy_states_distribution * fermi_function
            ax.plot(energy_states_distribution, energy_range, linewidth=2, color='black', linestyle='-')
            # ax.plot(trapped_carriers_distribution, energy_range, linewidth=2, color='red', linestyle='-')
            ax.plot(fermi_function * max(energy_states_distribution), energy_range,
                    linewidth=2, color='black', linestyle='--')
            ax.fill_between(trapped_carriers_distribution, 0, energy_range, color='blue', alpha=0.5)
        ax.arrow(0, fermi_level, ax.get_xlim()[1], 0,
                 linewidth=2, head_width=0, head_length=0, fc='black', ec='black')
        ax.set_title('Traps distribution in the Band Gap')
        ax.set_ylabel('Energy, ' + energy_unit)
        ax.set_ylim([-np.float(band_gap), 0])
        x_label = 'Traps distribution, 1/' + energy_unit
        if trap_concentration_units != 1 and trap_concentration_units != '':
            x_label += ' * ' + trap_concentration_units
        ax.set_xlabel(x_label)
        ax.set_xlim([0, max([ax.get_xlim()[1], trap_concentration])])
        if fancy_labels:
            ticks = ax.get_yticks()
            labels = np.array([lbl.get_text() for lbl in ax.get_yticklabels()])
            # print ticks
            # print labels
            if u'Ev' in labels:
                ticks = np.append(ticks, band_gap - trap_energy_level)
                labels = np.append(labels, 'Ec-%2.2g' % trap_energy_level)
            else:
                #ticks = np.array([0, band_gap, band_gap - trap_energy_level])
                ticks = np.array([-band_gap, -trap_energy_level, 0])
                labels = np.array(['Ev', 'Ec-%2.2g' % trap_energy_level, 'Ec'])
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)

    def __str__(self, *args, **kwargs):
        description = self.name
        for state in self.charge_states:
            description += ' Charge state: ' + str(state[0]) + '@' + str(state[1])
        return description
