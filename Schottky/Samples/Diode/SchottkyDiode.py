#!-*- coding: utf-8 -*-

'''
Created on 02 дек. 2014 г.

@author: anton
'''

from __future__ import division

import datetime
import warnings

import mpmath as mp
import numpy as np
from scipy.constants.codata import epsilon0
from scipy.optimize import fmin

from Schottky import constants
from Schottky.Helpers import smooth_dd, Psi_zero, to_numeric  # , interp_Fn#, fermihalf, gen_interpolated_P_E_function
from Schottky.Notation import q, k
from Schottky.Samples.Metal import Metal as MetalElectrode
from Schottky.Samples.Semiconductor import Semiconductor  # , Trap, Dopant, Dislocation, BondingInterface


class SchottkyDiode(object):
    '''
    Schottky Diode class
    '''

    equipment = 'Schottky Diode Simulator'
    description = 'Simulates properties of a Schottky diode'

    def __init__(self, Project, label, Metal=None, Semiconductor=None, Area=None, Rs=None, T=None, Va=None,
                 DeadLayer=None, L=None):
        '''
        Constructor
        '''
        self.Project = Project
        self.Project.log_record('Project database opened', 'Information')
        self.tool_id = self.Project.add_tool('1D_SchottkySym', '1D Schottky Diode Simulator')
        self.measurements_types = self.register_measurement_types()

        self.label = self.Project.add_project_info_if_changed('Project name', label)

        self.set_electrode(Metal)
        self.set_semiconductor(Semiconductor)
        self.set_contact_area(Area)
        self.set_Rs(Rs)
        self.set_DeadLayer(DeadLayer)
        self.FieldInversionPoint = 0.0
        self.FieldInversionPsi = 0.0
        self.set_L(L)

        self.rho_z_Psi_memo = {}
        self.Psi_memo = {}

        self.set_T(T)
        self.set_Va(Va)

    def __str__(self, *args, **kwargs):
        record = 'Schottky diode: %s' % (self.label)
        return record

    def set_electrode(self, Metal=None):
        measurement_type = self.Project.add_measurement_type('Metal work function',
                                                             measurement_type_description='Metal work function measurement')
        try:
            metal_label = self.Project.get_project_info('Schottky electrode')
            measurement_id = self.Project.add_measurement(measurement_type, 'Electrode work function',
                                                          'Work function of metal electrode',
                                                          measurement_datetime=None, create_new=False)
            WF_arr, _, _, _ = self.Project.get_data_points(measurement_id, 'Metal electrode work function')
            if len(WF_arr) > 0:
                WF = WF_arr[-1]
        except:
            metal_label = None
            WF = 0
        if Metal is None:
            if metal_label is None:
                raise Exception('Metal electriode is needed for Schottky diode')
            else:
                self.Metal = MetalElectrode(metal_label, WF)
        else:
            if Metal.label != metal_label:
                self.Project.add_project_info_if_changed('Schottky electrode', Metal.label)
                _, WF = self.Project.add_monitored_datapoint(measurement_type, 'Electrode work function',
                                                             'Work function of metal electrode',
                                                             self.tool_id, 'Metal electrode work function', 'Metal WF',
                                                             'J', None, to_numeric(Metal.WF))
            self.Metal = Metal

    def set_semiconductor(self, Semicond=None):
        # measurement_type = self.Project.add_measurement_type('Metal work function', measurement_type_description='Metal work function measurement')
        try:
            semiconductor_label = self.Project.get_project_info('Semiconductor')
        except:
            semiconductor_label = None
        if Semicond is None:
            if semiconductor_label is None:
                raise Exception('Semiconductor is needed for Schottky diode')
            else:
                self.Semiconductor = Semiconductor(semiconductor_label, lookup=True)
        else:
            if Semicond.reference['name'] != semiconductor_label:
                self.Project.add_project_info_if_changed('Semiconductor', Semicond.reference['name'])
            self.Semiconductor = Semicond

    def register_measurement_types(self):
        measurements_dict = {
            'Length measurement': self.Project.add_measurement_type('Length measurement',
                                                                    measurement_type_description='Length measurement'),
            'Area measurement': self.Project.add_measurement_type('Area measurement',
                                                                  measurement_type_description='Area measurement'),
            'Resistance measurement': self.Project.add_measurement_type('Resistance measurement',
                                                                        measurement_type_description='Resistance measurement'),
            'Temperature measurement': self.Project.add_measurement_type('Temperature measurement',
                                                                         measurement_type_description='Temperature measurement'),
            'Built-in voltage measurement': self.Project.add_measurement_type('Built-in voltage measurement',
                                                                              measurement_type_description='Built-in voltage measurement'),
            'Bias voltage': self.Project.add_measurement_type('Bias voltage monitor',
                                                              measurement_type_description='Applied bias voltage measurement'),
            'Fields and Electric Field measurement': self.Project.add_measurement_type(
                'Fields and Electric Field',
                measurement_type_description='Fields and Electric Field in the diode measurement'),
            'Traps kinetics measurement': self.Project.add_measurement_type(
                'Transient experiment',
                measurement_type_description='Time-resolved measurement of traps kinetic'),
        }
        return measurements_dict

    def set_contact_area(self, Area):
        point_name_long = 'Area'
        point_name_short = 'A'
        point_unit_name = 'm^-2'
        default_d = 1.5e-3
        default = np.pi * (default_d ** 2) / 4.0
        _, self.Area = self.Project.add_monitored_datapoint(self.measurements_types['Area measurement'], 'Contact area',
                                                            'Contact area measurement',
                                                            self.tool_id, point_name_long, point_name_short,
                                                            point_unit_name, default, Area)

    def set_contact_diameter(self, d=1.5e-3):
        if not isinstance(d, (float, int)):
            Area = None
        else:
            Area = np.pi * (d ** 2) / 4.0
        self.set_contact_area(Area)

    def set_L(self, L=10e-6):
        point_name_long = 'Schottky diode semiconductor thickness'
        point_name_short = 'L'
        point_unit_name = 'm'
        _, self.L = self.Project.add_monitored_datapoint(self.measurements_types['Length measurement'],
                                                         'Diode thickness', 'Thickness of semiconductor',
                                                         self.tool_id, point_name_long, point_name_short,
                                                         point_unit_name, 10e-6, L)

    def set_Rs(self, Rs=10):
        point_name_long = 'Serial Resistance'
        point_name_short = 'Rs'
        point_unit_name = 'Ohm'
        _, self.Rs = self.Project.add_monitored_datapoint(self.measurements_types['Resistance measurement'],
                                                          'Serial Resistance', 'Serial resistance of diode base',
                                                          self.tool_id, point_name_long, point_name_short,
                                                          point_unit_name, 10, Rs)

    def set_DeadLayer(self, DeadLayer=0.0):
        point_name_long = 'Dead Layer thickness'
        point_name_short = 'DL'
        point_unit_name = 'm'
        _, self.DeadLayer = self.Project.add_monitored_datapoint(self.measurements_types['Length measurement'],
                                                                 'Dead Layer', 'Dead Layer thickness',
                                                                 self.tool_id, point_name_long, point_name_short,
                                                                 point_unit_name, 0, DeadLayer)

    def set_T(self, T=300):
        point_name_long = 'Temperature'
        point_name_short = 'T'
        point_unit_name = 'K'
        self.T_measurement_id, self.T = self.Project.add_monitored_datapoint(
            self.measurements_types['Temperature measurement'], 'Temperature', 'Diode temperature',
            self.tool_id, point_name_long, point_name_short, point_unit_name, 300, T)
        try:
            self.Project.log_record('Temperature: %2.2f K' % (T), 'Information')
        except:
            pass

    def set_Va(self, V=0):
        point_name_long = 'Bias'
        point_name_short = 'Va'
        point_unit_name = 'V'
        self.Va_measurement_id, self.Va = self.Project.add_monitored_datapoint(self.measurements_types['Bias voltage'],
                                                                               'Va', 'Diode bias voltage monitor',
                                                                               self.tool_id, point_name_long,
                                                                               point_name_short, point_unit_name, 0, V)
        try:
            self.Project.log_record('Va: %2.2f V' % (V), 'Information')
        except:
            pass

    def V_bi(self, eV=False):
        measurement_id = self.Project.add_measurement(self.measurements_types['Built-in voltage measurement'], 'Vbi',
                                                      'Built-in voltage',
                                                      measurement_datetime=None, create_new=False)
        point_name_long = 'Vbi'
        Vbi_array, _, measurement_date_array, unit_name_array = self.Project.get_data_points(measurement_id,
                                                                                             point_name_long)
        for i in range(len(Vbi_array)):
            T, _ = self.Project.get_data_point_at_datetime(self.T_measurement_id, 'Temperature',
                                                           measurement_date_array[i], interpolation_type='step')
            if T == self.T:
                Vbi = Vbi_array[i]
                if eV:
                    if unit_name_array[i] == 'J':
                        Vbi /= to_numeric(q)
                else:
                    if unit_name_array[i] == 'V':
                        Vbi *= to_numeric(q)
                return np.float(Vbi)
        point_name_short = 'Vbi'
        point_unit_name = 'V' if eV else 'J'
        point_order = self.Project.get_next_data_point_order(measurement_id, point_name_long)
        if eV:
            Vbi = to_numeric(self.Metal.WF / q) - self.Semiconductor.WF(self.T, z=1e3, eV=eV)
        else:
            Vbi = to_numeric(self.Metal.WF) - self.Semiconductor.WF(self.T, z=1e3, eV=eV)
        self.Project.add_datapoint(point_order, measurement_id, self.tool_id,
                                   point_name_long, point_name_short, point_unit_name,
                                   Vbi, datetime.datetime.now())
        return np.float(Vbi)

    def EfEc(self, Psi=Psi_zero, z=0, eV=False):
        coeff = 1 if eV else to_numeric(q)
        Psi_nodes = Psi(z)
        xi = np.float(self.Semiconductor.Ech_pot(T=self.T, z=1e3, eV=eV, debug=False))
        return -coeff * Psi_nodes + xi

    def dopants_set_eq_F(self, z, Psi, debug=False):
        '''
        Sets dopants filling level to equilibrium
        '''
        for dopant in self.Semiconductor.dopants:
            Ef = self.EfEc(Psi, z, eV=False)
            if isinstance(z, (mp.mpf, float, int)):
                F = dopant.equilibrium_f(self.T, self.Semiconductor, Ef, electron_volts=False, debug=debug)
                dF = dopant.d_equilibrium_f_d_fermi_energy(self.T, self.Semiconductor, Ef, electron_volts=False, debug=debug)
            elif isinstance(z, (np.ndarray, list, tuple)):
                F = dopant.equilibrium_f(self.T, self.Semiconductor, Ef, electron_volts=False, debug=debug)
                dF = dopant.d_equilibrium_f_d_fermi_energy(self.T, self.Semiconductor, Ef, electron_volts=False, debug=debug)
                # F = np.array([dopant.equilibrium_f(self.T, self.Semiconductor, Ef_z, eV=False, debug=debug) for Ef_z in Ef])
                # dF = np.array([dopant.d_equilibrium_f_d_fermi_energy(self.T, self.Semiconductor, Ef_z, eV=False, debug=debug) for Ef_z in Ef])
            else:
                raise ValueError('Wrong type of Z: ' + str(type(z)))
            dopant.set_F_interp(z, F)
            dopant.set_dF_interp(z, dF)

    def disl_eq_F(self, Ef, disl_traps, eV=True, debug=False):
        F = []
        dF = []
        for i, trap in enumerate(disl_traps):
            F.append(trap[0].equilibrium_f(self.T, self.Semiconductor, Ef, eV, debug))
            dF.append(trap[0].d_equilibrium_f_d_fermi_energy(self.T, self.Semiconductor, Ef, eV, debug))
            if debug: print 'F = %2.2g' % F[i]
            if debug: print 'dF = %2.2g' % dF[i]
        return np.array(F), np.array(dF)

    def BI_set_eq_F(self, BI, Psi, eV=False, debug=False):
        Ef = self.EfEc(Psi, BI.depth, eV)
        BI.dsl_tilt_f, BI.dsl_tilt_df = self.disl_eq_F(Ef, BI.dsl_tilt.traps, eV, debug)
        BI.dsl_twist_f, BI.dsl_twist_df = self.disl_eq_F(Ef, BI.dsl_twist.traps, eV, debug)
        BI.set_traps_f(BI.dsl_tilt_f, BI.dsl_twist_f)
        BI.set_traps_df(BI.dsl_tilt_df, BI.dsl_twist_df)

    def build_rho_z_Psi(self, Psi=Psi_zero, carriers_charge=False):
        if not (Psi, self.T, carriers_charge) in self.rho_z_Psi_memo:
            def rho_z_Psi(z, Psi):
                rho = 0
                if carriers_charge:
                    if Psi != Psi_zero:
                        n, p = self.n_carriers_theory(Psi, z)
                        if self.Semiconductor.dop_type == 'n':
                            p = 0
                        elif self.Semiconductor.dop_type == 'p':
                            n = 0
                        rho += p - n
                for dopant in self.Semiconductor.dopants:
                    Nd = dopant.concentration(z)
                    rho += Nd * (
                        (dopant.charge_states[1][0] - dopant.charge_states[0][0]) * dopant.F(z) +
                        dopant.charge_states[0][
                            0])
                for BI in self.Semiconductor.bonding_interfaces:
                    rho += smooth_dd(z - BI.depth, BI.smooth_dirac_epsilon) * BI.density_of_charge
                return to_numeric(q) * rho

            self.rho_z_Psi_memo[(Psi, self.T, carriers_charge)] = rho_z_Psi
        return self.rho_z_Psi_memo[(Psi, self.T, carriers_charge)]

    def get_phi_bn(self, Psi=Psi_zero, Va=0, SchottkyEffect=False):
        chg_sign = 1 if self.Semiconductor.dop_type == 'n' else -1
        F = to_numeric(q / (16 * np.pi * epsilon0 * self.Semiconductor.reference['epsilon']))
        if not SchottkyEffect:
            F = 0

        def f(z):
            if F != 0:
                if not isinstance(z, np.ndarray):
                    return chg_sign * (Psi(z) + chg_sign * F / z)
                else:
                    return chg_sign * np.array([Psi(zz) + chg_sign * F / zz for zz in z])
            else:
                if not isinstance(z, np.ndarray):
                    return chg_sign * Psi(z)
                else:
                    #return chg_sign * np.array([Psi(zz) for zz in z])
                    return chg_sign * Psi(z)
        z_0 = np.linspace(1e-10, np.float(self.L), num=50)
        extrema_z = np.zeros_like(z_0)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        for i in range(z_0.size):
            extrema_z[i] = fmin(f, z_0[i], disp=False)
            # print extrema_z[i]*1e6
        warnings.resetwarnings()
        extrema = f(extrema_z)
        # _, (ax1, ax2, ax3) = plt.subplots(3)
        # ax1.plot(z_0*1e6, f(z_0))
        # ax2.plot(z_0*1e6, extrema_z*1e6)
        # ax3.plot(extrema_z*1e6, extrema)
        # plt.show()

        z_max = extrema_z[np.argmin(extrema)]
        if z_max < 0:
            z_max = 1e-15
        xi = self.Semiconductor.Ech_pot(T=self.T, z=1e3, eV=True, debug=False)
        Eg = self.Semiconductor.band_gap(self.T, symbolic=False, electron_volts=True)
        if self.Semiconductor.dop_type == 'n':
            phi_bn = abs(-f(z_max) + xi + Va)
            delta_phi = abs(self.Ec(Psi, Va, z_max, eV=True) - phi_bn)
            phi_b = phi_bn - xi
        else:
            phi_bn = abs(f(z_max) + xi - Va - Eg)
            delta_phi = abs(self.Ev(Psi, Va, z_max, eV=True) - phi_bn)
            phi_b = phi_bn - Eg + xi
        # print 'Z_max =', z_max*1e6, 'um'
        return np.float(z_max), np.float(phi_bn), np.float(delta_phi), np.float(phi_b)

    def ThermionicEmissionCurrent(self, Va, phi_bn, debug=False):
        kT = to_numeric(k * self.T)
        q_n = to_numeric(q)
        A = self.Area
        Rs = self.Rs
        if self.Semiconductor.dop_type == 'n':
            Ar = self.Semiconductor.reference['A_R_coeff_n'] * constants['A_R']
        else:
            Ar = self.Semiconductor.reference['A_R_coeff_p'] * constants['A_R']
        Js = Ar * (self.T ** 2) * mp.exp(-q_n * phi_bn / kT)
        if debug: print 'Js, Is =', Js, A * Js
        J = -Js + kT / (q_n * A * Rs) * mp.lambertw((q_n * A * Rs * Js / kT) * mp.exp(q_n * (Va + A * Js * Rs) / kT))
        if debug: print 'J, I =', J, A * J
        Vd = Va - A * J * Rs
        return np.float(Vd), np.float(J)

    def Ec(self, Psi=Psi_zero, Va=0, z=0, eV=False):
        coeff = 1.0 if eV else to_numeric(q)
        Psi_nodes = Psi(z)
        xi = np.float(self.Semiconductor.Ech_pot(T=self.T, z=1e3, eV=eV, debug=False))
        type_sign = -1 if self.Semiconductor.dop_type == 'p' else 1
        Ec = -coeff * Psi_nodes + xi + coeff * type_sign * Va
        return Ec

    def Ev(self, Psi=Psi_zero, Va=0, z=0, eV=False):
        Eg = self.Semiconductor.band_gap(self.T, symbolic=False, electron_volts=eV)
        Ec = self.Ec(Psi, Va, z, eV=eV)
        Ev = Ec - Eg
        return Ev

    def n_carriers_theory(self, Psi=Psi_zero, z=0):
        '''
        Charge carrier concentration in the main band
        '''
        k_n = to_numeric(k)
        Ef = self.EfEc(Psi, z, eV=False)
        Eg = self.Semiconductor.band_gap(self.T, symbolic=False, electron_volts=False)
        F_n = np.exp(-Ef / (k_n * self.T))
        F_p = np.exp((Ef - Eg) / (k_n * self.T))
        n = np.float(self.Semiconductor.Nc(self.T, symbolic=False)) * F_n
        p = np.float(self.Semiconductor.Nv(self.T, symbolic=False)) * F_p
        if isinstance(z, np.ndarray):
            idx = np.where(z < self.DeadLayer)[0]
            # print idx
            if idx.size > 0:
                n[idx] = 0.0
                p[idx] = 0.0
            idx = np.where(z < self.FieldInversionPoint)[0]
            if idx.size > 0:
                n[idx] = 0.0
                p[idx] = 0.0
        elif z < self.DeadLayer or z < self.FieldInversionPoint:
            n = 0
            p = 0
        return n, p

    def dn_carriers_dEf_theory(self, Psi=Psi_zero, z=0):
        k_n = to_numeric(k)
        Ef = self.EfEc(Psi, z, eV=False)
        Eg = self.Semiconductor.band_gap(self.T, symbolic=False, electron_volts=False)
        F_n = np.exp(-Ef / (k_n * self.T))
        F_p = np.exp((Ef - Eg) / (k_n * self.T))
        dn = np.float(self.Semiconductor.Nc(self.T, symbolic=False)) / (k_n * self.T) * F_n
        dp = np.float(-self.Semiconductor.Nv(self.T, symbolic=False)) / to_numeric(k * self.T) * F_p
        if isinstance(z, np.ndarray):
            idx = np.where(z < self.DeadLayer)[0]
            if idx.size > 0:
                dn[idx] = 0.0
                dp[idx] = 0.0
            idx = np.where(z < self.FieldInversionPoint)[0]
            if idx.size > 0:
                dn[idx] = 0.0
                dp[idx] = 0.0
        elif z < self.DeadLayer or z < self.FieldInversionPoint:
            dn = 0
            dp = 0
        return dn, dp
