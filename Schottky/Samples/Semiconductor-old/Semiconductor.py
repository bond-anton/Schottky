# -*- coding: utf-8 -*-

'''
Created on 30 окт. 2014 г.

@author: anton
'''

from __future__ import division

import multiprocessing
from functools import partial

import mpmath as mp
import numpy as np
import sympy as sym
from joblib import Parallel, delayed
from matplotlib import pyplot as plt

from Schottky import constants
from Schottky.Helpers import to_pos_pwr, symplify_factor, solve_polynomial_eq
from Schottky.Notation import q, k
from Schottky.Reference import database
from Schottky.Samples.Semiconductor import BondingInterface
from Schottky.Samples.Semiconductor import Dopant


def EchPotEc(Semi, T):
    return mp.mpf(Semi.band_gap(T, symbolic=False, electron_volts=True)) - mp.mpf(Semi.Ech_pot(T, eV=True, debug=True))
    # return 2*T


Fermi = lambda E, Ef, T: 1 / (1 + sym.exp((E - Ef) / (k * T)))  # fermi distribution function


class Semiconductor(object):
    '''
    All properties of Semiconductor-old are realized in this class as methods
    '''

    equipment = 'Bulk Semiconductor-old Simulator'
    description = 'Simulates properties of a bulk semiconductor'

    # Total_charge = lambda coords: 0

    def __init__(self, name, lookup=False):
        '''
        Constructor
        '''

        self.reference = {}
        self.dopants = []
        self.bonding_interfaces = []
        self.dop_type = 'i'
        self.dop_level = 0

        self.reference['name'] = name
        if lookup:
            self.find_references()

        self.EchPot_memo = {}
        self.Eg_memo = {}
        self.WF_memo = {}

        self.Total_charge = self.__neutr_eq()

    def find_references(self):
        for rec in database:
            if rec['name'] == self.reference['name']:
                self.reference = rec
                break

    def add_dopant(self, dopant):
        assert isinstance(dopant, Dopant), "dopant must be an object of Dopant class"
        self.dopants.append(dopant)
        self.Total_charge = self.__neutr_eq()
        # print 'Calculating total coefficient done'
        dop_level = dopant.concentration(1e3)
        if dopant.energy_level(300, self, 0, electron_volts=False) < self.band_gap(300, symbolic=False, electron_volts=False) / 2:
            dop_type = 'n'
        else:
            dop_type = 'p'
        if self.dop_type == dop_type:
            self.dop_level += dop_level
        elif self.dop_type == 'i':
            self.dop_type = dop_type
            self.dop_level = dop_level
        else:
            self.dop_level -= dop_level
            if self.dop_level < 0:
                self.dop_type = dop_type
                self.dop_level = abs(self.dop_level)

    def add_bonding_interface(self, bonding_interface):
        assert isinstance(bonding_interface, BondingInterface), "BI must be an object of BondingInterface class"
        self.bonding_interfaces.append(bonding_interface)

    def Nc(self, T=sym.symbols('T'), symbolic=True):
        if isinstance(T, (tuple, list, np.ndarray)):
            Nc_Wrap = partial(self.Nc, symbolic=symbolic)
            return np.array(map(Nc_Wrap, T))
        Nc = sym.symbols('N_c0') * (T ** (sym.Rational(3 / 2)))
        if symbolic:
            return Nc
        else:
            return sym.N(self.__to_numeric(Nc))

    def Nv(self, T=sym.symbols('T'), symbolic=True):
        if isinstance(T, (tuple, list, np.ndarray)):
            Nv_Wrap = partial(self.Nv, symbolic=symbolic)
            return np.array(map(Nv_Wrap, T))
        Nv = sym.symbols('N_v0') * (T ** (sym.Rational(3 / 2)))
        if symbolic:
            return Nv
        else:
            return sym.N(self.__to_numeric(Nv))

    def v_T(self, carrier='e', T=sym.symbols('T'), symbolic=True):
        m = sym.symbols('m_' + carrier + '_coeff') * sym.symbols('m_e')
        v = sym.sqrt(3 * k * T / m)
        if symbolic:
            return v
        else:
            return sym.N(self.__to_numeric(v))

    def band_gap(self, temperature=sym.symbols('temperature'), symbolic=True, electron_volts=False):
        '''
        Energy gap temperature dependence approximation 
        '''
        if isinstance(temperature, (tuple, list, np.ndarray)):
            E_Gap_Wrap = partial(self.band_gap, symbolic=symbolic, eV=electron_volts)
            return np.array(map(E_Gap_Wrap, temperature))
        if (temperature, symbolic, electron_volts) in self.Eg_memo:
            return self.Eg_memo[(temperature, symbolic, electron_volts)]
        alpha, beta, E_g0 = sym.symbols('alpha beta E_g0')
        if symbolic:
            if electron_volts:
                self.Eg_memo[(temperature, symbolic, electron_volts)] = (E_g0 - alpha * temperature ** 2 / (temperature + beta)) / q
            else:
                self.Eg_memo[(temperature, symbolic, electron_volts)] = (E_g0 - alpha * temperature ** 2 / (temperature + beta))
        else:
            self.Eg_memo[(temperature, symbolic, electron_volts)] = np.float(self.__to_numeric(self.band_gap(temperature, symbolic=True, electron_volts=electron_volts)))
        return self.Eg_memo[(temperature, symbolic, electron_volts)]

    def Ech_pot(self, T=sym.symbols('T'), z=sym.symbols('z'), eV=False, debug=False):
        if isinstance(T, (tuple, list, np.ndarray)):
            Ech_pot_Wrap = partial(self.Ech_pot, z=z, eV=eV, debug=debug)
            return np.array(map(Ech_pot_Wrap, T))
            # num_cores = multiprocessing.cpu_count()
            # return np.array(Parallel(n_jobs=num_cores)(delayed(Ech_pot_Wrap)(i) for i in T))
        if (T, z, eV) in self.EchPot_memo:
            return self.EchPot_memo[(T, z, eV)]
        if debug: print 'T = %2.2f K' % T
        if T < 0:
            raise Exception('T must be >= 0')
        elif T < 1:
            if debug: print 'T < 1: Using derivative linear extrapolation'
            dT = 0.5
            dE = ((self.Ech_pot(1 + dT, z, eV) - self.Ech_pot(1, z, eV))) / dT
            self.EchPot_memo[(T, z, eV)] = self.Ech_pot(1, z, eV) + dE * (T - dT)
            return self.Ech_pot(1, z, eV) + dE * (T - dT)
        else:
            Eunit = 'J'
            Escale = k * T
            if eV:
                Escale /= q
                Eunit = 'eV'
            Escale = self.__to_numeric(Escale)
            if debug: print 'kT = %2.2g %s' % (Escale * 1000, 'm' + Eunit)
            T1, Y = sym.symbols('T Y')
            # x, y, z = sym.symbols('x y z')
            Total_charge = self.__to_numeric(self.Total_charge.subs([('z', z), (T1, T)]))
            if debug: print 'Total coefficient:', Total_charge
            coeffs = {p: Total_charge.collect(Y).coeff(Y ** p) for p in range(1, 5)}
            coeffs = dict((k, v) for k, v in coeffs.iteritems() if v)
            if debug: print coeffs
            poly_tmp = 0
            for i in coeffs.keys():
                poly_tmp += coeffs[i] * Y ** i
                if debug: print poly_tmp
            coeffs[0] = Total_charge - poly_tmp
            # for i in coeffs.keys():
            #    if abs(coeffs[i]) < 1e-150:
            #        coeffs[i] = 0
            # coeffs = dict((k, v) for k, v in coeffs.iteritems() if v)
            if debug: print coeffs
            # coeffs = {}
            # coeffs_sym = sym.poly(Total_charge, Y).all_coeffs()
            # for i in range(len(coeffs_sym)):
            #    coeffs[i] = coeffs_sym[len(coeffs_sym)-1-i]
            if debug: print 'Polynomial equation, order:', len(coeffs) - 1
            for i in coeffs.keys():
                if debug: print '  Y**%d' % i, '*', coeffs[i]
            if len(coeffs) < 4:
                if len(coeffs) == 2:
                    if debug: print 'We have linear equation with analytical solution'
                    sol = np.array([-coeffs[1] / coeffs[0]])
                elif len(coeffs) == 3:
                    if debug: print 'We have quadratic equation with analytical solution'
                    Discr = coeffs[1] ** 2 - 4 * coeffs[0] * coeffs[2]
                    sol = np.array([(-coeffs[1] - sym.sqrt(Discr)) / (2 * coeffs[0]),
                                    (-coeffs[1] + sym.sqrt(Discr)) / (2 * coeffs[0])])
                else:
                    raise Exception('We have no single solution please check everything')
                if debug: print 'Solution:', sol
                for s in sol:
                    s = self.__to_numeric(s.subs(T1, T))
                    if s > 0:
                        Ef = sym.log(s)
                if debug: print 'Ef in kT units =', Ef
                if debug: print 'Ef = %2.4f %s' % (Ef * Escale, Eunit)
            if debug: print "no analytic solution"
            for i in coeffs.keys():
                coeffs[i] = self.__to_numeric(coeffs[i].subs(T1, T)).evalf()
            try:
                if debug: print 'Trying to solve polynomial equation numerically'
                if debug: print 'Coeffs:', coeffs
                p1 = mp.mpf(1.0)
                EgT = self.band_gap(T, symbolic=False, electron_volts=eV)
                if debug: print 'Eg:', EgT / Escale
                p2 = mp.exp(EgT / Escale)
                sol = solve_polynomial_eq(coeffs, (p1, p2), debug=debug)
            except Exception as ex:
                if debug: print 'Exception:', ex
                pass
            if debug: print 'Numerical solution', sol
            if isinstance(sol, np.ndarray):
                print sol
                if sol.size > 0:
                    sol = sol[sol > 0][0]
            if debug: print 'Numerical solution', sol
            Ef = mp.log(sol)
            if debug: print 'Ef in kT units =', Ef
            if debug: print 'Ef = %2.4f %s' % (Ef * Escale, Eunit)
            self.EchPot_memo[(T, z, eV)] = Ef * Escale
            return Ef * Escale

    def WF(self, T, z=sym.symbols('z'), eV=False):
        if (T, z, eV) in self.WF_memo:
            return self.WF_memo[(T, z, eV)]
        if eV:
            self.WF_memo[(T, z, eV)] = self.Ech_pot(T, z=z, eV=eV, debug=False) + self.__to_numeric(
                self.reference['affinity'] / q)
        else:
            self.WF_memo[(T, z, eV)] = self.Ech_pot(T, z=z, eV=eV, debug=False) + self.__to_numeric(
                self.reference['affinity'])
        return self.WF_memo[(T, z, eV)]

    def __neutr_eq(self, eV=False):
        z = sym.symbols('z')
        E, T1, mu, Eg, Ef, Y, Nc, Nv = sym.symbols('E T mu Eg Ef Y Nc Nv')
        Escale = k * T1
        if eV:
            Escale /= q
        CB_charge = -Nc * sym.exp(-mu / Escale)
        VB_charge = Nv * sym.exp((mu - Eg) / Escale)
        # print '-----'
        Total_charge = CB_charge + VB_charge
        for dopant in self.dopants:
            dop_e = dopant.charge_states[0][1]
            if eV:
                F = Fermi(E, Ef, T1).subs([(E, mu), (Ef, dop_e / q), (T1, T1 / q)])
            else:
                F = Fermi(E, Ef, T1).subs([(E, mu), (Ef, dop_e)])
            Charge = dopant.concentration(z) * (
            (dopant.charge_states[1][0] - dopant.charge_states[0][0]) * F + dopant.charge_states[0][0])
            # print Charge
            Total_charge += Charge
        Total_charge = Total_charge.expand().subs(sym.exp(mu / Escale), Y)  # .simplify()
        factor = to_pos_pwr(Total_charge, debug=False)
        Total_charge = symplify_factor(Total_charge, factor).collect(Y)
        Total_charge = Total_charge.subs(Eg, self.band_gap(T1, symbolic=True, electron_volts=eV))
        Total_charge = Total_charge.subs(Nc, self.Nc(T1))
        Total_charge = Total_charge.subs(Nv, self.Nv(T1))
        # print Total_charge
        return Total_charge

    def mobility(self, z=1000, E=0, T=300, pn=None):
        if pn is None:
            Eg = self.band_gap(T, symbolic=False, electron_volts=False)
            # print Eg, self.__to_numeric(-Eg/(k*T)), mp.exp(self.__to_numeric(-Eg/(k*T)))
            pn = self.Nc(T, symbolic=False) * self.Nv(T, symbolic=False) * mp.exp(
                self.__to_numeric(-Eg / (k * T))) * 1e-12
            # print pn
        N = 0
        for dopant in self.dopants:
            N += dopant.concentration(z)
        N *= 1e-6
        # print N
        mobility = {'mobility_e': {'mu_L': 0, 'mu_I': 0, 'mu_ccs': 0, 'mu_tot': 0},
                    'mobility_h': {'mu_L': 0, 'mu_I': 0, 'mu_ccs': 0, 'mu_tot': 0}}
        for key in mobility.keys():
            mu_L = self.reference[key]['mu_L0'] * (T / 300.0) ** (-self.reference[key]['alpha'])
            mu_I = (self.reference[key]['A'] * (T ** (3 / 2)) / N) / (
            mp.log(1 + self.reference[key]['B'] * (T ** 2) / N) - self.reference[key]['B'] * (T ** 2) / (
            self.reference[key]['B'] * (T ** 2) + N))
            try:
                mu_ccs = (2e17 * (T ** (3 / 2)) / mp.sqrt(pn)) / (mp.log(1 + 8.28e8 * (T ** 2) * (pn ** (-1 / 3))))
                X = mp.sqrt(6 * mu_L * (mu_I + mu_ccs) / (mu_I * mu_ccs))
            except:
                mu_ccs = np.nan
                X = 0
            # print X
            mu_tot = mu_L * (1.025 / (1 + ((X / 1.68) ** (1.43))) - 0.025)
            Field_coeff = (1 + (mu_tot * E * 1e-2 / self.reference[key]['v_s']) ** self.reference[key]['beta']) ** (
            -1 / self.reference[key]['beta'])
            mobility[key]['mu_L'] = mu_L * 1e-4
            mobility[key]['mu_I'] = mu_I * 1e-4
            mobility[key]['mu_ccs'] = mu_ccs * 1e-4
            mobility[key]['mu_tot'] = mu_tot * 1e-4 * Field_coeff
        return mobility

    def __to_numeric(self, expr):
        output = expr
        for key in self.reference.keys():
            output = output.subs(key, self.reference[key])
        for key in constants.keys():
            output = output.subs(key, constants[key])
        return output

    def __str__(self, *args, **kwargs):
        return 'Semiconductor-old: ' + self.reference['name'] + '\n type: ' + self.dop_type + ' (%2.2g cm^-3)' % (
        self.dop_level / 1e6)

    def plot_bulk_diagram(self, T_start=0, T_stop=700, T_step=25, parallel=True):
        T = np.arange(T_start, T_stop + T_step, T_step)
        ax = plt.subplot(111)
        ax.plot(T, self.band_gap(T, symbolic=False, electron_volts=True),
                linewidth=2, color='black', linestyle='-')
        for dopant in self.dopants:
            ax.plot(T, self.band_gap(T, symbolic=False, electron_volts=True) - dopant.energy_level(T, self, 0, electron_volts=True),
                    linewidth=1.5, color='black', linestyle='--')
        # ax.plot(T, self.band_gap(T, symbolic=False, eV=True) - self.Ech_pot1(T, symbolic=False, eV=True),
        #        linewidth=1, color='green', linestyle=':')
        if parallel:
            num_cores = multiprocessing.cpu_count()
            EchPot = np.array(Parallel(n_jobs=num_cores)(delayed(EchPotEc)(self, i) for i in T))
            ax.plot(T, EchPot, linewidth=2, color='red', linestyle=':')
        else:
            ax.plot(T, self.band_gap(T, symbolic=False, electron_volts=True) - self.Ech_pot(T, eV=True),
                    linewidth=2, color='red', linestyle=':')
        ax.set_xlabel("T, K")
        ax.set_ylabel("E, eV")
        return ax
