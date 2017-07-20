# -*- coding: utf-8 -*-

'''
Created on 04 нояб. 2014 г.

@author: anton
'''
from __future__ import division, print_function

import sympy as sym
import numpy as np
import math as m

import mpmath as mp
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from Schottky import constants
from Schottky.Notation import k

def fermi(energy, fermi_energy, temperature, g_ratio):
    """
    fermi distribution function
    :param energy: Energy level
    :param fermi_energy: Fermi level enrgy
    :param temperature: Temperature in K
    :param g_ratio: Degeneracy ratio
    :return: sympy symbolic expression
    """
    return 1 / (1 + g_ratio * sym.exp((energy - fermi_energy) / (k * temperature)))


def d_fermi_d_delta_fermi_energy(energy, fermi_energy, temperature, g_ratio):
    """
    Derivative of fermi distribution function by delta_fermi_energy at delta_fermi_energy = 0
    :param energy: Energy level
    :param fermi_energy: Fermi level enrgy
    :param temperature: Temperature in K
    :param g_ratio: Degeneracy ratio
    :return: sympy symbolic expression
    """
    numerator = g_ratio * sym.exp((energy - fermi_energy) / (k * temperature)) / (k * temperature)
    denominator = (1 + g_ratio * sym.exp((energy - fermi_energy) / (k * temperature))) ** 2
    return numerator / denominator


def fermihalf(x, sgn):
    """ Series approximation to the F_{1/2}(x) or F_{-1/2}(x) 
        Fermi-Dirac integral """

    f = lambda k: mp.sqrt(x ** 2 + np.pi ** 2 * (2 * k - 1) ** 2)

    # if x < -100:
    #    return 0.0
    if x < -9 or True:
        if sgn > 0:
            return mp.exp(x)
        else:
            return mp.exp(x)

    if sgn > 0:  # F_{1/2}(x)
        a = np.array((1.0 / 770751818298, -1.0 / 3574503105, -13.0 / 184757992,
                      85.0 / 3603084, 3923.0 / 220484, 74141.0 / 8289, -5990294.0 / 7995))
        g = lambda k: mp.sqrt(f(k) - x)

    else:  # F_{-1/2}(x)
        a = np.array((-1.0 / 128458636383, -1.0 / 714900621, -1.0 / 3553038,
                      27.0 / 381503, 3923.0 / 110242, 8220.0 / 919))
        g = lambda k: -0.5 * mp.sqrt(f(k) - x) / f(k)

    F = np.polyval(a, x) + 2 * np.sqrt(2 * np.pi) * sum(map(g, range(1, 21)))
    return F

def centered_linspace(center, radius, step):
    start = center - radius
    stop = center + radius
    points_number = np.ceil(abs(stop - start) / step)
    return np.linspace(start, stop, num=points_number, endpoint=True)

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def synchronuos_detection(f, U0, Ym):
    T = 1 / f  # s
    w = 2 * np.pi * f
    t = np.linspace(0, T, Ym.size)
    y1 = U0 * np.sin(w * t)
    y2 = U0 * np.cos(w * t)
    Is = 2 / (U0 * T) * np.trapz(y1 * Ym, t)
    Ic = 2 / (U0 * T) * np.trapz(y2 * Ym, t)
    phi_m = np.arctan2(Ic, Is)
    lag = phi_m / w
    U_m = np.sqrt(Is ** 2 + Ic ** 2)
    return U_m, phi_m, lag


def exp_fit(x, y):
    def func(x, a, b):
        return a * np.exp(-b * x)

    popt, _ = curve_fit(func, x, y)
    print popt

    def fitted(x):
        return func(x, *popt)

    return fitted


def to_numeric(expr):
    if isinstance(expr, np.ndarray):
        return np.array(map(to_numeric, expr))
    output = expr
    for key in constants.keys():
        try:
            output = output.subs(key, constants[key])
        except AttributeError:
            pass
    return np.float(output)


def to_pos_pwr(expr, debug=False):
    factor = 1
    if debug: print expr, expr.func
    if expr.func == sym.Pow:
        if expr.args[1] < 0:
            factor *= sym.Pow(expr.args[0], -expr.args[1])
            if debug: print 'Factor', factor
    elif expr.func == sym.Add or expr.func == sym.Mul:
        if debug: print 'Function', expr.func
        for arg in expr.args:
            factor *= to_pos_pwr(arg, debug=debug)
            if debug: print 'Factor', factor
    return factor


def symplify_factor(expr, factor):
    output = 0
    if expr.func == sym.Add:
        for arg in expr.args:
            output += symplify_factor(arg, factor)
            # print output
    else:
        output += (expr * factor).expand()  # .simplify().expand()
        # print output
    return output


def to_positive_power(expr, debug=False):
    fc = 1
    if expr.func == sym.Add:
        if debug: print 'Summation'
        for arg in expr.args:
            if debug: print '--> Argument:', arg
            if arg.func == sym.Mul:
                if debug: print '---> Multiplication'
                for arg1 in arg.args:
                    if debug: print arg1, arg1.func
                    if debug: print arg1.args
                    if arg1.func == sym.Pow:
                        if arg1.args[1] < 0:
                            fc *= sym.Pow(arg1.args[0], -arg1.args[1])
    elif expr.func == sym.Mul:
        for arg1 in expr.args:
            if debug: print arg1.func
            if debug: print arg1.args
            if arg1.func == sym.Pow:
                if arg1.args[1] < 0:
                    fc *= sym.Pow(arg1.args[0], -arg1.args[1])
    if debug: print fc
    return fc


def build_mpf_fn(x, coeffs, debug=False):
    if debug: print 'Polynomial equation, order:', len(coeffs) - 1
    # coeffs = [mp.mpf(c) for c in coeffs]
    eq = 0
    for i in coeffs.keys():
        if debug: print mp.mpf(x) ** i * mp.mpf(coeffs[i])
        eq += mp.mpf(coeffs[i]) * (mp.mpf(x) ** i)
    return mp.mpf(eq)


def my_secant(eq, p1, p2, debug=False):
    tol = mp.mpf(0.1)
    max_count = 10000
    sol = 0
    for count in range(max_count):
        if debug: print (count + 1), ''
        if p1 == p2:
            sol = p1
            break
        y1 = eq(p1)
        y2 = eq(p2)
        if debug: print '-->', p1, '->', y1
        if debug: print '-->', p2, '->', y2
        if abs(y1) < abs(y2):
            sol = p1
            err = abs(y1)
        else:
            sol = p2
            err = abs(y2)
        if err < tol:
            break
        if mp.sign(y1) * mp.sign(y2) < 0:
            # p3 = (p1+p2)/mpf(2)
            x1 = mp.log(p1)
            x2 = mp.log(p2)
            # x1 = p1
            # x2 = p2
            # x3 = (x2*y1 - x1*y2)/(y1-y2)
            # if x3 == x1 or x3 == x2:
            x3 = (x1 + x2) / mp.mpf(2)
            p3 = mp.exp(x3)
            # p3 = x3
            if p3 == p1 or p3 == p2:
                break
            y3 = eq(p3)
            if debug: print '--->', x1, x2, x3, p3, '->', y3
            if mp.sign(y3) == mp.sign(y1):
                p1 = p3
            else:
                p2 = p3
        elif mp.sign(y1) * mp.sign(y2) == 0:
            if y1 == 0:
                sol = p1
            elif y2 == 0:
                sol = p2
            else:
                raise Exception('Strange: sign returns zero! without zeros $)')
            break
        else:
            raise Exception('Functin has same sign on both ends')
    if debug: print 'Solution:', sol
    return sol


def solve_polynomial_eq(coeffs, (p1, p2), debug=False):
    if debug: print '*******'
    eq = lambda x: build_mpf_fn(x, coeffs, debug=False)
    # xxx = np.linspace(np.float(p1), np.float(p2), num=100)
    # yyy = np.array([np.float(eq(x)) for x in xxx])
    # plt.plot(xxx, yyy)
    # plt.show()
    if len(coeffs) < 4:
        if len(coeffs) == 2:
            print 'We have linear equation with analytical solution'
            sol = np.array([-coeffs[1] / coeffs[0]])
        elif len(coeffs) == 3:
            print 'We have quadratic equation with analytical solution'
            Discr = coeffs[1] ** 2 - 4 * coeffs[0] * coeffs[2]
            sol = np.array(
                [(-coeffs[1] - sym.sqrt(Discr)) / (2 * coeffs[0]), (-coeffs[1] + sym.sqrt(Discr)) / (2 * coeffs[0])])
        else:
            raise Exception('We have no single solution please check everything')
    else:
        # p1 = mpf('1e550')
        # p2 = mpf('1e590')
        sol = my_secant(eq, p1, p2, debug=debug)
        # sol = findroot(eq, (p1, p2), solver='secant', verify=False)
        # sol = findroot(eq, (p1, p2), solver=Bisection, verify=False)
        if debug: print sol
    if debug: print '*******'
    return sol


def smooth_dd(x, eps=None):
    if eps is None:
        eps = 1.0e-8
    if isinstance(x, (float, int)):
        return 1 / (2 * eps) * (1 + np.cos(np.pi * x / eps)) if (x < eps and x > -eps) else 0
    elif isinstance(x, np.ndarray):
        return np.array(map(smooth_dd, x))
    else:
        return 1 / (2 * eps) * (1 + sym.cos(np.pi * x / eps)) * sym.Heaviside(x + eps) * sym.Heaviside(-x + eps)
        # raise TypeError('%s x is wrong' % type(x))


def Psi_approx(L, bc1, bc2):
    def Psi(z):
        # n = 1.5
        # if isinstance(z, np.ndarray):

        # print type(L), type(bc1), type(bc2), type(z)
        # psi = (bc2-bc1)/(L/n) * z + bc1
        # psi[np.where(z > L/n)] = 0
        # return psi
        return bc1 * np.exp(-20 * z / 1e-6)

    return Psi


def Psi_zero(z):
    if isinstance(z, (mp.mpf, float, int)):
        return 0.0
    elif isinstance(z, np.ndarray):
        return np.zeros(z.size)
    else:
        raise Exception('Wrong argument type z:%s' % str(type(z)))


def Psi_E_zero(z):
    if isinstance(z, (mp.mpf, float, int)):
        return z, 0, 0
    elif isinstance(z, np.ndarray):
        return z, np.zeros(z.shape), np.zeros(z.shape)
    else:
        raise Exception('Wrong argument type z:%s' % str(type(z)))


def interp_Fn(Z, F, interp_type='linear'):
    '''
    z and F must be 1D arrays of equal size
    interp_type could be one of
    'linear'
    'last'
    'zero'
    '''
    # print 'type:', interp_type
    def interp(z):
        interpolator = interp1d(Z, F, bounds_error=True)
        xs = interpolator.x
        ys = interpolator.y

        def pointwise(x):
            if x < xs[0]:
                if interp_type == 'linear':
                    return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
                elif interp_type == 'last':
                    # print 'here'
                    return ys[0]
                elif interp_type == 'zero':
                    return 0.0
            elif x > xs[-1]:
                if interp_type == 'linear':
                    return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
                elif interp_type == 'last':
                    # print 'here-here'
                    return ys[-1]
                elif interp_type == 'zero':
                    return 0.0
            else:
                # sss = interpolator(x)
                # print 'I am here', type(np.float(sss)), sss.size, sss.shape
                return np.float(interpolator(x))

        if isinstance(z, (np.ndarray, list, tuple)):
            return np.array(map(pointwise, z), dtype=np.float)
        else:
            return pointwise(z)

    return interp


def check_if_integer(x, threshold):
    L = m.floor(x)
    U = m.ceil(x)
    C = L if abs(L - x) < abs(U - x) else U
    # print x, L, U, C, C - x
    if abs(C - x) < threshold:
        return True
    else:
        return False


def gen_interpolated_P_E_function(z_nodes, Psi, Field):
    # P = interp1d(z_nodes, Psi, bounds_error=False, fill_value=0)
    # E = interp1d(z_nodes, Field, bounds_error=False, fill_value=0)
    P = interp_Fn(z_nodes, Psi, interp_type='last')
    E = interp_Fn(z_nodes, Field, interp_type='last')

    def interpolated_f(z):
        if not isinstance(z, np.ndarray):
            return z, P(z), E(z)
        elif isinstance(z, np.ndarray):
            Z = np.union1d(z, z_nodes)
            Pv = P(Z)
            Ev = E(Z)
            return Z, Pv, Ev
        else:
            print z, 'Wrong argument type z:%s' % str(type(z))
            raise

    return interpolated_f
