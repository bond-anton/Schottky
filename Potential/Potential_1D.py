__author__ = 'anton'

import warnings
import numpy as np
from scipy.constants import elementary_charge, epsilon_0
from scipy.optimize import root
from matplotlib import pyplot as plt


class GenericPotential(object):

    def __init__(self, label='Unknown Potential'):
        super(GenericPotential, self).__init__()
        self.label = label

    def potential(self, x):
        return np.empty_like(x)

    def field(self, x):
        return np.empty_like(x)

    def plot_potential(self, ax=None, x=None):
        if ax is None:
            _, ax = plt.subplots()

        if x is None:
            x = np.linspace(-20e-9, 20e-9, num=500, endpoint=True)
        ax.plot(x, self.potential(x), label=self.label)
        ax.set_xlim([x[0], x[-1]])
        ax.grid(True)

    def plot_field(self, ax=None, x=None):
        if ax is None:
            _, ax = plt.subplots()
        if x is None:
            x = np.linspace(-20e-9, 20e-9, num=500, endpoint=True)
        ax.plot(x, self.field(x), label=self.label)
        ax.set_xlim([x[0], x[-1]])


class SuperposedPotential(GenericPotential):

    def __init__(self, label='Unknown Superposed Potential', potentials=None):
        super(SuperposedPotential, self).__init__(label)
        if potentials is None:
            potentials = []
        self.potentials = potentials

    def potential(self, x):
        superposition = np.zeros_like(x)
        for potential in self.potentials:
            superposition += potential.potential(x)
        return superposition

    def field(self, x):
        superposition = np.zeros_like(x)
        for potential in self.potentials:
            superposition += potential.field(x)
        return superposition

    def get_potential_by_name(self, name):
        for potential in self.potentials:
            if potential.label == name:
                return potential
        return None

    def barrier_lowering(self):
        guess_r = -1e-9
        warnings.filterwarnings('ignore')
        solution_l = root(self.field, guess_r)
        warnings.resetwarnings()
        if solution_l.success:
            r0_l = solution_l.x[0]
            delta_phi_l = abs(self.potential(r0_l))
            print 'left:', r0_l, delta_phi_l
        else:
            r0_l = 0
            delta_phi_l = 0
        guess_r = 1e-9
        warnings.filterwarnings('ignore')
        solution_r = root(self.field, guess_r)
        warnings.resetwarnings()
        if solution_r.success:
            r0_r = solution_r.x[0]
            delta_phi_r = abs(self.potential(r0_r))
            print 'right:', r0_r, delta_phi_r
        else:
            r0_r = 0
            delta_phi_r = 0
        if r0_l == r0_r:
            r0 = np.array([r0_l])
            delta_phi = np.array([delta_phi_l])
        else:
            r0 = np.array([r0_l, r0_r])
            delta_phi = np.array([delta_phi_l, delta_phi_r])
        return delta_phi, r0

    def plot_potential(self, ax=None, x=None):
        if ax is None:
            _, ax = plt.subplots()
        delta_phi_arr, r0_arr = self.barrier_lowering()
        for r0, delta_phi in zip(r0_arr, delta_phi_arr):
            if r0 == 0:
                if x is None:
                    x = np.linspace(-20e-9, 20e-9, num=500, endpoint=True)
            else:
                if x is None:
                    x = np.linspace(-10 * abs(r0), 10 * abs(r0), num=500, endpoint=True)
                ax.annotate('', xy=(r0, self.potential(r0)), xycoords='data', xytext=(r0, 0),
                            textcoords='data', arrowprops={'arrowstyle': '<->'})
                label_pos = 'bottom'
                ax.annotate('% 2.2g V' % delta_phi, xy=(r0, 0), xycoords='data',
                            xytext=(0, 5), textcoords='offset points',
                            horizontalalignment='center', verticalalignment=label_pos)
                y0 = -2 * delta_phi
                y1 = 0
                ax.set_ylim([y0, y1])
        ax.plot(x, self.potential(x), linewidth=2, color='k', label=self.label)
        ax.set_xlim([x[0], x[-1]])
        ax.grid(True)


class ConstantPotential(GenericPotential):

    def __init__(self, label='Unknown Constant Potential', phi=0):
        super(ConstantPotential, self).__init__(label)
        self.phi = phi

    def potential(self, x):
        return np.ones_like(x) * self.phi

    def field(self, x):
        return np.zeros_like(x)


class LinearPotential(GenericPotential):

    def __init__(self, label='Unknown Linear Potential', external_field=0):
        super(LinearPotential, self).__init__(label)
        self.external_field = external_field

    def potential(self, x):
        return self.external_field * x

    def field(self, x):
        return np.ones_like(x) * self.external_field


class HyperbolicPotential(GenericPotential):

    def __init__(self, label='Unknown Hyperbolic Potential', a=1, attractive=True):
        super(HyperbolicPotential, self).__init__(label)
        self.a = a
        self.attractive = attractive

    def potential(self, x):
        potential_sign = -1 if self.attractive else 1
        return potential_sign * self.a / np.abs(x)

    def field(self, x):
        potential_sign = -1 if self.attractive else 1
        return -potential_sign * self.a / (np.abs(x) * x)


class CoulombPotential(HyperbolicPotential):

    def __init__(self, label='Unknown Coulomb Potential', epsilon=1, charge=elementary_charge, attractive=True):
        self.k = 1 / (4 * np.pi * epsilon_0 * epsilon)
        self.charge = charge
        super(CoulombPotential, self).__init__(label, self.k * self.charge, attractive)


class DislocationDeformationPotential(HyperbolicPotential):

    def __init__(self, label='Deformation Potential', deformation_modulus=0, burgers_vector=1):
        self.deformation_modulus = deformation_modulus
        self.burgers_vector = burgers_vector
        a = deformation_modulus * burgers_vector / (4 * np.pi)
        super(DislocationDeformationPotential, self).__init__(label, a)


class InfiniteCylinderPotential(GenericPotential):

    def __init__(self, label='Unknown Potential', a=1, radius=1):
        super(InfiniteCylinderPotential, self).__init__(label)
        self.radius = radius
        self.a = a

    def set_a(self, a):
        self.a = a

    def potential(self, x):
        return self.a * np.log(abs(x) / self.radius)

    def field(self, x):
        return self.a / x


class ChargedCylinderPotential(InfiniteCylinderPotential):

    def __init__(self, label='Unknown Potential', epsilon=1, charge_sign=1, linear_charge_density=1, radius=1):
        self.epsilon = epsilon
        self.charge_sign = charge_sign
        self.linear_charge_density = linear_charge_density
        a = charge_sign * elementary_charge * linear_charge_density / (2 * np.pi * epsilon_0 * epsilon)
        super(ChargedCylinderPotential, self).__init__(label, a, radius)

    def set_linear_charge_density(self, linear_charge_density):
        self.linear_charge_density = linear_charge_density
        a = self.charge_sign * elementary_charge * self.linear_charge_density / (2 * np.pi * epsilon_0 * self.epsilon)
        self.set_a(a)
