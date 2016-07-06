#__author__ = 'anton'

import numpy as np
from matplotlib import pyplot as plt
try:
    from mayavi import mlab
    USE_MAYAVI = True
except ImportError:
    USE_MAYAVI = False
USE_MAYAVI = False
from scipy.constants import elementary_charge, epsilon_0, Boltzmann

from Potential.Potential_3D import ConstantPotential, ConstantFieldPotential, HyperbolicPotential, CoulombPotential, \
    DislocationDeformationPotential,\
    SuperposedPotential, ChargedCylinderPotential

if USE_MAYAVI:
    figU = mlab.figure('Potential', bgcolor=(0, 0, 0))
    #figF = mlab.figure('Electric field', bgcolor=(0, 0, 0))


cp = CoulombPotential('Coulomb')
dp = DislocationDeformationPotential('Deformation', -5, 4e-10 * 1.0)
#print dp.a * 100
ef = ConstantFieldPotential('External Field', (5e5, 0, 0))
#ef2 = ConstantFieldPotential('External Field', (5e5, 0, 0))
cc = ChargedCylinderPotential('Charged Dislocation', charge_sign=1, linear_charge_density=1*1e8, radius=1e-9, epsilon=11.8)
zp = ConstantPotential('Zero', 0)


sp = SuperposedPotential('Superposed', [dp, cc, ef])
#sp = SuperposedPotential('Superposed', [zp, ef])
#print sp.barrier_lowering()


#ef.external_field = 1e6
sp.get_potential_by_name('External Field').external_field = (5*5e5, np.pi, 0)
#sp.get_potential_by_name('Charged Dislocation').set_linear_charge_density(0.5e7)
#print sp.barrier_lowering()


r = np.linspace(2e-9, 5e-8, num=100, endpoint=True)
theta = np.linspace(0, 2 * np.pi, num=100, endpoint=True)
phi = np.linspace(0, 2 * np.pi, num=100, endpoint=True)

x = np.linspace(-5e-9, 5e-9, num=25)
y = np.linspace(-5e-9, 5e-9, num=25)
z = np.linspace(-5e-9, 5e-9, num=25)


if USE_MAYAVI:
    figU = sp.plot_potential(figU, r, theta, phi, lowering=True)
    #sp.plot_potential_cartesian()
    #mlab.outline()
    #sp.plot_field(figF, x, y, z)
    mlab.show()
#print sp.barrier_lowering(np.pi)

T_range = np.linspace(30, 300, 28, endpoint=True)
poole_frenkel = []
for T in T_range:
    kT = Boltzmann * T / elementary_charge

    theta = np.linspace(0, np.pi, num=100, endpoint=True)
    barrier_lowering = np.array([sp.barrier_lowering(theta_i) for theta_i in theta])
    poole_frenkel.append(0.5 * np.trapz(np.sin(theta) * np.exp(abs(barrier_lowering[:, 0]) / kT), theta))
    print('emission boost @T=%2.2f K: %2.4g' % (T, poole_frenkel[-1]))
plt.semilogy(T_range, np.array(poole_frenkel), 'r-o')
plt.show()