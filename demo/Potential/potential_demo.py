from __future__ import division, print_function
import numpy as np
from matplotlib import pyplot as plt

from Potential.Potential_1D import CoulombPotential, LinearPotential, SuperposedPotential, ChargedCylinderPotential, \
    DislocationDeformationPotential, ConstantPotential

cp = CoulombPotential('Coulomb')
dp = DislocationDeformationPotential('Deformation', 5, 4e-10 * 1.0)
print(dp.a * 100)
ef = LinearPotential('External Field', -1e6 * 1)
cc = ChargedCylinderPotential('Charged Dislocation', charge_sign=-1, linear_charge_density=1e7, radius=1e-9)
zp = ConstantPotential('Zero', 0)

sp = SuperposedPotential('Superposed', [dp, cc, ef])
#sp = SuperposedPotential('Superposed', [zp, ef])
print(sp.barrier_lowering())
_, ax = plt.subplots()

#ef.external_field = 1e6
sp.get_potential_by_name('External Field').external_field = 1e6
#sp.get_potential_by_name('Charged Dislocation').set_linear_charge_density(0.5e7)
print(sp.barrier_lowering())
x = np.linspace(-5e-7, 5e-7, num=500)

cc.plot_potential(ax=ax, x=x)
dp.plot_potential(ax=ax, x=x)
ef.plot_potential(ax=ax, x=x)
sp.plot_potential(ax=ax)

#sp.plot_field(ax=ax)
ax.legend()
ax.grid(True)
plt.show()

'''

lin_dens = []
lowering = []
radius = []
for l in np.linspace(0.0, 10.0, num=100, endpoint=True):
    lin_dens.append(l * 1.0e7)
    sp.get_potential_by_name('Charged Dislocation').set_linear_charge_density(lin_dens[-1])
    poole_frenkel = sp.barrier_lowering()
    print poole_frenkel
    lowering.append(poole_frenkel[0][0])
    radius.append(abs(poole_frenkel[1][0]))

_, (ax1, ax2) = plt.subplots(2)
ax1.plot(np.array(lin_dens) * 1e-2, lowering, 'b-o')
ax2.plot(np.array(lin_dens) * 1e-2, np.array(radius) * 1e9, 'b-o')
plt.show()

'''