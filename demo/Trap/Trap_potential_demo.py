import numpy as np
from Schottky.Potential import ExternalField, PointChargeCoulombPotential, PointLikeInExternalField
from Schottky.Potential import SphericallySymmetricInExternalField, HyperbolicInExternalField
from Schottky.Visual.Potential import draw_1d_profile_superposition_scalar_polar, draw_r_min_polar
from Schottky.Visual.Potential import draw_energy_lowering_polar_phi, draw_energy_lowering_polar_theta

from Schottky import constant
from Schottky.Trap import Trap

from matplotlib import pyplot as plt

from time import time


t = Trap('My trap', True,
         0.3 * constant.q, 0.8 * constant.q,
         1e-15, 1e-15)

ext_field_direction = np.array([0.0, 0.0, 1.0])
ext_field_magnitude = 1.0e7
ext_field = ExternalField('External field', ext_field_direction, ext_field_magnitude)

r = 0.5e-9
q = 1.6e-19
epsilon = 11.0
point_charge = PointChargeCoulombPotential('Point charge', q, r, epsilon)

t_pot = PointLikeInExternalField('Coulomb potential', t, point_charge, ext_field,
                                 r_min=1.0e-10, r_max=1.0e-7,
                                 phi_resolution=10, theta_resolution=50)

t_pot_sym = SphericallySymmetricInExternalField('Coulomb potential', t, point_charge, ext_field,
                                                r_min=1.0e-10, r_max=1.0e-7,
                                                phi_resolution=10, theta_resolution=50)

t_pot_hyp = HyperbolicInExternalField('Coulomb potential', t, point_charge, ext_field,
                                      r_min=1.0e-10, r_max=1.0e-7,
                                      phi_resolution=10, theta_resolution=50)

t0 = time()
emission_rate_enhancement = t_pot.emission_rate_enhancement()
print('E3/E0 gen = %2.2f, elapsed %2.6f s' % (emission_rate_enhancement, time() - t0))
t0 = time()
emission_rate_enhancement = t_pot_sym.emission_rate_enhancement()
print('E3/E0 sym = %2.2f, elapsed %2.6f s' % (emission_rate_enhancement, time() - t0))
t0 = time()
emission_rate_enhancement = t_pot_hyp.emission_rate_enhancement()
print('E3/E0 hyp = %2.2f, elapsed %2.6f s' % (emission_rate_enhancement, time() - t0))


theta = np.pi * 0.5
phi = 0.0

r_min = t_pot.max_energy_r_point(theta, phi)
delta_e = t_pot.energy_lowering_point(theta, phi)

ax = draw_1d_profile_superposition_scalar_polar(t_pot, theta, phi, 0, 20, num=1000, scale=1.0e-9,
                                                bw=False, draw_sum=True)

ax.plot([r_min * 1e9], [-delta_e * 1000], 'ro')
ax.legend()
ax.grid(True)

plt.show()

ax = draw_r_min_polar(t_pot, theta, scale=1.0e-9)
ax.legend()
plt.show()

ax = draw_energy_lowering_polar_phi(t_pot, theta, meV=True)
ax.legend()
plt.show()

ax = draw_energy_lowering_polar_theta(t_pot, phi, meV=True)
ax.legend()
plt.show()
