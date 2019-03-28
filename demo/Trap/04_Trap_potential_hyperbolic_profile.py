import numpy as np
from Schottky.Potential import ExternalField, PointChargeCoulombPotential
from Schottky.Potential import HyperbolicInExternalField

from Schottky import constant
from Schottky.Trap import Trap

import pstats, cProfile

ext_field_direction = np.array([0.0, 0.0, 1.0])
ext_field_magnitude = 1.0e7
ext_field = ExternalField('External field', ext_field_direction, ext_field_magnitude)

r = 0.5e-9
q = 1.6e-19
epsilon = 11.0
point_charge = PointChargeCoulombPotential('Point charge', q, r, epsilon)

t_pot = HyperbolicInExternalField('Coulomb potential', point_charge, ext_field,
                                  r_min=1.0e-10, r_max=1.0e-7,
                                  phi_resolution=360, theta_resolution=50)

t = Trap('My trap', True,
         0.3 * constant.q, 0.8 * constant.q,
         1e-15, 1e-15,
         e_potential=t_pot)

cProfile.runctx()

emission_rate_enhancement = t_pot.emission_rate_enhancement()

print('E3/E0 =', emission_rate_enhancement)
