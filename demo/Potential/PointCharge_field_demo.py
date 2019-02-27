import numpy as np
from Schottky.Potential import PointChargeCoulombPotential

import BDSpaceVis as Visual

from mayavi import mlab


r = 1
q = 1.6e-19
epsilon = 1.0

my_field = PointChargeCoulombPotential('External field', q, r, epsilon)

my_field.epsilon = -11.0  # Should give a warning

my_field.epsilon = 11.0
my_field.r = 2
my_field.charge *= -2

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure


field_vis = Visual.FieldView(fig, my_field, scalar_field_visible=False)

grid = np.mgrid[-10:10:40j, -10:10:40j, -5:5:10j]

field_vis.set_grid(grid)
field_vis.set_cs_visible(False)
field_vis.draw()

mlab.show()
