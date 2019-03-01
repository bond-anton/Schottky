import numpy as np
from Schottky.Potential import PointChargeCoulombPotential
from Schottky.Visual.Potential import draw_1d_profile_scalar

import BDSpaceVis as Visual

from mayavi import mlab
from matplotlib import pyplot as plt


r = 3e-9
q = 1.6e-19
epsilon = 1.0

my_field = PointChargeCoulombPotential('External field', q, r, epsilon)

my_field.epsilon = -11.0  # Should give a warning

my_field.epsilon = 11.0
my_field.r = 5e-9
my_field.charge *= -0.1

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure


field_vis = Visual.FieldView(fig, my_field, scalar_field_visible=True)

x_start = -1e-7
x_stop = 1e-7
y_start = -1e-7
y_stop = 1e-7
z_start = -5e-8
z_stop = 5e-8

grid = np.mgrid[x_start:x_stop:40j, y_start:y_stop:40j, z_start:z_stop:10j]

field_vis.set_grid(grid)
field_vis.set_cs_visible(False)
field_vis.draw()
field_vis.set_scale_factor(my_field.r)

fig.scene.isometric_view()
mlab.show()

ax = draw_1d_profile_scalar(my_field, np.array([0, 0, 1.0]), -200, 200, num=1000, scale=1.0e-9)
plt.show()