import numpy as np
from Schottky.Potential import ExternalField, PointChargeCoulombPotential
from Schottky.Visual.Potential import draw_1d_profile_scalar, draw_1d_profile_superposition_scalar

from BDSpace.Field import SuperposedField
import BDSpaceVis as Visual

from mayavi import mlab
from matplotlib import pyplot as plt


ext_field_direction = np.array([0.0, 0.0, 1.0])
ext_field_magnitude = 5.0e8
ext_field = ExternalField('External field', ext_field_direction, ext_field_magnitude)

r = 0.5e-9
q = 1.6e-19
epsilon = 11.0
point_charge = PointChargeCoulombPotential('Point charge', q, r, epsilon)

superposed_field = SuperposedField('Point charge in external field', [ext_field, point_charge])

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure


field_vis = Visual.FieldView(fig, superposed_field, scalar_field_visible=True)

x_start = -5e-9
x_stop = 5e-9
y_start = -5e-9
y_stop = 5e-9
z_start = -5e-9
z_stop = 5e-9

grid = np.mgrid[x_start:x_stop:40j, y_start:y_stop:40j, z_start:z_stop:20j]

field_vis.set_grid(grid)
field_vis.set_cs_visible(False)
field_vis.draw()
field_vis.set_scale_factor(point_charge.r)
# field_vis.set_scale_factor(1e-6)

fig.scene.isometric_view()
mlab.show()

ax = draw_1d_profile_superposition_scalar(superposed_field, np.array([0, 0, 1.0]), -20, 20, num=1000, scale=1.0e-9,
                                          bw=False, draw_sum=True)
ax.legend()
plt.show()
