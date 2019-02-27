import numpy as np
from Schottky.Potential import ExternalField

import BDSpaceVis as Visual

from mayavi import mlab


direction = np.array([0.0, 0.0, 1.0])
magnitude = 1.0e6

my_field = ExternalField('External field', direction, magnitude)
my_field.magnitude *= 2
my_field.direction = np.array([1.0, 1.0, 1.0])

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure


field_vis = Visual.FieldView(fig, my_field, scalar_field_visible=False)

grid = np.mgrid[-10:10:40j, -10:10:40j, -5:5:10j]

field_vis.set_grid(grid)
field_vis.set_cs_visible(False)
field_vis.draw()

mlab.show()
