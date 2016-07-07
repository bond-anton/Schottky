from __future__ import division, print_function
import numpy as np
from mayavi import mlab

from Schottky.Potential import UniformField, ChargedSphere, ChargedCylinder
from Space.Field import SuperposedField
import Space_visualization as Visual


pos_charged_sphere_field = ChargedSphere(q=1, r=1)
pos_charged_cylinder_field = ChargedCylinder(l=1, r=1)
external_field = UniformField(strength=0.25, direction=[0, 1, 0])


fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure

grid = np.mgrid[-10:10:20j, -10:10:20j, -5:5:10j]

superposed_field = SuperposedField('Superposed Field', [
    pos_charged_cylinder_field,
    #pos_charged_sphere_field,
    #external_field
])

superposed_field_vis = Visual.FieldView(fig, superposed_field,
                                        scalar_field_visible=True,
                                        vector_field_visible=True)
superposed_field_vis.set_grid(grid)
superposed_field_vis.set_cs_visible(False)
superposed_field_vis.draw()

mlab.show()
