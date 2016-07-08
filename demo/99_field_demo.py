from __future__ import division, print_function
import numpy as np
from mayavi import mlab

from Schottky.Potential import UniformField, ChargedSphere, ChargedCylinder, HyperbolicCylinder
from Space.Coordinates import transforms as gt
from Space.Field import SuperposedField
import Space_visualization as Visual

radius = 0.2
# pos_charged_sphere_field = ChargedSphere(q=1, r=1)
pos_charged_cylinder_field = ChargedCylinder(l=-3e7, r=radius)
deformation_cylinder_field = HyperbolicCylinder(a=-0.1 / 1.6e-19, r=radius)
external_field = UniformField(strength=0.05, direction=[0, 1, 0])


fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure

superposed_field = SuperposedField('Superposed Field', [
    pos_charged_cylinder_field,
    deformation_cylinder_field,
    external_field
])

'''
superposed_field_vis = Visual.FieldView(fig, superposed_field,
                                        scalar_field_visible=True,
                                        vector_field_visible=True)
grid = np.mgrid[-10:10:20j, -10:10:20j, -5:5:10j]
superposed_field_vis.set_grid(grid)
superposed_field_vis.set_cs_visible(False)
superposed_field_vis.draw()
'''

num_r = 500
num_phi = 720
r = np.linspace(0, 10, num=num_r, endpoint=True)
phi = np.linspace(0, 2 * np.pi, num=num_phi, endpoint=True)
r_grid, p_grid, z_grid = np.meshgrid(r, phi, np.array([0]))
positions = np.vstack([r_grid.ravel(), p_grid.ravel(), z_grid.ravel()]).T
xyz = gt.cylindrical_to_cartesian(positions)

scalar_field = superposed_field.scalar_field(xyz).reshape((num_phi, num_r, 1))

print(scalar_field[:, :, 0].shape)
r_max_arg = np.argmax(scalar_field[:, :, 0], axis=1)
r_max = r_grid[0, r_max_arg, 0]
field_max = np.diagonal(scalar_field[:, r_max_arg, 0])

mlab.mesh(
    # r_grid[:, :, 0],
    # p_grid[:, :, 0],
    r_grid[:, :, 0] * np.cos(p_grid[:, :, 0]),
    r_grid[:, :, 0] * np.sin(p_grid[:, :, 0]),
    scalar_field[:, :, 0] * 50,
    # np.gradient(scalar_field[:, :, 0])[1] * 50,
    # np.gradient(scalar_field[:, :, 0])[1],
    colormap='RdBu',
    # representation='wireframe'
)


mlab.plot3d(
    # r_grid[0, r_max_arg, 0],
    # p_grid[:, 0, 0],
    r_max * np.cos(p_grid[:, 0, 0]),
    r_max * np.sin(p_grid[:, 0, 0]),
    field_max * 50,
    tube_radius=50 * 5e-3,
    color=(1, 0, 0))


mlab.show()
