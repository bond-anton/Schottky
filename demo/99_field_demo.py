from __future__ import division, print_function

import numpy as np
from Space.Coordinates import transforms as gt
from Space.Field import SuperposedField
from mayavi import mlab

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, HyperbolicCylinder


def find_first_local_extremum(surface, r_grid, radius):
    surface_gradient = np.gradient(surface)[1]
    start = min(np.where(r_grid[0, :, 0] > radius)[0])
    step = 10
    stop = start + step
    while stop <= r_grid[0, :, 0].size:
        min_arg = np.argmin(abs(surface_gradient[:, start:stop]), axis=1)
        if (stop - start - min_arg > 3).all() and (min_arg > 0).all():
            return min_arg + start
        stop += step
    return np.argmax(surface, axis=1)

radius = 0.5
# pos_charged_sphere_field = ChargedSphere(coefficient=1, radius=1)
charged_cylinder_field = ChargedCylinder(l=3e7, r=radius)
deformation_cylinder_field = HyperbolicCylinder(a=-0.1 / 1.6e-19, r=radius)
external_field = UniformElectrostaticField(strength=0.1, direction=[0, 1, 0])


fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure

superposed_field = SuperposedField('Superposed Field', [
    charged_cylinder_field,
    deformation_cylinder_field,
    #external_field
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

r_max_arg = find_first_local_extremum(scalar_field[:, :, 0], r_grid, radius)

r_max = r_grid[0, r_max_arg, 0]
field_max = np.diagonal(scalar_field[:, r_max_arg, 0])

energy_scale = 50

'''
mlab.mesh(
    # r_grid[:, :, 0],
    # p_grid[:, :, 0],
    r_grid[:, :, 0] * np.cos(p_grid[:, :, 0]),
    r_grid[:, :, 0] * np.sin(p_grid[:, :, 0]),
    np.zeros_like(scalar_field[:, :, 0]),
    #external_sacalar_field[:, :, 0] * energy_scale,
    # np.gradient(scalar_field[:, :, 0])[1] * energy_scale,
    # np.gradient(scalar_field[:, :, 0])[1],
    colormap='RdBu',
    representation='wireframe'
)
'''
mlab.mesh(
    # r_grid[:, :, 0],
    # p_grid[:, :, 0],
    r_grid[:, :, 0] * np.cos(p_grid[:, :, 0]),
    r_grid[:, :, 0] * np.sin(p_grid[:, :, 0]),
    scalar_field[:, :, 0] * energy_scale,
    #abs(scalar_field_gradient) * energy_scale,
    colormap='RdBu',
    # representation='wireframe'
)


mlab.plot3d(
    # r_grid[0, r_max_arg, 0],
    # p_grid[:, 0, 0],
    r_max * np.cos(p_grid[:, 0, 0]),
    r_max * np.sin(p_grid[:, 0, 0]),
    field_max * energy_scale,
    tube_radius=energy_scale * 2e-3,
    color=(1, 0, 0))

mlab.show()
