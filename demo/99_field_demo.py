from __future__ import division, print_function

import Space_visualization as Visual
from ScientificProjects.Client import Client

import numpy as np

from mayavi import mlab

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, HyperbolicCylinder, SuperpositionField
from Schottky.Simulators.Field import FieldSimulator


def find_first_local_extrema(surface, r_grid, radius):
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


client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


radius = 0.5
charged_cylinder_field = ChargedCylinder(client=client, name='Charged cylinder field',
                                         charge_density=3e7, radius=radius, epsilon=1)
deformation_cylinder_field = HyperbolicCylinder(client=client, name='Deformation potential',
                                                coefficient=-0.1, radius=radius)
external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                           strength=0.01, direction=[1, 0, 0])


superposed_field = SuperpositionField(client=client, name='Superposed Field',
                                      fields=[charged_cylinder_field,
                                              deformation_cylinder_field,
                                              external_field])

field_simulator = FieldSimulator(client=client, field=superposed_field)

fig = mlab.figure(bgcolor=(0.2, 0.2, 0.2))

superposed_field_vis = Visual.FieldView(fig, field_simulator,
                                        grid=np.mgrid[-5:5:10j, -5:5:10j, -5:5:10j],
                                        scalar_field_visible=True,
                                        vector_field_visible=True)
superposed_field_vis.set_cs_visible(False)
superposed_field_vis.draw()
external_field.rotate_euler_angles([0, 0, np.pi / 4])

'''
x = np.linspace(-5, 5, num=20, endpoint=True)
y = np.linspace(-5, 5, num=20, endpoint=True)
r = np.linspace(0, 5, num=20, endpoint=True)
theta = np.linspace(0, np.pi, num=10, endpoint=True)
phi = np.linspace(0, np.pi, num=20, endpoint=True)
z = np.linspace(0, 5, num=12, endpoint=True)
#z = np.array([0])

r_grid, p_grid, z_grid, scalar_field, vector_field = field_simulator.measure_field_cylindrical_coordinates(
    r_range=r, phi_range=phi, z_range=z, length_unit='nm', force_recalculate=False
)

r_max_arg = find_first_local_extrema(scalar_field[:, :, 0], r_grid, radius)

r_max = r_grid[0, r_max_arg, 0]
field_max = np.diagonal(scalar_field[:, r_max_arg, 0])

energy_scale = 50
'''
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


'''
vector_field = mlab.quiver3d(r_grid * np.cos(p_grid),
                             r_grid * np.sin(p_grid),
                             z_grid,
                             vector_field[:, :, :, 0],
                             vector_field[:, :, :, 1],
                             vector_field[:, :, :, 2],
                             name=superposed_field.name)

x_grid, y_grid, z_grid, scalar_field, vector_field = field_simulator.measure_field_cartesian_coordinates(
    x_range=x, y_range=y, z_range=z, length_unit='nm', force_recalculate=False
)

vector_field = mlab.quiver3d(x_grid,
                             y_grid,
                             z_grid,
                             vector_field[:, :, :, 0],
                             vector_field[:, :, :, 1],
                             vector_field[:, :, :, 2],
                             name=superposed_field.name)

r_grid, theta_grid, phi_grid, scalar_field, vector_field = field_simulator.measure_field_spherical_coordinates(
    r_range=r, theta_range=theta, phi_range=phi, length_unit='nm', force_recalculate=False
)

vector_field = mlab.quiver3d(r_grid * np.sin(theta_grid) * np.cos(phi_grid),
                             r_grid * np.sin(theta_grid) * np.sin(phi_grid),
                             r_grid * np.cos(theta_grid),
                             vector_field[:, :, :, 0],
                             vector_field[:, :, :, 1],
                             vector_field[:, :, :, 2],
                             name=superposed_field.name)


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


client.user_manager.sign_out()


'''

client.user_manager.sign_out()

mlab.show()
