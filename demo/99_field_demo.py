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
charged_cylinder_field = ChargedCylinder(client=client, name='Charged cylinder field 2',
                                         charge_density=3e7, radius=radius, epsilon=1)
deformation_cylinder_field = HyperbolicCylinder(client=client, name='Deformation potential 2',
                                                coefficient=-0.1, radius=radius)
external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field 2',
                                           strength=0.01, direction=[1, 0, 0])


superposed_field = SuperpositionField(client=client, name='Superposed Field 2',
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

client.user_manager.sign_out()

mlab.show()
