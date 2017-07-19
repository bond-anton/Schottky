from __future__ import division, print_function
import numpy as np
from mayavi import mlab

from ScientificProjects.Client import Client

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, HyperbolicCylinder, SuperpositionField
from Schottky.Simulators.Field import FieldSimulator

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


radius = 0.5
pos_charged_cylinder_field = ChargedCylinder(client=client, name='Pos charged cylinder field 3',
                                             charge_density=3e6, radius=radius, epsilon=1)
pos_charged_cylinder_field.set_origin([0, 0, 0])

one_by_r_field = HyperbolicCylinder(client=client, name='the Trap',
                                    coefficient=-0.1, radius=radius)
one_by_r_field.set_origin([0, 0, 0])


external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                           strength=0.01, direction=[0, 0, 1])


superposed_field = SuperpositionField(client=client, name='Superposed Field',
                                      fields=[
                                          pos_charged_cylinder_field,
                                          one_by_r_field,
                                          external_field
                                      ])

field_simulator = FieldSimulator(client=client, field=superposed_field)

r = np.linspace(0, 25, num=500, endpoint=True)

r_grid, p_grid, scalar_field, field_max, r_max = field_simulator.barrier_lowering_azimuthal(r_range=r, radius=radius)

energy_scale = 50


fig = mlab.figure(bgcolor=(0.2, 0.2, 0.2))
mlab.mesh(
    r_grid * np.cos(p_grid),
    r_grid * np.sin(p_grid),
    scalar_field * energy_scale,
    colormap='RdBu'
)

mlab.plot3d(
    r_max * np.cos(p_grid[:, 0]),
    r_max * np.sin(p_grid[:, 0]),
    field_max * energy_scale,
    tube_radius=energy_scale * 2e-3,
    color=(1, 0, 0))

client.user_manager.sign_out()

mlab.show()
