from __future__ import division, print_function

import Space_visualization as Visual
from ScientificProjects.Client import Client

import numpy as np

from mayavi import mlab

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, SuperpositionField
from Schottky.Simulators.Field import FieldSimulator

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


radius = 0.5
pos_charged_cylinder_field = ChargedCylinder(client=client, name='Pos charged cylinder field',
                                             charge_density=3e7, radius=radius, epsilon=1)
#pos_charged_cylinder_field.set_origin([-2, 0, 0])
#pos_charged_cylinder_field.rotate_axis_angle([1, 0, 0], np.pi / 4)
neg_charged_cylinder_field = ChargedCylinder(client=client, name='Neg charged cylinder field',
                                             charge_density=-3e7, radius=radius, epsilon=1)
#neg_charged_cylinder_field.set_origin([2, 0, 0])
#neg_charged_cylinder_field.rotate_axis_angle([1, 0, 0], -np.pi / 4)

external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                           strength=0.01, direction=[1, 0, 0])


superposed_field = SuperpositionField(client=client, name='Superposed Field',
                                      #fields=[
                                      #        pos_charged_cylinder_field,
                                      #        neg_charged_cylinder_field,
                                      #        external_field
                                      #        ]
                                      )

field_simulator = FieldSimulator(client=client, field=superposed_field)

fig = mlab.figure(bgcolor=(0.2, 0.2, 0.2))

superposed_field_vis = Visual.FieldView(fig, field_simulator,
                                        grid=np.mgrid[-5:5:20j, -5:5:20j, -5:5:20j],
                                        scalar_field_visible=False,
                                        vector_field_visible=True,
                                        cs_visible=False)
superposed_field_vis.draw()


client.user_manager.sign_out()

mlab.show()
