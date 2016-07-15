from __future__ import division, print_function
import numpy as np
from matplotlib import pyplot as plt
#from mayavi import mlab

from ScientificProjects.Client import Client

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, HyperbolicCylinder, SuperpositionField
from Schottky.Samples.Trap import Trap
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


radius = 0.5
pos_charged_cylinder_field = ChargedCylinder(client=client, name='Pos charged cylinder field',
                                             charge_density=3e7, radius=radius, epsilon=1)
pos_charged_cylinder_field.set_origin([0, 0, 0])

one_by_r_field = HyperbolicCylinder(client=client, name='the Trap',
                                    coefficient=-0.1, radius=radius)
one_by_r_field.set_origin([0, 0, 0])


external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                           strength=0.01, direction=[1, 0, 0])


superposed_field = SuperpositionField(client=client, name='Superposed Field',
                                      fields=[
                                          pos_charged_cylinder_field,
                                          one_by_r_field,
                                          external_field
                                      ])

my_trap = Trap(client=client, name='Shallow Trap', description='Test Shallow Trap',
               charge_state={'empty': 0, 'full': 1}, activation_energy={'empty': 0.15, 'full': 0.15},
               band='Ec', energy_distribution_function='Single level', energy_spread=0.3,
               electron_capture_cross_section=1e-21, electron_capture_cross_section_activation_energy=0.0,
               hole_capture_cross_section=1e-21, hole_capture_cross_section_activation_energy=0.0,
               trap_potential=superposed_field)

my_trap_simulator = ChargeCarrierTrap(client=client, trap=my_trap)
print(my_trap_simulator.parts['Field Simulator'])

client.user_manager.sign_out()

#mlab.show()
