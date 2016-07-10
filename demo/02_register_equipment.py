from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky import constants
from Schottky import Simulators
from Schottky.Samples.Metal import Metal
from Schottky.Samples.Trap import Trap
from Schottky.Samples.Fields import UniformElectrostaticField

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_metal = Metal(client=client, name='Gold',
                 work_function=5.1)

my_metal = Metal(client=client, name='Gold2')

print(my_metal)
print(my_metal.work_function)
my_metal.set_work_function(work_function=5.1)
print(my_metal.work_function)

my_trap = Trap(client=client, name='Test trap')
my_trap.set_band('Ec')
print(my_trap)

my_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                     strength=1, direction=[1, 0, 0])
print(my_field.sample.id)
print(my_field.__class__.__module__)

#sim = Simulators.Simulator(client)

#print(sim)
#print(sim.equipment)

client.user_manager.sign_out()
