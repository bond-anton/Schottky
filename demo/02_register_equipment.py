from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky import constants
from Schottky import Simulators
from Schottky.Samples.Metal import Metal
from Schottky.Samples.Semiconductor.Trap import Trap

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

#sim = Simulators.Simulator(client)

#print(sim)
#print(sim.equipment)

client.user_manager.sign_out()
