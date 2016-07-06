from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky import constants
from Schottky import Simulators
from Schottky.Samples.Metal import Metal

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_metal = Metal(client=client, name='Gold',
                 work_function=5.2)

my_metal = Metal(client=client, name='Gold3')

print(my_metal)
print(my_metal.work_function)
my_metal.set_work_function(work_function=5.1)
print(my_metal.work_function)

#sim = Simulators.Simulator(client)

#print(sim)
#print(sim.equipment)

client.user_manager.sign_out()
