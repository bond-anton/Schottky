from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky.Samples.Semiconductor import Semiconductor

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_semiconductor = Semiconductor(client=client, name='n-type Silicon')
print(my_semiconductor)

client.user_manager.sign_out()
