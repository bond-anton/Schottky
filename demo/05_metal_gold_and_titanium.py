from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky.Samples.Metal import Metal

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

gold = Metal(client=client, name='Au', work_function=5.1, description='Gold')

print(gold)

titanium = Metal(client=client, name='Ti', work_function=4.33, description='Titanium')

print(titanium)

client.user_manager.sign_out()
