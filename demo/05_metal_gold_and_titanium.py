from __future__ import division, print_function

from BDProjects.Client import Client

from Schottky.Samples.Metal import Metal

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

gold = Metal(client=client, name='Gold', work_function=5.1, description='Gold')

print(gold)

titanium = Metal(client=client, name='Titanium', work_function=4.33, description='Titanium')

print(titanium)

client.user_manager.sign_out()
