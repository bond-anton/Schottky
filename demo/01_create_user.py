from __future__ import division, print_function

from ScientificProjects.Client import Client

client = Client(config_file_name='config.ini')

client.user_manager.create_user('Anton', 'Bondarenko', 'bond.anton@gmail.com', 'bond_anton', 'secret_password')
