from __future__ import division, print_function

from BDProjects.Client import Connector, Client

client = Client(Connector(config_file_name='config.ini'))
client.user_manager.sign_in('administrator', 'admin')
client.user_manager.create_user('bond_anton', 'secret_password', 'bond.anton@gmail.com', 'Anton', 'Bondarenko')
client.user_manager.sign_out()
client.user_manager.sign_in('bond_anton', 'secret_password')

client.user_manager.project_manager.create_project(name='Schottky diode',
                                                   description='Schottky diode simulations',
                                                   data_dir='./data/files')

client.user_manager.sign_out()
