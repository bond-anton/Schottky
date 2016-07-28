from __future__ import division, print_function

from ScientificProjects.Client import Client

from Schottky.Samples.Trap import Trap

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


boron = Trap(client=client, name='Boron in Si', description='Boron in Silicon',
             charge_state={'empty': 0, 'full': -1}, activation_energy={'empty': 0.044, 'full': 0.044},
             band='Ev', energy_distribution_function='Single level',
             electron_capture_cross_section=1e-15, electron_capture_cross_section_activation_energy=0.0,
             hole_capture_cross_section=1e-15, hole_capture_cross_section_activation_energy=0.0)

print(boron)

phosphorous = Trap(client=client, name='Phosphorous in Si', description='Phosphorous in Silicon',
                   charge_state={'empty': +1, 'full': 0}, activation_energy={'empty': 0.045, 'full': 0.045},
                   band='Ec', energy_distribution_function='Single level',
                   electron_capture_cross_section=1e-15, electron_capture_cross_section_activation_energy=0.0,
                   hole_capture_cross_section=1e-15, hole_capture_cross_section_activation_energy=0.0)

print(phosphorous)

client.user_manager.sign_out()

