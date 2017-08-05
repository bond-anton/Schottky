from __future__ import division, print_function

from BDProjects.Client import Client

from Schottky import constants
from Schottky.Samples.Trap import Trap
from Schottky.Samples.Semiconductor import Semiconductor, Dopant

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

epsilon = 11.8
affinity = 4.05

effective_mass = {
    'electron': 0.19 * constants['m_e'],
    'hole': 0.16 * constants['m_e']
}

bands_density_of_states = {
    'Nc': 6.2e21,
    'Nv': 3.5e21
}

band_gap_parameters = {
    'Eg_0': 1.17,
    'alpha': 4.73e-4,
    'beta': 636.0
}

electron_mobility_parameters = {
    'mu_L0': 1430,  # cm2/(V*s)
    'v_s': 1.0e7,  # cm/s
    'A': 4.61e17,  # cm-1 * V-1 * s-1 * K-3/2
    'B': 1.52e15,  # cm-3 * K-2
    'alpha': 2.2,
    'beta': 2
}

hole_mobility_parameters = {
    'mu_L0': 495,  # cm2/(V*s)
    'v_s': 1.0e7,  # cm/s
    'A': 1.0e17,  # cm-1 * V-1 * s-1 * K-3/2
    'B': 6.25e14,  # cm-3 * K-2
    'alpha': 2.2,
    'beta': 1,
}

thermo_emission_coefficient = {
    'electron': 2.1,
    'hole': 0.6
}

phosphorous = Trap(client=client, name='Phosphorous in Si')
phosphorous_dopant = Dopant(client=client, name='Phosphorous', trap=phosphorous,
                            concentration=1e21, description='Phosphorous dopant')
n_type_silicon = Semiconductor(client=client, name='n-type Silicon', epsilon=epsilon,
                               affinity=affinity, effective_mass=effective_mass,
                               effective_bands_density_of_states=bands_density_of_states,
                               band_gap_parameters=band_gap_parameters,
                               electron_mobility_parameters=electron_mobility_parameters,
                               hole_mobility_parameters=hole_mobility_parameters,
                               thermo_emission_coefficient=thermo_emission_coefficient,
                               dopants=[phosphorous_dopant],
                               description='n-type Silicon sample doped with Phosphorous'
                               )

print(n_type_silicon)

boron = Trap(client=client, name='Boron in Si')
boron_dopant = Dopant(client=client, name='Boron', trap=boron,
                      concentration=1e21, description='Boron dopant')
p_type_silicon = Semiconductor(client=client, name='p-type Silicon', epsilon=epsilon,
                               affinity=affinity, effective_mass=effective_mass,
                               effective_bands_density_of_states=bands_density_of_states,
                               band_gap_parameters=band_gap_parameters,
                               electron_mobility_parameters=electron_mobility_parameters,
                               hole_mobility_parameters=hole_mobility_parameters,
                               thermo_emission_coefficient=thermo_emission_coefficient,
                               dopants=[boron_dopant],
                               description='p-type Silicon sample doped with Boron'
                               )

print(p_type_silicon)

client.user_manager.sign_out()
