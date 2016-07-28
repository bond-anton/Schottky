from __future__ import division, print_function
import numpy as np
from matplotlib import pyplot as plt

from ScientificProjects.Client import Client

from Schottky.Samples.Semiconductor import Semiconductor
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_semiconductor = Semiconductor(client=client, name='n-type Silicon')
print(my_semiconductor)

my_simulator = BulkSemiconductor(client=client, semiconductor=my_semiconductor,
                                 description='Bulk %s properties measurement' % my_semiconductor.name)

print(my_simulator)

temperature_range = np.linspace(10, 700, num=1001, endpoint=True)
band_gap = my_simulator.band_gap(temperature=temperature_range)
plt.plot(temperature_range, band_gap, color='k', linewidth=2, linestyle='-')
for part in my_simulator.parts.values():
    if isinstance(part, ChargeCarrierTrap):
        plt.plot(temperature_range, band_gap - part.energy_level(band_gap)['empty'],
                 color='k', linewidth=1, linestyle='--')
plt.show()

bands_density_of_states = my_simulator.bands_density_of_states(temperature_range)
plt.plot(temperature_range, bands_density_of_states['Nc'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, bands_density_of_states['Nv'],
         color='r', linewidth=2, linestyle='-')
plt.show()

thermal_velocity = my_simulator.carriers_thermal_velocity(temperature_range)
plt.plot(temperature_range, thermal_velocity['electron'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, thermal_velocity['hole'],
         color='r', linewidth=2, linestyle='-')
plt.show()


mobility = my_simulator.mobility(temperature=temperature_range, field=None)

print(mobility['electron']['lattice'])

plt.plot(temperature_range, mobility['electron']['total'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, mobility['electron']['lattice'],
         color='b', linewidth=1, linestyle='--')
plt.plot(temperature_range, mobility['electron']['impurities'],
         color='b', linewidth=1, linestyle=':')
plt.show()

mobility = my_simulator.mobility(temperature=10, field=None)
print('electrons: %2.2f cm^2 / (V*s)' % mobility['electron']['total'])
print('electrons: %2.2f cm^2 / (V*s)' % mobility['electron']['lattice'])
print('electrons: %2.2f cm^2 / (V*s)' % mobility['electron']['impurities'])
print('electrons: %2.2f cm^2 / (V*s)' % mobility['electron']['carrier-carrier'])
print('holes: %2.2f cm^2 / (V*s)' % mobility['hole']['total'])

client.user_manager.sign_out()
