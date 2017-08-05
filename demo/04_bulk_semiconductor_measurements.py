from __future__ import division, print_function
import numpy as np
from matplotlib import pyplot as plt

from BDProjects.Client import Client

from Schottky.Samples.Semiconductor import Semiconductor
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_semiconductor = Semiconductor(client=client, name='p-type Silicon')
print(my_semiconductor)

my_simulator = BulkSemiconductor(client=client, semiconductor=my_semiconductor,
                                 description='Bulk %s properties measurement' % my_semiconductor.name)

print(my_simulator)

temperature_range = np.linspace(0, 750, num=1001, endpoint=True)

band_gap = my_simulator.band_gap(temperature=5.0)
print(band_gap)

band_gap = my_simulator.band_gap(temperature=temperature_range)

plt.plot(temperature_range, band_gap, color='k', linewidth=2, linestyle='-')
for part in my_simulator.parts.values():
    if isinstance(part, ChargeCarrierTrap):
        plt.plot(temperature_range, band_gap - part.energy_level(band_gap),
                 color='k', linewidth=1, linestyle='--')
plt.show()

bands_density_of_states = my_simulator.effective_bands_density_of_states(temperature=temperature_range)
plt.plot(temperature_range, bands_density_of_states['DOS C.band'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, bands_density_of_states['DOS V.band'],
         color='r', linewidth=2, linestyle='-')
plt.show()

exit(0)

thermal_velocity = my_simulator.carriers_thermal_velocity(temperature_range, use_storage=True)
plt.plot(temperature_range, thermal_velocity['electron'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, thermal_velocity['hole'],
         color='r', linewidth=2, linestyle='-')
plt.show()


mobility = my_simulator.mobility(temperature=temperature_range, field=None, use_storage=True)

plt.plot(temperature_range, mobility['electron']['total'],
         color='b', linewidth=2, linestyle='-')
plt.plot(temperature_range, mobility['hole']['total'],
         color='r', linewidth=2, linestyle='-')
plt.show()

e_f = my_simulator.electrochemical_potential(temperature_range, use_storage=True)
plt.plot(temperature_range, band_gap, color='k', linewidth=2, linestyle='-')
plt.plot(temperature_range, band_gap / 2, color='k', linewidth=1, linestyle=':')
for part in my_simulator.parts.values():
    if isinstance(part, ChargeCarrierTrap):
        plt.plot(temperature_range, band_gap - part.energy_level(band_gap),
                 color='k', linewidth=1, linestyle='--')
plt.plot(temperature_range, band_gap - e_f, 'k-o')
plt.show()


e_f = my_simulator.electrochemical_potential(temperature=0, use_storage=True)
print(e_f)

client.user_manager.sign_out()
