from __future__ import division, print_function
import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt

from ScientificProjects.Client import Client

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators.SchottkyDiodeSimulator import SchottkyDiodeSimulator

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_diode = SchottkyDiode(client=client, name='Au-n_Si')
print(my_diode)

diode_simulator = SchottkyDiodeSimulator(client=client, diode=my_diode, description='Schottky diode simulator')
print(diode_simulator)

temperature_range = np.linspace(0, 700, num=1001, endpoint=True)
v_bi = diode_simulator.v_bi(temperature=temperature_range)
plt.plot(temperature_range, v_bi, color='k', linewidth=2, linestyle='-')
plt.show()

client.user_manager.sign_out()
