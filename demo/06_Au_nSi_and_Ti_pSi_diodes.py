from __future__ import division, print_function
import numpy as np

from BDProjects.Client import Client

from Schottky.Samples.Metal import Metal
from Schottky.Samples.Semiconductor import Semiconductor
from Schottky.Samples.Diode import SchottkyDiode

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

area = np.pi * (0.15 ** 2) / 4  # cm2
thickness = 10 * 1e-4  # cm
resistance = 50  # Ohm

gold = Metal(client=client, name='Gold')
n_silicon = Semiconductor(client=client, name='n-type Silicon')
gold_silicon = SchottkyDiode(client=client, name='Au-n_Si',
                             area=area, thickness=thickness, serial_resistance=resistance,
                             metal=gold, semiconductor=n_silicon, description='Au - n-Si Schottky diode')
print(gold_silicon)

titanium = Metal(client=client, name='Titanium')
p_silicon = Semiconductor(client=client, name='p-type Silicon')
titanium_silicon = SchottkyDiode(client=client, name='Ti-p_Si',
                                 area=area, thickness=thickness, serial_resistance=resistance,
                                 metal=titanium, semiconductor=p_silicon, description='Ti - p-Si Schottky diode')
print(titanium_silicon)

client.user_manager.sign_out()
