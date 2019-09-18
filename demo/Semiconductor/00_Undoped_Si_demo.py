import numpy as np
from matplotlib import pyplot as plt

from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database

from Schottky.Visual.Semiconductor import draw_bands_diagram, draw_bands_diagram_t
# from Schottky.Visual.Semiconductor import draw_dopants_profile, draw_dopants_occupation_diagram_t


reference = database[0]

silicon = Semiconductor('Si', reference)

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
# ax3 = draw_dopants_profile(silicon)
# ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)

e_f = silicon.el_chem_pot(temperature)
fig, ax = plt.subplots()
ax.plot(temperature, np.asarray(e_f) / 1.6e-19)

plt.show()

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)
band_gap_boltzmann = silicon.band_gap_boltzmann(temperature)
e_i_boltzmann = silicon.e_i_boltzmann(temperature)
fig, ax = plt.subplots()
ax.plot(temperature, np.asarray(e_i_boltzmann) / np.asarray(band_gap_boltzmann), '-r')

plt.show()

n_i = silicon.n_i(temperature)
fig, ax = plt.subplots()
ax.plot(temperature, np.asarray(n_i), '-k')
plt.show()
