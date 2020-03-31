import numpy as np
from matplotlib import pyplot as plt

from BDFunction1D.Standard import Constant

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky import constant

from Schottky.Visual.Semiconductor import draw_bands_diagram, draw_bands_diagram_t
from Schottky.Visual.Semiconductor import draw_dopants_profile, draw_dopants_occupation_diagram_t


reference = database[0]

silicon = Semiconductor('Si', reference)

phosphorus = Dopant('P', True, Constant(1e21),
                    0.045 * constant.q, silicon.band_gap_t(0.0) - 0.045 * constant.q,
                    1e-15, 1e-15)
phosphorus.charge_state = {0: +1, 1: 0}
phosphorus.color = 'b'
phosphorus.linestyle = '--'

silicon.dopants = [phosphorus]

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)
plt.show()

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)
n_i = silicon.n_i(temperature)
e_f = silicon.el_chem_pot(temperature)
n_e = silicon.n_e(e_f, temperature)
n_h = silicon.n_h(e_f, temperature)
fig, ax = plt.subplots()
ax.semilogy(temperature, np.asarray(n_e)*1e-6, '-b')
ax.semilogy(temperature, np.asarray(n_h)*1e-6, '-r')
ax.semilogy(temperature, np.asarray(n_i)*1e-6, '-k')
plt.show()
n_e = silicon.n_e_n_i(e_f, temperature)
n_h = silicon.n_h_n_i(e_f, temperature)
fig, ax = plt.subplots()
ax.semilogy(temperature, np.asarray(n_h)*1e-10, '-r')
ax.semilogy(temperature, np.asarray(n_e)*1e-10, '-b')
ax.semilogy(temperature, np.ones_like(temperature)*1e-10, '-k')
plt.show()
