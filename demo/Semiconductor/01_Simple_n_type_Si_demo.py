import numpy as np
from matplotlib import pyplot as plt

from BDMesh.Mesh1DUniform import Mesh1DUniform
from BDMesh.TreeMesh1DUniform import TreeMesh1DUniform

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky import constant

from Schottky.Visual.Semiconductor import draw_bands_diagram, draw_bands_diagram_t
from Schottky.Visual.Semiconductor import draw_dopants_profile, draw_dopants_occupation_diagram_t


reference = database[0]

silicon = Semiconductor('Si', reference)

c_p = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f_p = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c_p.solution = np.ones(c_p.num) * 1e21  # 1e21 m^-3 = 1e15 cm^-3
f_p.solution = np.zeros(f_p.num)
phosphorus = Dopant('P', True, TreeMesh1DUniform(c_p, aligned=True), TreeMesh1DUniform(f_p, aligned=True),
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
n_e = silicon.n_e_n_i(e_f, temperature)
n_h = silicon.n_h_n_i(e_f, temperature)
fig, ax = plt.subplots()
ax.semilogy(temperature, np.asarray(n_e)*1e-6, '-b')
ax.semilogy(temperature, np.asarray(n_h)*1e-6, '-r')
ax.semilogy(temperature, np.asarray(n_i)*1e-6, '-k')
plt.show()
fig, ax = plt.subplots()
ax.semilogy(temperature, np.asarray(n_h), '-r')
# ax.semilogy(temperature, np.asarray(n_e), '-b')
plt.show()