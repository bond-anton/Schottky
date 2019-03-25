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
c_p.solution = np.ones(c_p.num) * 1e21
f_p.solution = np.zeros(f_p.num)
phosphorus = Dopant('P', True, TreeMesh1DUniform(c_p, aligned=True), TreeMesh1DUniform(f_p, aligned=True),
                    0.045 * constant.q, silicon.band_gap_t(0.0) - 0.045 * constant.q,
                    1e-15, 1e-15)
phosphorus.charge_state = {0: +1, 1: 0}
phosphorus.color = 'b'
phosphorus.linestyle = '--'

c_t = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f_t = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c_t.solution = np.ones(c_t.num) * 1e17
f_t.solution = np.zeros(f_t.num)
deep_trap = Dopant('Et', True, TreeMesh1DUniform(c_t, aligned=True), TreeMesh1DUniform(f_t, aligned=True),
                   0.3 * constant.q,  silicon.band_gap_t(0.0) - 0.3 * constant.q,
                   1e-15, 1e-15)
deep_trap.charge_state = {0: -1, 1: -2}
deep_trap.capture_barrier_ev = {0: 0.6, 1: 0.0}
deep_trap.color = 'g'
deep_trap.linestyle = '--'
deep_trap.marker = ''

silicon.dopants = [phosphorus, deep_trap]

temperature = np.linspace(1.0, 600.0, num=600, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)
plt.show()
