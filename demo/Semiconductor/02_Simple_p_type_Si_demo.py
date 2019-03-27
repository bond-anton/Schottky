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

c_b = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f_b = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c_b.solution = np.ones(c_b.num) * 1e21
f_b.solution = np.zeros(f_b.num)
boron = Dopant('B', False, TreeMesh1DUniform(c_b, aligned=True), TreeMesh1DUniform(f_b, aligned=True),
               silicon.band_gap_t(0.0) - 0.045 * constant.q,  0.045 * constant.q,
               1e-15, 1e-15)
boron.charge_state = {0: 0, 1: -1}
boron.color = 'r'
boron.linestyle = '--'

silicon.dopants = [boron]

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)
plt.show()
