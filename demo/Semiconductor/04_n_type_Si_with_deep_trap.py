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

deep_trap = Dopant('Et', True, Constant(1e17),
                   0.3 * constant.q,  silicon.band_gap_t(0.0) - 0.3 * constant.q,
                   1e-15, 1e-15)
deep_trap.charge_state = {0: 0, 1: -1}
deep_trap.color = 'g'
deep_trap.linestyle = '--'

silicon.dopants = [phosphorus, deep_trap]

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)
plt.show()
