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
deep_trap.charge_state = {0: -1, 1: -2}
deep_trap.capture_barrier_ev = {0: 0.6, 1: 0.0}
deep_trap.color = 'g'
deep_trap.linestyle = '--'
deep_trap.marker = ''

silicon.dopants = [phosphorus, deep_trap]

temperature = np.linspace(1.0, 600.0, num=600, endpoint=True)

f_threshold = 1.0e-28
max_iter = 100
verbose = True

ax1 = draw_bands_diagram(silicon, 300, f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
ax2 = draw_bands_diagram_t(silicon, temperature, f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature,
                                        f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
plt.show()
