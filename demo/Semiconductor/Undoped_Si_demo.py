import numpy as np
from matplotlib import pyplot as plt

from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database

from Schottky.Visual.Semiconductor import draw_bands_diagram, draw_bands_diagram_t
from Schottky.Visual.Semiconductor import draw_dopants_profile, draw_dopants_occupation_diagram_t


reference = database[0]

silicon = Semiconductor('Si', reference)

temperature = np.linspace(1.0, 600.0, num=200, endpoint=True)

ax1 = draw_bands_diagram(silicon, 300)
ax2 = draw_bands_diagram_t(silicon, temperature)
ax3 = draw_dopants_profile(silicon)
ax4 = draw_dopants_occupation_diagram_t(silicon, temperature)
plt.show()
