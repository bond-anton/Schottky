from matplotlib import pyplot as plt

from BDFunction1D.Standard import Constant

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky.Metal import Metal
from Schottky.SchottkyDiode import SchottkyDiode
from Schottky import constant

from Schottky.Visual.SchottkyDiode import plot_n_eh, plot_generation_recombination, \
    plot_ep, plot_ef, plot_qfe, plot_qfh, plot_band_diagram


reference = database[0]

silicon = Semiconductor('Si', reference)

boron = Dopant('B', False, Constant(1e21),
               silicon.band_gap_t(300.0) - 0.045 * constant.q,  0.045 * constant.q,
               1e-15, 1e-15)
boron.charge_state = {0: 0, 1: -1}
boron.color = 'b'
boron.linestyle = '--'

silicon.dopants = [boron]

electrode = Metal('Ti', 4.1 * constant.q)
diode = SchottkyDiode('Ti-p-Si', electrode, silicon,
                      length=1.0e-5)
diode.contact_diameter = 1.0e-3
print('Diode area is %2.2f mm^2' % (diode.area * 1e6))

temperature = 300.0
bias = 1.0

diode.temperature = temperature
diode.bias = bias

plot_band_diagram(diode)
plt.show()


plot_n_eh(diode)
plt.show()

# plot_generation_recombination(measurement)
# plt.show()
#
#
#
plot_ep(diode)
plt.show()

plot_ef(diode)
plt.show()
#
# plot_qfe(measurement)
# plt.show()
#
# plot_qfh(measurement)
# plt.show()
