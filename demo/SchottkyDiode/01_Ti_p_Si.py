import numpy as np
from matplotlib import pyplot as plt

from BDMesh.Mesh1DUniform import Mesh1DUniform
from BDMesh.TreeMesh1DUniform import TreeMesh1DUniform

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky.Metal import Metal
from Schottky.SchottkyDiode import SchottkyDiode
from Schottky import constant

from Schottky.DCMeasurement import DCMeasurement

from Schottky.Visual.DCMeasurement import plot_n_eh, plot_generation_recombination, plot_ep


reference = database[0]

silicon = Semiconductor('Si', reference)

c_p = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
f_p = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
c_p.solution = np.ones(c_p.num) * 1e21  # 1e21 m^-3 = 1e15 cm^-3
f_p.solution = np.zeros(f_p.num)
boron = Dopant('B', True, TreeMesh1DUniform(c_p, aligned=True), TreeMesh1DUniform(f_p, aligned=True),
               silicon.band_gap_t(0.0) - 0.045 * constant.q,  0.045 * constant.q,
               1e-15, 1e-15)
boron.charge_state = {0: +1, 1: 0}
boron.color = 'b'
boron.linestyle = '--'

silicon.dopants = [boron]

electrode = Metal('Ti', 4.1 * constant.q)
diode = SchottkyDiode('Ti-p-Si', electrode, silicon,
                      length=1.0e-5)
diode.contact_diameter = 1.0e-3
print('Diode area is %2.2f mm^2' % (diode.area * 1e6))

measurement = DCMeasurement('The measurement', diode, initial_step=1e-7)
plot_n_eh(measurement)
plt.show()

plot_generation_recombination(measurement)
plt.show()

measurement.prepare_psi0(1.0, 300)

plot_ep(measurement)
plt.show()
