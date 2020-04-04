from matplotlib import pyplot as plt

from BDFunction1D.Standard import Constant

from Schottky.Dopant import Dopant
from Schottky.Semiconductor import Semiconductor
from Schottky.Reference import database
from Schottky.Metal import Metal
from Schottky.SchottkyDiode import SchottkyDiode
from Schottky import constant

from Schottky.Visual.SchottkyDiode import plot_n_eh, plot_generation_recombination, \
    plot_ep, plot_ef, plot_any_fun, plot_qfe, plot_qfh, plot_qfeh, plot_band_diagram, plot_charge_density


reference = database[0]

silicon = Semiconductor('Si', reference)

phosphorus = Dopant('P', True, Constant(1e21),
                    0.045 * constant.q, silicon.band_gap_t(0.0) - 0.045 * constant.q,
                    1e-15, 1e-15)
phosphorus.charge_state = {0: +1, 1: 0}
phosphorus.color = 'b'
phosphorus.linestyle = '--'

silicon.dopants = [phosphorus]

electrode = Metal('Au', 5.1 * constant.q)
diode = SchottkyDiode('Au-n-Si', electrode, silicon,
                      length=1.0e-5)
diode.contact_diameter = 1.0e-3
print('Diode area is %2.2f mm^2' % (diode.area * 1e6))

temperature = 300.0
bias = -0.5

print(diode.phi_b_n_ev_t(temperature), diode.phi_b_p_ev_t(temperature))
print(diode.n0_t(temperature), diode.p0_t(temperature))


diode.temperature = temperature
diode.bias = bias

for i in range(1):
    diode.stationary_grad_qf_e_solver()
    diode.stationary_grad_qf_h_solver()
    rho = diode.poisson_eq_solver()


print(diode.thermionic_emission_current_e())
print(diode.thermionic_emission_current_h())

plot_band_diagram(diode)
plt.show()

plot_n_eh(diode)
plt.show()

plot_charge_density(diode)
plt.show()

#
# plot_generation_recombination(diode)
# plt.show()
#
plot_ep(diode)
plt.show()

plot_ef(diode)
plt.show()

plot_any_fun(diode, rho)
plt.show()

#
# plot_qfe(diode)
# plt.show()

#
# plot_qfh(diode)
# plt.show()


plot_qfeh(diode)
plt.show()
