from __future__ import division, print_function
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

from BDProjects.Client import Client

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators.SchottkyDiodeSimulator import SchottkyDiodeSimulator

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')

my_diode = SchottkyDiode(client=client, name='Au-n_Si')
#my_diode = SchottkyDiode(client=client, name='Ti-p_Si')
print(my_diode)

diode_simulator = SchottkyDiodeSimulator(client=client, diode=my_diode, description='Schottky diode simulator')
print(diode_simulator)

#temperature_range = np.linspace(0, 700, num=1001, endpoint=True)
#v_bi = diode_simulator.v_bi(temperature=temperature_range)
#plt.plot(temperature_range, v_bi, color='k', linewidth=2, linestyle='-')
#plt.show()

temperature = 300
#def psi(x):
#    v_bi = diode_simulator.v_bi(temperature=temperature)
#    return v_bi - x * v_bi / my_diode.thickness

#z = np.linspace(0, my_diode.thickness, num=100, endpoint=True)

#carriers_concentration = diode_simulator._free_carriers_concentration(z=z, psi=psi, temperature=temperature)
#dopants_occupation = diode_simulator._dopants_equilibrium_occupation(z=z, psi=psi, temperature=temperature)
#charge = diode_simulator._semiconductor_charge(z=z, psi=psi, temperature=temperature)
#_, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
#ax1.plot(z, psi(z), '-k')
#for dopant_name in dopants_occupation:
#    ax3.plot(z, dopants_occupation[dopant_name])
#ax2.plot(z, charge)
#ax4.plot(z, carriers_concentration['electrons'], '-b')
#ax4.plot(z, carriers_concentration['holes'], '-r')
#plt.show()

flat_grid = diode_simulator.potential(bias=1.5, temperature=temperature)
psi = interp1d(flat_grid.physical_nodes, flat_grid.solution, bounds_error=False, fill_value=0.0)
rho = diode_simulator._semiconductor_charge(z=flat_grid.physical_nodes, psi=psi, temperature=temperature)
_, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
ax1.plot(flat_grid.physical_nodes, flat_grid.solution, 'r-o')
ax2.plot(flat_grid.physical_nodes, flat_grid.residual, 'r-o')
ax3.plot(flat_grid.physical_nodes, rho, 'r-o')
ax4.plot(flat_grid.physical_nodes, flat_grid.residual / max(rho))
plt.show()

#bias = np.linspace(-1, 0.5, num=101)
#j = diode_simulator.thermionic_emission_current(bias=bias, temperature=300)
#plt.plot(bias, j, color='k', linewidth=2, linestyle='-')
#plt.show()


client.user_manager.sign_out()
