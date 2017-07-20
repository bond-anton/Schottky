from __future__ import division, print_function
import numpy as np
#import matplotlib
#matplotlib.use('Qt4Agg')
from matplotlib import pyplot as plt
#from mayavi import mlab

from BDProjects.Client import Client

from Schottky.Samples.Fields import UniformElectrostaticField, ChargedCylinder, HyperbolicCylinder, SuperpositionField
from Schottky.Samples.Trap import Trap
from Schottky.Simulators.ChargeCarrierTrap import ChargeCarrierTrap

client = Client(config_file_name='config.ini')

client.user_manager.sign_in('bond_anton', 'secret_password')
client.user_manager.project_manager.open_project('Schottky diode')


radius = 0.5
pos_charged_cylinder_field = ChargedCylinder(client=client, name='Pos charged cylinder field',
                                             charge_density=3e7, radius=radius, epsilon=1)
pos_charged_cylinder_field.set_origin([0, 0, 0])

one_by_r_field = HyperbolicCylinder(client=client, name='the Trap',
                                    coefficient=-0.1, radius=radius)
one_by_r_field.set_origin([0, 0, 0])


external_field = UniformElectrostaticField(client=client, name='Uniform electrostatic field',
                                           strength=0.01, direction=[1, 0, 0])


superposed_field = SuperpositionField(client=client, name='Superposed Field',
                                      fields=[
                                          pos_charged_cylinder_field,
                                          one_by_r_field,
                                          external_field
                                      ])

my_trap = Trap(client=client, name='Shallow Trap', description='Test Shallow Trap',
               charge_state={'empty': 0, 'full': 1}, activation_energy=0.15,
               band='Ec', energy_distribution_function='Single level', energy_spread=0.3,
               electron_capture_cross_section=1e-21, electron_capture_cross_section_activation_energy=0.0,
               hole_capture_cross_section=1e-21, hole_capture_cross_section_activation_energy=0.0,
               trap_potential=superposed_field)

my_trap_simulator = ChargeCarrierTrap(client=client, trap=my_trap)
print(my_trap_simulator.parts['Field Simulator'])

f = [1.0]
t = [0]
temperature = 300
band_gap = 1
threshold = 1e-4
v_n = 1e7
v_p = 1e7
n_c = 1e18
n_v = 1e18
#for i in range(50):
while t[-1] < 0.2:
    #print('t: %g s, F: %2.2f' % (t[-1], f[-1]))
    e_n, e_p = my_trap_simulator.emission_rate(temperature, band_gap, v_n, v_p, n_c, n_v,
                                               poole_frenkel_n=1.0, poole_frenkel_p=1.0)
    #print(e_n, e_p)
    d_t = min(threshold / e_n, threshold / e_p)
    d_f_n = d_t * e_n * f[-1]
    d_f_p = d_t * e_p * (1 - f[-1])
    #print('d_t: %g s, d_f_n: %g, d_f_p: %g' % (d_t, d_f_n, d_f_p))
    f.append(f[-1] - (d_f_n - d_f_p))
    t.append(t[-1] + d_t)

client.user_manager.sign_out()

plt.plot(t, f, '-bo')
plt.show()

#mlab.show()
