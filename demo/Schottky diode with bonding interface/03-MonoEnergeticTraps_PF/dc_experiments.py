# author: anton

from __future__ import division

import matplotlib

print matplotlib.rcsetup.all_backends
#matplotlib.use("Qt4Agg")
from diode_description import *
from pathos.pools import ProcessPool as Pool

from matplotlib import pyplot as plt
# plt.xkcd()
import pandas as pd
# import numpy as np
import mpmath as mp

from ProjectManager import Project
from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Samples.Diode import Poisson, Visual
from Schottky.Helpers import Psi_approx

temperature_start = 40.0
temperature_stop = 40.0
temperature_step = 20.0
temperature_numpoints = np.floor(abs((temperature_start - temperature_stop) / temperature_step)) + 1
temperature_range = np.linspace(temperature_start, temperature_stop, temperature_numpoints, endpoint=True)

voltage_start = -2.0
voltage_stop = 0.0
voltage_step = 0.25
voltage_numpoints = np.floor(abs((voltage_start - voltage_stop) / voltage_step)) + 1
voltage_range = np.linspace(voltage_start, voltage_stop, voltage_numpoints, endpoint=True)
voltage_start = 0.0
voltage_stop = 3.0
voltage_step = 0.05
voltage_numpoints = np.floor(abs((voltage_start - voltage_stop) / voltage_step)) + 1
voltage_range = np.union1d(voltage_range, np.linspace(voltage_start, voltage_stop, voltage_numpoints, endpoint=True))

poisson_relative_error = 1e-1


def calculate_for_temperature(args):
    silicon, electrode, temperature, voltage_range = args
    print temperature
    db_name = prefix + '_dc_%03.2fK.db' % temperature
    db_name = join(data_dir, db_name)
    project = Project(db_name=db_name, backend='sqlite', hostname='', overwrite=False)
    diode = SchottkyDiode(project, 'Au-nSi_BW', electrode, silicon, DeadLayer=1.5e-7, L=5e-6)
    diode.set_T(temperature)
    potential = []
    field = []
    diode_voltage = []
    diode_voltage_error = []
    current_density = []
    current_density_error = []
    bonding_interface_f = {}
    dopants_f_sum = {}
    z = []
    dopants_f = []
    for voltage in voltage_range:
        diode.set_Va(voltage)
        type_sign = -1 if silicon.dop_type == 'n' else 1
        potential_i = Psi_approx(diode.L, -(diode.V_bi(eV=True) + type_sign * voltage), 0.0)
        potential_i, field_i, z_nodes, rho_err_points, \
            diode_voltage_i, diode_voltage_error_i, \
            current_density_i, current_density_error_i, \
            bonding_interface_f_i, dopants_f_i, \
            measurement_id = Poisson.Reccurent_Poisson_solver(diode, potential_i, Vd_error=1e-6,
                                                              equilibrium_filling=True, t=mp.inf,
                                                              initial_condition_id=-1,
                                                              rho_rel_err=poisson_relative_error, max_iter=100,
                                                              debug=False)
        potential.append(potential_i)
        field.append(field_i)
        diode_voltage.append(diode_voltage_i)
        diode_voltage_error.append(diode_voltage_error_i)
        current_density.append(current_density_i)
        current_density_error.append(current_density_i)
        z.append(z_nodes)
        dopants_f.append(dopants_f_i)
        for trap_key in bonding_interface_f_i.keys():
            try:
                bonding_interface_f[trap_key].append(bonding_interface_f_i[trap_key])
            except KeyError:
                bonding_interface_f[trap_key] = [bonding_interface_f_i[trap_key]]
        for dopant_key in dopants_f_i.keys():
            try:
                dopants_f_sum[dopant_key].append(np.sum(dopants_f_i[dopant_key]))
            except KeyError:
                dopants_f_sum[dopant_key] = [np.sum(dopants_f_i[dopant_key])]

    data = {temperature: {'applied_voltage': voltage_range,
                          'diode_voltage': diode_voltage, 'diode_voltage_error': diode_voltage_error,
                          'current_density': current_density, 'current_density_error': current_density_error}}
    for trap_key in bonding_interface_f.keys():
        data[temperature][trap_key] = bonding_interface_f[trap_key]
    for dopant_key in dopants_f_sum.keys():
        data[temperature][dopant_key + '_sum'] = dopants_f_sum[dopant_key]
    return potential, field, data, z, dopants_f

args = np.empty((len(temperature_range), 4), dtype=object)
args[:, 0] = silicon
args[:, 1] = electrode
args[:, 2] = temperature_range
args[:, 3] = [voltage_range for _ in range(len(temperature_range))]

pool = Pool()
results = np.array(pool.map(calculate_for_temperature, args))
potential = results[:, 0]
field = results[:, 1]
data = results[:, 2]
dopants_f = results[:, 4]

_, ax_iv = plt.subplots()
ax_iv.set_title('I-V temperature dependence')
ax_iv.set_xlabel('Voltage drop on diode, V')
ax_iv.set_ylabel('Current density, A / m^2')
for record in data:
    temperature = record.keys()[0]
    print 'T = %03.2f K' % temperature
    df = pd.DataFrame(record[temperature])
    df.plot(x='diode_voltage', y='current_density', ax=ax_iv, style='-o', label='%2.2f K' % temperature)

ax_iv.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()


def bb_animation(f, ax, temperature, z, data, dopants_f, potential):
    db_name = prefix + '_dc_%03.2fK.db' % temperature
    db_name = join(data_dir, db_name)
    project = Project(db_name=db_name, backend='sqlite', hostname='', overwrite=False)
    diode = SchottkyDiode(project, 'Au-nSi_BW', electrode, silicon, DeadLayer=1.5e-7, L=5e-6)
    diode.set_T(temperature)
    for i, record in enumerate(data):
        current_temperature = record.keys()[0]
        if current_temperature == temperature:
            df = pd.DataFrame(record[temperature])
            break
    print i
    print potential[i]
    bi_trap_labels = []
    bi_f = []
    for bi in diode.Semiconductor.bonding_interfaces:
        for trap in bi.dsl_tilt.traps:
            bi_trap_labels.append(bi.label + '_tilt_' + trap[0].name + '_F')
        for trap in bi.dsl_tilt.traps:
            bi_trap_labels.append(bi.label + '_twist_' + trap[0].name + '_F')
    for voltage in df['applied_voltage'].values:
        bi_f_i = {}
        for trap_label in bi_trap_labels:
            try:
                df_i = df.loc[df['applied_voltage'] == voltage]
                bi_f_i[trap_label] = df_i[trap_label].values[0]
                print bi_f_i[trap_label]
            except KeyError:
                pass
        bi_f.append(bi_f_i)
    ani = Visual.BandsBendingAnimation(f, ax, diode, potential[i],
                                       df['applied_voltage'].values,
                                       df['diode_voltage'].values,
                                       z * 1e-6,
                                       bi_f, dopants_f[i],
                                       np.ones(len(potential[i])) * temperature,
                                       eV=True, draw_metal=True, label_bands=True,
                                       fancy_ticks=False, SchottkyEffect=False, draw_energies=False,
                                       interval=250)
    return ani


z = np.linspace(0, 3, num=200, endpoint=True)
f, ax = plt.subplots()
ani = bb_animation(f, ax, 40, z, data, dopants_f, potential)
plt.show()

