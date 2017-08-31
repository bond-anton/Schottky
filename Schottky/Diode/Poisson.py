# -*- coding: utf-8 -*-

'''
Created on 18 aug 2015 г.

@author: anton
'''
from __future__ import division
import time
import datetime

import numpy as np
import mpmath as mp

from matplotlib import pyplot as plt

from scipy.constants.codata import epsilon0
from scipy.interpolate import interp1d

from Schottky.Notation import q, k
from Schottky.Helpers import smooth_dd, Psi_zero, to_numeric, interp_Fn

from NumericalDE.FiniteDifference1D import dirichlet_non_linear_poisson_solver_mesh, points_for_refinement, adjust_range
from NumericalDE.Mesh import UniformMesh1D, Uniform1DMeshesTree


def add_Electric_Field_Measurement(SchDiode, Va, equilibrium_filling=True, t=mp.inf, initial_condition_id=-1):
    '''
    '''
    if equilibrium_filling:
        measurement_name = 'Potential and Field t=%s [T=%2.2f, Va=%2.2f]' % (t, SchDiode.T, Va)
        measurement_description = 'DC Potential and Electric Field at t=%s, T=%2.2f, Va=%2.2f' % (t, SchDiode.T, Va)
    else:
        if initial_condition_id < 0:
            raise Exception('Nonequilibrium experiments require initial condition')
        measurement_name = 'Potential and Field ic_id=%d t=%s [T=%2.2f, Va=%2.2f]' \
                           % (initial_condition_id, t, SchDiode.T, Va)
        measurement_description = 'Nonequilibrium Potential and Electric Field at t=%s, T=%2.2f, Va=%2.2f' \
                                  % (t, SchDiode.T, Va)
    measurement_id = SchDiode.Project.add_measurement(
        SchDiode.measurements_types['Potential and Electric Field measurement'],
        measurement_name, measurement_description,
        measurement_datetime=None, create_new=False)
    return measurement_id


def gen_point_names_Electric_Field_measurement(SchDiode):
    point_names_long = {
        'ic_id': ['initial condition id', 'ic_id', ''],
        'z': ['z coordinate', 'z', 'm'],
        'Psi': ['Potential', 'Psi', 'V'],
        'E': ['Electric Field', 'E', 'V/m'],
        'rho_rel_err': ['Charge density relative error', 'rho_rel_err', ''],
        'J': ['Current density', 'J', 'A/m'],
        'J_err': ['Current density error', 'J err', 'A/m'],
        'Vd': ['Voltage bias', 'Vd', 'V'],
        'Vd_err': ['Voltage bias error', 'Vd err', 'V'],
        'n': ['Free electrons concentration', 'n', 'm^(-3)'],
        'p': ['Free holes concentration', 'p', 'm^(-3)'],
    }
    for BI in SchDiode.Semiconductor.bonding_interfaces:
        for trap in BI.dsl_tilt.traps:
            point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'] = [BI.label + '_tilt_' + trap[0].name + '_F',
                                                                           'F', '']
        for trap in BI.dsl_twist.traps:
            point_names_long[BI.label + '_twist_' + trap[0].name + '_F'] = [BI.label + '_twist_' + trap[0].name + '_F',
                                                                            'F', '']
    for dopant in SchDiode.Semiconductor.dopants:
        point_names_long[dopant.name + '_F'] = [dopant.name + '_F', 'F', '']
    return point_names_long


def load_diode_state(SchDiode, measurement_id, initial_condition_id=-1, debug=False):
    point_names_long = gen_point_names_Electric_Field_measurement(SchDiode)
    if debug: print '==> Looking for saved solution in the database'
    points_names = [point_names_long['z'][0], point_names_long['Psi'][0], point_names_long['E'][0],
                    point_names_long['Vd'][0], point_names_long['Vd_err'][0], point_names_long['J'][0],
                    point_names_long['J_err'][0], point_names_long['rho_rel_err'][0], point_names_long['ic_id'][0]]
    for BI in SchDiode.Semiconductor.bonding_interfaces:
        for trap in BI.dsl_tilt.traps:
            points_names.append(point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'][0])
        for trap in BI.dsl_twist.traps:
            points_names.append(point_names_long[BI.label + '_twist_' + trap[0].name + '_F'][0])
    for dopant in SchDiode.Semiconductor.dopants:
        points_names.append(point_names_long[dopant.name + '_F'][0])
    try:
        data = SchDiode.Project.get_data_points_by_names(measurement_id, points_names)
    except:
        if debug: print '==> No solutions found'
        return False, Psi_zero, Psi_zero, np.zeros(1), np.zeros(1), 0, 0, 0, 0, {}, {}, measurement_id

    z_nodes = data[points_names[0]][:, 0]
    Psi_points = data[points_names[1]][:, 0]
    E_points = data[points_names[2]][:, 0]
    Vd = data[points_names[3]][:, 0]
    Vd_err = data[points_names[4]][:, 0]
    J = data[points_names[5]][:, 0]
    J_err = data[points_names[6]][:, 0]
    rho_rel_error_array = data[points_names[7]][:, 0]
    ic_id = data[points_names[8]][:, 0]
    if len(rho_rel_error_array) > 0 and (ic_id[0] == initial_condition_id or initial_condition_id == -1):
        if debug: print '==> Solution found'
        Psi = interp_Fn(z_nodes, Psi_points, interp_type='last')
        E = interp_Fn(z_nodes, E_points, interp_type='last')
        SchDiode.FieldInversionPoint, _, _, _ = SchDiode.get_phi_bn(Psi, Vd, SchottkyEffect=False)
        print '*** !!! Inversion pt (load) !!!', SchDiode.FieldInversionPoint
        BI_F = {}
        for BI in SchDiode.Semiconductor.bonding_interfaces:
            F_tilt = []
            F_twist = []
            for trap in BI.dsl_tilt.traps:
                F_i = data[point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'][0]][:, 0]
                BI_F[BI.label + '_tilt_' + trap[0].name + '_F'] = F_i[0]
                F_tilt.append(F_i[0])
            for trap in BI.dsl_twist.traps:
                F_i = data[point_names_long[BI.label + '_twist_' + trap[0].name + '_F'][0]][:, 0]
                BI_F[BI.label + '_twist_' + trap[0].name + '_F'] = F_i[0]
                F_twist.append(F_i[0])
            BI.set_traps_f(np.array(F_tilt), np.array(F_twist))
        dopants_F = {}
        for dopant in SchDiode.Semiconductor.dopants:
            F_i = data[point_names_long[dopant.name + '_F'][0]][:, 0]
            dopants_F[dopant.name + '_F'] = F_i
            # print 'l1'
            dopant.set_F_interp(z_nodes, F_i)
            # print 'l2'
        return True, Psi, E, z_nodes, rho_rel_error_array, Vd, Vd_err, J, J_err, BI_F, dopants_F, measurement_id
    else:
        if debug: print '==> No solutions found'
        return False, Psi_zero, Psi_zero, np.zeros(1), np.zeros(1), 0, 0, 0, 0, {}, {}, measurement_id


def save_diode_state(SchDiode, measurement_id, initial_condition_id,
                     z_nodes, Psi_points, E_points, rho_rel_err,
                     Vd, Vd_err, J, J_err, debug=False):
    point_names_long = gen_point_names_Electric_Field_measurement(SchDiode)
    if debug: print 'Writing data to database'
    db_write_start_time = time.time()
    measurement_time = datetime.datetime.now()
    saved_rho_rel_error_array, _, _, _ = SchDiode.Project.get_data_points(measurement_id,
                                                                          point_names_long['rho_rel_err'][0])
    if len(saved_rho_rel_error_array) > 0:
        SchDiode.Project.delete_measurement_datapoints(measurement_id)
    SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                   point_names_long['ic_id'][0], point_names_long['ic_id'][1],
                                   point_names_long['ic_id'][2],
                                   initial_condition_id, measurement_time)
    SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                   point_names_long['Vd'][0], point_names_long['Vd'][1], point_names_long['Vd'][2],
                                   Vd, measurement_time)
    SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                   point_names_long['Vd_err'][0], point_names_long['Vd_err'][1],
                                   point_names_long['Vd_err'][2],
                                   Vd_err, measurement_time)
    SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                   point_names_long['J'][0], point_names_long['J'][1], point_names_long['J'][2],
                                   J, measurement_time)
    SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                   point_names_long['J_err'][0], point_names_long['J_err'][1],
                                   point_names_long['J_err'][2],
                                   J_err, measurement_time)
    for BI in SchDiode.Semiconductor.bonding_interfaces:
        for i, trap in enumerate(BI.dsl_tilt.traps):
            SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                           point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'][0],
                                           point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'][1],
                                           point_names_long[BI.label + '_tilt_' + trap[0].name + '_F'][2],
                                           BI.dsl_tilt_f[i],
                                           measurement_time)
        for i, trap in enumerate(BI.dsl_twist.traps):
            SchDiode.Project.add_datapoint(0, measurement_id, SchDiode.tool_id,
                                           point_names_long[BI.label + '_twist_' + trap[0].name + '_F'][0],
                                           point_names_long[BI.label + '_twist_' + trap[0].name + '_F'][1],
                                           point_names_long[BI.label + '_twist_' + trap[0].name + '_F'][2],
                                           BI.dsl_twist_f[i],
                                           measurement_time)
    point_order_arr = np.arange(len(z_nodes))
    SchDiode.Project.add_datapoint_array(point_order_arr, measurement_id, SchDiode.tool_id,
                                         point_names_long['z'][0], point_names_long['z'][1], point_names_long['z'][2],
                                         z_nodes, measurement_time)
    SchDiode.Project.add_datapoint_array(point_order_arr, measurement_id, SchDiode.tool_id,
                                         point_names_long['Psi'][0], point_names_long['Psi'][1],
                                         point_names_long['Psi'][2],
                                         Psi_points, measurement_time)
    SchDiode.Project.add_datapoint_array(point_order_arr, measurement_id, SchDiode.tool_id,
                                         point_names_long['E'][0], point_names_long['E'][1], point_names_long['E'][2],
                                         E_points, measurement_time)
    SchDiode.Project.add_datapoint_array(point_order_arr, measurement_id, SchDiode.tool_id,
                                         point_names_long['rho_rel_err'][0], point_names_long['rho_rel_err'][1],
                                         point_names_long['rho_rel_err'][2],
                                         rho_rel_err, measurement_time)
    for dopant in SchDiode.Semiconductor.dopants:
        SchDiode.Project.add_datapoint_array(point_order_arr, measurement_id, SchDiode.tool_id,
                                             point_names_long[dopant.name + '_F'][0],
                                             point_names_long[dopant.name + '_F'][1],
                                             point_names_long[dopant.name + '_F'][2],
                                             dopant.F(z_nodes),
                                             measurement_time)
    db_write_elapsed_time = time.time() - db_write_start_time
    if debug: print 'Writing to database took ', db_write_elapsed_time, 's'


def Reccurent_Poisson_solver(SchDiode, Psi=Psi_zero, Vd_guess=None, Vd_error=1e-6,
                             equilibrium_filling=True, fast_traps=None,
                             t=mp.inf, initial_condition_id=-1,
                             rho_rel_err=1e-3, max_iter=100,
                             save_to_db=True, debug=False):
    recurrent_solver_start_time = time.time()
    Va = SchDiode.Va
    kT_eV = to_numeric(k * SchDiode.T / q)
    if save_to_db:
        measurement_id = add_Electric_Field_Measurement(SchDiode, Va, equilibrium_filling, t, initial_condition_id)
        Solution_found, Psi_found, E, z_nodes, rho_rel_err_points, Vd, Vd_err, J, J_err, BI_F, dopants_F, measurement_id = load_diode_state(
            SchDiode, measurement_id, initial_condition_id, debug)
    else:
        Solution_found = False
        measurement_id = -1
    if Solution_found:
        Psi = Psi_found
        Vd_guess = Vd[-1]
        J_tmp = J[-1]
        if max(abs(rho_rel_err_points)) <= rho_rel_err:
            if debug: print '==> Solution satisfy rel_error condition'
            return Psi, E, z_nodes, rho_rel_err_points, Vd[0], Vd_err[0], J[0], J_err[
                0], BI_F, dopants_F, measurement_id
        if debug or 1:
            print '==> Solution found does not satisfy rel_error condition'
            print '==> Recalculating...'
            SchDiode.FieldInversionPoint, PHI_bn, _, PHI_b = SchDiode.get_phi_bn(Psi, Vd_guess, SchottkyEffect=False)
    else:
        if debug: print '==> No solutions found'
        if Vd_guess is None:
            Vd_guess = Va
            PHI_b = abs(SchDiode.V_bi(eV=True))
            SchDiode.FieldInversionPoint = 0
        else:
            SchDiode.FieldInversionPoint, PHI_bn, _, PHI_b = SchDiode.get_phi_bn(Psi, Vd_guess, SchottkyEffect=False)
            if Va < PHI_b:
                Vd_guess = Va
        J_tmp = 0

    converged = False
    current_iter = 0
    Vd_tmp_corr = 0
    Vd_guess_monitor_length = 5
    Vd_guess_monitor = np.arange(Vd_guess_monitor_length, dtype=np.float)
    Vd_guess_monitor_count = 0
    stalled = False
    while not converged and not stalled and current_iter < max_iter:
        if debug: print '\nIteration:', current_iter, '**'
        if debug or 1: print 'Vd_tmp =', Vd_guess
        if debug: print 'PHI_b =', PHI_b
        if Vd_guess >= PHI_b:
            Vd_guess = PHI_b - 1 * kT_eV
        if debug: print 'Vd_tmp =', Vd_guess
        Vd_guess_monitor[Vd_guess_monitor_count] = Vd_guess
        Vd_guess_monitor_count += 1
        if Vd_guess_monitor_count == Vd_guess_monitor_length:
            Vd_guess_monitor_count = 0
        if np.unique(Vd_guess_monitor).size == 1:
            stalled = True
            print 'Solution process stalled on iteration', current_iter, '!!!'
            continue
        nodes_num = int(np.floor(SchDiode.L * 1e6 + 1) * 100)
        nodes, _ = np.linspace(0, SchDiode.L, num=nodes_num + 1, endpoint=True, retstep=True)
        Meshes = Poisson_eq_num_solver_amr(SchDiode, nodes, Psi, Vd_guess, equilibrium_filling, fast_traps,
                                           max_iterations=50, residual_threshold=rho_rel_err,
                                           int_residual_threshold=5e-14,
                                           max_level=5, mesh_refinement_threshold=1e-19, debug=debug)
        z_nodes, Psi_points, rho_rel_err_points = Meshes.flatten(debug=False)
        E_points = -np.gradient(Psi_points, z_nodes, edge_order=2)
        Psi = interp_Fn(z_nodes, Psi_points, interp_type='last')
        E = interp_Fn(z_nodes, E_points, interp_type='last')
        SchDiode.FieldInversionPoint, PHI_bn, _, PHI_b = SchDiode.get_phi_bn(Psi, Vd_guess, SchottkyEffect=False)
        print '*** !!! Inversion pt (calc) !!!', SchDiode.FieldInversionPoint
        if debug: print 'PHI_b, PHI_bn =', PHI_b, PHI_bn
        Vd_tmp_corr, J_tmp_corr = SchDiode.ThermionicEmissionCurrent(Va, PHI_bn, debug=True)
        if debug: print 'V_corr, J =', Vd_tmp_corr, J_tmp_corr
        Vd_err = abs(Vd_guess - Vd_tmp_corr)
        J_err = abs(J_tmp - J_tmp_corr)
        if debug or 1: print 'Vd err =', Vd_err
        if debug or 1: print 'J err =', J_err, '\n'
        if Vd_err < Vd_error:
            converged = True
        else:
            Vd_guess = Vd_tmp_corr
            J_tmp = J_tmp_corr
        current_iter += 1
    Vd = Vd_tmp_corr
    J = J_tmp
    if debug: print 'Calculation converged'
    recurrent_solver_elapsed_time = time.time() - recurrent_solver_start_time
    if debug: print 'Total recurrent solver execution time =', recurrent_solver_elapsed_time, 's\n'
    if save_to_db:
        save_diode_state(SchDiode, measurement_id, initial_condition_id,
                         z_nodes, Psi_points, E_points, rho_rel_err_points,
                         Vd, Vd_err, J, J_err, debug)
    BI_F = {}
    for BI in SchDiode.Semiconductor.bonding_interfaces:
        for i, trap in enumerate(BI.dsl_tilt.traps):
            BI_F[BI.label + '_tilt_' + trap[0].name + '_F'] = BI.dsl_tilt_f[i]
        for i, trap in enumerate(BI.dsl_twist.traps):
            BI_F[BI.label + '_twist_' + trap[0].name + '_F'] = BI.dsl_twist_f[i]
    dopants_F = {}
    for dopant in SchDiode.Semiconductor.dopants:
        dopants_F[dopant.name + '_F'] = dopant.F(z_nodes)
    recurrent_solver_elapsed_time = time.time() - recurrent_solver_start_time
    if debug: print 'Total recurrent solver execution time =', recurrent_solver_elapsed_time, 's'
    return Psi, E, z_nodes, rho_rel_err_points, Vd, Vd_err, J, J_err, BI_F, dopants_F, measurement_id


def Poisson_eq_num_solver(SchDiode, nodes, Psi=Psi_zero, Vd=0,
                          equilibrium_filling=True, fast_traps=None,
                          bc1=None, bc2=0, max_iterations=1000, threshold=1e-7,
                          debug=False):
    if fast_traps is None:
        fast_traps = []
    type_sign = -1 if SchDiode.Semiconductor.dop_type == 'n' else 1
    if bc1 is None:
        bc1 = -(SchDiode.V_bi(eV=True) + type_sign * Vd)
    Eps0Eps = epsilon0 * SchDiode.Semiconductor.reference['epsilon']
    mesh = UniformMesh1D(nodes[0], nodes[-1], nodes[1] - nodes[0], bc1, bc2)
    # print nodes[0], nodes[-1], (nodes[1] - nodes[0])*1e6
    mesh.int_residual = threshold + 1
    int_residual_array = []
    DPsi = np.ones_like(mesh.local_nodes)
    if debug:
        plt.ion()
        _, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5)
        ax1.set_autoscaley_on(True)
        ax2.set_autoscaley_on(True)
        ax3.set_autoscaley_on(True)
        ax4.set_autoscaley_on(True)
        ax4.grid(1)
        ax5.set_autoscalex_on(True)
        ax5.set_autoscaley_on(True)
        # print 'here 3'
        # print mesh.phys_nodes()[[0,-1]]
        Psi_line, = ax1.plot(mesh.local_nodes, Psi(mesh.phys_nodes()))
        DPsi_line, = ax2.plot(mesh.local_nodes, DPsi)
        f_line, = ax3.plot(mesh.local_nodes, DPsi, 'b-')
        d2Psi_line, = ax3.plot(mesh.local_nodes, DPsi, 'r-')
        E_line, = ax4.plot(mesh.local_nodes, DPsi)
        R_line, = ax5.plot(int_residual_array)
        plt.draw()
    i = 1
    int_res_diff_counter = 1
    W = 1
    while abs(mesh.int_residual) > threshold and i < max_iterations and np.max(abs(DPsi)) > 4 * np.finfo(np.float).eps:
        if debug: print 'Iteration:', i
        # time.sleep(1)
        if equilibrium_filling:
            SchDiode.dopants_set_eq_F(nodes, Psi, debug=False)
            for BI in SchDiode.Semiconductor.bonding_interfaces:
                SchDiode.BI_set_eq_F(BI, Psi, eV=False, debug=False)
        elif len(fast_traps) > 0:
            fermi_level = SchDiode.EfEc(Psi, nodes, eV=False)
            for dopant in SchDiode.Semiconductor.dopants:
                if dopant.name in fast_traps:
                    #print 'Setting dopant', dopant.name, 'occupation to equilibrium value'
                    dopant_f = dopant.equilibrium_f(SchDiode.T, SchDiode.Semiconductor, fermi_level,
                                                    electron_volts=False, debug=False)
                    dopant_df = dopant.d_equilibrium_f_d_fermi_energy(SchDiode.T, SchDiode.Semiconductor, fermi_level,
                                                                      electron_volts=False, debug=debug)
                    dopant.set_F_interp(nodes, dopant_f)
                    dopant.set_dF_interp(nodes, dopant_df)
            for bi in SchDiode.Semiconductor.bonding_interfaces:
                # continue
                fermi_level_at_bi = SchDiode.EfEc(Psi, bi.depth, eV=False)
                bi_tilt_f = bi.dsl_tilt_f
                bi_tilt_df = bi.dsl_tilt_df
                for i, trap in enumerate(bi.dsl_tilt.traps):
                    trap_name = bi.label + '_tilt_' + trap[0].name
                    if trap_name in fast_traps:
                        #print 'Setting trap', trap_name, 'occupation to equilibrium value'
                        bi_tilt_f[i] = trap[0].equilibrium_f(SchDiode.T, SchDiode.Semiconductor, fermi_level_at_bi,
                                                             electron_volts=False, debug=False)
                        bi_tilt_df[i] = trap[0].d_equilibrium_f_d_fermi_energy(SchDiode.T, SchDiode.Semiconductor,
                                                                               fermi_level_at_bi,
                                                                               electron_volts=False, debug=False)
                bi_twist_f = bi.dsl_twist_f
                bi_twist_df = bi.dsl_twist_df
                for i, trap in enumerate(bi.dsl_twist.traps):
                    trap_name = bi.label + '_twist_' + trap[0].name
                    if trap_name in fast_traps:
                        print 'Setting trap', trap_name, 'occupation to equilibrium value'
                        bi_twist_f[i] = trap[0].equilibrium_f(SchDiode.T, SchDiode.Semiconductor, fermi_level_at_bi,
                                                              electron_volts=False, debug=False)
                        bi_twist_df[i] = trap[0].d_equilibrium_f_d_fermi_energy(SchDiode.T, SchDiode.Semiconductor,
                                                                                fermi_level_at_bi,
                                                                                electron_volts=False, debug=False)
                bi.set_traps_f(bi_tilt_f, bi_twist_f)
                bi.set_traps_df(bi_tilt_df, bi_twist_df)

        rho_z_Psi = SchDiode.build_rho_z_Psi(Psi, carriers_charge=True)

        def rho_z_Psi_eps(z, Psi):
            return -rho_z_Psi(z, Psi) / Eps0Eps + 1e-200

        def d_rho_dDPsi(z, Psi):
            dn_nodes, dp_nodes = SchDiode.dn_carriers_dEf_theory(Psi, z)
            if SchDiode.Semiconductor.dop_type == 'n':
                dp_nodes = np.zeros_like(z)
            elif SchDiode.Semiconductor.dop_type == 'p':
                dn_nodes = np.zeros_like(z)
            d_carriers_nodes = dn_nodes - dp_nodes
            dN_dopants_nodes = np.zeros_like(z)
            for dopant in SchDiode.Semiconductor.dopants:
                dN_dopants_nodes += dopant.concentration(z) * (
                    dopant.charge_states[1][0] - dopant.charge_states[0][0]) * dopant.dF(z)
            dN_BI_nodes = np.zeros_like(z)
            for BI in SchDiode.Semiconductor.bonding_interfaces:
                dN_BI_nodes += smooth_dd(z - BI.depth, BI.smooth_dirac_epsilon) * BI.d_density_of_charge
            return (d_carriers_nodes - dN_dopants_nodes - dN_BI_nodes) * to_numeric(q ** 2) / Eps0Eps

        mesh, Psi, DPsi = dirichlet_non_linear_poisson_solver_mesh(mesh, Psi, rho_z_Psi_eps, d_rho_dDPsi, rel=True, W=W,
                                                                   debug=False)
        int_residual_array.append(mesh.int_residual)
        if int_res_diff_counter > 7:
            int_res_diff = np.array(int_residual_array[-5:]) - np.array(int_residual_array[-6:-1])
            int_res_diff_pos = np.where(int_res_diff >= 0)[0]
            # int_res_diff_neg = np.where(int_res_diff < 0)[0]
            if int_res_diff_pos.size > 1:
                W /= 2
                print 'GENERATION: adding relaxation factor:', W
                int_res_diff_counter = 1
        if debug: print 'Integrated residual:', mesh.int_residual
        if debug:
            Psi_line.set_ydata(mesh.solution)
            f_line.set_ydata(rho_z_Psi_eps(mesh.phys_nodes(), Psi))
            dPsi = np.gradient(mesh.solution, mesh.phys_nodes(), edge_order=2)
            d2Psi = np.gradient(dPsi, mesh.phys_nodes(), edge_order=2)
            DPsi_line.set_ydata(DPsi)
            d2Psi_line.set_ydata(d2Psi)
            E_line.set_ydata(-dPsi)
            R_line.set_data(np.arange(i), int_residual_array)
            R_line.set_ydata(int_residual_array)
            ax1.relim()
            ax1.autoscale_view()
            ax2.relim()
            ax2.autoscale_view()
            ax3.relim()
            ax3.autoscale_view()
            ax4.relim()
            ax4.autoscale_view()
            ax5.relim()
            ax5.autoscale_view()
            plt.draw()
            #time.sleep(10)
        i += 1
        int_res_diff_counter += 1
    if debug: plt.ioff()
    return mesh, Psi


def Poisson_eq_num_solver_mesh(SchDiode, mesh, Psi=Psi_zero, Vd=0,
                               equilibrium_filling=True, fast_traps=None,
                               max_iterations=1000, threshold=1e-7,
                               debug=False):
    mesh, Psi = Poisson_eq_num_solver(SchDiode, mesh.phys_nodes(), Psi, Vd, equilibrium_filling, fast_traps,
                                      mesh.bc1, mesh.bc2,
                                      max_iterations, threshold, debug)
    return mesh, Psi


def Poisson_eq_num_solver_amr(SchDiode, nodes, Psi=Psi_zero, Vd=0,
                              equilibrium_filling=True, fast_traps=None,
                              max_iterations=1000, residual_threshold=1e-3, int_residual_threshold=1e-6,
                              max_level=20, mesh_refinement_threshold=1e-12, debug=False):
    '''
    The recurrent NL Poisson solver with the Adaptive Mesh Refinement
    nodes is the initial physical mesh
    '''
    type_sign = -1 if SchDiode.Semiconductor.dop_type == 'n' else 1
    bc1 = -(SchDiode.V_bi(eV=True) + type_sign * Vd)
    bc2 = 0
    additional_length = SchDiode.L / 5
    additional_nodes = int(additional_length / (nodes[1] - nodes[0]))
    root_mesh_stop = nodes[-1] + additional_nodes * (nodes[1] - nodes[0])
    # print nodes*1e6
    root_mesh = UniformMesh1D(nodes[0], root_mesh_stop, nodes[1] - nodes[0], bc1, bc2)
    # print root_mesh.phys_nodes()[-1]*1e6
    # _, ax = plt.subplots()
    # ax.plot(root_mesh.phys_nodes()*1e6, Psi(root_mesh.phys_nodes()))
    # plt.show()
    # return
    Meshes = Uniform1DMeshesTree(root_mesh, refinement_coefficient=2, aligned=True, crop=[0, additional_nodes])
    converged = np.zeros(1)
    level = 0
    while (not converged.all() or level < Meshes.levels[-1]) and level <= max_level:
        if debug or 1: print '\nSolving for Meshes of level:', level
        converged = np.zeros(len(Meshes.Tree[level]))
        refinement_meshes = []
        for mesh_id, mesh in enumerate(Meshes.Tree[level]):
            # print 'here 2'
            mesh, Psi = Poisson_eq_num_solver_mesh(SchDiode, mesh, Psi, Vd, equilibrium_filling, fast_traps,
                                                   max_iterations, int_residual_threshold, debug)
            mesh.trim()
            # mesh.residual[-1] = 0
            Meshes.Tree[level][mesh_id] = mesh

            if max(abs(mesh.residual)) < residual_threshold:
                if debug: print 'CONVERGED!'
                converged[mesh_id] = True
                continue
            refinement_points_chunks = points_for_refinement(mesh, mesh_refinement_threshold)
            if len(refinement_points_chunks) == 0 or np.all(
                    np.array([block.size == 0 for block in refinement_points_chunks])):
                if debug: print 'CONVERGED!'
                converged[mesh_id] = True
                continue
            if level < max_level:
                if debug: print 'nodes for refinement:', refinement_points_chunks
                for block in refinement_points_chunks:
                    idx1, idx2, crop = adjust_range(block, mesh.num - 1, crop=[50, 10], step_scale=2)
                    # print idx1, idx2, type(idx1), type(idx2)
                    start_point = mesh.to_phys(mesh.local_nodes[idx1])
                    stop_point = mesh.to_phys(mesh.local_nodes[idx2])
                    ref_bc1 = mesh.solution[idx1]
                    ref_bc2 = mesh.solution[idx2]
                    refinement_mesh = UniformMesh1D(start_point, stop_point,
                                                    mesh.phys_step / Meshes.refinement_coefficient, ref_bc1, ref_bc2,
                                                    crop=crop)
                    if debug: print 'NEW Refinement mesh from', start_point * 1e6, 'to', stop_point * 1e6, 'phys_step', mesh.phys_step / Meshes.refinement_coefficient * 1e6
                    if debug: print 'NEW Refinement mesh:', refinement_mesh.phys_nodes()[[0, -1]]
                    if debug: print 'BC:', ref_bc1, ref_bc2
                    refinement_meshes.append(refinement_mesh)
                if debug:
                    _, ax = plt.subplots(1)
                    ax.plot(mesh.phys_nodes(), mesh.solution, 'b-o')
                    # ax.plot(refinement_mesh.phys_nodes(), Psi(refinement_mesh.phys_nodes()), 'r-o')
                    plt.show()
        flat_grid, flat_sol, flat_res = Meshes.flatten(debug=False)
        Psi = interp1d(flat_grid, flat_sol)
        if len(refinement_meshes) > 0:
            print 'adding refinement meshes for level', level
            for refinement_mesh in refinement_meshes:
                Meshes.add_mesh(refinement_mesh)
        # if debug: Meshes.plot_tree()
        if debug:
            _, (ax1, ax2, ax3) = plt.subplots(3)
            ax1.plot(flat_grid * 1e6, flat_sol, 'b-')
            ax2.plot(flat_grid * 1e6, flat_res, 'b-')
            # ax.plot(refinement_mesh.phys_nodes(), Psi(refinement_mesh.phys_nodes()), 'r-o')
            Meshes.plot_tree(ax3)
            plt.show()

        level += 1
    if debug: print 'Mesh tree has ', Meshes.levels[-1], 'refinement levels'
    Meshes.trim(debug=True)
    return Meshes
