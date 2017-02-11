__author__ = 'anton'

import numpy as np
from matplotlib import pyplot as plt

from Schottky.Notation import q, k
from Schottky.Helpers import to_numeric
from Schottky.Diode import Poisson, Visual


def add_kinetics_measurement(schottky_diode, initial_condition, fast_traps,
                             delta_t_min, delta_t_max, rho_rel_error, df_threshold):
    if initial_condition < 0:
        raise Exception('Kinetics experiments require initial condition')
    measurement_name = 'Kinetics for ic %d ' % initial_condition
    measurement_name += '[dt_min=%2.2g, dt_max=%2.2g, ' % (delta_t_min, delta_t_max)
    measurement_name += 'rho_rel_err=%2.2f, dF_threshold=%2.2g], fast_traps=[' % (rho_rel_error, df_threshold)
    measurement_description = 'Traps kinetics measurement for ic_id %d, fast traps=[' % initial_condition
    for i, fast_trap in enumerate(fast_traps):
        measurement_name += fast_trap
        measurement_description += fast_trap
        if i < len(fast_traps) - 1:
            measurement_name += ', '
            measurement_description += ', '
    measurement_name += ']'
    measurement_description += ']'
    measurement_type = schottky_diode.measurements_types['Traps kinetics measurement']
    measurement_id = schottky_diode.Project.add_measurement(measurement_type, measurement_name, measurement_description,
                                                            measurement_datetime=None, create_new=False)
    return measurement_id


def gen_point_names_kinetics_measurement(schottky_diode):
    point_names_long = {
        'dt_min': ['delta t min', 'dt_min', 's'],
        'dt_max': ['delta t max', 'dt_max', 's'],
        'rho_rel_err': ['Charge density relative error', 'rho_rel_err', ''],
        'dF_threshold': ['delta F threshold', 'dF_threshold', ''],
        't': ['time', 't', 's'],
        'ic_id': ['initial condition id', 'ic_id', ''],
        'Vd': ['Voltage bias', 'Vd', 'V'],
        'J': ['Current density', 'J', 'A/m'],
    }

    for dopant in schottky_diode.Semiconductor.dopants:
        dopant_key = dopant.name + '_F'
        point_names_long[dopant_key] = [dopant_key, 'F', '']

    for bi in schottky_diode.Semiconductor.bonding_interfaces:
        localized_trap_key = bi.label + '_tilt_' + trap[0].name + '_F'
        for trap in bi.dsl_tilt.traps:
            point_names_long[localized_trap_key] = [localized_trap_key, 'F', '']
        for trap in bi.dsl_twist.traps:
            localized_trap_key = bi.label + '_twist_' + trap[0].name + '_F'
            point_names_long[localized_trap_key] = [localized_trap_key, 'F', '']

    return point_names_long


def load_kinetics_data(schottky_diode, measurement_id, debug=False):
    results = {}
    point_names_long = gen_point_names_kinetics_measurement(schottky_diode)
    if debug:
        print '==> Looking for saved solution in the database'
    points_names = [point_names_long['dt_min'][0], point_names_long['dt_max'][0],
                    point_names_long['rho_rel_err'][0], point_names_long['dF_threshold'][0],
                    point_names_long['t'][0], point_names_long['ic_id'][0],
                    point_names_long['Vd'][0], point_names_long['J'][0]]
    for dopant in schottky_diode.Semiconductor.dopants:
        points_names.append(point_names_long[dopant.name + '_F'][0])
    for bi in schottky_diode.Semiconductor.bonding_interfaces:
        for trap in bi.dsl_tilt.traps:
            points_names.append(point_names_long[bi.label + '_tilt_' + trap[0].name + '_F'][0])
        for trap in bi.dsl_twist.traps:
            points_names.append(point_names_long[bi.label + '_twist_' + trap[0].name + '_F'][0])
    try:
        data = schottky_diode.Project.get_data_points_by_names(measurement_id, points_names)
    except:
        if debug: print '==> No solutions found'
        return False, 0, 0, 0, 0, results
    dt_min = data[points_names[0]][:, 0]
    dt_max = data[points_names[1]][:, 0]
    rho_rel_error = data[points_names[2]][:, 0]
    df_threshold = data[points_names[3]][:, 0]
    results['t_points'] = data[points_names[4]][:, 0]
    results['ic_id'] = data[points_names[5]][:, 0]
    results['Vd'] = data[points_names[6]][:, 0]
    results['J'] = data[points_names[7]][:, 0]
    if debug:
        print '==> Solution found'
    for dopant in schottky_diode.Semiconductor.dopants:
        dopant_key = dopant.name + '_F'
        results[dopant_key] = data[point_names_long[dopant_key][0]][:, 0]
    for bi in schottky_diode.Semiconductor.bonding_interfaces:
        for trap in bi.dsl_tilt.traps:
            localized_trap_key = bi.label + '_tilt_' + trap[0].name + '_F'
            results[localized_trap_key] = data[point_names_long[localized_trap_key][0]][:, 0]
        for trap in bi.dsl_twist.traps:
            localized_trap_key = bi.label + '_twist_' + trap[0].name + '_F'
            results[localized_trap_key] = data[point_names_long[localized_trap_key][0]][:, 0]
    return True, dt_min, dt_max, rho_rel_error, df_threshold, results


def dopants_df_dt(schottky_diode, initial_condition_id):
    ic_found, potential, field, z_nodes, _, \
        diode_voltage_drop, _, current_density, _, bonding_interfaces_f, dopants_f, \
        measurement_id = Poisson.load_diode_state(schottky_diode, initial_condition_id, debug=True)
    if not ic_found:
        raise(Exception('Initial condition not found. Please check everything'))
    #schottky_diode.dopants_set_eq_F(z_nodes, potential, debug=False)
    _, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True)
    ax1.plot(z_nodes, -potential(z_nodes))
    for dopant_name in dopants_f.keys():
        ax2.plot(z_nodes, dopants_f[dopant_name])
    fermi_level = schottky_diode.EfEc(potential, z_nodes, eV=False)
    n, p = schottky_diode.n_carriers_theory(potential, z_nodes)
    if schottky_diode.Semiconductor.dop_type == 'n':
        p = np.zeros_like(p)
    elif schottky_diode.Semiconductor.dop_type == 'p':
        n = np.zeros_like(n)
    ax4.plot(z_nodes, n)
    for dopant in schottky_diode.Semiconductor.dopants:
        dopant_f = dopant.equilibrium_f(schottky_diode.T, schottky_diode.Semiconductor, fermi_level,
                                        electron_volts=False, debug=False)
        ax2.plot(z_nodes, dopant_f)
        df_dt, tau = dopant.df_dt(schottky_diode.T, schottky_diode.Semiconductor, dopants_f[dopant_name], n, p,
                                  barrier_lowering_e=None, barrier_lowering_h=None, use_mpmath=False, debug=False)
        capture_e, capture_h, capture_tau_e, capture_tau_h = dopant.capture_rate(schottky_diode.T,
                                                                                 schottky_diode.Semiconductor,
                                                                                 dopants_f[dopant_name],
                                                                                 n, p, debug=False)
        emission_e, emission_h, emission_tau_e, emission_tau_h = dopant.emission_rate(schottky_diode.T,
                                                                                      schottky_diode.Semiconductor,
                                                                                      dopants_f[dopant_name],
                                                                                      barrier_lowering_e=None,
                                                                                      barrier_lowering_h=None,
                                                                                      use_mpmath=False, debug=False)
        print 'electrons: capture tau = %2.2g s, emission tau = %2.2g s' % (min(capture_tau_e), min(emission_tau_e))
        print 'holes: capture tau = %2.2g s, emission tau = %2.2g s' % (min(capture_tau_h), min(emission_tau_h))
        print 'tau = %2.2g s' % tau
        #ax3.plot(z_nodes, emission_h - emission_e)
        ax3.plot(z_nodes, df_dt)
        #ax5.plot(z_nodes, capture_tau_e)
    plt.show()


def traps_kinetics(schottky_diode, initial_condition_id, delta_t_min, delta_t_max, t_stop, fast_traps=None,
                   rho_rel_err=1e-1, df_threshold=1e-3, debug=False, debug_plot=False, z_limit=None):
    if z_limit is None:
        z_limit = 1e8
    t_points = []
    potential_t = []
    field_d = []
    z_t = []
    diode_voltage_drop_t = []
    current_density_t = []
    bonding_interfaces_f_t = []
    dopants_f_t = []

    ic_found, potential, z_nodes, _, _, _, _, _, _, _, _, _ = Poisson.load_diode_state(schottky_diode,
                                                                                       initial_condition_id,
                                                                                       debug=True)
    if not ic_found:
        raise(Exception('Initial condition not found. Please check everything'))

    if debug_plot:
        plt.ion()
        axes = Visual.prepare_debug_axes(['Potential', 'Dopants', 'Localized traps'])
        potential_lines = Visual.create_lines(axes['Potential'], ['Start potential', 'Current potential'])
        dopants_lines_names = [dopant.name + '_F' for dopant in schottky_diode.Semiconductor.dopants]
        localized_traps_lines_names = []
        for bi in schottky_diode.Semiconductor.bonding_interfaces:
            localized_traps_lines_names += [bi.label + '_tilt_' + trap[0].name + '_F' for trap in bi.dsl_tilt.traps]
            localized_traps_lines_names += [bi.label + '_twist_' + trap[0].name + '_F' for trap in bi.dsl_twist.traps]
        dopants_lines = Visual.create_lines(axes['Dopants'], dopants_lines_names)
        localized_traps_lines = Visual.create_lines(axes['Localized traps'], localized_traps_lines_names)

    t = 0
    last_state_id = initial_condition_id
    if fast_traps is None:
        fast_traps = []
    while t <= t_stop:
        print '\n\nt =', t
        potential, field, z_nodes, _, \
            diode_voltage_drop, _, current_density, _, \
            bonding_interface_f, dopants_f, \
            last_state_id = Poisson.Reccurent_Poisson_solver(schottky_diode, potential,
                                                             equilibrium_filling=False,
                                                             fast_traps=fast_traps,
                                                             t=t, initial_condition_id=last_state_id,
                                                             rho_rel_err=rho_rel_err, max_iter=100, debug=False)

        t_points.append(t)
        potential_t.append(potential)
        field_d.append(field)
        z_t.append(z_nodes)
        diode_voltage_drop_t.append(diode_voltage_drop)
        current_density_t.append(current_density)
        bonding_interfaces_f_t.append(bonding_interface_f.copy())
        dopants_f_t.append(dopants_f.copy())

        fermi_level = schottky_diode.EfEc(potential, z_nodes, eV=False)

        z_limit_idx = np.where(z_nodes < z_limit)

        if debug_plot:
            if t == 0:
                potential_lines['Start potential'].set_data(z_nodes * 1e6, -potential(z_nodes))
            potential_lines['Current potential'].set_data(z_nodes * 1e6, -potential(z_nodes))

            for dopant_line_name in dopants_lines_names:
                dopants_lines[dopant_line_name].set_data(z_nodes * 1e6, dopants_f[dopant_line_name])

            for localized_trap_line_name in localized_traps_lines_names:
                localized_trap_f_t = []
                for bonding_interfaces_f_t_i in bonding_interfaces_f_t:
                    localized_trap_f_t.append(bonding_interfaces_f_t_i[localized_trap_line_name])
                localized_traps_lines[localized_trap_line_name].set_data(t_points, localized_trap_f_t)
            Visual.autoscale_axes(axes)
            plt.draw()

        dt = delta_t_max if t_stop - t > delta_t_max else t_stop - t
        if dt == 0:
            break
        if debug:
            print '\n\nT = %2.2f K, t = %2.2g s, dt = %2.2g s' % (schottky_diode.T, t, dt)

        n, p = schottky_diode.n_carriers_theory(potential, z_nodes)
        if schottky_diode.Semiconductor.dop_type == 'n':
            p = np.zeros_like(p)
        elif schottky_diode.Semiconductor.dop_type == 'p':
            n = np.zeros_like(n)

        #fast_traps = []

        dopants_skip_list = []
        df_dopants = {}
        for dopant in schottky_diode.Semiconductor.dopants:
            dopant_key = dopant.name + '_F'
            if dopant.name in fast_traps:
                if debug:
                    print '\nDopant:', dopant.name
                    print 'This dopant is in a fast-traps list, skipping.'
                dopants_skip_list.append(dopant_key)
                continue
            poole_frenkel_e = np.ones_like(z_nodes, dtype=np.float)
            poole_frenkel_h = np.ones_like(z_nodes, dtype=np.float)
            barrier_lowering_e = np.zeros_like(n, dtype=np.float)
            barrier_lowering_h = np.zeros_like(p, dtype=np.float)
            field_z = field(z_nodes)
            if dopant.trap_potential is not None:
                kT = to_numeric(k * schottky_diode.T / q)
                theta_points = 100
                theta = np.linspace(0, np.pi, num=theta_points, endpoint=True)
                barrier_lowering = np.zeros((theta_points, len(z_nodes)), dtype=np.float)
                max_N_l = dopant.trap_potential.get_potential_by_name('Charged Dislocation')\
                        .max_linear_charge_density
                dsl_charge_density = max_N_l * dopant.F(z_nodes)
                for z_num, local_electric_field in enumerate(field_z):
                    #print z_num, 'of', len(z_nodes)
                    #dsl_charge_density = max_N_l * dopant.F(z_nodes[z_num])
                    local_electric_field_r = abs(local_electric_field)
                    local_electric_field_theta = 0 if local_electric_field >= 0 else np.pi
                    local_electric_field_3d = (local_electric_field_r, local_electric_field_theta, 0.0)
                    dopant.trap_potential.get_potential_by_name('Charged Dislocation')\
                        .set_linear_charge_density(dsl_charge_density[z_num])
                    dopant.trap_potential.get_potential_by_name('External Field')\
                        .external_field = local_electric_field_3d

                    loc_f = local_electric_field_r
                    loc_a = dopant.trap_potential.get_potential_by_name('Charged Dislocation').a
                    loc_b = dopant.trap_potential.get_potential_by_name('Deformation').a
                    r0 = np.zeros_like(theta)
                    if np.allclose(loc_f, 0.0):
                        print 'here', loc_f
                        r0 = loc_b / loc_a
                    else:
                        idx = np.where(loc_a**2 + 4* loc_b * loc_f * np.cos(theta) >= 0)
                        r0[idx] = (np.sqrt(loc_a**2 + 4* loc_b * loc_f * np.cos(theta[idx])) - loc_a) / (2 * loc_f * np.cos(theta[idx]))
                        zero_theta_idx = np.where(np.allclose(np.cos(theta), 0))
                        r0[zero_theta_idx] = loc_b / loc_a
                    non_zero_r_idx = np.where(not np.allclose(r0, 0.0))
                    bl_grid = dopant.trap_potential.potential(r0[non_zero_r_idx], theta[non_zero_r_idx], 0)
                    #print np.rad2deg(theta[non_zero_r_idx])
                    #print bl_grid[0,:,0].shape, theta.shape
                    #print bl_grid[0,:,0]
                    bl_flat = np.zeros_like(theta)
                    if loc_f > -1e3:
                        bl_flat[idx] = bl_grid[:,0,0]
                    #barrier_lowering[:,z_num] = np.array([dopant.trap_potential.barrier_lowering(theta_i)[0] for theta_i in theta])
                    barrier_lowering[:, z_num] = bl_flat
                    #print barrier_lowering[:, z_num]
                    #print bl_flat - barrier_lowering[:, z_num]
                    #poole_frenkel = 0.5 * np.trapz(np.sin(theta) * np.exp(abs(barrier_lowering[:, 0]) / kT), theta)
                poole_frenkel = 0.5 * np.trapz(np.exp(abs(barrier_lowering) / kT), theta, axis=0)
                if np.sum(barrier_lowering[:, 0]) < 0:
                    poole_frenkel_e = poole_frenkel
                    #print 'emission boost e:', poole_frenkel
                elif np.sum(barrier_lowering[:, 0]) > 0:
                    poole_frenkel_h = poole_frenkel
                    #print 'emission boost h:', poole_frenkel



            df_dt, tau = dopant.df_dt(schottky_diode.T, schottky_diode.Semiconductor, dopants_f[dopant_key], n, p,
                                      poole_frenkel_e=poole_frenkel_e,
                                      poole_frenkel_h=poole_frenkel_h,
                                      barrier_lowering_e=barrier_lowering_e,
                                      barrier_lowering_h=barrier_lowering_h,
                                      use_mpmath=False, debug=False)
            df_dopants[dopant_key] = df_dt
            df_total = np.sum(df_dt[z_limit_idx]) / np.sum(dopants_f_t[0][dopant_key][z_limit_idx])
            #max_dt = df_threshold / np.max(np.abs(df_dt))
            max_dt = df_threshold / np.max(np.abs(df_total))
            if debug:
                print '\nDopant:', dopant.name
                print 'Z limit of %2.2g m: left %d points of %d' % (z_limit, len(z_nodes[z_limit_idx]), len(z_nodes))
                print 'Min time constant %2.2g s' % tau
                #print 'Max dF:', np.max(np.abs(df_dt)), 'th:', df_threshold
                print 'Max dF:', np.max(np.abs(df_total)), 'th:', df_threshold
                print 'Max dt:', max_dt, 'dt:', dt
            if dt > max_dt > delta_t_min:
                if debug:
                    print 'Setting dt to', max_dt
                dt = max_dt
            elif max_dt < delta_t_min:
                if debug:
                    print 'Traps are too fast. Setting dopant occupation to equilibrium value'
                dopant_f = dopant.equilibrium_f(schottky_diode.T, schottky_diode.Semiconductor, fermi_level,
                                                electron_volts=False, debug=False)
                dopant.set_F_interp(z_nodes, dopant_f)
                dopant.set_dF_interp(z_nodes, np.zeros_like(dopant_f))
                fast_traps.append(dopant.name)
                dopants_skip_list.append(dopant_key)

        localized_traps_skip_list = []
        df_bonding_interfaces = {}
        for bi in schottky_diode.Semiconductor.bonding_interfaces:
            fermi_level_at_bi = schottky_diode.EfEc(potential, bi.depth, eV=False)
            electric_field_at_bi = field(bi.depth)
            n_bi, p_bi = schottky_diode.n_carriers_theory(potential, bi.depth)
            if schottky_diode.Semiconductor.dop_type == 'n':
                p_bi = 0.0
            elif schottky_diode.Semiconductor.dop_type == 'p':
                n_bi = 0.0
            if debug:
                print '\nBI:', bi.label
                print 'n_BI = %2.4g' % n_bi
                print 'p_BI = %2.4g' % p_bi

            bi_tilt_f = []
            for trap_idx, trap in enumerate(bi.dsl_tilt.traps):
                localized_trap_key = bi.label + '_tilt_' + trap[0].name + '_F'
                print 'BI Density of Charge = %2.2g cm-2' % (bi.density_of_charge / 1e4)
                print 'TILT F = %2.2f' % bi.dsl_tilt_f[trap_idx]
                dsl_charge_density = bi.dsl_tilt_f[trap_idx] * trap[1]
                print 'TILT Density of Charge = %2.2g cm-1' % (dsl_charge_density / 1e2)
                print 'EXT FIELD = %2.2g V*cm' % (electric_field_at_bi / 1e2)
                electric_field_at_bi_r = abs(electric_field_at_bi)
                electric_field_at_bi_theta = 0 if electric_field_at_bi >= 0 else np.pi
                electric_field_at_bi_3d = (electric_field_at_bi_r, electric_field_at_bi_theta, 0.0)
                trap[0].trap_potential.get_potential_by_name('Charged Dislocation').set_linear_charge_density(dsl_charge_density)
                trap[0].trap_potential.get_potential_by_name('External Field').external_field = electric_field_at_bi_3d

                kT = to_numeric(k * schottky_diode.T / q)
                theta = np.linspace(0, np.pi, num=100, endpoint=True)
                barrier_lowering = np.array([trap[0].trap_potential.barrier_lowering(theta_i) for theta_i in theta])
                poole_frenkel = 0.5 * np.trapz(np.sin(theta) * np.exp(abs(barrier_lowering[:, 0]) / kT), theta)
                poole_frenkel_e = 1.0
                poole_frenkel_h = 1.0
                if np.sum(barrier_lowering[:, 0]) < 0:
                    poole_frenkel_e = poole_frenkel
                    print 'emission boost e: %2.4g' % poole_frenkel
                elif np.sum(barrier_lowering[:, 0]) > 0:
                    poole_frenkel_h = poole_frenkel
                    print 'emission boost h: %2.4g' % poole_frenkel
                barrier_lowering_e = 0.0
                barrier_lowering_h = 0.0
                df_dt, tau = trap[0].df_dt(schottky_diode.T, schottky_diode.Semiconductor,
                                           bonding_interface_f[localized_trap_key], n_bi, p_bi,
                                           poole_frenkel_e=poole_frenkel_e,
                                           poole_frenkel_h=poole_frenkel_h,
                                           barrier_lowering_e=barrier_lowering_e,
                                           barrier_lowering_h=barrier_lowering_h,
                                           use_mpmath=False, debug=False)
                df_bonding_interfaces[localized_trap_key] = df_dt
                try:
                    max_dt = df_threshold / abs(df_dt)
                except ZeroDivisionError:
                    max_dt = 1e250
                if debug:
                    print '\nTrap:', trap[0].name
                    print 'time constant %2.2g s' % tau
                    print 'dF: %2.4g, th: %2.2f'% (df_dt, df_threshold)
                    print 'Max dt:', max_dt, 'dt:', dt
                if dt > max_dt > delta_t_min:
                    if debug:
                        print 'Setting dt to', max_dt
                    dt = max_dt
                elif max_dt < delta_t_min:
                    if dt > 10 * max_dt:
                        if debug:
                            print 'Trap is too fast. Setting trap occupation to equilibrium value'
                        bonding_interface_f[localized_trap_key] = trap[0].equilibrium_f(schottky_diode.T,
                                                                                        schottky_diode.Semiconductor,
                                                                                        fermi_level_at_bi,
                                                                                        electron_volts=False, debug=False)
                        localized_traps_skip_list.append(localized_trap_key)
                    else:
                        if debug:
                            print 'Setting dt to', max_dt
                        dt = max_dt
                        # fast_traps.append(bi.label + '_tilt_' + trap[0].name)
                bi_tilt_f.append(bonding_interface_f[localized_trap_key])
            bi_twist_f = []
            for trap in bi.dsl_twist.traps:
                localized_trap_key = bi.label + '_twist_' + trap[0].name + '_F'
                trap[0].trap_potential.get_potential_by_name('External Field').external_field = electric_field_at_bi
                print 'EXT FIELD = %2.2g V*cm' % (electric_field_at_bi / 1e2)
                print trap[0].trap_potential.barrier_lowering()
                barrier_lowering_e = 0.0
                barrier_lowering_h = 0.0
                df_dt, tau = trap[0].df_dt(schottky_diode.T, schottky_diode.Semiconductor,
                                           bonding_interface_f[localized_trap_key], n_bi, p_bi,
                                           barrier_lowering_e=barrier_lowering_e,
                                           barrier_lowering_h=barrier_lowering_h,
                                           use_mpmath=False, debug=False)
                df_bonding_interfaces[localized_trap_key] = df_dt
                try:
                    max_dt = df_threshold / abs(df_dt)
                except ZeroDivisionError:
                    max_dt = 1e250
                if debug:
                    print '\nTrap:', trap[0].name
                    print 'time constant %2.2g s' % tau
                    print 'dF: %2.4g, th: %2.2f'% (df_dt, df_threshold)
                    print 'Max dt:', max_dt, 'dt:', dt
                if dt > max_dt > delta_t_min:
                    if debug:
                        print 'Setting dt to', max_dt
                    dt = max_dt
                elif max_dt < delta_t_min:
                    if dt > 10 * max_dt:
                        if debug:
                            print 'Trap is too fast. Setting trap occupation to equilibrium value'
                        bonding_interface_f[localized_trap_key] = trap[0].equilibrium_f(schottky_diode.T,
                                                                                        schottky_diode.Semiconductor,
                                                                                        fermi_level,
                                                                                        electron_volts=False, debug=False)
                        localized_traps_skip_list.append(localized_trap_key)
                    else:
                        if debug:
                            print 'Setting dt to', max_dt
                        dt = max_dt
                    # fast_traps.append(bi.label + '_twist_' + trap[0].name)
                bi_twist_f.append(bonding_interface_f[localized_trap_key])
            bi.set_traps_f(np.array(bi_tilt_f), np.array(bi_twist_f))
            bi.set_traps_df(np.zeros(len(bi_tilt_f)), np.zeros(len(bi_twist_f)))
        dt = np.float(dt)

        for dopant in schottky_diode.Semiconductor.dopants:
            dopant_key = dopant.name + '_F'
            if dopant_key in dopants_skip_list:
                if debug:
                    print '\nDopant', dopant.name, 'does not need an update'
                continue
            dopants_f_corr = dopants_f[dopant_key] + df_dopants[dopant_key] * dt
            dopants_f_corr[np.where(dopants_f_corr > 1.0)] = 1.0
            dopants_f_corr[np.where(dopants_f_corr < 0.0)] = 0.0
            #dopants_f_corr[np.where(z_nodes > z_limit)] = 1.0 # dopants_f_corr[z_limit_idx][-1]
            dopant.set_F_interp(z_nodes, dopants_f_corr)
            dopant.set_dF_interp(z_nodes, np.zeros_like(dopants_f_corr))

        for bi in schottky_diode.Semiconductor.bonding_interfaces:
            bi_tilt_f = []
            for trap in bi.dsl_tilt.traps:
                localized_trap_key = bi.label + '_tilt_' + trap[0].name + '_F'
                if localized_trap_key not in localized_traps_skip_list:
                    bonding_interface_f[localized_trap_key] += df_bonding_interfaces[localized_trap_key] * dt
                    if bonding_interface_f[localized_trap_key] > 1:
                        bonding_interface_f[localized_trap_key] = 1
                    elif bonding_interface_f[localized_trap_key] < 0:
                        bonding_interface_f[localized_trap_key] = 0
                if debug:
                    print 'F:', bonding_interface_f[localized_trap_key]
                bi_tilt_f.append(bonding_interface_f[localized_trap_key])
            bi_twist_f = []
            for trap in bi.dsl_twist.traps:
                localized_trap_key = bi.label + '_twist_' + trap[0].name + '_F'
                if localized_trap_key not in localized_traps_skip_list:
                    bonding_interface_f[localized_trap_key] += df_bonding_interfaces[localized_trap_key] * dt
                    if bonding_interface_f[localized_trap_key] > 1:
                        bonding_interface_f[localized_trap_key] = 1
                    elif bonding_interface_f[localized_trap_key] < 0:
                        bonding_interface_f[localized_trap_key] = 0
                if debug:
                    print 'F:', bonding_interface_f[localized_trap_key]
                bi_twist_f.append(bonding_interface_f[localized_trap_key])
            bi.set_traps_f(np.array(bi_tilt_f), np.array(bi_twist_f))
            bi.set_traps_df(np.zeros(len(bi_tilt_f)), np.zeros(len(bi_twist_f)))
        t += dt

    if debug_plot:
        plt.ioff()

    return t_points, potential_t, field_d, z_t, diode_voltage_drop_t, current_density_t, \
        bonding_interfaces_f_t, dopants_f_t, last_state_id
