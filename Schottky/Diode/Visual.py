from __future__ import division

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import transforms, ticker, cm
import matplotlib.animation as animation

from scipy.constants.codata import epsilon0

from Schottky.Notation import q  # , k
from Schottky.Helpers import Psi_zero, to_numeric


def V_bi_diagram(ax, SchDiode, T_start=0, T_stop=700, T_step=25):
    T_range = np.arange(T_start, T_stop + T_step, T_step)
    V_bi = np.zeros_like(T_range)
    T_now = SchDiode.T
    for i, T in enumerate(T_range):
        SchDiode.set_T(T)
        V_bi[i] = SchDiode.V_bi(eV=True)
    SchDiode.set_T(T_now)
    ax.plot(T, V_bi, linewidth=2, color='black', linestyle='-')
    ax.set_xlabel("T, K")
    ax.set_ylabel("V_bi, V")


def BandsBending_prepare_axes(ax, eV=True):
    E_units = 'eV' if eV else 'J'
    ax.set_title('Bands Bending Diagram')
    ax.set_ylabel('Energy, ' + E_units)
    ax.set_xlabel(r'Coordinate (z), um')


def Z_coordinate_ticks(ax, SchDiode, fancy_ticks=True):
    if fancy_ticks:
        ticks = ax.get_xticks()
        labels = np.array([lbl.get_text() for lbl in ax.get_xticklabels()])
        # ticks = np.append(ticks, 0)
        # labels = np.append(labels, '0')
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        for BI in SchDiode.Semiconductor.bonding_interfaces:
            ticks = np.append(ticks, BI.depth * 1e6)
            labels = np.append(labels, 'depth=%2.2f' % (BI.depth * 1e6))
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)


def annotate_energies(ax, SchDiode, arrow_z_pos, Psi, Vd, SchottkyEffect, metal_size, draw_metal, eV, draw_energies):
    annot_fig = {}
    if draw_energies:
        barrier_style = [1, 'black', '--']
        z_barrier, phi_bn, _, phi_b = SchDiode.get_phi_bn(Psi, Vd, SchottkyEffect=SchottkyEffect)
        if not eV:
            phi_bn *= to_numeric(q)
            phi_b *= to_numeric(q)
        z_barrier *= 1e6
        phi_bn_z = -metal_size[0] / 2
        phi_b_z = arrow_z_pos * 1e6
        allowed_band = SchDiode.Ec(Psi, Vd, np.float(phi_b_z * 1e-6), eV=eV)
        if SchDiode.Semiconductor.dop_type == 'p':
            phi_bn *= -1
            allowed_band = SchDiode.Ev(Psi, Vd, np.float(phi_b_z * 1e-6), eV=eV)
            if draw_metal:
                phi_bn_z = -metal_size[0] - 0.1
        annot_fig['zero_line'] = ax.hlines(0, phi_bn_z, 0,
                                           linewidth=barrier_style[0], color=barrier_style[1],
                                           linestyle=barrier_style[2])
        annot_fig['phi_bn_line'] = ax.hlines(phi_bn, phi_bn_z, z_barrier,
                                             linewidth=barrier_style[0], color=barrier_style[1],
                                             linestyle=barrier_style[2])
        annot_fig['phi_bn_arrows'] = ax.annotate('', xy=(phi_bn_z, phi_bn), xycoords='data', xytext=(phi_bn_z, 0),
                                                 textcoords='data', arrowprops={'arrowstyle': '<->'})
        annot_fig['phi_bn_annot'] = ax.annotate('$q\phi_{bn}$', xy=(phi_bn_z, phi_bn / 2), xycoords='data',
                                                xytext=(-5, 0), textcoords='offset points',
                                                horizontalalignment='right', verticalalignment='top')
        annot_fig['phi_b_line'] = ax.hlines(phi_bn, z_barrier, phi_b_z,
                                            linewidth=barrier_style[0], color=barrier_style[1],
                                            linestyle=barrier_style[2])
        annot_fig['phi_b_arrows'] = ax.annotate('', xy=(phi_b_z, phi_bn), xycoords='data',
                                                xytext=(phi_b_z, allowed_band), textcoords='data',
                                                arrowprops={'arrowstyle': '<->'})
        annot_fig['phi_b_annot'] = ax.annotate('$q\phi_{b}$', xy=(phi_b_z, (phi_bn + allowed_band) / 2),
                                               xycoords='data', xytext=(5, 0), textcoords='offset points',
                                               horizontalalignment='left', verticalalignment='top')
        return annot_fig


def BandsBending_draw_metal(ax, SchDiode, metal_size, bands_style, lbl_style, draw_energies, draw_metal):
    if draw_metal:
        # metal_size = [0.2, 1]  # um x eV
        metal_color = 'lightgrey'
        ax.add_patch(Rectangle((0, 0), -metal_size[0], -metal_size[1],
                               facecolor=metal_color, edgecolor=bands_style[1], linewidth=bands_style[0]))
        if SchDiode.Semiconductor.dop_type == 'p' and draw_energies:
            metal_label_y = -metal_size[1] / 4
        else:
            metal_label_y = -metal_size[1] / 2
        ax.text(-metal_size[0] / 2, metal_label_y, SchDiode.Metal.label,
                verticalalignment='bottom', horizontalalignment='center', color=lbl_style[1], fontsize=lbl_style[0])


def BandsBending_prepare_data(SchDiode, Psi, Vd, z, eV, SchottkyEffect, draw_metal):
    dz = np.gradient(z)
    Psi_points = Psi(z)
    E_points = -np.gradient(Psi_points, dz, edge_order=2)
    '''
    if SchDiode.Semiconductor.dop_type == 'n':
        idx = np.where(E_points > 0)[0]
    else:
        idx = np.where(E_points < 0)[0]
    if idx.size > 0:
        split = np.where(idx[1:]-idx[:-1] > 1)[0]+1
        if split.size > 0:
            idx = np.split(idx, split)[0]
        #print 'Inversion pt', z[idx][-1]*1e6
        FieldInversionPoint = z[idx[-1]]
        FieldInversionPsi = Psi_points[idx[-1]]
    else:
        FieldInversionPoint = 0
        FieldInversionPsi = 0
    '''
    FieldInversionPoint, PHI_bn, _, PHI_b = SchDiode.get_phi_bn(Psi, Vd, SchottkyEffect=False)
    FieldInversionPsi = Psi(FieldInversionPoint)
    idx = np.where(z <= FieldInversionPoint)[0]
    coeff = 1 if eV else to_numeric(q)
    Eg = SchDiode.Semiconductor.band_gap(SchDiode.T, symbolic=False, electron_volts=eV)
    Ec = SchDiode.Ec(Psi, Vd, z, eV=eV)
    type_sign = 1 if SchDiode.Semiconductor.dop_type == 'n' else -1
    Ef = np.zeros_like(z) + coeff * type_sign * Vd
    if FieldInversionPoint > 0:
        # xi = np.float(SchDiode.Semiconductor.Ech_pot(T=SchDiode.T, z=1e3, eV=eV, debug=False))
        Ef[0:idx[-1]] = coeff * FieldInversionPsi * type_sign
        # print idx[-1]
        # Ef[0:idx[-1]] = Ec[0:idx[-1]] - Ef[idx[-1]]
    Ev = Ec - Eg
    if draw_metal:
        Z = np.append(0, z)
        Ec = np.append(0, Ec)
        Ev = np.append(0, Ev)
        Ef = np.append(Ef[0], Ef)
    MirrEnergy = np.zeros_like(Ec)
    if SchottkyEffect:
        for ii, zz in enumerate(Z):
            if zz != 0:
                MirrEnergy[ii] = to_numeric(
                    q ** 2 / (16 * np.pi * epsilon0 * SchDiode.Semiconductor.reference['epsilon'] * zz))
                if eV:
                    MirrEnergy[ii] /= to_numeric(q)
            else:
                MirrEnergy[ii] = Ec[ii]
    return Z, Ec, Ev, Ef, Eg, MirrEnergy, FieldInversionPoint


def BandsBendingDiagram(ax, SchDiode, Psi=Psi_zero, Vd=0, z=0, BI_F={}, dopants_F={}, eV=True,
                        draw_metal=True, label_bands=True, fancy_ticks=True,
                        SchottkyEffect=False, draw_energies=False):
    print 'BB Diagram Start'
    Z, Ec, Ev, Ef, Eg, MirrEnergy, FieldInversionPoint = BandsBending_prepare_data(SchDiode, Psi, Vd, z, eV,
                                                                                   SchottkyEffect, draw_metal)
    bands_style = [2, 'black', '-']
    Ef_style = [1, 'black', ':']
    lbl_style = [12, 'black']
    inv = ax.transData.inverted()
    dxdy = inv.transform(ax.transData.transform((0, 0)) + np.array([bands_style[0], bands_style[0]]))
    metal_size = [0.2, 1]  # um x eV
    BandsBending_draw_metal(ax, SchDiode, metal_size, bands_style, lbl_style, draw_energies, draw_metal)
    arrow_z_pos = np.max(Z) / 2
    annot_fig = annotate_energies(ax, SchDiode, arrow_z_pos, Psi, Vd, SchottkyEffect, metal_size, draw_metal, eV,
                                  draw_energies)
    bands_lines = {}
    if SchDiode.Semiconductor.dop_type == 'n':
        bands_lines['Ec'], = ax.plot(Z * 1e6, Ec - MirrEnergy, linewidth=bands_style[0], color=bands_style[1],
                                     linestyle=bands_style[2])
        bands_lines['Ev'], = ax.plot(Z * 1e6, Ev, linewidth=bands_style[0], color=bands_style[1],
                                     linestyle=bands_style[2])
    else:
        bands_lines['Ec'], = ax.plot(Z * 1e6, Ec, linewidth=bands_style[0], color=bands_style[1],
                                     linestyle=bands_style[2])
        bands_lines['Ev'], = ax.plot(Z * 1e6, Ev + MirrEnergy, linewidth=bands_style[0], color=bands_style[1],
                                     linestyle=bands_style[2])
    Ef_idx = np.where(Z > FieldInversionPoint)
    # Ef_idx = np.where(Z >= 0)
    bands_lines['Ef'], = ax.plot(Z[Ef_idx] * 1e6, Ef[Ef_idx], linewidth=Ef_style[0], color=Ef_style[1],
                                 linestyle=Ef_style[2])
    bands_labels = {}
    if label_bands:
        trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
        if SchDiode.Semiconductor.dop_type == 'n':
            lbl_pos_y = [Ec[-1] + 3 * dxdy[1], Ev[-1] - 3 * dxdy[1], Ef[-1] - 3 * dxdy[1]]
            lbl_algn = ['bottom', 'top', 'top']
        else:
            lbl_pos_y = [Ec[-1] + 3 * dxdy[1], Ev[-1] - 3 * dxdy[1], Ef[-1] + 3 * dxdy[1]]
            lbl_algn = ['bottom', 'top', 'bottom']
        bands_labels['Ec'] = ax.text(0.9, lbl_pos_y[0], 'Ec', verticalalignment=lbl_algn[0], horizontalalignment='left',
                                     transform=trans, color=lbl_style[1], fontsize=lbl_style[0])
        bands_labels['Ev'] = ax.text(0.9, lbl_pos_y[1], 'Ev', verticalalignment=lbl_algn[1], horizontalalignment='left',
                                     transform=trans, color=lbl_style[1], fontsize=lbl_style[0])
        bands_labels['Ef'] = ax.text(0.9, lbl_pos_y[2], 'Ef', verticalalignment=lbl_algn[2], horizontalalignment='left',
                                     transform=trans, color=lbl_style[1], fontsize=lbl_style[0])
    dopant_lines = {}
    for dopant in SchDiode.Semiconductor.dopants:
        dopant_style = [1, 'black', '--']
        threshold_dopant = 0.05 if eV else to_numeric(q) * 0.05
        level = dopant.energy_level(SchDiode.T, SchDiode.Semiconductor, electron_volts=eV)
        if level > threshold_dopant and level < (Eg - threshold_dopant):
            if draw_metal:
                dopant_lines[dopant.name], = ax.plot(Z * 1e6, np.append(Ec[1], Ec[1:]) - level,
                                                     linewidth=dopant_style[0], color=dopant_style[1],
                                                     linestyle=dopant_style[2])
            else:
                dopant_lines[dopant.name], = ax.plot(Z * 1e6, Ec - level, linewidth=dopant_style[0],
                                                     color=dopant_style[1], linestyle=dopant_style[2])
    BI_levels = {}
    for j, BI in enumerate(SchDiode.Semiconductor.bonding_interfaces):
        for i, trap in enumerate(BI.dsl_tilt.traps):
            charge_state_idx = 0
            level = trap[0].energy_level(SchDiode.T, SchDiode.Semiconductor, charge_state_idx, electron_volts=eV)
            print 'Et =', level, 'eV'
            F_i = BI_F[BI.label + '_tilt_' + trap[0].name + '_F']
            gray = (1 - np.float(F_i), 1 - np.float(F_i), 1 - np.float(F_i))
            BI_levels[str(j) + '-tilt-' + trap[0].name] = ax.scatter(BI.depth * 1e6,
                                                                     SchDiode.Ec(Psi, Vd, BI.depth, eV=True) - level, s=40,
                                                                     edgecolors='black', c=gray[0], vmin=0, vmax=1,
                                                                     cmap=cm.gray)
        for i, trap in enumerate(BI.dsl_twist.traps):
            charge_state_idx = 0
            level = trap[0].energy_level(SchDiode.T, SchDiode.Semiconductor, charge_state_idx, electron_volts=eV)
            print 'Et =', level, 'eV'
            F_i = BI_F[BI.label + '_twist_' + trap[0].name + '_F']
            gray = (1 - np.float(F_i), 1 - np.float(F_i), 1 - np.float(F_i))
            BI_levels[str(j) + '-twist-' + trap[0].name] = ax.scatter(BI.depth * 1e6,
                                                                      SchDiode.Ec(Psi, Vd, BI.depth, eV=True) - level,
                                                                      s=40, edgecolors='black', c=gray[0], vmin=0,
                                                                      vmax=1, cmap=cm.gray)

    BandsBending_prepare_axes(ax, eV)
    Z_coordinate_ticks(ax, SchDiode, fancy_ticks)
    print 'BB Diagram Stop'
    return bands_lines, bands_labels, annot_fig, dopant_lines, BI_levels


def BandsBendingAnimation(fig, ax, SchDiode, Psi, Va, Vd, z, BI_F, dopants_F, T, eV=True,
                          draw_metal=True, label_bands=True, fancy_ticks=False,
                          SchottkyEffect=False, draw_energies=False, interval=100):
    # init Diagram
    bands_lines, bands_labels, annot_fig, dopant_lines, BI_levels = BandsBendingDiagram(ax, SchDiode, Psi[0], Vd[0], z,
                                                                                        BI_F[0], dopants_F[0], eV,
                                                                                        draw_metal, label_bands,
                                                                                        fancy_ticks,
                                                                                        SchottkyEffect, draw_energies)

    V_label = ax.text(0.05, 0.95, '', transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='left', color='black', fontsize=12)
    T_label = ax.text(0.55, 0.95, '', transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='left', color='black', fontsize=12)

    def init():
        Eg = SchDiode.Semiconductor.band_gap(SchDiode.T, symbolic=False, electron_volts=eV)
        bands_lines['Ec'].set_data([], [])
        bands_lines['Ev'].set_data([], [])
        bands_lines['Ef'].set_data([], [])
        for dopant in SchDiode.Semiconductor.dopants:
            threshold_dopant = 0.05 if eV else to_numeric(q) * 0.05
            level = dopant.energy_level(SchDiode.T, SchDiode.Semiconductor, electron_volts=eV)
            if level > threshold_dopant and level < (Eg - threshold_dopant):
                dopant_lines[dopant.name].set_data([], [])
        if label_bands:
            bands_labels['Ec'].set_text('')
            bands_labels['Ev'].set_text('')
            bands_labels['Ef'].set_text('')
        V_label.set_text('')
        T_label.set_text('')
        for key in BI_levels.keys():
            BI_levels[key].set_offsets(np.array([]))
            BI_levels[key].set_array(np.array([]))
        fig_list = bands_lines.values()
        fig_list.extend(bands_labels.values())
        fig_list.extend(dopant_lines.values())
        fig_list.extend([V_label, T_label])
        fig_list.extend(BI_levels.values())
        return fig_list

    def update_frame(i):
        try:
            V_label.set_text('Va = %2.2f V, Vd = %2.2f V' % (Va[i], Vd[i]))
            T_label.set_text('T = %3.2f K' % T[i])
            Z_coordinate_ticks(ax, SchDiode, fancy_ticks)
            Z, Ec, Ev, Ef, Eg, MirrEnergy, FieldInversionPoint = BandsBending_prepare_data(SchDiode, Psi[i], Vd[i], z,
                                                                                           eV, SchottkyEffect,
                                                                                           draw_metal)
            if SchDiode.Semiconductor.dop_type == 'n':
                bands_lines['Ec'].set_data(Z * 1e6, (Ec - MirrEnergy))
                bands_lines['Ev'].set_data(Z * 1e6, Ev)
            else:
                bands_lines['Ec'].set_data(Z * 1e6, Ec)
                bands_lines['Ev'].set_data(Z * 1e6, Ev + MirrEnergy)
            Ef_idx = np.where(Z > FieldInversionPoint)
            # Ef_idx = np.where(Z >= 0)
            bands_lines['Ef'].set_data(Z[Ef_idx] * 1e6, Ef[Ef_idx])
            if label_bands:
                inv = ax.transData.inverted()
                dxdy = inv.transform(ax.transData.transform((0, 0)) + np.array([2, 2]))
                # trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
                if SchDiode.Semiconductor.dop_type == 'n':
                    lbl_pos_y = [Ec[-1] + 3 * dxdy[1], Ev[-1] - 3 * dxdy[1], Ef[-1] - 3 * dxdy[1]]
                else:
                    lbl_pos_y = [Ec[-1] + 3 * dxdy[1], Ev[-1] - 3 * dxdy[1], Ef[-1] + 3 * dxdy[1]]
                bands_labels['Ec'].set_y(lbl_pos_y[0])
                bands_labels['Ev'].set_y(lbl_pos_y[1])
                bands_labels['Ef'].set_y(lbl_pos_y[2])
                bands_labels['Ec'].set_text('Ec')
                bands_labels['Ev'].set_text('Ev')
                bands_labels['Ef'].set_text('Ef')

            for dopant in SchDiode.Semiconductor.dopants:
                threshold_dopant = 0.05 if eV else to_numeric(q) * 0.05
                level = dopant.energy_level(SchDiode.T, SchDiode.Semiconductor, electron_volts=eV)
                if level > threshold_dopant and level < (Eg - threshold_dopant):
                    if draw_metal:
                        dopant_lines[dopant.name].set_data(Z * 1e6, np.append(Ec[1], Ec[1:]) - level)
                    else:
                        dopant_lines[dopant.name].set_data(Z * 1e6, Ec - level)

            for j, BI in enumerate(SchDiode.Semiconductor.bonding_interfaces):
                for l, trap in enumerate(BI.dsl_tilt.traps):
                    charge_state_idx = 0
                    level = trap[0].energy_level(SchDiode.T, SchDiode.Semiconductor, charge_state_idx, electron_volts=eV)
                    # print 'Et =', level, 'eV'
                    F_l = BI_F[i][BI.label + '_tilt_' + trap[0].name + '_F']
                    gray = np.array([1 - np.float(F_l), 1 - np.float(F_l), 1 - np.float(F_l), 1.0])
                    # print gray
                    # gray = np.array([(1 - np.float(BI.dsl_tilt_f[l]), 1 - np.float(BI.dsl_tilt_f[l]), 1 - np.float(BI.dsl_tilt_f[l]))])
                    BI_levels[str(j) + '-tilt-' + trap[0].name].set_offsets(
                        np.array([BI.depth * 1e6, SchDiode.Ec(Psi[i], Vd[i], BI.depth, eV=True) - level]))
                    BI_levels[str(j) + '-tilt-' + trap[0].name].set_array(gray)
                for l, trap in enumerate(BI.dsl_twist.traps):
                    charge_state_idx = 0
                    level = trap[0].energy_level(SchDiode.T, SchDiode.Semiconductor, charge_state_idx, electron_volts=eV)
                    # print 'Et =', level, 'eV'
                    F_l = BI_F[i][BI.label + '_twist_' + trap[0].name + '_F']
                    gray = (1 - np.float(F_l), 1 - np.float(F_l), 1 - np.float(F_l))
                    # gray = np.array([(1 - np.float(BI.dsl_twist_f[l]), 1 - np.float(BI.dsl_twist_f[l]), 1 - np.float(BI.dsl_twist_f[l]))])
                    BI_levels[str(j) + '-twist-' + trap[0].name].set_offsets(
                        np.array([BI.depth * 1e6, SchDiode.Ec(Psi[i], Vd[i], BI.depth, eV=True) - level]))
                    # BI_levels[str(j)+'-twist-'+trap[0].name].set_array(gray)


        except Exception as e:
            print '!!! ====>', e
            pass
        fig_list = bands_lines.values()
        fig_list.extend(dopant_lines.values())
        fig_list.extend(bands_labels.values())
        fig_list.extend([V_label, T_label])
        fig_list.extend(BI_levels.values())
        return fig_list

    ani = animation.FuncAnimation(fig, update_frame, frames=len(Psi), init_func=init, blit=True, interval=interval,
                                  repeat=True)

    return ani


def ElectricFieldDiagram(ax, SchDiode, E, z=0, dot=False):
    print 'Electric Field Diagram Start'
    linestyle = '-'
    if dot:
        linestyle = ':'
    idx = np.where(z < SchDiode.FieldInversionPoint)[0]
    ax.plot(z * 1e6, E(z) * 1e-2, linewidth=2, color='black', linestyle=linestyle)
    if idx.size > 0:
        ax.plot(z[idx] * 1e6, E(z[idx]) * 1e-2, linewidth=2, color='red', linestyle=linestyle)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_title('Electric Field Diagram')
    ax.set_ylabel('Electric Field, V/cm')
    ax.set_xlabel(r'Coordinate (z), um')
    ax.grid(True)
    print 'Electric Field Diagram Stop'


def RhoDiagram(ax, SchDiode, Psi=Psi_zero, z=0):
    print 'Rho Diagram Start'
    rho_semi_num = SchDiode.build_rho_z_Psi(Psi, carriers_charge=True)
    rho = rho_semi_num(z, Psi) / to_numeric(q)
    n, p = SchDiode.n_carriers_theory(Psi, z)
    dopants_rho = []
    for dopant in SchDiode.Semiconductor.dopants:
        Nd = dopant.concentration(z)
        dopants_rho.append(
            Nd * ((dopant.charge_states[1][0] - dopant.charge_states[0][0]) * dopant.F(z) + dopant.charge_states[0][0]))

    ax.plot(z * 1e6, abs(rho * 1e-6), linewidth=2, color='black', linestyle='-')

    if SchDiode.Semiconductor.dop_type == 'n':
        ax.plot(z * 1e6, n * 1e-6, linewidth=2, color='blue', linestyle='-')
    elif SchDiode.Semiconductor.dop_type == 'p':
        ax.plot(z * 1e6, p * 1e-6, linewidth=2, color='red', linestyle='-')
    for dopant_rho in dopants_rho:
        ax.plot(z * 1e6, dopant_rho * 1e-6, linewidth=2, color='green', linestyle='-')

    ax.set_yscale('log')
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(-7,7))
    ax.set_title('Net density of charge')
    ax.set_ylabel('Concentration, cm$^{-3}$')
    ax.set_xlabel(r'Coordinate (z), um')
    ax.grid(True)
    print 'Rho Diagram Stop'


def CarriersConcentrationDiagram(ax, SchDiode, Psi=Psi_zero, z=0):
    # q_n = to_numeric(q)
    # print 'Charge Carriers Diagram Start'
    n, p = SchDiode.n_carriers_theory(Psi, z)
    if SchDiode.Semiconductor.dop_type == 'nn':
        ax.plot(z * 1e6, n * 1e-6, 'b-')
    elif SchDiode.Semiconductor.dop_type == 'pp':
        ax.plot(z * 1e6, p * 1e-6, 'r-')
    else:
        ax.plot(z * 1e6, n * 1e-6, 'b-')
        ax.plot(z * 1e6, p * 1e-6, 'r-')
        # print 'Charge Carriers Diagram Stop'


def prepare_debug_axes(axes_names, share_x=False, autoscale=True):
    axes = {}
    axes_number = len(axes_names)
    _, axes_tuple = plt.subplots(axes_number, sharex=share_x)
    for i, name in enumerate(axes_names):
        axes[name] = axes_tuple[i]
        axes[name].set_autoscalex_on(autoscale)
        axes[name].set_autoscaley_on(autoscale)
    return axes


def autoscale_axes(axes):
    if isinstance(axes, dict):
        axes_list = axes.values()
    else:
        axes_list = axes
    for ax in axes_list:
        ax.relim()
        ax.autoscale_view()


def create_lines(ax, names):
    lines = {}
    for line_name in names:
        lines[line_name], = ax.plot([], [])
    return lines
