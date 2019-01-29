from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from Schottky import constant


def draw_bands_diagram(semiconductor, temperature, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    band_gap = semiconductor.band_gap_t(temperature, electron_volts=True)
    max_x = 0
    for dopant in semiconductor.dopants:
        mesh = dopant.concentration.flatten()
        x = mesh.physical_nodes * 1e6
        max_x = max(max_x, max(x))
    x = np.linspace(0.0, max_x, num=100, endpoint=True, dtype=np.double)
    ax.plot(x, np.zeros_like(x), 'k-', linewidth=2)
    ax.plot(x, np.ones_like(x) * band_gap, 'k-', linewidth=2)
    e_f = semiconductor.el_chem_pot_t(temperature) / constant.q
    ax.plot(x, np.ones_like(x) * (band_gap - e_f),
            color='k', linestyle='-.', linewidth=1)
    ax.text(max_x * 1.05, band_gap - e_f, 'Ef',
            horizontalalignment='center', verticalalignment='center')
    ax.text(max_x * 1.05, band_gap, 'Ec',
            horizontalalignment='center', verticalalignment='center')
    ax.text(max_x * 1.05, 0.0, 'Ev',
            horizontalalignment='center', verticalalignment='center')
    for dopant in semiconductor.dopants:
        if dopant.cb_bound:
            ax.plot(x, np.ones_like(x) * (-dopant.energy_c_ev + band_gap),
                    color=dopant.color, linestyle=dopant.linestyle, linewidth=1)
            ax.text(max_x * 1.05, -dopant.energy_c_ev + band_gap,
                    dopant.label, color=dopant.color,
                    horizontalalignment='center', verticalalignment='center')
        else:
            ax.plot(x, np.ones_like(x) * dopant.energy_v_ev,
                    color=dopant.color, linestyle=dopant.linestyle, linewidth=1)
            ax.text(max_x * 1.05, dopant.energy_v_ev,
                    dopant.label, color=dopant.color,
                    horizontalalignment='center', verticalalignment='center')
    ax.fill_between(x, -band_gap * 0.1, 0.0, facecolor='silver', alpha=0.5)
    ax.fill_between(x, band_gap, band_gap * 1.1, facecolor='silver', alpha=0.5)
    ax.set_xlim([0, max_x * 1.1])
    ax.set_xticks([])
    ax.set_ylim([-band_gap * 0.1, band_gap * 1.1])
    ax.set_xlabel('z, bulk')
    ax.set_ylabel('E, eV')
    return ax


def draw_dopants_profile(semiconductor, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    max_x = 0
    for dopant in semiconductor.dopants:
        mesh = dopant.concentration.flatten()
        x = mesh.physical_nodes * 1e6
        n = np.asarray(dopant.n_t(x)) * 1e-6
        max_x = max(max_x, max(x))
        ax.plot(x, n, color=dopant.color, linestyle=dopant.linestyle, linewidth=1, label=dopant.label)
    ax.set_xlim([0, max_x * 1.1])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    ax.legend()
    return ax


def draw_bands_diagram_t(semiconductor, temperature, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    i_max = np.argmax(temperature)
    i_min = np.argmin(temperature)
    t_max = temperature[i_max]
    t_min = temperature[i_min]
    band_gap = np.asarray(semiconductor.band_gap(temperature, electron_volts=True))
    ax.plot(temperature, np.zeros_like(temperature), 'k-', linewidth=2)
    ax.plot(temperature, band_gap, 'k-', linewidth=2)
    e_f = np.asarray(semiconductor.el_chem_pot(temperature)) / constant.q
    ax.plot(temperature, band_gap - e_f,
            color='k', linestyle='-.', linewidth=1)
    ax.text(t_max * 1.05, band_gap[i_max] - e_f[i_max], 'Ef',
            horizontalalignment='center', verticalalignment='center')
    ax.text(t_max * 1.05, band_gap[i_max], 'Ec',
            horizontalalignment='center', verticalalignment='center')
    ax.text(t_max * 1.05, 0.0, 'Ev',
            horizontalalignment='center', verticalalignment='center')
    for dopant in semiconductor.dopants:
        if dopant.cb_bound:
            ax.plot(temperature, -dopant.energy_c_ev + band_gap,
                    color=dopant.color, linestyle=dopant.linestyle, linewidth=1)
            ax.text(t_max * 1.05, -dopant.energy_c_ev + band_gap[i_max],
                    dopant.label, color=dopant.color,
                    horizontalalignment='center', verticalalignment='center')
        else:
            ax.plot(temperature, np.ones_like(temperature) * dopant.energy_v_ev,
                    color=dopant.color, linestyle=dopant.linestyle, linewidth=1)
            ax.text(t_max * 1.05, dopant.energy_v_ev,
                    dopant.label, color=dopant.color,
                    horizontalalignment='center', verticalalignment='center')
    ax.fill_between(temperature, -max(band_gap) * 0.1, 0.0, facecolor='silver', alpha=0.5)
    ax.fill_between(temperature, band_gap, max(band_gap) * 1.1, facecolor='silver', alpha=0.5)
    ax.set_xlim([t_min, t_max * 1.1])
    ax.set_ylim([-max(band_gap) * 0.1, max(band_gap) * 1.1])
    ax.set_xlabel('T, K')
    ax.set_ylabel('E, eV')
    return ax