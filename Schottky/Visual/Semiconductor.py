import numpy as np
from matplotlib import pyplot as plt

from Schottky import constant


def draw_bands_diagram(semiconductor, temperature, f_threshold=1.0e-23, max_iter=100, max_x=5, verbose=False, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    band_gap = semiconductor.band_gap_ev_t(temperature)
    x = np.linspace(0.0, max_x, num=100, endpoint=True, dtype=np.double)
    ax.plot(x, np.zeros_like(x), 'k-', linewidth=2)
    ax.plot(x, np.ones_like(x) * band_gap, 'k-', linewidth=2)
    e_f = semiconductor.el_chem_pot_t(temperature, f_threshold=f_threshold,
                                      max_iter=max_iter, verbose=verbose) / constant.q
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


def draw_dopants_profile(semiconductor, max_x=5, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    x = np.linspace(0.0, max_x, num=500) * 1e6
    for dopant in semiconductor.dopants:
        n = np.asarray(dopant.concentration.evaluate(x)) * 1e-6
        ax.semilogy(x, n, color=dopant.color, linestyle=dopant.linestyle, linewidth=1, label=dopant.label)
    ax.set_xlim([0, max_x])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    if semiconductor.dopants:
        ax.legend()
    return ax


def draw_bands_diagram_t(semiconductor, temperature, f_threshold=1.0e-23, max_iter=100, verbose=False, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    i_max = np.argmax(temperature)
    i_min = np.argmin(temperature)
    t_max = temperature[i_max]
    t_min = temperature[i_min]
    band_gap = np.asarray(semiconductor.band_gap_ev(temperature))
    e_i = np.asarray(semiconductor.e_i_ev(temperature))
    ax.plot(temperature, np.zeros_like(temperature), 'k-', linewidth=2)
    ax.plot(temperature, band_gap, 'k-', linewidth=2)
    e_f = np.asarray(semiconductor.el_chem_pot(temperature, f_threshold=f_threshold,
                                               max_iter=max_iter, verbose=verbose)) / constant.q
    ax.plot(temperature, band_gap - e_f,
            color='k', linestyle='-.', linewidth=1)
    ax.plot(temperature, band_gap - e_i,
            color='r', linestyle=':', linewidth=1)
    if abs(e_i[i_max] - e_f[i_max]) < 0.1 * band_gap[i_max]:
        ax.text(t_max * 1.05, band_gap[i_max] - e_f[i_max], 'Ef\nEi',
                horizontalalignment='center', verticalalignment='center')
    else:
        ax.text(t_max * 1.05, band_gap[i_max] - e_f[i_max], 'Ef',
                horizontalalignment='center', verticalalignment='center')
        ax.text(t_max * 1.05, band_gap[i_max] - e_i[i_max], 'Ei',
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


def draw_dopants_occupation_diagram_t(semiconductor, temperature,
                                      f_threshold=1.0e-23, max_iter=100, verbose=False, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    i_max = np.argmax(temperature)
    i_min = np.argmin(temperature)
    t_max = temperature[i_max]
    t_min = temperature[i_min]
    e_f = semiconductor.el_chem_pot(temperature)
    f = np.empty(temperature.shape[0], dtype=np.double)
    for dopant in semiconductor.dopants:
        for i in range(temperature.shape[0]):
            f[i] = semiconductor.trap_eq_occupation(dopant, e_f[i], temperature[i],
                                                    f_threshold=f_threshold, max_iter=max_iter, verbose=verbose)
        ax.plot(temperature, f, marker=dopant.marker, color=dopant.color,
                linestyle=dopant.linestyle, linewidth=1, label=dopant.label)
    ax.set_xlim([t_min, t_max])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('T, K')
    ax.set_ylabel('Occupation')
    ax.grid(True)
    if semiconductor.dopants:
        ax.legend()
    return ax
