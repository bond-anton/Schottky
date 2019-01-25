from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from Schottky import constant


def draw_bands_diagram(trap, f=0.0, e_f=None, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    max_x = 1.0
    x = np.linspace(-max_x, max_x, num=201, endpoint=True, dtype=np.double)
    band_gap = trap.energy_c_ev + trap.energy_v_ev
    phi_e = trap.capture_barrier_ev[1] * f
    phi_h = trap.capture_barrier_ev[0] * (1 - f)
    e_c = band_gap + np.exp(-abs(x) / (max_x / 3)) * phi_e - np.exp(-abs(x) / (max_x / 3)) * phi_h
    e_v = np.exp(-abs(x) / (max_x / 3)) * phi_e - np.exp(-abs(x) / (max_x / 3)) * phi_h
    if trap.cb_bound:
        e_t = trap.energy_v_ev + phi_e
    else:
        e_t = trap.energy_v_ev - phi_h
    ax.plot(x, e_c, 'k-', linewidth=2)
    ax.plot(x, e_v, 'k-', linewidth=2)
    ax.plot([0.0], [e_t], 'ro')
    ax.fill_between(x, e_v, min(e_v) - 0.1 * band_gap, facecolor='silver', alpha=0.5)
    ax.fill_between(x, e_c, max(e_c) + 0.1 * band_gap, facecolor='silver', alpha=0.5)
    if e_f is not None:
        ax.plot(x, np.ones_like(x) * (band_gap - e_f),
                color='k', linestyle='-.', linewidth=1)
        ax.text(max_x * 1.05, band_gap - e_f, 'Ef',
                horizontalalignment='center', verticalalignment='center')
    ax.text(max_x * 1.05, e_c[-1], 'Ec',
            horizontalalignment='center', verticalalignment='center')
    ax.text(max_x * 1.05, e_v[-1], 'Ev',
            horizontalalignment='center', verticalalignment='center')
    ax.set_xlim([-max_x, max_x * 1.1])
    ax.set_xticks([])
    ax.set_ylim([min(e_v) - 0.1 * band_gap, max(e_c) + 0.1 * band_gap])
    ax.set_xlabel('z, a.u.')
    ax.set_ylabel('E, eV')
    return ax
