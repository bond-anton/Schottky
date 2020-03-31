import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


def plot_band_diagram(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    ep = np.asarray(diode.electric_potential.evaluate(z)) - diode.v_bi + diode.bias
    qfe = np.asarray(diode.quasi_fermi_e.evaluate(z))
    qfh = np.asarray(diode.quasi_fermi_h.evaluate(z))
    phi_b_n = diode.phi_b_n_ev_t(diode.temperature)
    phi_b_p = diode.phi_b_p_ev_t(diode.temperature)
    ec = ep + phi_b_n
    ev = ep - phi_b_p
    length = diode.length * 1e6
    if phi_b_p > 0:
        metal_energy_width = phi_b_p * 1.25
    else:
        metal_energy_width = (phi_b_n + phi_b_p) * 0.75
    metal_energy_width = (phi_b_n + phi_b_p) * 0.75
    metal_length = length / 5.0
    metal_patch = Rectangle((-metal_length, -metal_energy_width), metal_length, metal_energy_width,
                            facecolor='lightgray', edgecolor='k', linewidth=1)
    ax.plot(z * 1e6, ec, color='k', linewidth=1)
    ax.plot(z * 1e6, ev, color='k', linewidth=1)
    ax.plot(z * 1e6, ec - qfe, color='b', linewidth=1, linestyle=':')
    ax.plot(z * 1e6, ec - qfh, color='r', linewidth=1, linestyle=':')
    ax.plot(z * 1e6, np.ones_like(z) * diode.bias, color='k', linewidth=1, linestyle='-.')
    ax.plot([0., 0.], [-phi_b_p, phi_b_n], color='k', linewidth=1)
    ax.add_patch(metal_patch)
    ax.text(-metal_length / 2.0, -metal_energy_width / 2, diode.metal.label,
            horizontalalignment='center', verticalalignment='center')

    ax.set_xlim([-metal_length * 1.1, length])
    ax.set_ylim([min(min(ev), -metal_energy_width) - 0.1, max(max(ec), 0) + 0.1])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('E, eV')
    return ax


def plot_ep(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    ep = np.asarray(diode.electric_potential.evaluate(z))
    ax.plot(z * 1e6, ep, color='k', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, ep, 0, facecolor='silver', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Electric Potential, V')
    return ax


def plot_qfe(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    qfe = np.asarray(diode.quasi_fermi_e.evaluate(z))
    ax.plot(z * 1e6, qfe, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z, qfe, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_qfh(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    qfh = np.asarray(diode.quasi_fermi_h.evaluate(z))
    ax.plot(z * 1e6, qfh, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, qfh, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_qfeh(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    qfh = np.asarray(diode.quasi_fermi_h.evaluate(z))
    ax.plot(z * 1e6, qfh, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, qfh, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_n_e(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    n = np.asarray(diode.n_e.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, n, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, n, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    return ax


def plot_n_h(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    n = np.asarray(diode.n_h.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, n, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, n, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    return ax


def plot_n_eh(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    n_e = np.asarray(diode.n_e.evaluate(z)) * 1e-6
    n_h = np.asarray(diode.n_h.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, n_e, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, n_e, 0, facecolor='blue', alpha=0.5)
    ax.plot(z * 1e6, n_h, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, n_h, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    ax.set_yscale('log')
    return ax

def plot_generation(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    g = np.asarray(diode.generation.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, g, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, g, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax


def plot_recombination(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    r = np.asarray(diode.recombination.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, r, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, r, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax

def plot_generation_recombination(diode, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    z = np.linspace(0.0, diode.length, num=500, dtype=np.float)
    g = np.asarray(diode.generation.evaluate(z)) * 1e-6
    r = np.asarray(diode.recombination.evaluate(z)) * 1e-6
    ax.plot(z * 1e6, g, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, g, 0, facecolor='red', alpha=0.5)
    ax.plot(z * 1e6, r, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z * 1e6, r, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax
