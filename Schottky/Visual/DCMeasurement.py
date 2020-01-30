import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


def plot_band_diagram(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    ep_flat_mesh = measurement.electric_potential.flatten()
    z_ep = np.asarray(ep_flat_mesh.physical_nodes) * 1e6
    ep = np.asarray(ep_flat_mesh.solution) - measurement.v_bi + measurement.bias
    flat_mesh = measurement.quasi_fermi_e.flatten()
    z_qfe = np.asarray(flat_mesh.physical_nodes) * 1e6
    qfe = np.asarray(flat_mesh.solution)
    flat_mesh = measurement.quasi_fermi_h.flatten()
    z_qfh = np.asarray(flat_mesh.physical_nodes) * 1e6
    qfh = np.asarray(flat_mesh.solution)
    phi_b_n = measurement.diode.phi_b_n_ev_t(measurement.temperature)
    phi_b_p = measurement.diode.phi_b_p_ev_t(measurement.temperature)
    ec = ep + phi_b_n
    ev = ep - phi_b_p
    length = measurement.diode.length * 1e6
    if phi_b_p > 0:
        metal_energy_width = phi_b_p * 1.25
    else:
        metal_energy_width = (phi_b_n + phi_b_p) * 0.75
    metal_energy_width = (phi_b_n + phi_b_p) * 0.75
    metal_length = length / 5.0
    metal_patch = Rectangle((-metal_length, -metal_energy_width), metal_length, metal_energy_width,
                            facecolor='lightgray', edgecolor='k', linewidth=1)
    ax.plot(z_ep, ec, color='k', linewidth=1)
    ax.plot(z_ep, ev, color='k', linewidth=1)
    ax.plot(z_qfe, ec-qfe, color='b', linewidth=1, linestyle=':')
    ax.plot(z_qfh, ec - qfh, color='r', linewidth=1, linestyle=':')
    ax.plot(z_ep, np.ones_like(z_ep) * measurement.bias, color='k', linewidth=1, linestyle='-.')
    ax.plot([0., 0.], [-phi_b_p, phi_b_n], color='k', linewidth=1)
    ax.add_patch(metal_patch)
    ax.text(-metal_length / 2.0, -metal_energy_width / 2, measurement.diode.metal.label,
            horizontalalignment='center', verticalalignment='center')

    ax.set_xlim([-metal_length * 1.1, length])
    ax.set_ylim([min(min(ev), -metal_energy_width) - 0.1, max(max(ec), 0) + 0.1])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('E, eV')
    return ax


def plot_ep(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.electric_potential.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    ep = np.asarray(flat_mesh.solution)
    ax.plot(z, ep, color='k', linewidth=2, linestyle='-')
    ax.fill_between(z, ep, 0, facecolor='silver', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Electric Potential, V')
    return ax


def plot_qfe(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.quasi_fermi_e.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    qfe = np.asarray(flat_mesh.solution)
    ax.plot(z, qfe, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z, qfe, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_qfh(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.quasi_fermi_h.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    qfh = np.asarray(flat_mesh.solution)
    ax.plot(z, qfh, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z, qfh, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_qfeh(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.quasi_fermi_h.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    qfh = np.asarray(flat_mesh.solution)
    ax.plot(z, qfh, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z, qfh, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('Energy, J')
    return ax


def plot_n_e(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.n_e.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    n = np.asarray(flat_mesh.solution) * 1e-6
    ax.plot(z, n, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z, n, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    return ax


def plot_n_h(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.n_h.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    n = np.asarray(flat_mesh.solution) * 1e-6
    ax.plot(z, n, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z, n, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    return ax


def plot_n_eh(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh_e = measurement.n_e.flatten()
    z_e = np.asarray(flat_mesh_e.physical_nodes) * 1e6
    n_e = np.asarray(flat_mesh_e.solution) * 1e-6
    flat_mesh_h = measurement.n_h.flatten()
    z_h = np.asarray(flat_mesh_h.physical_nodes) * 1e6
    n_h = np.asarray(flat_mesh_h.solution) * 1e-6
    ax.plot(z_e, n_e, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z_e, n_e, 0, facecolor='blue', alpha=0.5)
    ax.plot(z_h, n_h, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z_h, n_h, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('N, cm$^{-3}$')
    ax.set_yscale('log')
    return ax

def plot_generation(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.generation.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    g = np.asarray(flat_mesh.solution) * 1e-6
    ax.plot(z, g, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z, g, 0, facecolor='red', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax


def plot_recombination(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh = measurement.recombination.flatten()
    z = np.asarray(flat_mesh.physical_nodes) * 1e6
    r = np.asarray(flat_mesh.solution) * 1e-6
    ax.plot(z, r, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z, r, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax

def plot_generation_recombination(measurement, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    flat_mesh_g = measurement.generation.flatten()
    z_g = np.asarray(flat_mesh_g.physical_nodes) * 1e6
    g = np.asarray(flat_mesh_g.solution) * 1e-6
    flat_mesh_r = measurement.recombination.flatten()
    z_r = np.asarray(flat_mesh_r.physical_nodes) * 1e6
    r = np.asarray(flat_mesh_r.solution) * 1e-6
    ax.plot(z_g, g, color='red', linewidth=2, linestyle='-')
    ax.fill_between(z_g, g, 0, facecolor='red', alpha=0.5)
    ax.plot(z_r, r, color='blue', linewidth=2, linestyle='-')
    ax.fill_between(z_r, r, 0, facecolor='blue', alpha=0.5)
    ax.set_xlim([0, measurement.diode.length * 1e6])
    ax.set_xlabel('z, $\mu$m')
    ax.set_ylabel('dN/dt, cm$^{-3}$s$^{-1}$')
    return ax
