import numpy as np
from matplotlib import pyplot as plt


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
