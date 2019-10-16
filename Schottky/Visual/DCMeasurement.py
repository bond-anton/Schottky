import numpy as np
from matplotlib import pyplot as plt


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
