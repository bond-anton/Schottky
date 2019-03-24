from itertools import cycle
import numpy as np
from matplotlib import pyplot as plt


def draw_1d_profile_scalar(potential, direction, start, stop, num=100,
                           ax=None, scale=1.0, meV=True, electron=True,
                           color='k', linewidth=2, linestyle='-'):
    if scale == 1.0:
        units = 'm'
    elif scale == 1.0e-2:
        units = 'cm'
    elif scale == 1.0e-3:
        units = 'mm'
    elif scale == 1.0e-6:
        units = '$\mu$m'
    elif scale == 1.0e-9:
        units = 'nm'
    elif scale == 1.0e-10:
        units = '$\AA$'
    elif scale == 1.0e-12:
        units = 'pm'
    else:
        units = 'x%2.2g m'
    if ax is None:
        _, ax = plt.subplots()
        ax.set_xlabel('r along [%d, %d, %d], %s' % (direction[0], direction[1], direction[2], units))
        if electron:
            y_label = 'electron '
        else:
            y_label = 'hole '
        if meV:
            y_label += 'E, meV'
        else:
            y_label += 'E, eV'
        ax.set_ylabel(y_label)
    l = np.linspace(start, stop, num=num, dtype=np.double)
    r = direction / np.linalg.norm(direction)
    xyz = np.empty((l.shape[0], 3), dtype=np.double)
    xyz[:, 0] = l * scale * r[0]
    xyz[:, 1] = l * scale * r[1]
    xyz[:, 2] = l * scale * r[2]
    pot = np.asarray(potential.scalar_field(xyz))
    if meV:
        pot *= 1000
    if electron:
        pot *= -1
    ax.plot(l, pot, color=color, linewidth=linewidth, linestyle=linestyle, label=potential.name)
    return ax


def draw_1d_profile_scalar_polar(potential, theta, phi, r_start, r_stop, num=100,
                                 ax=None, scale=1.0, meV=True, electron=True,
                                 color='k', linewidth=2, linestyle='-'):
    if scale == 1.0:
        units = 'm'
    elif scale == 1.0e-2:
        units = 'cm'
    elif scale == 1.0e-3:
        units = 'mm'
    elif scale == 1.0e-6:
        units = '$\mu$m'
    elif scale == 1.0e-9:
        units = 'nm'
    elif scale == 1.0e-10:
        units = '$\AA$'
    elif scale == 1.0e-12:
        units = 'pm'
    else:
        units = 'x%2.2g m'
    if ax is None:
        _, ax = plt.subplots()
        ax.set_xlabel(r'r along $\theta$: %2.2f, $\phi$: %2.2f, %s' % (theta, phi, units))
        if electron:
            y_label = 'electron '
        else:
            y_label = 'hole '
        if meV:
            y_label += 'E, meV'
        else:
            y_label += 'E, eV'
        ax.set_ylabel(y_label)
    r = np.linspace(r_start, r_stop, num=num, dtype=np.double)
    rtp = np.empty((r.shape[0], 3), dtype=np.double)
    rtp[:, 0] = r * scale
    rtp[:, 1] = theta
    rtp[:, 2] = phi
    pot = np.asarray(potential.scalar_field_polar(rtp))
    if meV:
        pot *= 1000
    if electron:
        pot *= -1
    ax.plot(r, pot, color=color, linewidth=linewidth, linestyle=linestyle, label=potential.name)
    ax.set_xlim(r_start, r_stop)
    return ax


def draw_1d_profile_superposition_scalar(potential, direction, start, stop, num=100, ax=None,
                                         scale=1.0, meV=True, electron=True,
                                         bw=True, draw_sum=True):
    linestyles = cycle(['-.', '--', ':', '-'])
    if bw:
        colors = cycle(['k', 'gray', 'lightgray'])
        linewidths = cycle(3 * [1.5] + 3 * [2])
    else:
        colors = cycle(['C%d' % i for i in range(10)])
        linewidths = cycle(10 * [1.5] + 10 * [2])
    if scale == 1.0:
        units = 'm'
    elif scale == 1.0e-2:
        units = 'cm'
    elif scale == 1.0e-3:
        units = 'mm'
    elif scale == 1.0e-6:
        units = '$\mu$m'
    elif scale == 1.0e-9:
        units = 'nm'
    elif scale == 1.0e-10:
        units = '$\AA$'
    elif scale == 1.0e-12:
        units = 'pm'
    else:
        units = 'x%2.2g m'
    if ax is None:
        _, ax = plt.subplots()
        ax.set_xlabel('r along [%d, %d, %d], %s' % (direction[0], direction[1], direction[2], units))
        if electron:
            y_label = 'electron '
        else:
            y_label = 'hole '
        if meV:
            y_label += 'E, meV'
        else:
            y_label += 'E, eV'
        ax.set_ylabel(y_label)
    for component in potential.fields:
        draw_1d_profile_scalar(component, direction, start, stop, num=num,
                               ax=ax, scale=scale, meV=meV, electron=electron,
                               color=next(colors), linewidth=next(linewidths), linestyle=next(linestyles))
    if draw_sum:
        draw_1d_profile_scalar(potential, direction, start, stop, num=num,
                               ax=ax, scale=scale, meV=meV, electron=electron,
                               color='k', linewidth=2, linestyle='-')
    return ax


def draw_1d_profile_superposition_scalar_polar(potential, theta, phi, r_start, r_stop, num=100, ax=None,
                                               scale=1.0, meV=True, electron=True,
                                               bw=True, draw_sum=True):
    linestyles = cycle(['-.', '--', ':', '-'])
    if bw:
        colors = cycle(['k', 'gray', 'lightgray'])
        linewidths = cycle(3 * [1.5] + 3 * [2])
    else:
        colors = cycle(['C%d' % i for i in range(10)])
        linewidths = cycle(10 * [1.5] + 10 * [2])
    if scale == 1.0:
        units = 'm'
    elif scale == 1.0e-2:
        units = 'cm'
    elif scale == 1.0e-3:
        units = 'mm'
    elif scale == 1.0e-6:
        units = '$\mu$m'
    elif scale == 1.0e-9:
        units = 'nm'
    elif scale == 1.0e-10:
        units = '$\AA$'
    elif scale == 1.0e-12:
        units = 'pm'
    else:
        units = 'x%2.2g m'
    if ax is None:
        _, ax = plt.subplots()
        ax.set_xlabel(r'r along $\theta$: %2.2f, $\phi$: %2.2f, %s' % (theta, phi, units))
        if electron:
            y_label = 'electron '
        else:
            y_label = 'hole '
        if meV:
            y_label += 'E, meV'
        else:
            y_label += 'E, eV'
        ax.set_ylabel(y_label)
    for component in potential.fields:
        draw_1d_profile_scalar_polar(component, theta, phi, r_start, r_stop, num=num,
                                     ax=ax, scale=scale, meV=meV, electron=electron,
                                     color=next(colors), linewidth=next(linewidths), linestyle=next(linestyles))
    if draw_sum:
        draw_1d_profile_scalar_polar(potential, theta, phi, r_start, r_stop, num=num,
                                     ax=ax, scale=scale, meV=meV, electron=electron,
                                     color='k', linewidth=2, linestyle='-')
    return ax


def draw_r_min_polar(potential, theta, ax=None, scale=1.0):
    if scale == 1.0:
        units = 'm'
    elif scale == 1.0e-2:
        units = 'cm'
    elif scale == 1.0e-3:
        units = 'mm'
    elif scale == 1.0e-6:
        units = '$\mu$m'
    elif scale == 1.0e-9:
        units = 'nm'
    elif scale == 1.0e-10:
        units = '$\AA$'
    elif scale == 1.0e-12:
        units = 'pm'
    else:
        units = 'x%2.2g m'
    if ax is None:
        _, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    phi = np.linspace(0.0, np.pi * 2, num=potential.phi_resolution, endpoint=True)
    r = np.asarray(potential.max_energy_r(theta, phi))
    ax.plot(phi, r / scale, label=r'r$_{min}$ (%s), $\theta$=%2.2f' % (units, theta))
    return ax


def draw_energy_lowering_polar_phi(potential, theta, ax=None, meV=True):
    if ax is None:
        _, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    phi = np.linspace(0.0, np.pi * 2, num=potential.phi_resolution, endpoint=True)
    delta_e = np.asarray(potential.energy_lowering_phi_range(theta, phi))
    if meV:
        delta_e *= 1000
        units = 'meV'
    else:
        units = 'eV'
    ax.plot(phi, delta_e, label=r'$\Delta$E (%s), $\theta$=%2.2f' % (units, theta))
    return ax


def draw_energy_lowering_polar_theta(potential, phi, ax=None, meV=True):
    if ax is None:
        _, ax = plt.subplots()
    theta = np.linspace(0.0, np.pi, num=potential.theta_resolution, endpoint=True)
    delta_e = np.asarray(potential.energy_lowering_theta_range(phi, theta))
    if meV:
        delta_e *= 1000
        units = 'meV'
    else:
        units = 'eV'
    ax.plot(theta, delta_e, label=r'$\phi$=%2.2f' % phi)
    ax.set_xlabel(r'$\theta$, rad')
    ax.set_ylabel(r'$\Delta$E, %s' % units)
    ax.set_xlim(0.0, np.pi)
    return ax
