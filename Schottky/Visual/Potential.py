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
