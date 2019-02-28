import numpy as np
from matplotlib import pyplot as plt


def draw_1d_profile_scalar(potential, direction, start, stop, num=100, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    l = np.linspace(start, stop, num=num, dtype=np.double)
    r = direction / np.linalg.norm(direction)
    xyz = np.empty((l.shape[0], 3), dtype=np.double)
    xyz[:, 0] = l * r[0]
    xyz[:, 1] = l * r[1]
    xyz[:, 2] = l * r[2]
    pot = potential.scalar_field(xyz)
    ax.plot(l, pot, 'r-o')
    return ax
