from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from Schottky.Semiconductor import Semiconductor


def draw_bands_diagram(semiconductor, temperature, ax=None):
    if ax is None:
        _, ax = plt.subplots()
    return ax
