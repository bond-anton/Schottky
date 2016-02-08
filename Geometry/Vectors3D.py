import numpy as np


def cartesian_to_spherical(x, y, z):
    xy_square = x**2 + y**2
    r = np.sqrt(xy_square + z**2)
    theta = np.arctan2(np.sqrt(xy_square), z) # for elevation angle defined from Z-axis down
    phi = np.arctan2(y, x)
    return r, theta, phi


def spherical_to_cartesian(r, theta, phi):
    z = r * np.cos(theta)
    xy = r * np.sin(theta)
    x = xy * np.cos(phi)
    y = xy * np.sin(phi)
    return x, y, z


def angle_between_vectors((x1, y1, z1), (x2, y2, z2), debug=False):
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    r2 = np.sqrt(x2**2 + y2**2 + z2**2)
    dot_product = x1 * x2 + y1 * y2 + z1 * z2
    if debug:
        print 'a dot b:', dot_product
        print 'a, b:', r1, r2
        print 'ab:', r1 * r2
        print 'cos theta:', dot_product / (r1 * r2)
        print 'theta:', np.arccos(dot_product / (r1 * r2))
    return np.arccos(dot_product / (r1 * r2))
