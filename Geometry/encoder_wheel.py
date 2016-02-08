import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

def arc_patch(center, radius1, radius2, theta1, theta2, ax=None, resolution=50, **kwargs):
    # make sure ax is not empty
    if ax is None:
        ax = plt.gca()
    # generate the points
    theta = np.linspace(theta1, theta2, resolution)
    points1 = np.vstack((radius1 * np.cos(theta) + center[0],
                         radius1 * np.sin(theta) + center[1]))
    points2 = np.vstack((radius2 * np.cos(theta[::-1]) + center[0],
                         radius2 * np.sin(theta[::-1]) + center[1]))
    points = np.hstack((points1, points2))
    # build the polygon and add it to the axes
    poly = mpatches.Polygon(points.T, closed=True, **kwargs)
    ax.add_patch(poly)
    return poly

R = 150  # radius of the wheel in mm
R_inner = 75 # radius of the hole
N = 8   # maximal number of bits

dR = (R - R_inner) / N
print 'dR = %2.2f mm' % dR
_, ax = plt.subplots()
center = (0.0, 0.0)
for bit in range(1, N+1):
    print 'BIT:', bit
    dTheta = 2 * np.pi / 2**bit
    print 'dTheta = %2.4f deg' % np.degrees(dTheta)
    r1 = R_inner + dR * (bit - 1)
    r2 = R_inner + dR * bit
    dw1 = dTheta * r1
    dw2 = dTheta * r2
    print 'W = %2.4f mm, %2.4fmm' % (dw1, dw2)
    for segment in range(2**bit):
        th1 = (segment - 0.5) * dTheta
        th2 = th1 + dTheta
        fill = segment % 2 == 0
        arc_patch(center, r1, r2, th1, th2, ax=ax, resolution=100, fill=fill, color='black')


ax.set_title('Encoder Wheel')

# axis settings
ax.set_aspect('equal')
ax.set_xlim(-R, R)
ax.set_ylim(-R, R)

plt.show()
