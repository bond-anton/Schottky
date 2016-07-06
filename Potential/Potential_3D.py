__author__ = 'anton'

import warnings
import numpy as np
from scipy.constants import elementary_charge, epsilon_0, Boltzmann
from scipy.optimize import root
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm
#from Geometry import Vectors3D
try:
    from mayavi import mlab
    USE_MAYAVI = True
except ImportError:
    USE_MAYAVI = False

class GenericPotential(object):

    def __init__(self, label='Unknown Potential'):
        self.USE_MAYAVI = USE_MAYAVI
        super(GenericPotential, self).__init__()
        self.label = label

    def potential(self, r, theta, phi):
        return np.empty((np.array(r).size, np.array(theta).size, np.array(phi).size))

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        return np.empty_like(r_grid)

    def field(self, r, theta, phi):
        F_r = np.empty((np.array(r).size, np.array(theta).size, np.array(phi).size))
        F_theta = np.empty((np.array(r).size, np.array(theta).size, np.array(phi).size))
        F_phi = np.empty((np.array(r).size, np.array(theta).size, np.array(phi).size))
        return F_r, F_theta, F_phi

    def field_cartesian(self, x, y, z):
        F_x = np.empty((np.array(x).size, np.array(y).size, np.array(z).size))
        F_y = np.empty((np.array(x).size, np.array(y).size, np.array(z).size))
        F_z = np.empty((np.array(x).size, np.array(y).size, np.array(z).size))
        return F_x, F_y, F_z

    def barrier_lowering(self, theta):
        return 0, 0

    def plot_potential(self, fig=None, r=None, theta=None, phi=None, lowering=True):
        if not self.USE_MAYAVI:
            return None
        if fig is None:
            fig = mlab.figure('Potential')#, bgcolor=(0, 0, 0))
        else:
            mlab.figure(fig, bgcolor=fig.scene.background)

        if r is None:
            r = np.linspace(1e-9, 20e-9, num=50, endpoint=True)
        if theta is None:
            theta = np.linspace(0, np.pi, num=50, endpoint=True)
        if phi is None:
            phi = np.linspace(0, 2 * np.pi, num=50, endpoint=True)
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        volumetric_data = self.potential(r, theta, phi)
        energy_scale = np.max(np.abs(r)) / np.max(np.abs(volumetric_data[:, :, 0]))
        dim_idx = 2
        for i, dim in enumerate(r_grid.shape):
            if dim == 1:
                dim_idx = i
                break
        if dim_idx == 2:
            print('here')
            surf = mlab.mesh(r_grid[:, :, 0] * np.cos(theta_grid[:, :, 0]),
                             r_grid[:, :, 0] * np.sin(theta_grid[:, :, 0]),
                             volumetric_data[:, :, 0] * energy_scale,
                             #extent=(0, 1, 0, 1, 0, 1),
                             representation='wireframe', colormap='RdBu')
        elif dim_idx == 1:
            surf = mlab.mesh(r_grid[:, 0, :] * np.sin(theta_grid[0, 0, 0]) * np.cos(phi_grid[:, 0, :]),
                             r_grid[:, 0, :] * np.sin(theta_grid[0, 0, 0]) * np.sin(phi_grid[:, 0, :]),
                             volumetric_data[:, 0, :],
                             extent=(0, 1, 0, 1, 0, 1),
                             representation='wireframe', colormap='RdBu')
        elif dim_idx == 0:
            surf = mlab.mesh(theta_grid[0, :, :],
                             phi_grid[0, :, :],
                             volumetric_data[0, :, :],
                             extent=(0, 1, 0, 1, 0, 1),
                             representation='wireframe', colormap='RdBu')

        mlab.outline()
        if lowering:
            theta2 = np.linspace(0, 2*np.pi, num=50)
            barrier_lowering = np.array([self.barrier_lowering(theta_i) for theta_i in theta2])
            idx = np.where(barrier_lowering[:, 1] > 0)
            print(barrier_lowering[idx])

            mlab.plot3d(barrier_lowering[idx][:, 1] * np.cos(theta2[idx]),
                        barrier_lowering[idx][:, 1] * np.sin(theta2[idx]),
                        barrier_lowering[idx][:, 0] * energy_scale,
                        tube_radius=energy_scale * 5e-3,
                        color=(1, 0, 0))
            #mlab.outline()
            '''
            mlab.points3d(barrier_lowering[idx][:, 1] * np.cos(theta2[idx]),
                          barrier_lowering[idx][:, 1] * np.sin(theta2[idx]),
                          barrier_lowering[idx][:, 0] * energy_scale)
            mlab.outline()
            '''
        return fig

    def plot_potential_cartesian(self, fig=None, x=None, y=None, z=None):
        if not self.USE_MAYAVI:
            return None
        if fig is None:
            fig = mlab.figure('Potential')#, bgcolor=(0, 0, 0))
            #fig = mlab.figure('Potential', bgcolor=(0, 0, 0))
        else:
            mlab.figure(fig, bgcolor=fig.scene.background)

        if x is None:
            x = np.linspace(1e-9, 5e-8, num=100, endpoint=True)
        if y is None:
            y = np.linspace(1e-9, 5e-8, num=100, endpoint=True)
        if z is None:
            z = np.linspace(1e-9, 5e-8, num=100, endpoint=True)
        x_grid, y_grid, z_grid = np.meshgrid(x, y, z)
        r_grid, theta_grid, phi_grid = Vectors3D.cartesian_to_spherical(x_grid, y_grid, z_grid)
        volumetric_data = self.potential_on_grid(r_grid, theta_grid, phi_grid)
        #volumetric_data = 1 / (x_grid**2 + y_grid**2 + z_grid**2)
        #mlab.contour3d(x_grid, y_grid, z_grid, volumetric_data, colormap='jet')
        #mlab.contour3d(volumetric_data, colormap='jet')

        source = mlab.pipeline.scalar_field(volumetric_data)
        min = volumetric_data.min()
        max = volumetric_data.max()
        vol = mlab.pipeline.volume(source)
        #vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
        #                           vmax=min + 0.9 * (max - min))


        mlab.outline()
        return fig

    def plot_field(self, fig=None, x=None, y=None, z=None):
        if not self.USE_MAYAVI:
            return None
        if fig is None:
            fig = mlab.figure('Field')#, bgcolor=(0, 0, 0))
            #fig = mlab.figure('Field', bgcolor=(0, 0, 0))
        else:
            mlab.figure(fig, bgcolor=fig.scene.background)

        if x is None:
            x = np.linspace(-20e-9, 20e-9, num=10, endpoint=True)
        if y is None:
            y = np.linspace(-20e-9, 20e-9, num=10, endpoint=True)
        if z is None:
            z = np.linspace(-20e-9, 20e-9, num=10, endpoint=True)
        x_grid, y_grid, z_grid = np.meshgrid(x, y, z)
        F_x, F_y, F_z = self.field_cartesian(x, y, z)
        mlab.quiver3d(x_grid * 1e9, y_grid * 1e9, z_grid * 1e9, F_x, F_y, F_z,
                      scale_factor=1e-7)
        mlab.outline()
        return fig


class SuperposedPotential(GenericPotential):

    def __init__(self, label='Unknown Superposed Potential', potentials=None):
        super(SuperposedPotential, self).__init__(label)
        if potentials is None:
            potentials = []
        self.potentials = potentials

    def potential(self, r, theta, phi):
        superposition = np.zeros((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float)
        for potential in self.potentials:
            superposition += potential.potential(r, theta, phi)
        return superposition

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        superposition = np.zeros_like(r_grid)
        for potential in self.potentials:
            superposition += potential.potential_on_grid(r_grid, theta_grid, phi_grid)
        return superposition

    def field(self, r, theta, phi):
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        x_grid, y_grid, z_grid = Vectors3D.spherical_to_cartesian(r_grid, theta_grid, phi_grid)
        super_x = np.zeros_like(x_grid)
        super_y = np.zeros_like(y_grid)
        super_z = np.zeros_like(z_grid)
        for potential in self.potentials:
            F_r, F_theta, F_phi = potential.field(r, theta, phi)
            F_x, F_y, F_z = Vectors3D.spherical_to_cartesian(F_r, F_theta, F_phi)
            super_x += F_x
            super_y += F_y
            super_z += F_z
        super_r, super_theta, super_phi = Vectors3D.cartesian_to_spherical(super_x, super_y, super_z)
        return super_r, super_theta, super_phi

    def field_cartesian(self, x, y, z):
        x_grid, y_grid, z_grid = np.meshgrid(x, y, z)
        super_x = np.zeros_like(x_grid)
        super_y = np.zeros_like(y_grid)
        super_z = np.zeros_like(z_grid)
        for potential in self.potentials:
            F_x, F_y, F_z = potential.field_cartesian(x, y, z)
            super_x += F_x
            super_y += F_y
            super_z += F_z
        return super_x, super_y, super_z

    def get_potential_by_name(self, name):
        for potential in self.potentials:
            if potential.label == name:
                return potential
        return None

    def barrier_lowering(self, theta):
        guess_r = 1e-9
        warnings.filterwarnings('ignore')
        def equation(r):
            field = self.field(r, theta, 0)
            field_x, field_y, field_z = Vectors3D.spherical_to_cartesian(field[0][0][0], field[1][0][0], field[2][0][0])
            r_x, r_y, r_z = Vectors3D.spherical_to_cartesian(r, theta, 0)
            F_r = field_x * r_x + field_y * r_y + field_z * r_z
            return F_r

        solution = root(equation, guess_r)
        warnings.resetwarnings()
        if solution.success:
            r0 = solution.x[0]
            delta_phi = self.potential(r0, theta, 0)[0][0][0]
            #print 'left:', r0, delta_phi
        else:
            r0 = 0
            delta_phi = 0
        return delta_phi, r0


class ConstantPotential(GenericPotential):

    def __init__(self, label='Unknown Constant Potential', U=0):
        super(ConstantPotential, self).__init__(label)
        self.U = U

    def potential(self, r, theta, phi):
        return np.ones((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float) * self.U

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        return np.ones_like(r_grid, dtype=np.float) * self.U

    def field(self, r, theta, phi):
        F_r = np.zeros((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float)
        F_theta = np.zeros((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float)
        F_phi = np.zeros((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float)
        return F_r, F_theta, F_phi

    def field_cartesian(self, x, y, z):
        F_x = np.zeros((np.array(x).size, np.array(y).size, np.array(z).size), dtype=np.float)
        F_y = np.zeros((np.array(x).size, np.array(y).size, np.array(z).size), dtype=np.float)
        F_z = np.zeros((np.array(x).size, np.array(y).size, np.array(z).size), dtype=np.float)
        return F_x, F_y, F_z


class ConstantFieldPotential(GenericPotential):

    def __init__(self, label='Unknown Linear Potential', external_field=(0, 0, 0)):
        super(ConstantFieldPotential, self).__init__(label)
        self.external_field = external_field

    def potential(self, r, theta, phi):
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        if self.external_field[0] == 0:
            return np.zeros_like(r_grid)
        F_r = np.ones_like(r_grid) * self.external_field[0]
        F_theta = np.ones_like(r_grid) * self.external_field[1]
        F_phi = np.ones_like(r_grid) * self.external_field[2]
        x, y, z = Vectors3D.spherical_to_cartesian(r_grid, theta_grid, phi_grid)
        F_x, F_y, F_z = Vectors3D.spherical_to_cartesian(F_r, F_theta, F_phi)
        angles_grid = Vectors3D.angle_between_vectors((x, y, z), (F_x, F_y, F_z))
        return F_r * r_grid * np.cos(angles_grid)

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        if self.external_field[0] == 0:
            return np.zeros_like(r_grid)
        F_r = np.ones_like(r_grid) * self.external_field[0]
        F_theta = np.ones_like(r_grid) * self.external_field[1]
        F_phi = np.ones_like(r_grid) * self.external_field[2]
        x, y, z = Vectors3D.spherical_to_cartesian(r_grid, theta_grid, phi_grid)
        F_x, F_y, F_z = Vectors3D.spherical_to_cartesian(F_r, F_theta, F_phi)
        angles_grid = Vectors3D.angle_between_vectors((x, y, z), (F_x, F_y, F_z))
        return F_r * r_grid * np.cos(angles_grid)

    def field(self, r, theta, phi):
        grid = np.ones((np.array(r).size, np.array(theta).size, np.array(phi).size), dtype=np.float)
        F_r = grid * self.external_field[0]
        F_theta = grid * self.external_field[1]
        F_phi = grid * self.external_field[2]
        return F_r, F_theta, F_phi

    def field_cartesian(self, x, y, z):
        grid = np.ones((np.array(x).size, np.array(y).size, np.array(z).size), dtype=np.float)
        ex, ey, ez = Vectors3D.spherical_to_cartesian(self.external_field[0],
                                                      self.external_field[1],
                                                      self.external_field[2])
        F_x = grid * ex
        F_y = grid * ey
        F_z = grid * ez
        return F_x, F_y, F_z


class HyperbolicPotential(GenericPotential):

    def __init__(self, label='Unknown Hyperbolic Potential', a=1, attractive=True):
        super(HyperbolicPotential, self).__init__(label)
        self.a = a
        self.attractive = attractive

    def potential(self, r, theta, phi):
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        potential_sign = -1 if self.attractive else 1
        return potential_sign * self.a / np.abs(r_grid)

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        potential_sign = -1 if self.attractive else 1
        return potential_sign * self.a / np.abs(r_grid)

    def field(self, r, theta, phi):
        potential_sign = -1 if self.attractive else 1
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        F_r = -potential_sign * self.a / (r_grid**2)
        F_theta = theta_grid
        F_phi = phi_grid
        return F_r, F_theta, F_phi

    def field_cartesian(self, x, y, z):
        potential_sign = -1 if self.attractive else 1
        x_grid, y_grid, z_grid = np.meshgrid(x, y, z)
        r_grid, theta_grid, phi_grid = Vectors3D.cartesian_to_spherical(x_grid, y_grid, z_grid)
        F_r = -potential_sign * self.a / r_grid**2
        F_x, F_y, F_z = Vectors3D.spherical_to_cartesian(F_r, theta_grid, phi_grid)
        return F_x, F_y, F_z


class CoulombPotential(HyperbolicPotential):

    def __init__(self, label='Unknown Coulomb Potential', epsilon=1, charge=elementary_charge, attractive=True):
        self.k = 1 / (4 * np.pi * epsilon_0 * epsilon)
        self.charge = charge
        super(CoulombPotential, self).__init__(label, self.k * self.charge, attractive)


class DislocationDeformationPotential(HyperbolicPotential):

    def __init__(self, label='Deformation Potential', deformation_modulus=0, burgers_vector=1):
        self.deformation_modulus = deformation_modulus
        self.burgers_vector = burgers_vector
        a = deformation_modulus * burgers_vector / (4 * np.pi)
        super(DislocationDeformationPotential, self).__init__(label, a)


class InfiniteCylinderPotential(GenericPotential):

    def __init__(self, label='Unknown Potential', a=1, radius=1):
        super(InfiniteCylinderPotential, self).__init__(label)
        self.radius = radius
        self.a = a

    def set_a(self, a):
        self.a = a

    def potential(self, r, theta, phi):
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        return self.a * np.log(np.abs(r_grid) / self.radius)

    def potential_on_grid(self, r_grid, theta_grid, phi_grid):
        return self.a * np.log(np.abs(r_grid) / self.radius)

    def field(self, r, theta, phi):
        r_grid, theta_grid, phi_grid = np.meshgrid(r, theta, phi)
        F_r = self.a / r_grid
        F_theta = theta_grid
        F_phi = phi_grid
        return F_r, F_theta, F_phi

    def field_cartesian(self, x, y, z):
        x_grid, y_grid, z_grid = np.meshgrid(x, y, z)
        r_grid, theta_grid, phi_grid = Vectors3D.cartesian_to_spherical(x_grid, y_grid, z_grid)
        F_r = self.a / r_grid
        F_x, F_y, F_z = Vectors3D.spherical_to_cartesian(F_r, theta_grid, phi_grid)
        return F_x, F_y, F_z


class ChargedCylinderPotential(InfiniteCylinderPotential):

    def __init__(self, label='Unknown Potential', epsilon=1, charge_sign=1, linear_charge_density=1, radius=1):
        self.epsilon = epsilon
        self.charge_sign = charge_sign
        self.linear_charge_density = linear_charge_density
        a = charge_sign * elementary_charge * linear_charge_density / (2 * np.pi * epsilon_0 * epsilon)
        super(ChargedCylinderPotential, self).__init__(label, a, radius)

    def set_linear_charge_density(self, linear_charge_density):
        self.linear_charge_density = linear_charge_density
        a = self.charge_sign * elementary_charge * self.linear_charge_density / (2 * np.pi * epsilon_0 * self.epsilon)
        self.set_a(a)
