import math
import numpy as np


class BondingInterface(object):
    """
    Bonded wafers of the same kind semiconductor
    """

    def __init__(self, depth, eps, twist, tilt, dsl_twist, dsl_tilt, label=None):
        """
        Constructor for Bonding Interface
        """
        if label is None:
            self.label = 'Bonding interface @%2.2f nm' % (depth * 1.0e9)
        else:
            self.label = label
        self.depth = depth
        self.smooth_dirac_epsilon = eps
        self.twist = twist
        self.tilt = tilt
        self.dsl_tilt = dsl_tilt
        self.dsl_twist = dsl_twist
        self.tilt_distance = dsl_tilt.b / math.sin(math.radians(tilt))
        self.screw_distance = dsl_twist.b / math.sin(math.radians(twist))
        self.dsl_twist_f = np.zeros(len(self.dsl_twist.traps))
        self.dsl_twist_df = np.zeros(len(self.dsl_twist.traps))
        self.dsl_tilt_f = np.zeros(len(self.dsl_tilt.traps))
        self.dsl_tilt_df = np.zeros(len(self.dsl_tilt.traps))
        self.dsl_twist_total_density_of_states = []
        self.dsl_tilt_total_density_of_states = []
        self.density_of_charge = 0
        self.d_density_of_charge = 0
        for trap in self.dsl_twist.traps:
            self.dsl_twist_total_density_of_states.append(2 * trap[1] / self.screw_distance)
        for trap in self.dsl_tilt.traps:
            self.dsl_tilt_total_density_of_states.append(trap[1] / self.tilt_distance)

    def set_traps_f(self, tilt_f, twist_f):
        if tilt_f.size == self.dsl_tilt_f.size:
            self.dsl_tilt_f = tilt_f
        if twist_f.size == self.dsl_twist_f.size:
            self.dsl_twist_f = twist_f
        density_of_charge = 0
        for i, trap in enumerate(self.dsl_tilt.traps):
            density_of_charge += self.dsl_tilt_total_density_of_states[i] \
                                 * ((trap[0].charge_states[1][0] - trap[0].charge_states[0][0])
                                    * self.dsl_tilt_f[i] + trap[0].charge_states[0][0])
        for i, trap in enumerate(self.dsl_twist.traps):
            density_of_charge += self.dsl_twist_total_density_of_states[i] \
                                 * ((trap[0].charge_states[1][0] - trap[0].charge_states[0][0])
                                    * self.dsl_twist_f[i] + trap[0].charge_states[0][0])
        self.density_of_charge = density_of_charge

    def set_traps_df(self, tilt_df, twist_df):
        if tilt_df.size == self.dsl_tilt_df.size:
            self.dsl_tilt_df = tilt_df
        if twist_df.size == self.dsl_twist_df.size:
            self.dsl_twist_df = twist_df
        d_density_of_charge = 0
        for i, trap in enumerate(self.dsl_tilt.traps):
            d_density_of_charge += self.dsl_tilt_total_density_of_states[i] * (
                (trap[0].charge_states[1][0] - trap[0].charge_states[0][0])
                * self.dsl_tilt_df[i] + 0 * trap[0].charge_states[0][0])
        for i, trap in enumerate(self.dsl_twist.traps):
            d_density_of_charge += self.dsl_twist_total_density_of_states[i] * (
                (trap[0].charge_states[1][0] - trap[0].charge_states[0][0])
                * self.dsl_twist_df[i] + 0 * trap[0].charge_states[0][0])
        self.d_density_of_charge = np.float(d_density_of_charge)

    def energy_distribution_diagram(self, ax, temperature, semiconductor, electron_volts=True, fancy_labels=False):
        dsl_distance = [self.screw_distance, self.tilt_distance]
        for i, dislocation in enumerate([self.dsl_twist, self.dsl_tilt]):
            for trap in dislocation.traps:
                trap[0].energy_distribution_diagram(ax, temperature, semiconductor,
                                                    trap_concentration=(2 - i) * trap[1] / dsl_distance[i] * 1e-4, trap_concentration_units='1/cm^2',
                                                    electron_volts=electron_volts, fancy_labels=fancy_labels)

    def __str__(self, *args, **kwargs):
        description = 'Bonding interface @%2.2f nm' % (self.depth * 1.0e9) + '\n'
        description += 'Twist misorientation angle: %2.2f, dist = %2.2f nm' % (
            self.twist, self.screw_distance * 1e9) + '\n'
        description += 'Tilt misorientation angle: %2.2f, dist = %2.2f nm' % (
            self.tilt, self.tilt_distance * 1e9) + '\n'
        return description
