from __future__ import division, print_function
# import timeit
import numbers
import numpy as np

from Schottky.Samples.Trap import Trap
from Schottky.Simulators import Simulator
from Schottky.Simulators.Field import FieldSimulator
from Schottky import constants


class ChargeCarrierTrap(Simulator):

    def __init__(self, client, trap, description=None):
        assert isinstance(trap, Trap), 'Valid Trap Sample object expected'
        self.trap = trap
        samples = [self.trap]
        name = 'Charge Carrier Trap Simulator'
        label = name + ' [%s]' % trap.name
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}
        measurement_types = [{'name': 'Capture and Emission Kinetics',
                              'description': 'Measurement of charge carrier trap capture and emission kinetics',
                              'children': []}]
        measurements = [
            {'name': 'Capture rate',
             'description': 'measure charge carriers capture rate',
             'type': 'Capture and Emission Kinetics'},
            {'name': 'Emission rate',
             'description': 'measure charge carriers emission rate',
             'type': 'Capture and Emission Kinetics'},
        ]
        parts = None
        if trap.trap_potential:
            parts = [FieldSimulator(client=client, field=trap.trap_potential)]
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category, label=label,
            measurement_types=measurement_types,
            measurements=measurements)

    def energy_level(self, band_gap):
        if isinstance(band_gap, (list, tuple, np.ndarray)):
            energy = np.ones_like(np.array(band_gap))
        else:
            assert isinstance(band_gap, numbers.Number), 'Provide a numeric value of band_gap'
            energy = 1
        if self.trap.band == 'Ec':
            return energy * self.trap.activation_energy
        else:
            return band_gap - energy * self.trap.activation_energy

    def capture_rate(self, temperature, n, p, v_n, v_p):
        sigma_n, sigma_p = self.trap.capture_cross_section(temperature)
        capture_rate_n = sigma_n * v_n * n
        capture_rate_p = sigma_p * v_p * p
        return capture_rate_n, capture_rate_p

    def emission_rate(self, temperature, band_gap, v_n, v_p, n_c, n_v,
                      poole_frenkel_n=1.0, poole_frenkel_p=1.0):
        energy_scale = constants['k'] * temperature
        energy_level = self.energy_level(band_gap=band_gap)
        activation_energy_n = energy_level
        activation_energy_p = band_gap - energy_level
        sigma_n, sigma_p = self.trap.capture_cross_section(temperature)
        g_ratio_n = 1.0
        g_ratio_p = 1.0
        factor_n = sigma_n * v_n * n_c * g_ratio_n * poole_frenkel_n
        factor_p = sigma_p * v_p * n_v * g_ratio_p * poole_frenkel_p
        emission_rate_n = factor_n * np.exp(-activation_energy_n / energy_scale)
        emission_rate_p = factor_p * np.exp(-activation_energy_p / energy_scale)
        return emission_rate_n, emission_rate_p
