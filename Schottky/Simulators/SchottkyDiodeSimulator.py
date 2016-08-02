from __future__ import division, print_function
# import timeit
import numbers
import numpy as np

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators import Simulator
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky import constants


class SchottkyDiodeSimulator(Simulator):

    def __init__(self, client, diode, description=None):
        assert isinstance(diode, SchottkyDiode), 'Valid SchottkyDiode Sample object expected'
        self.diode = diode
        samples = [self.diode]
        name = 'Schottky Diode Simulator'
        category = {
            'name': 'Software',
            'description': 'Measurement, automation, control, simulation, and other software tools',
            'subcategory': {'name': 'Simulation',
                            'description': 'Simulation software',
                            'subcategory': None}}
        measurement_types = [{'name': 'Current Measurement',
                              'description': 'Measurement of electrical current',
                              'children': []}]
        measurements = [
            {'name': 'IV',
             'description': 'measure I-V dependence',
             'type': 'Current Measurement'},
            {'name': 'Emission rate',
             'description': 'measure charge carriers emission rate',
             'type': 'Capture and Emission Kinetics'},
        ]
        parts = [BulkSemiconductor(client=client, semiconductor=diode.semiconductor)]
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=measurements)

    def v_bi(self, temperature):
        metal_wf = self.diode.metal.work_function
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        semiconductor_wf = xi + self.diode.semiconductor.affinity
        return metal_wf - semiconductor_wf
