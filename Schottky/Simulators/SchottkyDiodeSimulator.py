from __future__ import division, print_function
import numpy as np

from Schottky.Samples.Diode import SchottkyDiode
from Schottky.Simulators import Simulator
from Schottky.Simulators.BulkSemiconductor import BulkSemiconductor
from Schottky import constants
from ._helpers import storage_manager, prepare_array


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
        measurement_types = [
            {'name': 'Diode energetics',
             'description': 'Measurement of diode energetics',
             'children': []},
            {'name': 'Current Measurement',
             'description': 'Measurement of electrical current',
             'children': []}]
        self.measurement_details = {
            'v_bi': {
                'name': 'Built-in potential temperature dependence',
                'description': 'Measurement of Diode Built-in potential temperature dependence',
                'type': 'Diode energetics'},
            'iv': {
                'name': 'IV',
                'description': 'measure I-V dependence',
                'type': 'Current Measurement'},
        }

        self.measurement_specs = {
            'v_bi': {'parameters': None,
                     'input data': None,
                     'variables': [{
                         'name': 'Temperature',
                         'description': 'Sample temperature',
                         'units': 'K'
                     }],
                     'result': [{
                         'name': 'Built-in potential',
                         'description': 'Diode\'s Built-in potential',
                         'units': 'eV'
                         }]},
        }

        parts = [BulkSemiconductor(client=client, semiconductor=diode.semiconductor)]
        Simulator.__init__(
            self,
            client=client, name=name, description=description,
            samples=samples, parts=parts,
            category=category,
            measurement_types=measurement_types,
            measurements=list(self.measurement_details.values()))

    @storage_manager('v_bi', use_storage=True)
    def v_bi(self, temperature=0.0):
        temperature = prepare_array(temperature)
        metal_wf = self.diode.metal.work_function
        xi = self.parts['Bulk Semiconductor Simulator'].electrochemical_potential(temperature=temperature)
        semiconductor_wf = xi + self.diode.semiconductor.affinity
        bi_voltage = metal_wf - semiconductor_wf
        return bi_voltage

    def thermionic_emission_current(self, bias, temperature):
        pass
