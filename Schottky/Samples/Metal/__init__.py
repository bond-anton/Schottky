from __future__ import division
import numpy as np

from Schottky import constants
from Schottky.Samples import Sample


class Metal(Sample):

    def __init__(self, name, ):
        self.label = label
        self.WF = WF

    def __str__(self, *args, **kwargs):
        return 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (str(self.label), self.to_numeric(self.WF / q),
                                                                self.to_numeric(self.WF))
