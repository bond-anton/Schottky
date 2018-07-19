from __future__ import division, print_function

from Schottky.Notation import q
from Schottky import constants
from Schottky.Helpers import to_numeric


class Metal(object):
    '''
    Metal electrode class
    '''

    def __init__(self, label, WF):
        '''
        Constructor
        '''
        self.label = label
        self.WF = WF

    def __str__(self, *args, **kwargs):
        return 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (str(self.label), to_numeric(self.WF / q),
                                                                to_numeric(self.WF))
