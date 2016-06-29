# -*- coding: utf-8 -*-

'''
Created on 1 дек. 2014 г.

@author: anton
'''

from __future__ import division
import sympy as sym

from Schottky.Notation import q
from Schottky import constants
from Schottky.Helpers import to_numeric


class Metal(object):
    '''
    Metal electrode class
    '''

    equipment = 'Metal electrode Simulator'
    description = 'Simulates properties of a metal electrode'

    def __init__(self, label, WF):
        '''
        Constructor
        '''
        self.label = label
        self.WF = WF

    def __str__(self, *args, **kwargs):
        return 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (str(self.label), self.to_numeric(self.WF / q),
                                                                self.to_numeric(self.WF))
