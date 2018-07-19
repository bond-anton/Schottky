from __future__ import division, print_function
from scipy.constants import c, N_A, k, e, m_e, epsilon_0

from Schottky import Constants

import unittest


class TestHelpers(unittest.TestCase):

    def setUp(self):
        pass

    def test_constants(self):
        constants = Constants()
        print(c)
        print(constants.c)
        self.assertEqual(constants.c, c)
        self.assertEqual(constants.avogadro, N_A)
        self.assertEqual(constants.k, k)
        self.assertEqual(constants.q, e)
        self.assertEqual(constants.m_e, m_e)
        self.assertEqual(constants.epsilon_0, epsilon_0)
        self.assertEqual(constants.A_R, 1.20173e6)
