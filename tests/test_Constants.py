from __future__ import division, print_function
from scipy.constants import c, N_A, k, e, m_e, epsilon_0

from Schottky.Constants import Constants
from Schottky import constant

import unittest


class TestConstants(unittest.TestCase):

    def setUp(self):
        pass

    def test_constants(self):
        constants = Constants()
        self.assertEqual(constants.c, c)
        self.assertEqual(constants.avogadro, N_A)
        self.assertEqual(constants.k, k)
        self.assertEqual(constants.q, e)
        self.assertEqual(constants.m_e, m_e)
        self.assertEqual(constants.epsilon_0, epsilon_0)
        self.assertEqual(constants.A_R, 1.20173e6)

    def test_constant(self):
        self.assertEqual(constant.c, c)
        self.assertEqual(constant.avogadro, N_A)
        self.assertEqual(constant.k, k)
        self.assertEqual(constant.q, e)
        self.assertEqual(constant.m_e, m_e)
        self.assertEqual(constant.epsilon_0, epsilon_0)
        self.assertEqual(constant.A_R, 1.20173e6)
