from __future__ import division, print_function
import numpy as np

from Schottky import constant
from Schottky.Trap import Trap

import unittest


class TestTrap(unittest.TestCase):

    def setUp(self):
        pass

    def test_label(self):
        t = Trap('My trap', True,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)
        self.assertEqual(t.label, 'My trap')
        t.label = 'Another trap'
        self.assertEqual(t.label, 'Another trap')
        with self.assertRaises(TypeError):
            t.label = 3
        with self.assertRaises(TypeError):
            Trap(3, True,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)

    def test_energy(self):
        t = Trap('My trap', True,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)
        np.testing.assert_allclose(t.energy_c, 0.3 * constant.q, atol=1e-15)
        np.testing.assert_allclose(t.energy_v, 0.8 * constant.q, atol=1e-15)
        np.testing.assert_allclose(t.energy_c_ev, 0.3, atol=1e-15)
        np.testing.assert_allclose(t.energy_v_ev, 0.8, atol=1e-15)
        t.energy_c = 0.2 * constant.q
        np.testing.assert_allclose(t.energy_c, 0.2 * constant.q, atol=1e-15)
        t.energy_c_ev = 0.3 * constant.q
        np.testing.assert_allclose(t.energy_c, 0.3 * constant.q, atol=1e-15)
        t.energy_v = 0.7 * constant.q
        np.testing.assert_allclose(t.energy_v, 0.7 * constant.q, atol=1e-15)
        t.energy_v_ev = 0.75 * constant.q
        np.testing.assert_allclose(t.energy_v, 0.75 * constant.q, atol=1e-15)

    def test_cs(self):
        t = Trap('My trap', True,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)
        np.testing.assert_allclose(t.e_cs0, 1e-15)
        np.testing.assert_allclose(t.h_cs0, 1e-15)
        np.testing.assert_allclose(t.e_cs_activation, 0)
        np.testing.assert_allclose(t.h_cs_activation, 0)
        np.testing.assert_allclose(t.e_cs(300.0), t.e_cs0)
        np.testing.assert_allclose(t.h_cs(300.0), t.h_cs0)

    def test_str(self):
        t = Trap('My trap', True,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)
        s = 'Trap: %s\nEc-Et: %2.2f eV (%2.2g J)' % (t.label, t.energy_c_ev, t.energy_c)
        self.assertEqual(str(t), s)
