from __future__ import division, print_function

from Schottky import constant
from Schottky.Metal import Metal

import unittest


class TestMetal(unittest.TestCase):

    def setUp(self):
        pass

    def test_label(self):
        m = Metal('Au', 5.1 * constant.q)
        self.assertEqual(m.label, 'Au')
        m.label = 'Al'
        self.assertEqual(m.label, 'Al')
        with self.assertRaises(TypeError):
            m.label = 3
        with self.assertRaises(TypeError):
            Metal(3, 5.1 * constant.q)

    def test_work_function(self):
        m = Metal('Au', 5.1 * constant.q)
        self.assertEqual(m.work_function, 5.1 * constant.q)
        self.assertEqual(m.work_function_ev, 5.1)
        m.work_function = 4.8 * constant.q
        self.assertEqual(m.work_function, 4.8 * constant.q)
        self.assertEqual(m.work_function_ev, 4.8)
        m.work_function_ev = 5.1
        self.assertEqual(m.work_function, 5.1 * constant.q)
        self.assertEqual(m.work_function_ev, 5.1)
        with self.assertRaises(TypeError):
            m.work_function = 'xxx'
        with self.assertRaises(TypeError):
            Metal('Au', 'xxx')

    def test_str(self):
        m = Metal('Au', 5.1 * constant.q)
        s = 'Metal: %s, Workfunction: %2.2f eV (%2.2g J)' % (m.label,
                                                             m.work_function_ev,
                                                             m.work_function)
        self.assertEqual(str(m), s)
