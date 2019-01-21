from __future__ import division, print_function
import numpy as np


from BDMesh import Mesh1DUniform, TreeMesh1DUniform

from Schottky import constant
from Schottky.Dopant import Dopant

import unittest


class TestDopant(unittest.TestCase):

    def setUp(self):
        pass

    def test_label(self):

        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)

        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)
        self.assertEqual(t.label, 'My trap')
        t.label = 'Another trap'
        self.assertEqual(t.label, 'Another trap')
        with self.assertRaises(TypeError):
            t.label = 3
        with self.assertRaises(TypeError):
            Dopant(3,
                 0.3 * constant.q, 0.8 * constant.q,
                 1e-15, 1e-15)

    def test_energy(self):
        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)
        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
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
        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)

        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                   0.3 * constant.q, 0.8 * constant.q,
                   1e-15, 1e-15)
        np.testing.assert_allclose(t.e_cs0, 1e-15)
        np.testing.assert_allclose(t.h_cs0, 1e-15)
        np.testing.assert_allclose(t.e_cs_activation, 0)
        np.testing.assert_allclose(t.h_cs_activation, 0)
        np.testing.assert_allclose(t.e_cs(300.0), t.e_cs0)
        np.testing.assert_allclose(t.h_cs(300.0), t.h_cs0)

    def test_concentration(self):
        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)

        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                   0.3 * constant.q, 0.8 * constant.q,
                   1e-15, 1e-15)

        query_z = np.linspace(0, 5e-6, num=101)
        result = np.asarray(t.n_t(query_z))
        np.testing.assert_allclose(result, np.ones(101) * 1e15)

        c.solution = np.asarray(c.physical_nodes) * 1e15
        t.concentration = TreeMesh1DUniform(c, aligned=True)
        query_z = np.linspace(0, 5e-6, num=101)
        result = np.asarray(t.n_t(query_z))
        np.testing.assert_allclose(result, query_z * 1e15)

    def test_str(self):
        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)

        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                   0.3 * constant.q, 0.8 * constant.q,
                   1e-15, 1e-15)
        s = 'Dopant: %s\nEc-Et: %2.2f eV (%2.2g J)\nEt-Ev: %2.2f eV (%2.2g J)' % (t.label,
                                                                                  t.energy_c_ev,
                                                                                  t.energy_c,
                                                                                  t.energy_v_ev,
                                                                                  t.energy_v)
        self.assertEqual(str(t), s)

    def test_coerce_occupation_mesh(self):
        c = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        f = Mesh1DUniform(0.0, 5e-6, physical_step=1e-6)
        c.solution = np.ones(c.num) * 1e15
        f.solution = np.zeros(f.num)

        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                   0.3 * constant.q, 0.8 * constant.q,
                   1e-15, 1e-15)
        np.testing.assert_allclose(f.solution, np.zeros(f.num))
        f.solution = np.ones(f.num)
        np.testing.assert_allclose(f.solution, np.ones(f.num))
        f.solution = np.ones(f.num) * 2
        np.testing.assert_allclose(f.solution, np.ones(f.num) * 2)
        t.coerce_mesh_occupation(f)
        np.testing.assert_allclose(f.solution, np.ones(f.num))
        f.solution = -2 * np.ones(f.num)
        np.testing.assert_allclose(f.solution, -2 * np.ones(f.num))
        t.coerce_mesh_occupation(f)
        np.testing.assert_allclose(f.solution, np.zeros(f.num))
        f.solution = np.ones(f.num) * 2
        t = Dopant('My trap', TreeMesh1DUniform(c, aligned=True), TreeMesh1DUniform(f, aligned=True),
                   0.3 * constant.q, 0.8 * constant.q,
                   1e-15, 1e-15)
        np.testing.assert_allclose(t.occupation.tree[0][0].solution, 2 * np.ones(f.num))
        t.coerce_mesh_tree_occupation(t.occupation)
        np.testing.assert_allclose(t.occupation.tree[0][0].solution, np.ones(f.num))
