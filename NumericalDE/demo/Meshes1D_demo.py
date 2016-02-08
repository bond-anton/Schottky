#-*- coding: utf-8 -*-
"""
Created on 31 мая 2015 г.

@author: anton
"""

import numpy as np
from NumericalDE.Mesh import Uniform1DMeshesTree, UniformMesh1D

root_mesh = UniformMesh1D(0.0, 10.0, 1, bc1=1, bc2=0)
child_mesh_1_1 = UniformMesh1D(0.0, 9.0, 0.5, bc1=1, bc2=0)
#child_mesh_1_2 = UniformMesh1D(3.0, 17.0, 0.5, bc1=1, bc2=0)
child_mesh_1_3 = UniformMesh1D(6.0, 8.0, 0.5, bc1=1, bc2=0)
child_mesh_2_1 = UniformMesh1D(1.0, 2.0, 0.25, bc1=1, bc2=0)
child_mesh_2_2 = UniformMesh1D(2.0, 5.0, 0.25, bc1=1, bc2=0)
child_mesh_3_1 = UniformMesh1D(1.0, 1.5, 0.125, bc1=1, bc2=0)
Meshes = Uniform1DMeshesTree(root_mesh)
Meshes.add_mesh(child_mesh_1_1)
print Meshes.Tree
#Meshes.add_mesh(child_mesh_1_2)
Meshes.add_mesh(child_mesh_1_3)
print Meshes.Tree
Meshes.add_mesh(child_mesh_2_1)
Meshes.add_mesh(child_mesh_2_2)
Meshes.add_mesh(child_mesh_3_1)

print Meshes.Tree
Meshes.plot_tree()
