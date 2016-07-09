# -*- coding: utf-8 -*-
"""
Created on 30 мая 2015 г.

@author: anton
"""

from __future__ import division
import math as m

import numpy as np
from matplotlib import pyplot as plt

from Schottky.Helpers import check_if_integer


class UniformMesh1D:
    """
    One dimensional grid
    """

    def __init__(self, phys_boundary1, phys_boundary2, phys_step, bc1, bc2, crop=None):
        if not crop:
            crop = [0, 0]
        if phys_boundary1 < phys_boundary2:
            self.phys_boundary1 = phys_boundary1
            self.phys_boundary2 = phys_boundary2
            self.bc1 = bc1
            self.bc2 = bc2
        else:
            self.phys_boundary1 = phys_boundary2
            self.phys_boundary2 = phys_boundary1
            self.bc1 = bc2
            self.bc2 = bc1
        self.crop = np.array(crop)
        self.phys_step = abs(phys_step)
        self.num = int(np.ceil((self.phys_boundary2 - self.phys_boundary1) / self.phys_step) + 1)
        if self.phys_boundary1 + (self.num - 1) * self.phys_step > self.phys_boundary2:
            self.num -= 1
        self.local_nodes, self.local_step = np.linspace(0.0, 1.0, num=self.num, endpoint=True, retstep=True)
        self.J = self.phys_step / self.local_step
        self.solution = np.zeros(self.num)
        self.residual = np.zeros(self.num)
        self.int_residual = 0

    def trim(self, debug=False):
        if debug:
            print 'Cropping', np.sum(self.crop), 'elements'
            print 'phys_boundary1', self.phys_boundary1
            print 'phys_boundary2', self.phys_boundary2
        self.phys_boundary1 += self.crop[0] * self.phys_step
        self.phys_boundary2 -= self.crop[1] * self.phys_step
        if debug:
            print 'phys_boundary1', self.phys_boundary1
            print 'phys_boundary2', self.phys_boundary2
        self.num = int(np.ceil(self.num - np.sum(self.crop)))
        self.local_nodes, self.local_step = np.linspace(0.0, 1.0, num=self.num, endpoint=True, retstep=True)
        self.J = self.phys_step / self.local_step
        self.solution = self.solution[self.crop[0]:self.solution.size - self.crop[1]]
        self.residual = self.residual[self.crop[0]:self.residual.size - self.crop[1]]
        self.bc1 = self.solution[0]
        self.bc2 = self.solution[-1]
        self.crop = np.array([0, 0])

    def to_phys(self, x):
        return self.phys_boundary1 + self.J * x

    def phys_nodes(self):
        return self.to_phys(self.local_nodes)

    def to_local(self, x):
        return (x - self.phys_boundary1) / self.J

    def local_f(self, f, args=None):
        def f_l(x, arguments=args):
            if arguments is not None:
                return f(self.to_phys(x), arguments)
            else:
                return f(self.to_phys(x))

        return f_l

    def is_inside_of(self, mesh):
        if mesh.phys_boundary1 <= self.phys_boundary1 and mesh.phys_boundary2 >= self.phys_boundary2:
            return True
        else:
            return False

    def inner_mesh_indexes(self, mesh):
        if mesh.is_inside_of(self):
            local_start = self.to_local(mesh.phys_boundary1)
            local_stop = self.to_local(mesh.phys_boundary2)
            idx1 = np.where(abs(self.local_nodes - local_start) < self.local_step / 2)[0][0]
            idx2 = np.where(abs(self.local_nodes - local_stop) < self.local_step / 2)[0][0]
            return [idx1, idx2]
        else:
            return [None, None]

    def overlap_with(self, mesh):
        if self.is_inside_of(mesh) or mesh.is_inside_of(self):
            return True
        elif mesh.phys_boundary1 <= self.phys_boundary1 <= mesh.phys_boundary2:
            return True
        elif mesh.phys_boundary2 >= self.phys_boundary2 >= mesh.phys_boundary1:
            return True
        else:
            return False

    # noinspection PyUnresolvedReferences
    def is_aligned_with(self, mesh):
        min_size = min(mesh.num, self.num)
        min_step = min(mesh.phys_step, self.phys_step)
        max_step = max(mesh.phys_step, self.phys_step)
        step_ratio = max_step / min_step
        # print step_ratio, m.rnd(step_ratio), abs(m.floor(step_ratio) - step_ratio)
        # if abs(m.floor(step_ratio) - step_ratio) < 1e-6:#step_ratio.is_integer():
        if check_if_integer(step_ratio, 1e-8):
            shift = (mesh.to_phys(mesh.local_nodes)[:min_size] - self.to_phys(self.local_nodes)[:min_size]) / min_step
            if shift[0].is_integer() or abs(round(shift[0]) - shift[0]) < 1e-6:
                return True
            else:
                print 'SHIFT!!! %3.16f' % shift[0]
                return False
        else:
            print abs(m.floor(step_ratio) - step_ratio)
            return False

    def merge_with(self, mesh):
        if self.overlap_with(mesh) and abs(self.phys_step - mesh.phys_step) < 0.1 * self.phys_step:
            if self.is_aligned_with(mesh):
                self.bc1 = self.bc1 if self.phys_boundary1 <= mesh.phys_boundary1 else mesh.bc1
                self.bc2 = self.bc2 if self.phys_boundary2 >= mesh.phys_boundary2 else mesh.bc2
                self.crop[0] = self.crop[0] if self.phys_boundary1 <= mesh.phys_boundary1 else mesh.crop[0]
                self.crop[1] = self.crop[1] if self.phys_boundary2 >= mesh.phys_boundary2 else mesh.crop[1]
                phys_nodes = np.union1d(self.phys_nodes(), mesh.phys_nodes())
                self.phys_boundary1 = phys_nodes[0]
                self.phys_boundary2 = phys_nodes[-1]
                self.num = int((self.phys_boundary2 - self.phys_boundary1) / self.phys_step + 1)
                self.solution = np.zeros(self.num)
                self.residual = np.zeros(self.num)
                self.local_nodes, self.local_step = np.linspace(0.0, 1.0, num=self.num, endpoint=True, retstep=True)
                self.J = self.phys_step / self.local_step
            else:
                print 'meshes are not aligned and could not be merged'
        else:
            print abs(self.phys_step - mesh.phys_step), self.phys_step
            print 'meshes do not overlap or have different step size'


class Uniform1DMeshesTree(object):
    """
    Manages the tree of meshes
    """

    def __init__(self, root_mesh, refinement_coefficient=2, aligned=True, crop=None):
        if not crop:
            crop = [0, 0]
        self.Tree = {0: [root_mesh]}
        self.refinement_coefficient = refinement_coefficient
        self.aligned = aligned
        self.levels = self.Tree.keys()
        self.crop = np.array(crop)

    def get_root_mesh(self):
        return self.Tree[0][0]

    def add_mesh(self, mesh):
        level = m.log(self.get_root_mesh().phys_step / mesh.phys_step, self.refinement_coefficient)
        # print level, np.floor(level), int(m.floor(level)), abs(m.floor(level) - level), level.is_integer()
        # if abs(m.floor(level) - level) > 1e-10:
        if not check_if_integer(level, 1e-10):
            raise Exception('all child meshes must have step multiple to the root mesh with refinement coefficient')
        if self.aligned and not self.get_root_mesh().is_aligned_with(mesh):
            raise Exception('all child meshes must be aligned with the root mesh')
        # level = int(m.floor(level))
        lower_estimation = m.floor(level)
        upper_estimation = m.ceil(level)
        if abs(lower_estimation - level) < abs(upper_estimation - level):
            level = int(lower_estimation)
        else:
            level = int(upper_estimation)
        try:
            self.Tree[level].append(mesh)
        except KeyError:
            self.Tree[level] = [mesh]
        self.cleanup()

    def get_mesh_level(self, mesh):
        for level in self.Tree.keys():
            for tree_mesh in self.Tree[level]:
                if tree_mesh == mesh:
                    return level
        print 'mesh not found in a tree'
        return -1

    def get_children(self, mesh):
        children = {}
        level = self.get_mesh_level(mesh)
        upper_levels = np.array(self.Tree.keys())
        upper_levels = upper_levels[np.where(upper_levels > level)]
        for upper_level in upper_levels:
            for tree_mesh in self.Tree[upper_level]:
                if tree_mesh.is_inside_of(mesh):
                    try:
                        children[upper_level].append(tree_mesh)
                    except KeyError:
                        children[upper_level] = [tree_mesh]
        return children

    def del_mesh(self, mesh, del_children=True):
        level = self.get_mesh_level(mesh)
        if level > 0:
            children = self.get_children(mesh)
            if children == {}:
                self.Tree[level].remove(mesh)
            elif del_children:
                for child_level in sorted(children.keys(), reverse=True):
                    print 'deleting children at level', child_level
                    for child in children[child_level]:
                        self.del_mesh(child, del_children=False)
                self.del_mesh(mesh, del_children=False)
            else:
                print 'mesh has children, use del_children=True flag'
        elif level == 0:
            print 'Can not delete root mesh'
        else:
            print 'mesh not found in a tree'
        self.cleanup()

    def trim(self, debug=False):
        if debug:
            print 'Cropping', np.sum(self.crop), 'elements (of root_mesh):', self.crop
        self.get_root_mesh().crop = self.crop
        self.get_root_mesh().trim()
        level = 1
        trimmed = True if level > self.levels[-1] else False
        while not trimmed:
        # for level in self.levels[1:]:
            if debug:
                print 'trimming level', level
            meshes_for_delete = []
            for mesh in self.Tree[level]:
                mesh.trim()
                crop = [0, 0]
                # print mesh.phys_boundary1*1e6, self.get_root_mesh().phys_boundary1*1e6
                # print mesh.phys_boundary2*1e6, self.get_root_mesh().phys_boundary2*1e6
                left_offset = (self.get_root_mesh().phys_boundary1 - mesh.phys_boundary1) / mesh.phys_step
                right_offset = (mesh.phys_boundary2 - self.get_root_mesh().phys_boundary2) / mesh.phys_step
                # print left_offset, right_offset
                crop[0] = int(m.ceil(left_offset)) if left_offset > 0 else 0
                crop[1] = int(m.ceil(right_offset)) if right_offset > 0 else 0
                if crop[0] == 0 and crop[1] > 0:
                    if crop[1] >= mesh.num:
                        if debug:
                            print 'Deleting mesh'
                        meshes_for_delete.append(mesh)
                        continue
                elif crop[1] == 0 and crop[0] > 0:
                    if crop[0] >= mesh.num:
                        if debug:
                            print 'Deleting mesh'
                        meshes_for_delete.append(mesh)
                        continue
                if debug:
                    print 'Cropping mesh by', crop
                mesh.crop = np.array(crop)
                mesh.trim()
            for mesh in meshes_for_delete:
                self.del_mesh(mesh, del_children=True)
            level += 1
            if level > self.levels[-1]:
                trimmed = True
        self.crop = [0, 0]

    def cleanup(self):
        self.merge_overlaps()
        # self.remove_coarse_duplicates()
        self.recalculate_levels()

    def remove_coarse_duplicates(self):
        for level in self.Tree.keys():
            for mesh in self.Tree[level]:
                upper_levels = np.array(self.Tree.keys())
                upper_levels = upper_levels[np.where(upper_levels > level)]
                for upper_level in upper_levels:
                    for tree_mesh in self.Tree[upper_level]:
                        if mesh.is_inside_of(tree_mesh):
                            self.Tree[level].remove(mesh)
                            print 'mesh overlaps with coarse mesh. DELETING COARSE.'

    def recalculate_levels(self):
        tidy = False
        while not tidy:
            for level in self.Tree.keys():
                if not self.Tree[level]:
                    self.Tree.pop(level)
                    break
            tidy = True
        offset = min(self.Tree.keys())
        if offset != 0:
            for level in self.Tree.keys():
                self.Tree[level - offset] = self.Tree.pop(level)
        self.levels = self.Tree.keys()

    def merge_overlaps(self, debug=False):
        # print 'checking overlaps'
        level = 0
        while level < len(self.Tree.keys()):
            overlap_found = False
            i_list = []
            for i in range(len(self.Tree[level])):
                i_list.append(i)
                for j in range(len(self.Tree[level])):
                    if j not in i_list:
                        # print i, j
                        if self.Tree[level][i].overlap_with(self.Tree[level][j]):
                            if debug:
                                print 'meshes overlap. MERGING.'
                            overlap_found = True
                            self.Tree[level][i].merge_with(self.Tree[level][j])
                            self.Tree[level].remove(self.Tree[level][j])
                            break
                if overlap_found:
                    break
            if not overlap_found:
                level += 1

    def plot_tree(self, ax=None):
        colors = ['radius', 'g', 'b', 'c', 'm', 'y', 'k', 'radius', 'g', 'b', 'c', 'm', 'y', 'k',
                  'radius', 'g', 'b', 'c', 'm', 'y', 'k', 'radius', 'g', 'b', 'c', 'm', 'y', 'k']
        styles = ['-', ':', '--', '-', ':', '--', '-', '-', ':', '--', '-', ':', '--', '-',
                  '-', ':', '--', '-', ':', '--', '-', '-', ':', '--', '-', ':', '--', '-']
        show = False
        if ax is None:
            _, ax = plt.subplots()
            show = True
        for level in self.Tree.keys():
            for i, mesh in enumerate(self.Tree[level]):
                ax.plot(mesh.to_phys(mesh.local_nodes), np.ones(mesh.num) * level, colors[level] + styles[i] + 'o')
        ax.set_ylim([-1, max(self.Tree.keys()) + 1])
        ax.grid()
        if show:
            plt.show()

    def flatten(self, debug=False):
        flat_grid = self.get_root_mesh().phys_nodes()
        flat_sol = self.get_root_mesh().solution
        flat_res = self.get_root_mesh().residual
        if debug:
            print 'root_mesh is from', flat_grid[0] * 1e6, 'to', flat_grid[-1] * 1e6
        for level in self.levels[1:]:
            if debug:
                print 'working with level', level
            for mesh in self.Tree[level]:
                if debug:
                    print 'flat_grid is from', flat_grid[0] * 1e6, 'to', flat_grid[-1] * 1e6
                    print 'merging mesh from', mesh.phys_boundary1 * 1e6,
                    print 'to', mesh.phys_boundary2 * 1e6, mesh.phys_step * 1e6
                ins_idx1 = np.where(flat_grid <= mesh.phys_boundary1 + mesh.phys_step / 10)[0][-1]
                ins_idx2 = np.where(flat_grid >= mesh.phys_boundary2 - mesh.phys_step / 10)[0][0]
                if ins_idx2 == flat_grid.size:
                    flat_grid = np.hstack((flat_grid[0:ins_idx1], mesh.phys_nodes()))
                    flat_sol = np.hstack((flat_sol[0:ins_idx1], mesh.solution))
                    flat_res = np.hstack((flat_res[0:ins_idx1], mesh.residual))
                else:
                    if ins_idx2 < flat_grid.size - 1:
                        ins_idx2 += 1
                    if len(flat_grid[ins_idx2:]) == 1 and flat_grid[ins_idx2] == mesh.phys_nodes()[-1]:
                        flat_grid = np.hstack((flat_grid[0:ins_idx1], mesh.phys_nodes()))
                        flat_sol = np.hstack((flat_sol[0:ins_idx1], mesh.solution))
                        flat_res = np.hstack((flat_res[0:ins_idx1], mesh.residual))
                    else:
                        flat_grid = np.hstack((flat_grid[0:ins_idx1], mesh.phys_nodes(), flat_grid[ins_idx2:]))
                        flat_sol = np.hstack((flat_sol[0:ins_idx1], mesh.solution, flat_sol[ins_idx2:]))
                        flat_res = np.hstack((flat_res[0:ins_idx1], mesh.residual, flat_res[ins_idx2:]))
                if debug:
                    print 'flat_grid is from', flat_grid[0] * 1e6, 'to', flat_grid[-1] * 1e6
        return flat_grid, flat_sol, flat_res


def list_merge_overlaps(mesh_list):
    tidy = False
    while not tidy:
        i_list = []
        overlap_found = False
        for i in range(len(mesh_list)):
            i_list.append(i)
            for j in range(len(mesh_list)):
                if j not in i_list:
                    # print i, j
                    if mesh_list[i].overlap_with(mesh_list[j]):
                        overlap_found = True
                        # print 'meshes overlap. MERGING.'
                        mesh_list[i].merge_with(mesh_list[j])
                        mesh_list.remove(mesh_list[j])
                        break
            if overlap_found:
                break
        tidy = True
    return mesh_list
