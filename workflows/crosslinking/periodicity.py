"""Periodic dimension of a bonded cluster: the rank of its lattice periods.

Each bond carries the integer cell offset it crosses (minimum image, in
fractional coordinates, so any triclinic cell works). Walking the bonds from
one atom assigns every atom a lift (cell offset). Meeting an atom again at a
different lift closes a cycle that winds through the periodic cell; the
difference of the two lifts is a lattice period of the cluster. A cycle that
closes in the same image (an ordinary ring) adds nothing.

The periodic dimension is the number of linearly independent periods: 0 for
a finite molecule, 1 for a string, 2 for a sheet, 3 for a network infinite in
all directions. This counts a string wound only along (1,1,1) as 1, where
checking each axis separately would call it percolating along x, y and z.
The same rank criterion defines the gel point of periodic systems in the
RWTH Aachen work behind puls-group/percolation-analyzer.
"""
import numpy as np


def add_period(basis, period):
    """Append period to basis (list of int 3-vectors) when it is linearly independent; returns basis."""
    period = np.asarray(period, dtype=int)
    if len(basis) < 3 and period.any() and \
            np.linalg.matrix_rank(np.array(basis + [period], dtype=float)) > len(basis):
        basis.append(period)
    return basis


def walk(nodes, edges, start):
    """Lifts of the nodes reached from start, and a basis of the periods met on the way.

    edges: node -> [(other, shift)] where shift is the integer cell offset of
    other seen from node across that bond. Only nodes in `nodes` are visited.
    """
    lift, basis, pending = {start: np.zeros(3, dtype=int)}, [], [start]
    while pending:
        i = pending.pop()
        for j, shift in edges[i]:
            if j not in nodes:
                continue
            expected = lift[i] + shift
            if j not in lift:
                lift[j] = expected
                pending.append(j)
            elif len(basis) < 3:
                add_period(basis, expected - lift[j])
    return lift, basis


def bond_edges(bonds, xyz, cell):
    """edges for walk() from (i, j) bonds, coordinates and cell rows (same length unit)."""
    inverse = np.linalg.inv(np.asarray(cell, dtype=float))
    xyz = np.asarray(xyz, dtype=float)
    edges = [[] for _ in xyz]
    for i, j in bonds:
        shift = -np.rint((xyz[j] - xyz[i]) @ inverse).astype(int)
        edges[i].append((j, shift))
        edges[j].append((i, -shift))
    return edges


def components(bonds, xyz, cell):
    """Connected components with their periodic dimension: [(atoms, dimension, periods)]."""
    edges = bond_edges(bonds, xyz, cell)
    everything, seen, found = set(range(len(edges))), set(), []
    for root in range(len(edges)):
        if root in seen:
            continue
        lift, basis = walk(everything, edges, root)
        seen.update(lift)
        found.append((sorted(lift), len(basis), [b.tolist() for b in basis]))
    return found
