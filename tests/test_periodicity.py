import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('periodicity', ROOT / 'workflows/crosslinking/periodicity.py')
periodicity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(periodicity)

BOX = np.diag([10.] * 3)


def ring(points):
    """Bonds closing a ring through consecutive points."""
    return [(i, (i + 1) % len(points)) for i in range(len(points))]


def dimension(bonds, xyz):
    (_, d, periods), = periodicity.components(bonds, np.asarray(xyz, float), BOX)
    return d, periods


class PeriodicDimension(unittest.TestCase):
    def test_finite_ring_has_no_period(self):
        self.assertEqual(dimension(ring(range(4)), [[1, 1, 1], [2.5, 1, 1], [2.5, 2.5, 1], [1, 2.5, 1]]), (0, []))

    def test_string_wound_along_x_is_one_dimensional(self):
        d, periods = dimension(ring(range(4)), [[1, 1, 1], [3.5, 1, 1], [6, 1, 1], [8.5, 1, 1]])
        self.assertEqual(d, 1)
        self.assertIn(periods, ([[1, 0, 0]], [[-1, 0, 0]]))

    def test_diagonal_string_is_one_dimensional_although_it_crosses_x_y_and_z(self):
        # Checking each axis would call this percolating along x, y and z.
        d, periods = dimension(ring(range(4)), [[1, 1, 1], [3.5, 3.5, 3.5], [6, 6, 6], [8.5, 8.5, 8.5]])
        self.assertEqual(d, 1)
        self.assertEqual(np.abs(periods).tolist(), [[1, 1, 1]])

    def test_sheet_with_a_component_on_every_axis_is_two_dimensional(self):
        # Periods (1,0,0) and (0,1,1) through a shared atom 0.
        xyz = [[1, 1, 1], [3.5, 1, 1], [6, 1, 1], [8.5, 1, 1], [1, 3.5, 3.5], [1, 6, 6], [1, 8.5, 8.5]]
        bonds = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (4, 5), (5, 6), (6, 0)]
        self.assertEqual(dimension(bonds, xyz)[0], 2)

    def test_three_independent_windings_are_three_dimensional(self):
        xyz = [[1, 1, 1], [3.5, 1, 1], [6, 1, 1], [8.5, 1, 1], [1, 3.5, 1], [1, 6, 1], [1, 8.5, 1],
               [1, 1, 3.5], [1, 1, 6], [1, 1, 8.5]]
        bonds = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (4, 5), (5, 6), (6, 0), (0, 7), (7, 8), (8, 9), (9, 0)]
        self.assertEqual(dimension(bonds, xyz)[0], 3)

    def test_diagonal_periods_can_span_three_dimensions(self):
        # (1,1,0), (0,1,1), (1,0,1): no axis-aligned winding, still rank 3; works in a triclinic cell too.
        xyz, bonds = [[1., 1, 1]], []
        for direction in ((1, 1, 0), (0, 1, 1), (1, 0, 1)):
            start = len(xyz)
            xyz += [[1 + 2.5 * k * c for c in direction] for k in (1, 2, 3)]
            bonds += [(0, start), (start, start + 1), (start + 1, start + 2), (start + 2, 0)]
        self.assertEqual(dimension(bonds, xyz)[0], 3)
        shear = np.array([[1., 0, 0], [.3, 1, 0], [.2, .1, 1]])
        (_, d, _), = periodicity.components(bonds, np.asarray(xyz) @ shear, BOX @ shear)
        self.assertEqual(d, 3)

    def test_network_topology_reports_the_dimension_and_mass_by_dimension(self):
        from workflows.network_properties.analysis import topology
        atoms = [{'mass': 12.}] * 4
        string = dict(atoms=atoms, bonds=[[0, 1, 1], [1, 2, 1], [2, 3, 1], [3, 0, 1]])
        result = topology([string], np.array([[1, 1, 1], [3.5, 3.5, 3.5], [6, 6, 6], [8.5, 8.5, 8.5]]) / 10,
                          BOX / 10)
        self.assertEqual(result['components'][0]['periodic_dimension'], 1)
        self.assertTrue(result['percolates_periodically'])
        self.assertFalse(result['percolates_in_three_dimensions'])
        self.assertEqual(result['mass_fraction_by_periodic_dimension'], {'0': 0., '1': 1., '2': 0., '3': 0.})


if __name__ == '__main__':
    unittest.main()
