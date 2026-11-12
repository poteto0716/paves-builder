"""PAVES structure JSON / OpenMM state -> LAMMPS data and extended XYZ for OVITO."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.to_ovito import load, main


def structure(positions, box=10.):
    atom = lambda z, e, m, q: dict(atomic_number=z, element=e, mass=m, partial_charge=q)
    return dict(units=dict(positions='angstrom', box='angstrom'),
                box_vectors=[[box, 0, 0], [0, box, 0], [0, 0, box]], positions=positions,
                graphs=[dict(atoms=[atom(8, 'O', 15.999, -.8), atom(1, 'H', 1.008, .4), atom(1, 'H', 1.008, .4)],
                             bonds=[[0, 1, 1], [0, 2, 1]])])


class ToOvitoTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_data_file_wraps_with_image_flags_and_keeps_bonds(self):
        # The second H lies outside the cell: wrapped back with image flag -1.
        path = self.root / 'water.structure.json'
        path.write_text(json.dumps(structure([[.2, 5, 5], [1.1, 5, 5], [-.5, 5, 5]])))
        self.assertEqual(main([str(path)]), 0)
        text = (self.root / 'water.data').read_text()
        self.assertIn('3 atoms', text)
        self.assertIn('2 bonds', text)
        self.assertIn('1 1.008000 # H', text)
        self.assertIn('2 15.999000 # O', text)
        atoms = text.split('Atoms # full\n\n')[1].split('\n\nBonds')[0].splitlines()
        self.assertEqual(atoms[2].split()[-3:], ['-1', '0', '0'])
        self.assertAlmostEqual(float(atoms[2].split()[4]), 9.5)
        self.assertIn('Bonds\n\n1 1 1 2\n2 1 1 3', text)

    def test_state_xml_supplies_coordinates_in_nm(self):
        path = self.root / 'final.structure.json'
        path.write_text(json.dumps(structure([[0, 0, 0]] * 3)))
        state = self.root / 'state.xml'
        state.write_text('<?xml version="1.0" ?><State><PeriodicBoxVectors>'
                         '<A x="2" y="0" z="0"/><B x="0" y="2" z="0"/><C x="0" y="0" z="2"/>'
                         '</PeriodicBoxVectors><Positions><Position x=".1" y=".2" z=".3"/>'
                         '<Position x=".2" y=".2" z=".3"/><Position x=".1" y=".3" z=".3"/></Positions></State>')
        loaded = load(path, state)
        np.testing.assert_allclose(loaded['cell'], np.eye(3) * 20)
        np.testing.assert_allclose(loaded['positions'][0], [1, 2, 3])

    def test_replica_directory_becomes_phase_files_and_trajectory(self):
        for phase, box in (('00_minimize', 10.), ('02_high_equilibration', 12.)):
            (self.root / 'replica_01' / phase).mkdir(parents=True)
            (self.root / 'replica_01' / phase / 'final.structure.json').write_text(
                json.dumps(structure([[1, 1, 1], [2, 1, 1], [1, 2, 1]], box)))
        trajectory = self.root / 'history.xyz'
        self.assertEqual(main([str(self.root / 'replica_01'), '--trajectory', str(trajectory)]), 0)
        self.assertTrue((self.root / 'replica_01/02_high_equilibration/02_high_equilibration.data').is_file())
        lines = trajectory.read_text().splitlines()
        self.assertEqual(len(lines), 10)
        self.assertIn('phase=02_high_equilibration', lines[6])
        self.assertIn('Lattice="12.00000000 0.00000000', lines[6])


if __name__ == '__main__':
    unittest.main()
