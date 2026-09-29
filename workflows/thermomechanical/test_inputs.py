import tempfile
import unittest
from pathlib import Path

from thermomechanical import chemistry_from_file, normalize_block_counts, parse_assignments


class ChemistryInputTest(unittest.TestCase):
    def test_bulk_pcff_enables_tail_correction_but_adhesion_fields_do_not(self):
        root = Path(__file__).resolve().parents[2]
        fields = root / 'examples' / 'forcefields'
        bulk = (fields / 'pcff.ff').read_text()
        self.assertRegex(bulk, r'(?m)^pair_modify\s+mix\s+sixthpower\s+tail\s+yes\s*$')
        for name in ('pcff_iff_long_bulk.ff', 'pcff_iff_long_bulk_mod.ff'):
            self.assertNotRegex((fields / name).read_text(), r'(?m)\btail\s+yes\b')

    def test_fixed_blocks_are_normalized_for_automatic_dp(self):
        lines = normalize_block_counts(['block A 70', 'block B 30'])
        self.assertEqual(lines, ['block         A fraction 0.7',
                                 'block         B fraction 0.3'])

    def test_file_keeps_chemistry_and_replaces_size_records(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'copolymer.paves'
            path.write_text("""polypaves-build 1
name old
sequence_seed 19
monomer A '*CC*'
monomer B '*CCCC*'
architecture random-exact
block A 3
block B 1
degree 40
chains 8
density 0.9
""")
            lines, seed = chemistry_from_file(path)
        self.assertEqual(seed, 19)
        self.assertIn("monomer A '*CC*'", lines)
        self.assertIn('block         A fraction 0.75', lines)
        self.assertNotIn('degree 40', lines)
        self.assertNotIn('chains 8', lines)

    def test_name_value_parser(self):
        self.assertEqual(parse_assignments(['A=0.7', 'B=0.3'], '--fraction'),
                         {'A': '0.7', 'B': '0.3'})
        with self.assertRaises(ValueError):
            parse_assignments(['missing_separator'], '--fraction')


if __name__ == '__main__':
    unittest.main()
