from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "smiles_parameter_check", ROOT / "scripts" / "smiles_parameter_check.py"
)
TOOL = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(TOOL)


class TableTests(unittest.TestCase):
    def test_csv_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            path.write_text("ID,SMILES,comment\nA,CC.O,test\n", encoding="utf-8")
            headers, rows = TOOL.read_table(path)
            self.assertEqual(headers, ["ID", "SMILES", "comment"])
            self.assertEqual(rows, [{"ID": "A", "SMILES": "CC.O", "comment": "test"}])

    def test_xlsx_round_trip(self):
        rows = [{field: "" for field in TOOL.OUTPUT_FIELDS}]
        rows[0].update(source_row=2, id="salt-α", input_smiles="[Na+].[Cl-]",
                       fragment_index=1, fragment_smiles="[Na+]", forcefield="sage",
                       parameter_assignment="○", charge_sum_e=1.0, atom_count=1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "results.xlsx"
            TOOL.write_table(path, rows)
            headers, actual = TOOL.read_table(path, "results")
            self.assertEqual(headers, list(TOOL.OUTPUT_FIELDS))
            self.assertEqual(actual[0]["id"], "salt-α")
            self.assertEqual(actual[0]["charge_sum_e"], "1.0")
            self.assertEqual(actual[0]["parameter_assignment"], "○")

    def test_unknown_extension_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "csv or .xlsx"):
            TOOL.read_table(Path("input.tsv"))


class CheckTests(unittest.TestCase):
    def test_dot_fragments_and_failures_are_reported_independently(self):
        def fake(fragment, forcefield, **_):
            if fragment == "bad" or forcefield == "pcff":
                raise RuntimeError("missing parameter\nwith detail")
            return ({"[Na+]": 1.0, "[Cl-]": -1.0}[fragment], 1)

        with patch.object(TOOL, "parameterize_fragment", side_effect=fake):
            results = TOOL.check_records(
                [{"ID": "salt", "SMILES": "[Na+].[Cl-]"}, {"ID": "broken", "SMILES": "bad"}],
                id_column="ID", smiles_column="SMILES", forcefields=("pcff", "sage"), progress=False,
            )
        self.assertEqual(len(results), 6)
        sage = [row for row in results if row["forcefield"] == "sage" and row["id"] == "salt"]
        self.assertEqual([row["charge_sum_e"] for row in sage], [1.0, -1.0])
        self.assertTrue(all(row["parameter_assignment"] == "○" for row in sage))
        failed = next(row for row in results if row["id"] == "broken" and row["forcefield"] == "sage")
        self.assertEqual(failed["parameter_assignment"], "×")
        self.assertEqual(failed["error_type"], "RuntimeError")
        self.assertEqual(failed["error_message"], "missing parameter with detail")

    def test_empty_fragment_is_a_failure(self):
        results = TOOL.check_records(
            [{"id": "empty", "smiles": "CC..O"}], id_column="id", smiles_column="smiles",
            forcefields=("sage",), progress=False,
        )
        self.assertEqual(len(results), 3)
        self.assertEqual(results[1]["parameter_assignment"], "×")
        self.assertEqual(results[1]["error_message"], "empty SMILES fragment")


if __name__ == "__main__":
    unittest.main()
