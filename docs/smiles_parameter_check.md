# Batch SMILES parameter check

`scripts/smiles_parameter_check.py` reads compound IDs and SMILES from a CSV or
XLSX table, splits disconnected SMILES at `.`, and checks every fragment with
PAVES. Each fragment is parameterized independently, so an electrically neutral
salt cannot hide incorrect neutralization of both ions.

Run the tool in an environment containing PAVES and the optional dependencies
for the selected force fields:

```bash
python scripts/smiles_parameter_check.py compounds.xlsx \
  --forcefields pcff,sage,gaff2 \
  --sage-version 2.3.0 \
  --ambertools-path /path/to/ambertools/environment \
  --output parameter_check.xlsx
```

CSV input and output are also supported:

```bash
python scripts/smiles_parameter_check.py compounds.csv -o parameter_check.csv
```

The first row must contain column names. By default, the tool recognizes `id`,
`name`, `compound_id`, `compound`, or `identifier` as the ID column, and
`smiles`, `canonical_smiles`, or `isomeric_smiles` as the SMILES column. Column
matching is case-insensitive. Select nonstandard columns explicitly:

```bash
python scripts/smiles_parameter_check.py compounds.xlsx \
  --id-column RegistryNumber --smiles-column Structure
```

For XLSX input, the first worksheet is read unless `--sheet NAME` is given.
XLSX support uses the Python standard library and does not require pandas or
openpyxl. Formula cells are read from their cached values; save the workbook in
a spreadsheet application before running the check if formulas were changed.

The output contains one row for each input row, dot-separated fragment, and
force field:

| Column | Meaning |
|---|---|
| `source_row` | Row number in the input table, including the header |
| `id` | Value from the selected ID column |
| `input_smiles` | Original disconnected SMILES |
| `fragment_index` | One-based position in the dot-separated SMILES |
| `fragment_smiles` | SMILES sent to PAVES |
| `forcefield` | `pcff`, `sage`, or `gaff2` |
| `parameter_assignment` | `○` on success or `×` on failure |
| `charge_sum_e` | Sum of assigned partial charges in elementary-charge units |
| `atom_count` | Number of atoms, including hydrogens added by PAVES |
| `error_type`, `error_message` | Failure details; blank after success |

Assignment failures are expected screening results and therefore do not stop
the batch or make the command fail. Invalid input/output tables and command-line
arguments return a nonzero exit status.
