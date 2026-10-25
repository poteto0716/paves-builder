#!/usr/bin/env python3
"""Batch-check PAVES parameter assignment and fragment charge sums.

The input is a CSV or XLSX table containing an ID column and a SMILES column.
Disconnected SMILES are split on ``.`` and every fragment is parameterized
independently with each requested force field.  One output row is written for
every input-row/fragment/force-field combination.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape
import zipfile


FORCEFIELDS = ("pcff", "sage", "gaff2")
OUTPUT_FIELDS = (
    "source_row",
    "id",
    "input_smiles",
    "fragment_index",
    "fragment_smiles",
    "forcefield",
    "parameter_assignment",
    "charge_sum_e",
    "atom_count",
    "error_type",
    "error_message",
)
ID_ALIASES = ("id", "name", "compound_id", "compound", "identifier")
SMILES_ALIASES = ("smiles", "canonical_smiles", "isomeric_smiles")

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
_CELL_REFERENCE = re.compile(r"([A-Z]+)([0-9]+)$")
_ILLEGAL_XML = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _column_number(reference: str) -> int:
    match = _CELL_REFERENCE.fullmatch(reference)
    if not match:
        raise ValueError(f"invalid XLSX cell reference: {reference!r}")
    number = 0
    for character in match.group(1):
        number = number * 26 + ord(character) - ord("A") + 1
    return number


def _column_name(number: int) -> str:
    result = ""
    while number:
        number, remainder = divmod(number - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def _xlsx_target(path: Path, sheet: str | None) -> str:
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        item.attrib["Id"]: item.attrib["Target"]
        for item in relationships.findall(f"{{{_PKG_REL_NS}}}Relationship")
    }
    sheets = workbook.find(f"{{{_MAIN_NS}}}sheets")
    entries = [] if sheets is None else list(sheets)
    if not entries:
        raise ValueError("XLSX workbook has no worksheets")
    selected = entries[0] if sheet is None else next(
        (entry for entry in entries if entry.attrib.get("name") == sheet), None
    )
    if selected is None:
        names = ", ".join(repr(entry.attrib.get("name")) for entry in entries)
        raise ValueError(f"worksheet {sheet!r} not found; available worksheets: {names}")
    relationship_id = selected.attrib[f"{{{_REL_NS}}}id"]
    target = targets[relationship_id].lstrip("/")
    if target.startswith("xl/"):
        return target
    return "xl/" + target.lstrip("/")


def _read_xlsx(path: Path, sheet: str | None = None) -> list[list[str]]:
    """Read the displayed values from one worksheet using only the stdlib."""
    target = _xlsx_target(path, sheet)
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{_MAIN_NS}}}si"):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{_MAIN_NS}}}t")))
        root = ET.fromstring(archive.read(target))

    rows: list[list[str]] = []
    sheet_data = root.find(f"{{{_MAIN_NS}}}sheetData")
    if sheet_data is None:
        return rows
    for row in sheet_data.findall(f"{{{_MAIN_NS}}}row"):
        values: dict[int, str] = {}
        for cell in row.findall(f"{{{_MAIN_NS}}}c"):
            column = _column_number(cell.attrib["r"])
            kind = cell.attrib.get("t")
            if kind == "inlineStr":
                value = "".join(node.text or "" for node in cell.iter(f"{{{_MAIN_NS}}}t"))
            else:
                node = cell.find(f"{{{_MAIN_NS}}}v")
                value = "" if node is None or node.text is None else node.text
                if kind == "s" and value:
                    value = shared[int(value)]
                elif kind == "b":
                    value = "TRUE" if value == "1" else "FALSE"
            values[column] = value
        if values:
            rows.append([values.get(i, "") for i in range(1, max(values) + 1)])
    return rows


def _xlsx_text(value: object) -> str:
    value = _ILLEGAL_XML.sub("", str(value))[:32767]
    preserve = ' xml:space="preserve"' if value[:1].isspace() or value[-1:].isspace() else ""
    return f"<is><t{preserve}>{escape(value)}</t></is>"


def _write_xlsx(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    all_rows = [list(fields), *[[row.get(field, "") for field in fields] for row in rows]]
    xml_rows = []
    for row_number, row in enumerate(all_rows, 1):
        cells = []
        for column_number, value in enumerate(row, 1):
            if value is None or value == "":
                continue
            reference = f"{_column_name(column_number)}{row_number}"
            if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
                cells.append(f'<c r="{reference}"><v>{value}</v></c>')
            else:
                cells.append(f'<c r="{reference}" t="inlineStr">{_xlsx_text(value)}</c>')
        xml_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    last_cell = f"{_column_name(len(fields))}{len(all_rows)}"
    worksheet = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<worksheet xmlns="{_MAIN_NS}"><dimension ref="A1:{last_cell}"/>'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="15"/><sheetData>'
        + "".join(xml_rows)
        + f'</sheetData><autoFilter ref="A1:{_column_name(len(fields))}{len(all_rows)}"/></worksheet>'
    )
    content_types = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>'''
    package_relationships = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="{_PKG_REL_NS}">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>'''
    workbook = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="{_MAIN_NS}" xmlns:r="{_REL_NS}"><sheets><sheet name="results" sheetId="1" r:id="rId1"/></sheets></workbook>'''
    workbook_relationships = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="{_PKG_REL_NS}">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>'''
    styles = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="{_MAIN_NS}"><fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>'''
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", package_relationships)
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_relationships)
        archive.writestr("xl/styles.xml", styles)
        archive.writestr("xl/worksheets/sheet1.xml", worksheet)


def read_table(path: Path, sheet: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            raw_rows = list(csv.reader(handle))
    elif suffix == ".xlsx":
        raw_rows = _read_xlsx(path, sheet)
    else:
        raise ValueError("input must have a .csv or .xlsx extension")
    raw_rows = [row for row in raw_rows if any(str(value).strip() for value in row)]
    if not raw_rows:
        raise ValueError("input table is empty")
    headers = [str(value).strip() for value in raw_rows[0]]
    if not all(headers):
        raise ValueError("input header contains a blank column name")
    if len(set(name.casefold() for name in headers)) != len(headers):
        raise ValueError("input header contains duplicate column names")
    records = []
    for row in raw_rows[1:]:
        padded = [*row, *([""] * (len(headers) - len(row)))]
        records.append(dict(zip(headers, [str(value) for value in padded[: len(headers)]])))
    return headers, records


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    elif suffix == ".xlsx":
        _write_xlsx(path, OUTPUT_FIELDS, rows)
    else:
        raise ValueError("output must have a .csv or .xlsx extension")


def _find_column(headers: list[str], requested: str | None, aliases: tuple[str, ...], label: str) -> str:
    by_folded = {header.casefold(): header for header in headers}
    if requested is not None:
        if requested.casefold() not in by_folded:
            raise ValueError(f"{label} column {requested!r} not found")
        return by_folded[requested.casefold()]
    matches = [by_folded[name] for name in aliases if name in by_folded]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"multiple possible {label} columns found: {', '.join(matches)}; select one explicitly")
    raise ValueError(f"no {label} column found; use --{label}-column")


def _charge_sum(parameterized: object, forcefield: str) -> float:
    if forcefield in ("pcff", "gaff2"):
        return math.fsum(float(value) for value in parameterized.charges)
    import openmm
    from openmm import unit

    nonbonded = next(
        force for force in parameterized.openmm_system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    )
    return math.fsum(
        nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
        for index in range(nonbonded.getNumParticles())
    )


def parameterize_fragment(
    smiles: str,
    forcefield: str,
    *,
    sage_version: str | None = None,
    ambertools_path: str | None = None,
) -> tuple[float, int]:
    import polypaves as paves

    options: dict[str, object] = {}
    if forcefield == "sage" and sage_version:
        options["forcefield_version"] = sage_version
    if forcefield == "gaff2" and ambertools_path:
        options["ambertools_path"] = ambertools_path
    config = paves.pack(
        [paves.Solvent(smiles, name="fragment")],
        counts=[1],
        cell=(40.0, 40.0, 40.0),
        forcefield=forcefield,
        trials=1,
        relax=0,
        postrelax=0,
        _return_config=True,
        **options,
    )
    parameterized = config.parameterize()
    structure = getattr(parameterized, "structure", parameterized)
    return _charge_sum(parameterized, forcefield), int(structure.n_atoms)


def check_records(
    records: list[dict[str, str]],
    *,
    id_column: str,
    smiles_column: str,
    forcefields: tuple[str, ...],
    sage_version: str | None = None,
    ambertools_path: str | None = None,
    progress: bool = True,
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    jobs = sum(max(1, len(record[smiles_column].strip().split("."))) for record in records) * len(forcefields)
    completed = 0
    for source_row, record in enumerate(records, 2):
        original = record[smiles_column].strip()
        fragments = original.split(".") if original else [""]
        for fragment_index, fragment in enumerate(fragments, 1):
            fragment = fragment.strip()
            for forcefield in forcefields:
                completed += 1
                result: dict[str, object] = {
                    "source_row": source_row,
                    "id": record[id_column],
                    "input_smiles": original,
                    "fragment_index": fragment_index,
                    "fragment_smiles": fragment,
                    "forcefield": forcefield,
                    "parameter_assignment": "×",
                    "charge_sum_e": None,
                    "atom_count": None,
                    "error_type": "",
                    "error_message": "",
                }
                if progress:
                    print(
                        f"[{completed}/{jobs}] row={source_row} fragment={fragment_index} "
                        f"forcefield={forcefield}",
                        file=sys.stderr,
                        flush=True,
                    )
                try:
                    if not fragment:
                        raise ValueError("empty SMILES fragment")
                    charge, atom_count = parameterize_fragment(
                        fragment,
                        forcefield,
                        sage_version=sage_version,
                        ambertools_path=ambertools_path,
                    )
                    if not math.isfinite(charge):
                        raise ValueError("PAVES returned a non-finite charge sum")
                    result.update(
                        parameter_assignment="○",
                        charge_sum_e=charge,
                        atom_count=atom_count,
                    )
                except Exception as error:  # Every failed assignment belongs in the output table.
                    result["error_type"] = type(error).__name__
                    result["error_message"] = " ".join(str(error).split())
                results.append(result)
    return results


def _forcefields(value: str) -> tuple[str, ...]:
    fields = tuple(dict.fromkeys(item.strip().lower() for item in value.split(",") if item.strip()))
    unknown = sorted(set(fields) - set(FORCEFIELDS))
    if not fields or unknown:
        message = "choose one or more of: " + ", ".join(FORCEFIELDS)
        if unknown:
            message += "; unknown: " + ", ".join(unknown)
        raise argparse.ArgumentTypeError(message)
    return fields


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Check PAVES parameter assignment and total charge for each dot-separated SMILES fragment."
    )
    result.add_argument("input", type=Path, help="input .csv or .xlsx table")
    result.add_argument("-o", "--output", type=Path, help="output .csv or .xlsx (default: <input>_paves_check.csv)")
    result.add_argument("--id-column", help="ID column name (auto-detected by default)")
    result.add_argument("--smiles-column", help="SMILES column name (auto-detected by default)")
    result.add_argument("--sheet", help="XLSX input worksheet name (default: first sheet)")
    result.add_argument("--forcefields", type=_forcefields, default=FORCEFIELDS,
                        help="comma-separated force fields (default: pcff,sage,gaff2)")
    result.add_argument("--sage-version", help="Sage version, for example 2.3.0")
    result.add_argument("--ambertools-path", help="AmberTools installation or environment path for GAFF2")
    result.add_argument("--quiet", action="store_true", help="do not print per-job progress")
    return result


def main(argv: list[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    try:
        if not arguments.input.is_file():
            raise FileNotFoundError(arguments.input)
        output = arguments.output or arguments.input.with_name(arguments.input.stem + "_paves_check.csv")
        if output.resolve() == arguments.input.resolve():
            raise ValueError("output must not overwrite the input table")
        headers, records = read_table(arguments.input, arguments.sheet)
        id_column = _find_column(headers, arguments.id_column, ID_ALIASES, "id")
        smiles_column = _find_column(headers, arguments.smiles_column, SMILES_ALIASES, "smiles")
        results = check_records(
            records,
            id_column=id_column,
            smiles_column=smiles_column,
            forcefields=arguments.forcefields,
            sage_version=arguments.sage_version,
            ambertools_path=arguments.ambertools_path,
            progress=not arguments.quiet,
        )
        write_table(output, results)
    except (FileNotFoundError, OSError, ValueError, KeyError, IndexError,
            ET.ParseError, zipfile.BadZipFile) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    successes = sum(row["parameter_assignment"] == "○" for row in results)
    failures = len(results) - successes
    print(f"Wrote {len(results)} rows to {output} ({successes} succeeded, {failures} failed).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
