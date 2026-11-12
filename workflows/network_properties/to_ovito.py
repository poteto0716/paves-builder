"""Convert PAVES structure JSON (optionally with OpenMM state.xml coordinates) for OVITO.

The OpenMM state.xml written for each phase has only coordinates, velocities
and the box (nm); elements, bonds and charges are in the PAVES
*.structure.json (angstrom) written next to it. This script reads the JSON
and writes:

- a LAMMPS data file (atom_style full) with element-named atom types, partial
  charges, molecule IDs and bonds, wrapped into the cell with image flags.
  OVITO reads it directly (File > Load File) including the bonds;
- optionally one multi-frame extended XYZ trajectory of several structures
  (e.g. all thermal-history phases), whose cell may change per frame.
  It has no bonds; add them in OVITO with Create bonds.

Examples:
    # one structure -> relaxed_300K.data
    python workflows/network_properties/to_ovito.py replica_01/relaxed_300K.structure.json
    # topology from the JSON, coordinates and box from an OpenMM state
    python workflows/network_properties/to_ovito.py replica_01/02_high_equilibration/final.structure.json \\
        --state replica_01/02_high_equilibration/state.xml -o high.data
    # every phase of a replica: one .data per phase and a trajectory in phase order
    python workflows/network_properties/to_ovito.py replica_01 --trajectory replica_01/history.xyz
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

SUFFIX = '.structure.json'


def load(path, state=None):
    """Atoms, bonds, positions (angstrom) and cell rows (angstrom) of one structure."""
    data = json.loads(Path(path).read_text())
    units = data.get('units', {})
    if units and (units.get('positions'), units.get('box')) != ('angstrom', 'angstrom'):
        raise ValueError(f'{path}: expected angstrom units, found {units}')
    atoms, bonds, molecule, offset = [], [], [], 0
    for index, graph in enumerate(data['graphs']):
        atoms.extend(graph['atoms'])
        molecule.extend([index + 1] * len(graph['atoms']))
        bonds.extend((int(i) + offset, int(j) + offset, order) for i, j, order in graph['bonds'])
        offset += len(graph['atoms'])
    positions = np.asarray(data['positions'], float)
    cell = np.asarray(data['box_vectors'], float)
    if state is not None:
        positions, cell = read_state(state)
    if cell.shape == (3,):
        cell = np.diag(cell)
    if positions.shape != (len(atoms), 3):
        raise ValueError(f'{path}: {len(positions)} positions for {len(atoms)} atoms')
    return dict(atoms=atoms, bonds=bonds, molecule=molecule, positions=positions, cell=cell,
                time_ps=data.get('time_ps'))


def read_state(path):
    """Positions and box of an OpenMM serialized State, converted from nm to angstrom."""
    root = ET.parse(path).getroot()
    box = root.find('PeriodicBoxVectors')
    if box is None:
        raise ValueError(f'{path}: no PeriodicBoxVectors')
    cell = np.array([[float(box.find(k).get(a)) for a in 'xyz'] for k in 'ABC']) * 10
    positions = np.array([[float(p.get(a)) for a in 'xyz'] for p in root.find('Positions')]) * 10
    return positions, cell


def wrap(positions, cell):
    fractional = positions @ np.linalg.inv(cell)
    image = np.floor(fractional).astype(int)
    return (fractional - image) @ cell, image


def write_data(structure, path, title=''):
    cell = structure['cell']
    # LAMMPS cells are lower triangular: a along x, b in the xy plane.
    if abs(cell[0, 1]) > 1e-8 or abs(cell[0, 2]) > 1e-8 or abs(cell[1, 2]) > 1e-8:
        raise ValueError('cell is not lower triangular; rotate it before writing LAMMPS data')
    positions, image = wrap(structure['positions'], cell)
    atoms = structure['atoms']
    elements = sorted({a['element'] for a in atoms}, key=lambda e: min(a['atomic_number'] for a in atoms
                                                                         if a['element'] == e))
    types = {e: i for i, e in enumerate(elements, 1)}
    masses = {e: next(a['mass'] for a in atoms if a['element'] == e) for e in elements}
    lines = [f'LAMMPS data from PAVES structure {title}'.rstrip(), '',
             f'{len(atoms)} atoms', f'{len(elements)} atom types',
             f'{len(structure["bonds"])} bonds', '1 bond types', '',
             f'0.0 {cell[0, 0]:.8f} xlo xhi', f'0.0 {cell[1, 1]:.8f} ylo yhi', f'0.0 {cell[2, 2]:.8f} zlo zhi']
    if np.any(np.abs([cell[1, 0], cell[2, 0], cell[2, 1]]) > 1e-8):
        lines.append(f'{cell[1, 0]:.8f} {cell[2, 0]:.8f} {cell[2, 1]:.8f} xy xz yz')
    # OVITO takes the atom type names from these comments.
    lines += ['', 'Masses', ''] + [f'{types[e]} {masses[e]:.6f} # {e}' for e in elements]
    lines += ['', 'Atoms # full', '']
    for i, (atom, xyz, flags) in enumerate(zip(atoms, positions, image)):
        lines.append(f'{i + 1} {structure["molecule"][i]} {types[atom["element"]]} '
                     f'{atom.get("partial_charge", 0.):.6f} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f} '
                     f'{flags[0]} {flags[1]} {flags[2]}')
    if structure['bonds']:
        lines += ['', 'Bonds', '']
        lines += [f'{k} 1 {i + 1} {j + 1}' for k, (i, j, _) in enumerate(structure['bonds'], 1)]
    Path(path).write_text('\n'.join(lines) + '\n')


def xyz_frame(structure, comment=''):
    cell, atoms = structure['cell'], structure['atoms']
    positions, _ = wrap(structure['positions'], cell)
    lattice = ' '.join(f'{v:.8f}' for v in cell.reshape(-1))
    header = (f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:charge:R:1:molecule:I:1 pbc="T T T"'
              + (f' {comment}' if comment else ''))
    rows = [f'{a["element"]} {x:.6f} {y:.6f} {z:.6f} {a.get("partial_charge", 0.):.6f} {m}'
            for a, (x, y, z), m in zip(atoms, positions, structure['molecule'])]
    return f'{len(atoms)}\n{header}\n' + '\n'.join(rows) + '\n'


def inputs(paths):
    """Expand a replica directory into its phase structures in run order."""
    found = []
    for path in map(Path, paths):
        if path.is_dir():
            phases = sorted(path.glob(f'0*/final{SUFFIX}'))
            if not phases:
                raise ValueError(f'{path}: no <phase>/final{SUFFIX} files')
            found += phases
        else:
            found.append(path)
    return found


def output_name(path):
    name = path.name[:-len(SUFFIX)] if path.name.endswith(SUFFIX) else path.stem
    # Phase files are all called final.structure.json; name them after the phase.
    return path.parent / f'{path.parent.name if name == "final" else name}.data'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputs', nargs='+', help='*.structure.json files or replica directories')
    parser.add_argument('-o', '--output', type=Path, help='LAMMPS data path (single input only)')
    parser.add_argument('--state', type=Path, help='OpenMM state.xml supplying positions and box (single input)')
    parser.add_argument('--trajectory', type=Path, help='also write all inputs as one extended XYZ trajectory')
    parser.add_argument('--no-data', action='store_true', help='write only the trajectory')
    args = parser.parse_args(argv)
    paths = inputs(args.inputs)
    if (args.output or args.state) and len(paths) != 1:
        parser.error('--output and --state take exactly one structure')
    frames = []
    for path in paths:
        structure = load(path, args.state)
        if not args.no_data:
            target = args.output or output_name(path)
            write_data(structure, target, title=str(path))
            print(f'{path} -> {target} ({len(structure["atoms"])} atoms, {len(structure["bonds"])} bonds)')
        if args.trajectory:
            label = path.parent.name if path.name == f'final{SUFFIX}' else path.name[:-len(SUFFIX)]
            time = structure['time_ps']
            frames.append(xyz_frame(structure, f'phase={label}' + (f' time_ps={time}' if time is not None else '')))
    if args.trajectory:
        args.trajectory.write_text(''.join(frames))
        print(f'{len(frames)} frames -> {args.trajectory}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
