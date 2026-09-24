"""DCD -> LAMMPS dump for OVITO: `id x y z` only, cell 0..L as in final.data.

    python dcd_to_ovito.py projects/NAME/runs/09_pull     # writes .../09_pull/trajectory_ovito.dump

In OVITO: open final.data (LAMMPS data, atom_style full), then add a
"Load trajectory" modifier with trajectory_ovito.dump. Types, charges, bonds
and molecules come from final.data; only positions change per frame.

Why not mdtraj's save_lammpstrj: it writes type 1 for every atom, which the
Load trajectory modifier copies over the data file's types, and a cell origin
that moves from frame to frame. TIMESTEP here is the frame index.
"""
import sys
from pathlib import Path


def convert(run_dir, force=False):
    """Write run_dir/trajectory_ovito.dump unless it is newer than the DCD; returns its path."""
    import mdtraj as md
    d = Path(run_dir)
    dump = d / 'trajectory_ovito.dump'
    if not force and dump.exists() and dump.stat().st_mtime >= (d / 'trajectory.dcd').stat().st_mtime:
        return dump
    t = md.load(str(d / 'trajectory.dcd'), top=str(d / 'topology.pdb'))
    xyz = t.xyz * 10.0                      # nm -> A
    box = t.unitcell_lengths * 10.0
    with open(dump, 'w') as fh:
        for f in range(t.n_frames):
            Lx, Ly, Lz = box[f]
            fh.write(f'ITEM: TIMESTEP\n{f}\nITEM: NUMBER OF ATOMS\n{t.n_atoms}\n'
                     f'ITEM: BOX BOUNDS pp pp pp\n0 {Lx:.6f}\n0 {Ly:.6f}\n0 {Lz:.6f}\n'
                     'ITEM: ATOMS id x y z\n')
            fh.writelines(f'{i + 1} {x:.4f} {y:.4f} {z:.4f}\n' for i, (x, y, z) in enumerate(xyz[f]))
    print(f'{t.n_frames} frames, {t.n_atoms} atoms -> {dump}')
    return dump


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    convert(sys.argv[1], force=True)
