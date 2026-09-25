#!/usr/bin/env python3
"""Convert structures between OpenMM and LAMMPS for the same polypaves-built system.

    # OpenMM State -> LAMMPS data (types, charges, topology from a template data file)
    python openmm_lammps.py to-lammps --state runs/07_relax/state.xml \\
        --template runs/04_assemble/out/interface.data --out relaxed.data \\
        [--styles runs/04_assemble/out/interface.in.styles]

    # LAMMPS data -> OpenMM State (positions, velocities if present, box)
    python openmm_lammps.py to-openmm --data relaxed.data \\
        --system runs/04_assemble/out/interface.openmm_system.xml --out state.xml

Both directions keep the atom order of the polypaves build. Coordinates written to
LAMMPS are wrapped into the cell with image flags, so molecules are whole on
reading; coordinates read from LAMMPS are unwrapped with those flags. Units:
LAMMPS real (A, A/fs), OpenMM (nm, nm/ps).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.unit as u

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdtools as T  # noqa: E402


def read_data(path):
    """Positions (A, unwrapped with image flags), velocities (A/fs or None), box (lo, L)."""
    lines = Path(path).read_text().splitlines()
    n = 0
    lo, L = np.zeros(3), np.zeros(3)
    for line in lines[:60]:
        w = line.split('#')[0].split()
        if len(w) == 2 and w[1] == 'atoms':
            n = int(w[0])
        if len(w) == 4 and w[2] in ('xlo', 'ylo', 'zlo'):
            a = 'xyz'.index(w[2][0])
            lo[a], L[a] = float(w[0]), float(w[1]) - float(w[0])
    pos = np.zeros((n, 3))
    vel = None

    def section(name):
        k = next((i for i, l in enumerate(lines) if l.split('#')[0].strip() == name), None)
        if k is None:
            return []
        out = []
        for line in lines[k + 2:]:
            w = line.split('#')[0].split()
            if not w:
                break
            out.append(w)
        return out
    atoms = section('Atoms')
    if not atoms:
        sys.exit(f'{path}: no Atoms section')
    for w in atoms:                                    # atom_style full: id mol type q x y z [ix iy iz]
        i = int(w[0]) - 1
        p = np.array(list(map(float, w[4:7])))
        if len(w) >= 10:
            p += np.array(list(map(int, w[7:10]))) * L
        pos[i] = p
    v = section('Velocities')
    if v:
        vel = np.zeros((n, 3))
        for w in v:
            vel[int(w[0]) - 1] = list(map(float, w[1:4]))
    return pos, vel, lo, L


def to_lammps(a):
    st = T.load_xml(a.state)
    vel = None
    try:
        vel = T.velocities_A_per_fs(st)
    except Exception:
        pass
    T.write_lammps_data(a.template, a.out, T.positions_A(st), T.box_A(st), velocities=vel,
                        title=f'converted from {Path(a.state).name} by openmm_lammps.py')
    if a.styles:
        out_styles = Path(a.out).with_suffix('.in.styles')
        T.write_styles(a.styles, out_styles, Path(a.out).name)
        print(f'wrote {a.out} and {out_styles}')
    else:
        print(f'wrote {a.out}')


def to_openmm(a):
    pos, vel, lo, L = read_data(a.data)
    system = T.load_xml(a.system)
    if system.getNumParticles() != len(pos):
        sys.exit(f'{a.data} has {len(pos)} atoms, the System {system.getNumParticles()}')
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(L * 0.1)])
    ctx.setPositions((pos - lo) * 0.1)
    if vel is not None:
        ctx.setVelocities(vel * 0.1 * 1000)            # A/fs -> nm/ps
    Path(a.out).write_text(mm.XmlSerializer.serialize(ctx.getState(getPositions=True, getVelocities=vel is not None)))
    print(f'wrote {a.out} ({len(pos)} atoms, box {" x ".join(f"{x:.4f}" for x in L)} A'
          + (', with velocities)' if vel is not None else ')'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    l = sub.add_parser('to-lammps')
    l.add_argument('--state', required=True)
    l.add_argument('--template', required=True, help='a data file polypaves wrote for the same system')
    l.add_argument('--out', required=True)
    l.add_argument('--styles', help='the in.styles polypaves wrote; copied and pointed at --out')
    o = sub.add_parser('to-openmm')
    o.add_argument('--data', required=True)
    o.add_argument('--system', required=True, help='the openmm_system.xml of the same system')
    o.add_argument('--out', required=True)
    a = ap.parse_args()
    {'to-lammps': to_lammps, 'to-openmm': to_openmm}[a.cmd](a)


if __name__ == '__main__':
    main()
