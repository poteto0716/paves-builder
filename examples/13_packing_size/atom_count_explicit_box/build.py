from pathlib import Path

import polyse


here = Path(__file__).resolve().parent
system = polyse.build(
    monomer="*CC*",
    dp=10,
    total_atoms=700,
    atom_count_mode="floor",
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    cell=(30.0, 30.0, 32.0),
    periodic=(True, True, True),
    temperature=413,
    terminator="*C",
    seed=1234,
)

print(system.n_atoms)
print(system.requested)
system.write_lammps(here / "output", prefix="atom_count_explicit_box")

