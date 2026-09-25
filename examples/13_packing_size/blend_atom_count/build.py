from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
short = polypaves.Polymer("*CC*", dp=10, terminator="*C", name="short")
long = polypaves.Polymer("*CCCC*", dp=10, terminator="*C", name="long")
system = polypaves.pack(
    [short, long],
    total_atoms=1200,
    mole_fractions={"short": 0.5, "long": 0.5},
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    density=0.785,
    temperature=413,
    tolerance=0.10,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="blend_atom_count")

