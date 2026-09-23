from pathlib import Path

import polyse


here = Path(__file__).resolve().parent
solvent = polyse.Solvent("CCCC", name="solvent")
polymer = polyse.Polymer("*CC*", dp=10, terminator="*C", name="polymer")
system = polyse.pack(
    [solvent, polymer],
    total_atoms=1000,
    mole_fractions={"solvent": 0.9, "polymer": 0.1},
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    cell=(34.0, 34.0, 34.0),
    periodic=(True, True, True),
    temperature=413,
    tolerance=0.10,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="solvent_and_polymer_atom_count")

