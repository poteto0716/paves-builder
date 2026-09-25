from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
copolymer = polypaves.Copolymer(
    monomers={"A": "*CC*", "B": "*CCCC*"},
    sequence=["A", "A", "B", "A", "B", "B", "A", "A"],
    terminator="*C",
    name="explicit",
)
system = polypaves.pack(
    [copolymer],
    counts={"explicit": 2},
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    density=0.785,
    temperature=413,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="explicit")

