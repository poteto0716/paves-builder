from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
polymer = polypaves.Polymer(
    "*CC*",
    dp=20,
    terminator=("*C", "*c1ccccc1"),
    name="asymmetric_end_groups",
)
system = polypaves.pack(
    [polymer],
    counts={"asymmetric_end_groups": 6},
    forcefield=here.parents[1] / "forcefields/pcff.ff",
    density=0.85,
    temperature=413,
    seed=1234,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="asymmetric_end_groups")

