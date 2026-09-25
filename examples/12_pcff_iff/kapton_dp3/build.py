from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
kapton = polypaves.Polymer(
    "*N1C(=O)c2cc3c(cc2C1=O)C(=O)N(C3=O)c4ccc(Oc5ccc(*)cc5)cc4",
    dp=3,
    terminator=("*c1ccccc1", "*c1ccccc1"),
    name="kapton",
)
system = polypaves.pack(
    [kapton],
    counts={"kapton": 4},
    forcefield=here.parents[1] / "forcefields/pcff_iff_long_bulk.ff",
    density=1.42,
    temperature=300,
    seed=1234,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="kapton")

