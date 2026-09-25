from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
system = polypaves.build(
    monomer="*CC(C)(C(=O)OC)*",
    forcefield=here.parents[1] / "forcefields/pcff.ff",
    chains=1,
    dp=10,
    terminator="*C",
    density=1.13,
    temperature=413,
    seed=12345,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="pmma10")

