from pathlib import Path

import polyse


here = Path(__file__).resolve().parent
system = polyse.build(
    monomer="*C(c1ccccc1)C*",
    forcefield=here.parents[1] / "forcefields/pcff.ff",
    chains=2,
    dp=10,
    terminator="*C",
    density=0.969,
    temperature=413,
    seed=1234,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="ps10")

