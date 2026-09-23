from pathlib import Path

import polyse


here = Path(__file__).resolve().parent
system = polyse.build(
    monomer="*CC*",
    dp=10,
    chains=8,
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    cell=(28.0, 34.0, 30.0),
    periodic=(True, True, True),
    temperature=413,
    terminator="*C",
    seed=1234,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="chain_count_explicit_box")

