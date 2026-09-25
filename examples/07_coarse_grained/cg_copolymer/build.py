from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
copolymer = polypaves.Copolymer(
    monomers={"A": "*[Aa]*", "B": "*[Bb]*"},
    dp=30,
    fractions={"A": 0.5, "B": 0.5},
    architecture="random_exact",
    terminator="*[Aa]",
    name="cg_random",
)
system = polypaves.pack(
    [copolymer],
    counts={"cg_random": 4},
    forcefield=here.parents[1] / "forcefields/cg_copolymer.ff",
    density=0.85,
    temperature=1.0,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="cg_random")

