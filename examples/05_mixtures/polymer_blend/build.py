from pathlib import Path

import polypaves


here = Path(__file__).resolve().parent
components = [
    polypaves.Solvent("CC", name="ethane"),
    polypaves.Polymer("*CC*", dp=20, terminator="*C", name="pe20"),
    polypaves.Polymer("*CCCC*", dp=10, terminator="*C", name="pb10"),
    polypaves.Copolymer(
        monomers={"A": "*CC*", "B": "*CCCC*"},
        dp=14,
        fractions={"A": 0.5, "B": 0.5},
        architecture="random_exact",
        terminator="*C",
        name="copolymer",
    ),
]
system = polypaves.pack(
    components,
    total_atoms=1600,
    weight_fractions={"ethane": 0.1, "pe20": 0.3, "pb10": 0.3, "copolymer": 0.3},
    forcefield=here.parents[1] / "forcefields/opls_alkane.ff",
    density=0.75,
    temperature=413,
    tolerance=0.12,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="polymer_blend")

