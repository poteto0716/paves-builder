from pathlib import Path

import polyse


here = Path(__file__).resolve().parent
components = [
    polyse.Solvent("Cc1ccccc1", name="toluene"),
    polyse.Polymer("*C(c1ccccc1)C*", dp=4, terminator="*C", name="polystyrene"),
    polyse.Polymer("*CC(C)(C(=O)OC)*", dp=4, terminator="*C", name="pmma"),
]
system = polyse.pack(
    components,
    total_atoms=800,
    weight_fractions={"toluene": 0.2, "polystyrene": 0.4, "pmma": 0.4},
    forcefield=here.parents[1] / "forcefields/pcff.ff",
    density=1.0,
    temperature=413,
    tolerance=0.12,
    seed=1234,
    sequence_seed=7,
)

print(system.n_atoms)
print(system.composition)
system.write_lammps(here / "output", prefix="blend_in_solvent")

