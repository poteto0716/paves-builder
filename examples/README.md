# Python API examples

These examples use the public Python API directly through `import polypaves`.
They do not invoke the native CLI or read `.paves` input files. Each script
locates its resources relative to its own path, so it can be run from any
working directory inside the repository.

| Category | Script | Description |
|---|---|---|
| PCFF | [PMMA DP10](01_pcff_polystyrene/pmma_dp10/build.py) | One PMMA chain with ester groups |
| PCFF | [PS DP10](01_pcff_polystyrene/ps_dp10/build.py) | Two PS chains with aromatic side groups |
| Copolymer | [Explicit sequence](03_copolymers/explicit_sequence/build.py) | `Copolymer(sequence=...)` |
| Mixture | [Four-component blend](05_mixtures/polymer_blend/build.py) | Solvent, two polymers, and a copolymer |
| Mixture | [PCFF solvent blend](05_mixtures/polymer_blend_in_solvent/build.py) | Weight fractions and a target atom count |
| Coarse-grained | [CG copolymer](07_coarse_grained/cg_copolymer/build.py) | Exact random sequence and junction typing |
| End groups | [Asymmetric end groups](10_sequences_and_ends/asymmetric_end_groups/build.py) | Different head and tail groups |
| PCFF-IFF | [Kapton DP3](12_pcff_iff/kapton_dp3/build.py) | Kapton with ring closures |
| System size | [Atom count and box](13_packing_size/atom_count_explicit_box/build.py) | Target atom count for one component |
| System size | [Blend atom count](13_packing_size/blend_atom_count/build.py) | Mole fractions and a target atom count |
| System size | [Chain count and box](13_packing_size/chain_count_explicit_box/build.py) | Explicit cell dimensions and chain count |
| System size | [Solvent and polymer](13_packing_size/solvent_and_polymer_atom_count/build.py) | Mole fractions, target atom count, and explicit cell dimensions |

Run an example with:

```bash
python examples/03_copolymers/explicit_sequence/build.py
```

The OPLS and coarse-grained descriptors in `examples/forcefields/` are part of
the examples. The PCFF and IFF descriptors reference the parameter and typing
template files under `external/`. Generated files are written to each example's
`output/` directory and are ignored by Git.
