# PAVES

**Polymer Automation for Virtual Evaluation and Simulation**

*Pave the way from polymer chemistry to simulation.*

PAVES is a Python package that constructs molecular-dynamics systems from
monomer SMILES and composition specifications. It calls a C++ engine directly
from Python and returns a `System` object that exposes coordinates, atom types,
charges, bonds, composition, and other simulation data. The project name is
**PAVES**; the Python package and command-line tools are named **polypaves**
because the `paves` name is already used by another package on PyPI. Releases
through 0.2.x used the name polyse.

This public repository distributes a Linux binary wheel, Python API examples,
documentation, and force-field data. It does not include the engine's C/C++
source code, headers, object files, or build tree. See [LICENSE](LICENSE) for
the terms of use.

## Supported platform and installation

The bundled wheel is validated on **Ubuntu 24.04 x86_64 with CPython 3.13**.
It cannot be installed with a different Python ABI or on another platform.

```bash
git clone <repository-url> polypaves
cd polypaves
python3 -m pip install dist/polypaves-0.3.0-cp313-cp313-linux_x86_64.whl
python3 -c 'import polypaves; print(polypaves.__version__)'
sha256sum -c SHA256SUMS
```

At runtime, Ubuntu 24.04 packages `libc6`, `libstdc++6`, and `libgcc-s1` are
required. A compiler, CMake, and development source code are not required.

## Build your first molecular system

```python
import polypaves

system = polypaves.build(
    monomer="*CC(C)(C(=O)OC)*",  # PMMA
    forcefield="pcff",
    chains=20,
    dp=100,
    density=1.18,
    temperature=413,
)

print(system.n_atoms)
system.write_lammps("pmma_system")
```

`dp` is the number of repeat units per chain, and `chains` is the number of
chains. `write_lammps()` writes `system.data`, `system.in.styles`, and
`system.identity` to the requested directory.

## Compose and pack multiple components

`Polymer`, `Copolymer`, `Solvent`, and `Slab` objects can be combined in one
system.

```python
import polypaves

pmma = polypaves.Polymer("*CC(C)(C(=O)OC)*", dp=20, name="pmma")
toluene = polypaves.Solvent("Cc1ccccc1", name="toluene")

# Specify the number of molecules of each component.
system = polypaves.pack(
    [pmma, toluene],
    counts={"pmma": 4, "toluene": 100},
    forcefield="pcff",
    density=1.0,
)

# Or specify a target atom count and weight fractions.
system = polypaves.pack(
    [pmma, toluene],
    total_atoms=10_000,
    weight_fractions={"pmma": 0.7, "toluene": 0.3},
    forcefield="pcff",
    density=1.0,
)

print(system.composition)  # Realized counts, mole fractions, and weight fractions.
```

A `System` exposes `positions`, `atom_types`, `masses`, `charges`, `bonds`,
`angles`, `dihedrals`, `impropers`, `identity`, `box`, `composition`, and
`report`. See the [Python API documentation](docs/python_api.md) for details.

## Examples

The [examples](examples/README.md) directory contains executable Python API
examples. Every input is a `.py` file and runs without a custom input format or
a CLI subprocess.

```bash
python examples/01_pcff_polystyrene/pmma_dp10/build.py
```

Generated files are written below each example's `output/` directory and are
ignored by Git. See [external data](docs/external_data.md) for force-field
provenance and [binary validation](docs/validation.md) for the validation scope.

## Binary distribution and implementation privacy

The wheel contains the simulation engine and high-level API as stripped native
extensions. Their original implementation source is not included in this public
repository. Software that runs on a user's computer cannot be made completely
immune to technical analysis, but symbols, debug information, and private build
paths are removed. [LICENSE](LICENSE) prohibits decompilation, disassembly, and
implementation recovery. See the [binary distribution policy](docs/binary_distribution.md)
for details.
