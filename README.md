# PAVES

**Pipeline for Automated Virtual Experiments and Screening**

*Pave the way from polymer structure to properties.*

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

The bundled wheels target **Ubuntu 24.04 x86_64** with **CPython 3.11 or 3.13**.
Use Python 3.11 for the OpenFF crosslinking and sorption workflows. Each wheel
requires its matching Python ABI and cannot be installed on another platform.

```bash
git clone <repository-url> polypaves
cd polypaves
python3 -m pip install dist/polypaves-0.4.0-cp313-cp313-linux_x86_64.whl
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

## Crosslinking and water sorption

The binary distribution includes reaction SMARTS, periodic network MD,
saved-charge network loading, automatic Sage solvent assignment, periodic solvent
packing, and fixed-loading NVT/NPT simulation with checkpoint resume.

- [Crosslinked epoxy construction](workflows/crosslinking/README.md): monomer
  SMILES, activation rules, reaction/MD cycles, charged MOL2 and structure JSON.
- [Water and other solvent simulations](workflows/sorption/README.md): retain the
  network's saved charges, specify `Solvent("SMILES")`, pack a molecule count,
  and relax/sample the wet network.
- [Crosslinked network properties](workflows/network_properties/README.md):
  independent epoxy construction per replica, 800→100 K relaxation and Tg/CTE,
  300 K thermodynamics, void fraction, elastic response and
  thermal transport, followed by arbitrary neutral-solvent insertion free energies
  and self diffusion, with sampling diagnostics, per-replica CSV summaries, and
  local GPU job submission with explicit GPU indices and independent seeds.

These workflows use the **0.4.0 CPython 3.11 wheel** and the
[OpenFF environment](workflows/crosslinking/environment.yml). The public workflow
scripts call compiled PAVES APIs; engine and API implementation source is absent.
Fixed solvent loading is an initial condition, not a verified equilibrium uptake.

## Docker Compose and OVITO Jupyter

The [Docker Compose setup](docs/docker.md) starts JupyterLab with PAVES, OpenFF,
OpenMM, MDTraj, and OVITO already installed. It opens the existing adhesion
visualization notebook and saves edits and simulation results to this checkout.
CPU execution is supported; an additional Compose file enables NVIDIA GPU access.

```bash
export PAVES_UID=$(id -u)
export PAVES_GID=$(id -g)
docker compose up -d --build
docker compose logs jupyter
```

Open the token-bearing localhost URL in the logs. The image installs the public
wheel and does not require the private PAVES implementation source.

## Screening tools

[`scripts/smiles_parameter_check.py`](docs/smiles_parameter_check.md) reads IDs
and SMILES from CSV or XLSX, checks parameter assignment with PCFF, Sage, and/or
GAFF2, and reports the charge sum of every dot-separated fragment. Failures are
recorded per fragment without stopping the remaining batch.

## Binary distribution and implementation privacy

The wheel contains the simulation engine and high-level API as stripped native
extensions. Their original implementation source is not included in this public
repository. Software that runs on a user's computer cannot be made completely
immune to technical analysis, but symbols, debug information, and private build
paths are removed. [LICENSE](LICENSE) prohibits decompilation, disassembly, and
implementation recovery. See the [binary distribution policy](docs/binary_distribution.md)
for details.
