"""Public Python API for polypaves."""

from importlib.metadata import version as _distribution_version

from .api import (Polymer, Copolymer, Solvent, Slab, BoxRegion, System,
                  BuildError, build, pack, register_forcefield)
from .config import (PavesConfig, ForceFieldConfig, BackendConfig, SimulationConfig,
                     InputParser, from_file, DEFAULT_SAGE_VERSION, ChargeConfig, FragmentConfig)
from .openff import SageForceField, SageCompatibilityError, OpenFFDependencyError
from .openmm import OpenMMBackend
from .structure import load_structure, write_topology, PCFFForceField

__version__ = _distribution_version("polypaves")
__all__ = [
    "Polymer", "Copolymer", "Solvent", "Slab", "BoxRegion", "System", "BuildError",
    "build", "pack", "register_forcefield",
    "PavesConfig", "ForceFieldConfig", "BackendConfig", "SimulationConfig", "InputParser",
    "from_file", "DEFAULT_SAGE_VERSION", "ChargeConfig", "FragmentConfig",
    "SageForceField", "SageCompatibilityError", "OpenFFDependencyError", "OpenMMBackend",
    "load_structure", "write_topology", "PCFFForceField", "GAFF2ForceField",
]


def __getattr__(name):
    if name == "GAFF2ForceField":
        from .gaff2 import GAFF2ForceField
        return GAFF2ForceField
    raise AttributeError(name)
