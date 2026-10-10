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
    "ReactionRule", "ReactionEngine", "ReactiveMDRunner", "ReactionGraph",
    "make_solvent", "pack_network", "PackedNetwork", "configure_logging",
]


def __getattr__(name):
    if name in ('make_solvent', 'pack_network', 'PackedNetwork'):
        from . import sorption
        return getattr(sorption, name)
    if name in ('ReactionRule', 'ReactionEngine', 'ReactiveMDRunner', 'ReactionGraph'):
        from . import reaction
        return getattr(reaction, name)
    if name == 'configure_logging':
        from .progress import configure_logging
        return configure_logging
    if name == "GAFF2ForceField":
        from .gaff2 import GAFF2ForceField
        return GAFF2ForceField
    raise AttributeError(name)
