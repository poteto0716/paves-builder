"""Monomer builds and crosslinking runs through the public crosslinking CLI.

A monomer project packs the monomer (no start sites, no reaction) and
relaxes it at its temperature T. A crosslink project starts every replica
from the relaxed T structure of the same-numbered monomer replica, creates
the requested start sites there and crosslinks at the same T. The ring
chemistry comes from the activation/reaction YAMLs (epoxide by default,
oxetane preset), copied into each project's definitions/ directory.
"""
from __future__ import annotations

from collections import deque
import hashlib
import inspect
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CROSSLINKING = ROOT / 'workflows/crosslinking'
DEFINITIONS = ('activation.yaml', 'reactions.yaml')
# Ring-opening presets: YAML directory and default monomer.
CHEMISTRIES = {
    'epoxy': dict(definitions=CROSSLINKING, smiles=None),  # crosslinking.py default, TMP triglycidyl ether
    # 3-ethyl-3-{[(3-ethyloxetan-3-yl)methoxy]methyl}oxetane (di-oxetane, "OXT-221").
    'oxetane': dict(definitions=CROSSLINKING / 'oxetane', smiles='CCC1(COCC2(CC)COC2)COC1'),
}
CHILD_ERROR = re.compile(r'^[\w.]+(?:Error|Exception|Interrupt)\b')
# crosslinking.run arguments the pipeline sets itself, never stored in project.json.
RESERVED = {'project', 'force', 'platform', 'seed', 'template_cache', 'start_structure'}
# Packing and pre-reaction options: used by the monomer build, meaningless when starting from a structure.
PACKING = ('smiles', 'total_atoms', 'atom_count_mode', 'density', 'compression_pressure_mpa',
           'compression_steps', 'pack_relax', 'pack_postrelax', 'reference_temperature',
           'reference_npt_steps', 'reference_frame_steps', 'pre_reaction_sampling_steps',
           'pre_reaction_frame_steps')
# A monomer build creates no start site and runs no cycle; the relaxed T structure is saved directly.
MONOMER_FIXED = dict(activation_count=0, activation_fraction=None, cycles=0, target_conversion=None,
                     reference_npt_steps=0, pre_reaction_sampling_steps=0, final_relaxation_steps=0,
                     final_md_steps=0, template_workers=1, save_trajectory=False, remove_sol_fraction=False,
                     max_reactions_per_cycle=None)
# NPT after the start sites are added to the relaxed monomer (1 fs steps): 100 ps.
CROSSLINK_EQUILIBRATION_STEPS = 100_000
MONOMER_RESULTS, NETWORK_RESULTS = 'monomer', 'network'


def defaults():
    from workflows.crosslinking.crosslinking import run
    return {name: parameter.default for name, parameter in inspect.signature(run).parameters.items()
            if name not in RESERVED}


def _given(args, **names):
    """CLI values that were given (not None), keyed by crosslinking.run argument names."""
    return {key: getattr(args, flag, None) for key, flag in names.items() if getattr(args, flag, None) is not None}


def choose_chemistry(chemistry, smiles):
    """Chemistry preset for a monomer: the given one, else the ring the SMILES carries (epoxy by default).

    A preset whose ring the SMILES lacks while it has the other one is refused: its
    activation rule would find no site and the monomer would be built without rings.
    """
    if not smiles:
        return chemistry or 'epoxy'
    found = [name for name, preset in CHEMISTRIES.items()
             if reactive_sites(smiles, (preset['definitions'] / 'activation.yaml').read_bytes())]
    if chemistry is None:
        if len(found) > 1:
            raise ValueError(f'the monomer has both {" and ".join(found)} rings; choose --chemistry')
        return found[0] if found else 'epoxy'
    if found and chemistry not in found:
        raise ValueError(f'the monomer has no {chemistry} ring but has {" / ".join(found)} rings; '
                         f'use --chemistry {found[0]} (or omit --chemistry)')
    return chemistry


def monomer_options(args, temperature):
    from workflows.crosslinking.crosslinking import DEFAULT_SMILES
    chemistry = getattr(args, 'chemistry', None) or 'epoxy'
    options = {k: v for k, v in defaults().items() if k in PACKING or k in MONOMER_FIXED
               or k in ('temperature', 'pressure_atm', 'equilibration_steps', 'equilibration_timestep',
                        'minimization_iterations', 'nvt_timestep')}
    options['smiles'] = CHEMISTRIES[chemistry]['smiles'] or DEFAULT_SMILES
    options.update(_given(args, smiles='smiles', total_atoms='total_atoms', atom_count_mode='atom_count_mode',
                          density='packing_density', compression_pressure_mpa='compression_pressure_mpa',
                          compression_steps='compression_steps', equilibration_steps='equilibration_steps'))
    options.update(MONOMER_FIXED, temperature=float(temperature))
    if getattr(args, 'smoke', False):
        options.update(compression_steps=2, equilibration_steps=2, minimization_iterations=100)
        if getattr(args, 'total_atoms', None) is None:
            options['total_atoms'] = 500
    return options


def crosslink_options(args, *, temperature, count=None, fraction=None, smoke=False):
    options = {k: v for k, v in defaults().items() if k not in PACKING}
    options.update(equilibration_steps=CROSSLINK_EQUILIBRATION_STEPS)
    options.update(_given(args, cycles='cycles', target_conversion='target_conversion',
        max_reactions_per_cycle='max_reactions_per_cycle', reaction_cutoff='reaction_cutoff',
        max_reaction_cutoff='max_reaction_cutoff', minimization_iterations='reaction_minimization_iterations',
        post_reaction_md_steps='post_reaction_md_steps', relaxation_stages='relaxation_stages',
        stop_after_idle_cycles='stop_after_idle_cycles', template_workers='template_workers',
        reaction_minimization_radius='reaction_minimization_radius',
        final_relaxation_steps='final_relaxation_steps', final_md_steps='final_md_steps',
        equilibration_steps='equilibration_steps'))
    for flag in ('save_trajectory', 'remove_sol_fraction', 'fixed_reaction_cutoff'):
        if getattr(args, flag, False):
            options[flag] = True
    if getattr(args, 'reaction_cutoff', None) is None:
        options['reaction_cutoff'] = min(options['reaction_cutoff'], options['max_reaction_cutoff'])
    options.update(activation_count=count, activation_fraction=fraction, temperature=float(temperature))
    if count == 0 or fraction == 0:
        # No active site, no reaction: skip the cycles; a conversion target does not apply.
        options.update(cycles=0, target_conversion=None)
    if smoke:
        options.update(equilibration_steps=2, cycles=min(options['cycles'], 1), nvt_steps=2, npt_steps=2,
                       template_workers=1, final_relaxation_steps=2, final_md_steps=2,
                       max_reactions_per_cycle=1, reaction_cutoff=10., max_reaction_cutoff=10.,
                       minimization_iterations=100)
    return options


def definitions_for(chemistry, activation_yaml=None, reactions_yaml=None):
    preset = CHEMISTRIES[chemistry]
    paths = dict(zip(DEFINITIONS, (activation_yaml, reactions_yaml)))
    return {name: (Path(path) if path else preset['definitions'] / name).read_bytes() for name, path in paths.items()}


def with_selection(definitions, *, count=None, fraction=None):
    """The copied activation.yaml states the start count/fraction actually used."""
    from workflows.crosslinking.crosslinking import activation_yaml_with_selection
    result = dict(definitions)
    result['activation.yaml'] = activation_yaml_with_selection(
        definitions['activation.yaml'].decode(), count=count, fraction=fraction).encode()
    return result


def write_definitions(project, definitions):
    directory = Path(project) / 'definitions'
    directory.mkdir(parents=True)
    for name, content in definitions.items():
        (directory / name).write_bytes(content)


def definition_digest(project):
    hashes = {name: hashlib.sha256((Path(project) / 'definitions' / name).read_bytes()).hexdigest()
              for name in DEFINITIONS}
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def reactive_sites(smiles, activation_yaml):
    """Activation target sites per monomer molecule (distinct atom sets)."""
    import yaml
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError(f'smiles must describe one connected molecule: {smiles}')
    rules = yaml.safe_load(activation_yaml)['activations']
    target = Chem.MolFromSmarts(rules[0]['target_smarts'])
    if target is None:
        raise ValueError(f"activation.yaml target_smarts is invalid: {rules[0]['target_smarts']}")
    return len({frozenset(m) for m in Chem.AddHs(molecule).GetSubstructMatches(target)})


def target_smarts(project):
    from polypaves.reaction import load_activation_definitions
    return load_activation_definitions(Path(project) / 'definitions/activation.yaml')[0].target_smarts


def _validate_common(options, label):
    known = defaults()
    if set(options) - set(known):
        raise ValueError(f'unknown {label} settings: {sorted(set(options) - set(known))}')
    for name in ('minimization_iterations', 'equilibration_steps', 'nvt_steps', 'npt_steps',
                 'final_relaxation_steps', 'final_md_steps', 'post_reaction_md_steps',
                 'stop_after_idle_cycles', 'template_workers', 'cycles', 'compression_steps'):
        if name in options and (type(options[name]) is not int or options[name] < (1 if name == 'template_workers' else 0)):
            raise ValueError(f'{label}.{name} must be a nonnegative integer')
    temperature = options.get('temperature')
    if isinstance(temperature, bool) or not isinstance(temperature, (int, float)) or not temperature > 0:
        raise ValueError(f'{label}.temperature must be positive')


def validate_monomer(cfg):
    options = cfg['build']
    _validate_common(options, 'build')
    if type(options['total_atoms']) is not int or options['total_atoms'] < 1:
        raise ValueError('build.total_atoms must be a positive integer')
    if options['atom_count_mode'] not in ('above', 'nearest', 'floor', 'exact'):
        raise ValueError('invalid build.atom_count_mode')
    if any(options.get(k) != v for k, v in MONOMER_FIXED.items()):
        raise ValueError('a monomer build has no start sites and no reaction cycles')
    if options['temperature'] != cfg['temperature_K']:
        raise ValueError('build.temperature must equal temperature_K')


def validate_crosslink(cfg):
    from polypaves.reaction import AdaptiveReactionCutoff
    from workflows.crosslinking.crosslinking import cycle_stages
    options = cfg['crosslinking']
    _validate_common(options, 'crosslinking')
    if set(options) & set(PACKING):
        raise ValueError(f'crosslinking starts from the monomer structure; remove {sorted(set(options) & set(PACKING))}')
    merged = {**defaults(), **options}
    AdaptiveReactionCutoff(initial_cutoff=merged['reaction_cutoff'], maximum_cutoff=merged['max_reaction_cutoff'],
        growth_factor=merged['cutoff_growth_factor'], zero_reaction_cycles=merged['cutoff_zero_cycles'])
    if merged['relaxation_stages'] is not None:
        cycle_stages(merged['relaxation_stages'], temperature=merged['temperature'],
                     pressure_atm=merged['pressure_atm'])
    target = merged['target_conversion']
    if target is not None and (isinstance(target, bool) or not 0 <= target <= 1):
        raise ValueError('crosslinking.target_conversion must lie between 0 and 1')
    for flag in ('save_trajectory', 'remove_sol_fraction', 'fixed_reaction_cutoff'):
        if type(merged[flag]) is not bool:
            raise ValueError(f'crosslinking.{flag} must be a boolean')
    count, fraction = merged['activation_count'], merged['activation_fraction']
    if (count is None) == (fraction is None):
        raise ValueError('a crosslink project needs exactly one activation count or fraction')
    if count is not None and (type(count) is not int or count < 0):
        raise ValueError('crosslinking.activation_count must be an integer >= 0')
    if fraction is not None and (isinstance(fraction, bool) or not 0 <= fraction < 1):
        raise ValueError('crosslinking.activation_fraction must satisfy 0 <= fraction < 1')
    if merged['temperature'] != cfg['temperature_K']:
        raise ValueError('crosslinking.temperature must equal the monomer temperature_K')


def parent_project(project, cfg):
    return (Path(project) / cfg['parent']).resolve()


def network_path(project, cfg, replica):
    base = Path(project) / f'replica_{replica:02d}'
    mode = cfg.get('input_mode', 'saved_network')
    if mode == 'monomer':
        return base / f'build/results/{MONOMER_RESULTS}/final.structure.json'
    if mode == 'crosslink':
        return base / f'crosslinking/results/{NETWORK_RESULTS}/final.structure.json'
    return Path(project) / cfg['network']


def build_command(work, name, options, *, seed, platform, start_structure=None):
    command = [sys.executable, '-u', str(CROSSLINKING / 'crosslinking.py'), '--project', name, '--seed', str(seed),
               '--platform', platform, '--template-cache', str(Path(work) / 'reaction_templates')]
    if start_structure is not None:
        command += ['--start-structure', str(start_structure)]
    for key, value in options.items():
        if value is None or value is False:
            continue
        flag = '--' + key.replace('_', '-')
        if isinstance(value, (list, dict)):
            value = json.dumps(value)
        command += [flag] if value is True else [flag, str(value)]
    # Only this pipeline-owned, incomplete construction is restarted. Completed
    # builds are retained and verified before reaching here.
    if (Path(work) / 'results' / name).exists():
        command.append('--force')
    return command


def run_crosslinking(project, work, name, options, *, seed, platform, lock=None, start_structure=None, label):
    """Run crosslinking.py in work (cwd) and return its results directory."""
    from workflows.network_properties.network_properties import write_json
    from workflows.network_properties.execution import die_with_parent
    work.mkdir(parents=True, exist_ok=True)
    for definition in DEFINITIONS:
        shutil.copy2(Path(project) / 'definitions' / definition, work / definition)
    command = build_command(work, name, options, seed=seed, platform=platform, start_structure=start_structure)
    write_json(work / 'invocation.json', dict(command=command, cwd=str(work), seed=seed))
    # Inherit the worker's single selected GPU; never share a working directory.
    # The child also holds the replica lock, so an orphaned build still blocks
    # a restart from deleting its output with --force.
    # stderr is relayed line by line to the job log and its tail kept, so a failure
    # names the crosslinking error itself rather than only the exit status.
    tail = deque(maxlen=200)
    with subprocess.Popen(command, cwd=work, stderr=subprocess.PIPE, text=True, errors='replace',
                          preexec_fn=die_with_parent,
                          pass_fds=() if lock is None else (lock.fileno(),)) as child:
        for line in child.stderr:
            sys.stderr.write(line)
            sys.stderr.flush()
            tail.append(line)
    if child.returncode:
        errors = [line.strip() for line in tail if CHILD_ERROR.match(line.strip())]
        cause = errors[-1] if errors else f'exit status {child.returncode}'
        failure = RuntimeError(f'{label}: {cause}')
        failure.add_note('--- crosslinking.py stderr (last lines) ---\n' + ''.join(tail).rstrip())
        raise failure
    output = work / 'results' / name
    settings = json.loads((output / 'settings.json').read_text())
    if settings['seed'] != seed:
        raise RuntimeError('crosslinking seed differs from the assigned replica seed')
    data = json.loads((output / 'final.structure.json').read_text())
    if not data.get('graphs') or not data.get('positions'):
        raise RuntimeError('crosslinking did not produce a saved structure')
    return output


def ensure_monomer_build(project, cfg, replica, lock=None):
    """Pack and relax the monomer at T (no start site); return the saved structure."""
    from workflows.network_properties.network_properties import digest, write_json
    from workflows.network_properties.execution import replica_seed
    project = Path(project)
    seed = replica_seed(project, cfg, replica)
    source = network_path(project, cfg, replica)
    work = source.parents[2]
    marker = work / 'done.json'
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved['seed'] != seed or digest(source) != saved['structure_sha256']:
            raise ValueError(f'Completed monomer build changed: {source}')
        return source
    print(f'Replica {replica}: pack and relax the monomer at {cfg["temperature_K"]:g} K, seed={seed}', flush=True)
    output = run_crosslinking(project, work, MONOMER_RESULTS, cfg['build'], seed=seed, platform=cfg['platform'],
                              lock=lock, label=f'Monomer build failed for replica {replica}')
    settings = json.loads((output / 'settings.json').read_text())
    record = dict(seed=seed, structure_sha256=digest(source), structure=str(source.relative_to(project)),
                  atoms=settings['actual_total_atoms'], molecules=settings['monomers'],
                  reactive_sites=settings['initial_reactive_sites'], smiles=settings['smiles'],
                  temperature_K=settings['temperature_K'])
    write_json(marker, record)
    return source


def ensure_crosslink(project, cfg, replica, lock=None):
    """Crosslink from the same-numbered monomer replica's relaxed T structure; return the network."""
    from workflows.network_properties.network_properties import digest, write_json
    from workflows.network_properties.execution import configuration_fingerprint, replica_seed
    from workflows.network_properties import monomer
    project = Path(project)
    parent = parent_project(project, cfg)
    parent_cfg = json.loads((parent / 'project.json').read_text())
    if configuration_fingerprint(parent_cfg) != cfg['parent_configuration_sha256']:
        raise ValueError(f'monomer project {parent} changed after this crosslink project was created')
    start = monomer.relaxed_path(parent, parent_cfg, replica)
    seed = replica_seed(project, cfg, replica)
    source = network_path(project, cfg, replica)
    work = source.parents[2]
    base = project / f'replica_{replica:02d}'
    marker = work / 'done.json'
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved['seed'] != seed or digest(source) != saved['network_sha256']:
            raise ValueError(f'Completed network changed: {source}')
        if not start.is_file() or digest(start) != saved['monomer']['structure_sha256']:
            raise ValueError(f'monomer replica_{replica:02d} changed after this network was built from it')
        return source
    if not monomer.replica_complete(parent, parent_cfg, replica):
        raise RuntimeError(f'monomer {parent.name}/replica_{replica:02d} has not finished; '
                           'run or submit the monomer project first')
    options = cfg['crosslinking']
    starts = (f"{options['activation_count']} starts" if options['activation_count'] is not None
              else f"start fraction {options['activation_fraction']:g}")
    print(f'Replica {replica}: crosslink from {start.relative_to(parent.parent)} ({starts}), seed={seed}', flush=True)
    output = run_crosslinking(project, work, NETWORK_RESULTS, options, seed=seed, platform=cfg['platform'],
                              lock=lock, start_structure=start,
                              label=f'Crosslinking failed for replica {replica}')
    settings = json.loads((output / 'settings.json').read_text())
    metrics = json.loads((output / 'metrics.json').read_text())
    if settings['start_structure_sha256'] != digest(start):
        raise RuntimeError('crosslinking started from another structure than the monomer replica')
    data = json.loads(source.read_text())
    target = options['target_conversion']
    record = dict(seed=seed, network_sha256=digest(source),
        network=str(source.relative_to(project)), atoms=len(data['positions']),
        monomer=dict(project=str(parent), replica=f'replica_{replica:02d}',
                     structure=str(start.relative_to(parent)), structure_sha256=digest(start)),
        conversion=metrics['conversion'], target_conversion=target,
        target_reached=target is None or metrics['conversion'] >= target,
        activation_count=settings.get('activation_count'),
        requested_activation_fraction=options['activation_fraction'],
        actual_activation_fraction=settings.get('actual_activation_fraction'),
        cycles_completed=metrics['cycles_completed'],
        epoxide_consumption_including_activation=metrics['epoxide_consumption_including_activation'],
        **{key: metrics[key] for key in ('kinetic_chains', 'volume_before_crosslinking_nm3',
                                         'activation_density_nm3', 'activation_density_mmol_cm3',
                                         'sol_mass_fraction', 'sol_monomers_removed', 'sol_oligomers_removed',
                                         'sol_atoms_removed', 'gel_initial_epoxide_sites',
                                         'gel_reacted_epoxide_sites', 'gel_epoxide_conversion') if key in metrics})
    write_json(base / 'network_origin.json', record)
    write_json(marker, record)
    return source
