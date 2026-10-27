"""Per-replica epoxy packing/crosslinking using the public crosslinking CLI."""
from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections import deque

ROOT = Path(__file__).resolve().parents[2]
CROSSLINKING = ROOT / 'workflows/crosslinking'
DEFINITIONS = ('activation.yaml', 'reactions.yaml')
CHILD_ERROR = re.compile(r'^[\w.]+(?:Error|Exception|Interrupt)\b')


def defaults():
    from workflows.crosslinking.crosslinking import run
    # Keep the standalone CLI and the full pipeline's defaults in one place.
    reserved = {'project', 'force', 'platform', 'seed', 'template_cache'}
    return {name: parameter.default for name, parameter in inspect.signature(run).parameters.items()
            if name not in reserved}


def configure(args):
    options = defaults()
    overrides = dict(smiles=getattr(args, 'epoxy_smiles', None),
        total_atoms=getattr(args, 'total_atoms', None), atom_count_mode=getattr(args, 'atom_count_mode', None),
        cycles=getattr(args, 'cycles', None), target_conversion=getattr(args, 'target_conversion', None),
        activation_fraction=getattr(args, 'activation_fraction', None),
        activation_count=getattr(args, 'activation_count', None),
        max_reactions_per_cycle=getattr(args, 'max_reactions_per_cycle', None),
        save_trajectory=getattr(args, 'save_trajectory', False),
        compression_pressure_mpa=getattr(args, 'compression_pressure_mpa', None),
        compression_steps=getattr(args, 'compression_steps', None),
        reaction_cutoff=getattr(args, 'reaction_cutoff', None),
        max_reaction_cutoff=getattr(args, 'max_reaction_cutoff', None),
        fixed_reaction_cutoff=getattr(args, 'fixed_reaction_cutoff', None) or None,
        minimization_iterations=getattr(args, 'reaction_minimization_iterations', None),
        post_reaction_md_steps=getattr(args, 'post_reaction_md_steps', None),
        relaxation_stages=getattr(args, 'relaxation_stages', None),
        stop_after_idle_cycles=getattr(args, 'stop_after_idle_cycles', None),
        template_workers=getattr(args, 'template_workers', None),
        reaction_minimization_radius=getattr(args, 'reaction_minimization_radius', None),
        density=getattr(args, 'packing_density', None),
        final_relaxation_steps=getattr(args, 'final_relaxation_steps', None),
        final_md_steps=getattr(args, 'final_md_steps', None))
    options.update({k:v for k,v in overrides.items() if v is not None})
    if overrides['reaction_cutoff'] is None:
        options['reaction_cutoff'] = min(options['reaction_cutoff'], options['max_reaction_cutoff'])
    if args.smoke:
        options.update(compression_steps=2, equilibration_steps=2, cycles=1, nvt_steps=2, npt_steps=2, template_workers=1,
            final_relaxation_steps=2, final_md_steps=2,
            max_reactions_per_cycle=1, reaction_cutoff=10., max_reaction_cutoff=10.,
            minimization_iterations=100)
        if overrides['total_atoms'] is None:
            options['total_atoms'] = 500
        if overrides['activation_count'] is None and overrides['activation_fraction'] is None:
            options['activation_count'] = 1
    paths = dict(zip(DEFINITIONS, (getattr(args, 'activation_yaml', None),
                                  getattr(args, 'reactions_yaml', None))))
    definitions = {name: (Path(path) if path else CROSSLINKING / name).read_bytes()
                   for name, path in paths.items()}
    from workflows.crosslinking.crosslinking import activation_yaml_with_selection
    # The project's copied activation.yaml states the start count/fraction actually used.
    definitions['activation.yaml'] = activation_yaml_with_selection(definitions['activation.yaml'].decode(),
        count=options['activation_count'], fraction=options['activation_fraction']).encode()
    hashes = {name: hashlib.sha256(content).hexdigest() for name,content in definitions.items()}
    fingerprint = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return dict(input_mode='epoxy', crosslinking=options, definition_sha256=hashes,
                input_sha256=fingerprint, schema_version=2), definitions


def definition_digest(project):
    hashes = {name: hashlib.sha256((Path(project) / 'crosslinking' / name).read_bytes()).hexdigest()
              for name in DEFINITIONS}
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def validate(cfg):
    if cfg.get('input_mode') != 'epoxy':
        return
    if cfg['platform'] not in ('CPU', 'CUDA'):
        raise ValueError('Epoxy construction supports CPU or CUDA')
    options = cfg['crosslinking']
    known = defaults()
    # Older projects may lack newer options; crosslinking.py then uses its defaults.
    if set(options) - set(known):
        raise ValueError(f'unknown crosslinking settings: {sorted(set(options) - set(known))}')
    merged = {**known, **options}
    from workflows.crosslinking.crosslinking import cycle_stages
    from polypaves.reaction import AdaptiveReactionCutoff
    AdaptiveReactionCutoff(initial_cutoff=merged['reaction_cutoff'], maximum_cutoff=merged['max_reaction_cutoff'],
        growth_factor=merged['cutoff_growth_factor'], zero_reaction_cycles=merged['cutoff_zero_cycles'])
    if merged['relaxation_stages'] is not None:
        cycle_stages(merged['relaxation_stages'], temperature=merged['temperature'],
                     pressure_atm=merged['pressure_atm'])
    for name, minimum in (('post_reaction_md_steps', 0), ('minimization_iterations', 0),
                          ('final_relaxation_steps', 0), ('final_md_steps', 0),
                          ('stop_after_idle_cycles', 0), ('template_workers', 1)):
        if type(merged[name]) is not int or merged[name] < minimum:
            raise ValueError(f'crosslinking.{name} must be an integer >= {minimum}')
    options = merged
    if type(options['total_atoms']) is not int or options['total_atoms'] < 1:
        raise ValueError('crosslinking.total_atoms must be a positive integer')
    if type(options['cycles']) is not int or options['cycles'] < 0:
        raise ValueError('crosslinking.cycles must be an integer >= 0')
    if options['atom_count_mode'] not in ('above', 'nearest', 'floor', 'exact'):
        raise ValueError('invalid crosslinking.atom_count_mode')
    target = options['target_conversion']
    if target is not None and (isinstance(target, bool) or not 0 <= target <= 1):
        raise ValueError('crosslinking.target_conversion must lie between 0 and 1')
    if type(options['save_trajectory']) is not bool:
        raise ValueError('crosslinking.save_trajectory must be a boolean')
    count, fraction = options['activation_count'], options['activation_fraction']
    if count is not None and fraction is not None:
        raise ValueError('use one activation count or fraction')


def network_path(project, cfg, replica):
    if cfg.get('input_mode') == 'epoxy':
        return Path(project) / f'replica_{replica:02d}/crosslinking/results/epoxy/final.structure.json'
    return Path(project) / cfg['network']


def build_command(project, cfg, replica):
    work = Path(project) / f'replica_{replica:02d}/crosslinking'
    command = [sys.executable, '-u', str(CROSSLINKING / 'crosslinking.py'),
               '--project', 'epoxy', '--seed', str(cfg['seed'] + replica),
               '--platform', cfg['platform'], '--template-cache', str(work / 'reaction_templates')]
    for key, value in cfg['crosslinking'].items():
        if value is None or value is False:
            continue
        flag = '--' + key.replace('_', '-')
        if isinstance(value, (list, dict)):
            value = json.dumps(value)
        command += [flag] if value is True else [flag, str(value)]
    # Only this pipeline-owned, incomplete construction is restarted. Completed
    # builds are retained and verified by ensure_network before reaching here.
    if (work / 'results/epoxy').exists():
        command.append('--force')
    return command


def ensure_network(project, cfg, replica, lock=None):
    from workflows.network_properties.network_properties import digest, write_json
    from workflows.network_properties.execution import die_with_parent
    source = network_path(project, cfg, replica)
    if cfg.get('input_mode') != 'epoxy':
        return source
    base = Path(project) / f'replica_{replica:02d}'
    work = base / 'crosslinking'
    work.mkdir(parents=True, exist_ok=True)
    marker = work / 'done.json'
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved['seed'] != cfg['seed'] + replica or digest(source) != saved['network_sha256']:
            raise ValueError(f'Completed network changed: {source}')
        return source
    for name in DEFINITIONS:
        shutil.copy2(Path(project) / 'crosslinking' / name, work / name)
    command = build_command(project, cfg, replica)
    write_json(work / 'invocation.json', dict(command=command, cwd=str(work), seed=cfg['seed'] + replica))
    print(f'Replica {replica}: build and crosslink epoxy, seed={cfg["seed"] + replica}', flush=True)
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
        failure = RuntimeError(f'Epoxy construction failed for replica {replica}: {cause}')
        failure.add_note('--- crosslinking.py stderr (last lines) ---\n' + ''.join(tail).rstrip())
        raise failure
    settings = json.loads((source.parent / 'settings.json').read_text())
    metrics = json.loads((source.parent / 'metrics.json').read_text())
    if settings['seed'] != cfg['seed'] + replica:
        raise RuntimeError('crosslinking seed differs from the assigned replica seed')
    data = json.loads(source.read_text())
    if not data.get('graphs') or not data.get('positions'):
        raise RuntimeError('crosslinking did not produce a saved network')
    target = cfg['crosslinking']['target_conversion']
    record = dict(seed=cfg['seed'] + replica, network_sha256=digest(source),
        network=str(source.relative_to(project)), atoms=settings['actual_total_atoms'],
        conversion=metrics['conversion'], target_conversion=target,
        target_reached=target is None or metrics['conversion'] >= target,
        cycles_completed=metrics['cycles_completed'],
        epoxide_consumption_including_activation=metrics['epoxide_consumption_including_activation'])
    write_json(base / 'network_origin.json', record)
    write_json(marker, record)
    return source
