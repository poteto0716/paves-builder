"""Compare project directories without pooling different materials or conditions."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
from urllib.parse import quote

from workflows.network_properties.execution import file_lock
from workflows.network_properties.replica_summary import collect as collect_replicas

MARKER = '.network_properties_collection.json'
META = ['project', 'project_directory', 'input_mode', 'epoxy_smiles', 'base_seed',
        'replicates_requested', 'replicates_completed', 'project_state', 'smoke',
        'target_total_atoms', 'target_conversion', 'cooling_rate_K_ns', 'pressure_atm',
        'input_sha256', 'configuration_sha256']
CONDITION = ['solvent_name', 'solvent_smiles', 'solvent_molecules']
STATS = ['n_expected', 'n_available', 'n_checks_passed', 'mean', 'std', 'sem',
         'checked_mean', 'checked_std', 'checked_sem']


def _json(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def _csv(path, fields, rows):
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', newline='', encoding='utf-8-sig') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def enable_collection(root):
    root = Path(root).resolve()
    with file_lock(root / '.projects_report.lock'):
        if not (root / MARKER).exists():
            _json(root / MARKER, dict(schema_version=1, kind='network_properties'))


def _discover(root):
    projects, errors = [], []
    for directory, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith('.') and not d.startswith('replica_')
                         and d not in ('jobs', 'frames', '__pycache__'))
        if 'project.json' not in files:
            continue
        path = Path(directory)
        try:
            cfg = json.loads((path / 'project.json').read_text())
        except (OSError, ValueError) as exc:
            errors.append(dict(project=str(path.relative_to(root)), error=str(exc)))
            dirs.clear()
            continue
        if (isinstance(cfg, dict) and 'replicates' in cfg and 'solvents' in cfg
                and ('network' in cfg or cfg.get('input_mode') == 'epoxy')):
            projects.append(path)
            # Never walk the potentially large trajectories inside a project.
            dirs.clear()
    return sorted(projects), errors


def _finite(value):
    return None if isinstance(value, float) and not math.isfinite(value) else value


def _project(root, project):
    from workflows.network_properties.gpu_jobs import replica_complete
    with file_lock(project / '.report.lock'):
        cfg = json.loads((project / 'project.json').read_text())
        if type(cfg['replicates']) is not int or cfg['replicates'] < 1:
            raise ValueError('replicates must be a positive integer')
        ids, properties = collect_replicas(project, cfg)
        completed = sum(replica_complete(project, cfg, i) for i in range(1, cfg['replicates'] + 1))
        available = any(p['n_available'] for p in properties.values())
        build = cfg.get('crosslinking', {})
        metadata = dict(project=project.relative_to(root).as_posix(), project_directory=str(project),
            input_mode=cfg.get('input_mode', 'saved_network'), epoxy_smiles=build.get('smiles'),
            base_seed=cfg.get('seed'), replicates_requested=cfg['replicates'],
            replicates_completed=completed,
            project_state='completed' if completed == cfg['replicates'] else 'partial' if available else 'pending',
            smoke=bool(cfg.get('smoke')), target_total_atoms=build.get('total_atoms'),
            target_conversion=build.get('target_conversion'), cooling_rate_K_ns=cfg.get('cooling_rate_K_ns'),
            pressure_atm=cfg.get('pressure_atm'), input_sha256=cfg.get('input_sha256'),
            configuration_sha256=hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest())
        metadata = {k:_finite(v) for k,v in metadata.items()}
        return dict(metadata=metadata, replica_ids=ids, properties=properties,
                    solvents=cfg.get('solvents', []), solvent_molecules=cfg.get('solvent_molecules'))


def _condition(record, key):
    name = key.split('.')[1] if key.startswith('solvents.') else ''
    solvent = next((s for s in record['solvents'] if s['name'] == name), {})
    return dict(solvent_name=name, solvent_smiles=solvent.get('smiles'),
                solvent_molecules=record['solvent_molecules'] if name else None)


def _markdown(root, records, errors):
    lines = ['# Project property comparison', '',
        '[Project statistics CSV](projects_summary.csv) · [All replica values CSV](replicas_summary.csv)', '',
        'Each project retains its own replica statistics. Different projects are never pooled.',
        'Cells show mean ± replica SD and [available/requested replicas]. An asterisk means some '
        'available estimates have not passed checks, including smoke runs. Check counts, SEM and '
        'checked-only statistics are in the CSV files.', '',
        '| Project | Input mode | Completed replicas | Smoke |', '|---|---|---:|---|']

    def escape(text):
        return str(text).replace('|', '\\|').replace('\n', ' ')

    def project_link(record):
        name = record['metadata']['project']
        report = Path(record['metadata']['project_directory']) / 'report.md'
        target = quote((Path(name) / 'report.md').as_posix(), safe='/')
        return f'[{escape(name)}]({target})' if report.is_file() else escape(name)

    for record in records:
        m = record['metadata']
        lines.append(f'| {project_link(record)} | {m["input_mode"]} | '
                     f'{m["replicates_completed"]}/{m["replicates_requested"]} | {m["smoke"]} |')
    groups = {}
    for record in records:
        for key, p in record['properties'].items():
            if p['display']:
                groups.setdefault(p['group'], {})[key,p['unit']] = p['label']
    for group, columns in groups.items():
        lines += ['', f'**{escape(group)}**', '',
            '| Project | ' + ' | '.join(f'{escape(label)} / {unit}' for (_,unit),label in columns.items()) + ' |',
            '|---|' + '---:|' * len(columns)]
        for record in records:
            cells = []
            for key,unit in columns:
                p = record['properties'].get(key)
                if not p or p['unit'] != unit or p['mean'] is None:
                    cells.append('—')
                    continue
                value = f'{p["mean"]:.6g}'
                if p['std'] is not None:
                    value += f' ± {p["std"]:.3g}'
                if p['n_checks_passed'] < p['n_available']:
                    value += '*'
                cells.append(f'{value} [{p["n_available"]}/{p["n_expected"]}]')
            lines.append('| ' + project_link(record) + ' | ' + ' | '.join(cells) + ' |')
    if errors:
        lines += ['', '**Directories requiring attention**', '']
        lines += [f'- {escape(e["project"])}: {escape(e["error"])}' for e in errors]
    return '\n'.join(lines) + '\n'


def collect_projects(root):
    root = Path(root).resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    with file_lock(root / '.projects_report.lock'):
        projects, errors = _discover(root)
        if not projects and not errors:
            raise ValueError(f'No network-property projects found under {root}')
        records = []
        for project in projects:
            try:
                records.append(_project(root, project))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                errors.append(dict(project=project.relative_to(root).as_posix(), error=str(exc)))
        aggregate, individuals = [], []
        for record in records:
            metadata = record['metadata']
            for key,p in record['properties'].items():
                common = dict(metadata, **_condition(record,key), property=key, label=p['label'], unit=p['unit'])
                aggregate.append(dict(common, **{k:p[k] for k in STATS}))
                for i,replica in enumerate(record['replica_ids'],1):
                    value = p['replicas'][replica]
                    seed = metadata['base_seed'] + i if isinstance(metadata['base_seed'],int) else None
                    individuals.append(dict(common, replica=replica, seed=seed, value=value['value'],
                        status=value['status'], checks_passed=value['checks_passed'],
                        source=value['source']))
        _csv(root / 'projects_summary.csv', META + CONDITION + ['property','label','unit'] + STATS, aggregate)
        _csv(root / 'replicas_summary.csv', META + CONDITION + ['property','label','unit',
            'replica','seed','value','status','checks_passed','source'], individuals)
        _json(root / 'projects_summary.json', dict(schema_version=1, root=str(root),
            projects=records, errors=errors))
        target = root / 'projects_report.md'
        temporary = target.with_suffix('.md.tmp')
        temporary.write_text(_markdown(root,records,errors))
        temporary.replace(target)
        if not (root / MARKER).exists():
            _json(root / MARKER, dict(schema_version=1, kind='network_properties'))
    print(f'Project comparison: {root / "projects_report.md"} ({len(records)} projects)')
    return dict(projects=records, errors=errors)
