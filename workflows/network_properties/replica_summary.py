"""Join physical-property estimates by replica ID, retaining missing/unchecked values."""
from __future__ import annotations

import csv
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np


NPT = (
    ('density_g_cm3', 'Density (300 K)', 'g/cm^3'),
    ('Cp_J_kg_K', 'Cp', 'J/(kg K)'),
    ('Cv_J_kg_K', 'Cv (NPT identity)', 'J/(kg K)'),
    ('compressibility_isothermal_1_Pa', 'Isothermal compressibility', '1/Pa'),
    ('bulk_modulus_isothermal_GPa', 'Isothermal bulk modulus', 'GPa'),
    ('compressibility_isentropic_1_Pa', 'Isentropic compressibility', '1/Pa'),
    ('bulk_modulus_isentropic_GPa', 'Isentropic bulk modulus', 'GPa'),
    ('alpha_volume_1_K', 'Volumetric CTE (300 K)', '1/K'),
    ('alpha_linear_1_K', 'Linear CTE (300 K)', '1/K'),
    ('bulk_sound_speed_m_s', 'Bulk sound speed', 'm/s'),
)
COOLING = (
    ('tg_K', 'Tg', 'K'),
    ('density_at_tg_g_cm3', 'Density at Tg', 'g/cm^3'),
    ('alpha_linear_glass_1_K', 'Linear CTE (glass)', '1/K'),
    ('alpha_linear_rubber_1_K', 'Linear CTE (rubber)', '1/K'),
    ('alpha_volume_glass_1_K', 'Volumetric CTE (glass)', '1/K'),
    ('alpha_volume_rubber_1_K', 'Volumetric CTE (rubber)', '1/K'),
)
ELASTIC = (
    ('bulk_Voigt_GPa', 'Bulk modulus (Voigt)', 'GPa'),
    ('bulk_Reuss_GPa', 'Bulk modulus (Reuss)', 'GPa'),
    ('bulk_Hill_GPa', 'Bulk modulus (Hill)', 'GPa'),
    ('shear_Voigt_GPa', 'Shear modulus (Voigt)', 'GPa'),
    ('shear_Reuss_GPa', 'Shear modulus (Reuss)', 'GPa'),
    ('shear_Hill_GPa', 'Shear modulus (Hill)', 'GPa'),
    ('young_Hill_GPa', 'Young modulus (Hill)', 'GPa'),
    ('poisson_Hill', 'Poisson ratio (Hill)', '1'),
)


def _read(path):
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        # A native stage may still be writing its result during manual reporting.
        return None


def _at(data, *keys):
    for key in keys:
        if isinstance(data, dict):
            data = data.get(key)
        elif isinstance(data, list) and isinstance(key, int) and key < len(data):
            data = data[key]
        else:
            return None
    return data


def _statistics(values):
    if not values:
        return dict(mean=None, std=None, sem=None)
    std = float(np.std(values, ddof=1)) if len(values) > 1 else None
    return dict(mean=float(np.mean(values)), std=std,
                sem=std / math.sqrt(len(values)) if std is not None else None)


def collect(project, cfg):
    """Collect only physical quantities; sample counts and fit diagnostics are not averaged.

    Means include all finite estimates, explicitly labeled by their check counts.
    A separate checked_* summary includes only estimates passing existing checks.
    """
    project = Path(project)
    replica_ids = [f'replica_{i:02d}' for i in range(1, cfg['replicates'] + 1)]
    properties = {}
    construction_checks = {replica: bool((_read(project / replica / 'network_origin.json') or {}).get('target_reached'))
                           for replica in replica_ids} if cfg.get('input_mode') == 'epoxy' else {}

    def add(replica, source, data, group, key, label, unit, value, checked, display=True):
        finite = isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
        value = float(value) if finite else None
        smoke = bool(cfg.get('smoke') or (data or {}).get('smoke'))
        checked = bool(value is not None and checked and not smoke and construction_checks.get(replica, True))
        status = ('pending' if data is None else 'unavailable') if value is None else (
            'smoke' if smoke else 'checks_passed' if checked else 'unchecked')
        entry = properties.setdefault(key, dict(label=label, unit=unit, group=group,
                                                 display=display, replicas={}))
        entry['replicas'][replica] = dict(value=value, checks_passed=checked, status=status,
                                         source=str(source.relative_to(project)))

    for replica in replica_ids:
        base = project / replica
        if cfg.get('input_mode') == 'epoxy':
            source = base / 'network_origin.json'
            build = _read(source)
            for key, label in (('conversion', 'Crosslink conversion'),
                    ('epoxide_consumption_including_activation', 'Epoxide consumption including activation')):
                add(replica, source, build, 'Epoxy construction', f'crosslinking.{key}', label, '1',
                    _at(build, key), _at(build, 'target_reached'))
        source = base / 'dry_result.json'
        dry = _read(source)
        for key, label, unit in COOLING:
            add(replica, source, dry, 'Cooling', f'cooling.{key}', label, unit,
                _at(dry, 'cooling', key), _at(dry, 'cooling', 'fit_checks_passed'))
        for key, label, unit in NPT:
            add(replica, source, dry, '300 K thermodynamics', f'npt.{key}', label, unit,
                _at(dry, 'npt', 'estimates', key), _at(dry, 'npt', 'sampling_checks_passed'))
        add(replica, source, dry, '300 K structure', 'structure.point_probe_void_fraction',
            'Point-probe void fraction (300 K)', '1',
            _at(dry, 'structure', 'point_probe_void_fraction'),
            _at(dry, 'npt', 'sampling_checks_passed'))
        for key, label, unit in (('Cv_direct_J_kg_K', 'Cv (direct NVT)', 'J/(kg K)'),
                ('dielectric_isotropic_fixed_charge', 'Dielectric response (fixed charge)', '1')):
            add(replica, source, dry, '300 K NVT', f'nvt.{key}', label, unit,
                _at(dry, 'nvt', key), _at(dry, 'nvt', 'sampling_checks_passed'))
        for i, axis in enumerate('xyz'):
            for j, other in enumerate('xyz'):
                add(replica, source, dry, 'Dielectric tensor', f'nvt.dielectric_{axis}{other}',
                    f'Dielectric {axis}{other} (fixed charge)', '1',
                    _at(dry, 'nvt', 'dielectric_tensor_fixed_charge', i, j),
                    _at(dry, 'nvt', 'sampling_checks_passed'), display=False)

        if cfg.get('mechanical', True):
            source = base / 'mechanical/result.json'
            mechanical = _read(source)
            checked = bool(_at(mechanical, 'sampling_checks_passed') and _at(mechanical, 'isotropic', 'stable'))
            for key, label, unit in ELASTIC:
                add(replica, source, mechanical, 'Elasticity', f'mechanical.{key}', label, unit,
                    _at(mechanical, 'isotropic', key), checked)
            voigt = ('xx', 'yy', 'zz', 'yz', 'xz', 'xy')
            for i, axis in enumerate(voigt):
                for j, other in enumerate(voigt):
                    add(replica, source, mechanical, 'Stiffness tensor', f'mechanical.C_{axis}_{other}_GPa',
                        f'C({axis},{other})', 'GPa', _at(mechanical, 'symmetric_stiffness_GPa', i, j),
                        checked, display=False)

        if cfg.get('thermal_conductivity', True):
            source = base / 'thermal_conductivity/result.json'
            transport = _read(source)
            for key, label, unit in (
                    ('conductivity_direction_mean_W_m_K', 'Thermal conductivity (direction mean)', 'W/(m K)'),
                    ('thermal_diffusivity_m2_s', 'Thermal diffusivity', 'm^2/s')):
                add(replica, source, transport, 'Thermal transport', f'transport.{key}', label, unit,
                    _at(transport, key), _at(transport, 'sampling_checks_passed'))
            axes = {a['axis']: a for a in (transport or {}).get('axes', [])}
            for i in cfg.get('transport_axes', [0, 1, 2]):
                add(replica, source, transport, 'Thermal transport', f'transport.conductivity_{"xyz"[i]}_W_m_K',
                    f'Thermal conductivity ({"xyz"[i]})', 'W/(m K)',
                    _at(axes.get(i), 'conductivity_W_m_K'),
                    _at(axes.get(i), 'sampling_checks_passed'), display=False)

        for solvent in cfg.get('solvents', []):
            name = solvent['name']
            source = base / 'solvents' / name / 'result.json'
            wet = _read(source)
            fe_source = source if wet is not None else source.parent / 'infinite_dilution_free_energy/result.json'
            md_source = source if wet is not None else source.parent / 'diffusion/result.json'
            free_energy = _at(wet, 'infinite_dilution_free_energy') if wet is not None else _read(fe_source)
            fixed_loading = _at(wet, 'fixed_loading') if wet is not None else _read(md_source)
            for key, label in (('coupling_free_energy_kJ_mol', 'Insertion free energy'),
                               ('chemical_potential_relative_kJ_mol', 'Relative chemical potential')):
                add(replica, fe_source, wet if wet is not None else free_energy,
                    f'Solvent: {name}', f'solvents.{name}.{key}', label, 'kJ/mol',
                    _at(free_energy, key), _at(free_energy, 'convergence_verified'))
            add(replica, md_source, wet if wet is not None else fixed_loading,
                f'Solvent: {name}', f'solvents.{name}.diffusion_coefficient_m2_s',
                'Solvent self diffusion', 'm^2/s',
                _at(fixed_loading, 'diffusion', 'diffusion_coefficient_m2_s'),
                bool(_at(fixed_loading, 'fixed_loading_converged') and
                     _at(fixed_loading, 'diffusion', 'converged')))

    for entry in properties.values():
        values = [v['value'] for v in entry['replicas'].values() if v['value'] is not None]
        checked = [v['value'] for v in entry['replicas'].values() if v['checks_passed']]
        entry.update(n_expected=len(replica_ids), n_available=len(values), n_checks_passed=len(checked),
                     **_statistics(values),
                     **{f'checked_{k}': v for k, v in _statistics(checked).items()})
    return replica_ids, properties


def write_tables(project, replica_ids, properties):
    project = Path(project)
    def write(name, fields, rows):
        path = project / name
        temporary = path.with_suffix(path.suffix + '.tmp')
        # A BOM allows spreadsheet programs to recognize UTF-8 labels/units.
        with temporary.open('w', newline='', encoding='utf-8-sig') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    individual, summary = [], []
    for key, entry in properties.items():
        row = dict(property=key, label=entry['label'], unit=entry['unit'])
        for replica in replica_ids:
            row[replica] = entry['replicas'][replica]['value']
            row[f'{replica}_status'] = entry['replicas'][replica]['status']
        individual.append(row)
        summary.append(dict(property=key, label=entry['label'], unit=entry['unit'],
            **{k: entry[k] for k in ('n_expected', 'n_available', 'n_checks_passed',
                                    'mean', 'std', 'sem', 'checked_mean', 'checked_std', 'checked_sem')}))
    write('replica_properties.csv', ['property', 'label', 'unit'] +
          [column for r in replica_ids for column in (r, f'{r}_status')], individual)
    write('property_summary.csv', ['property', 'label', 'unit', 'n_expected', 'n_available',
          'n_checks_passed', 'mean', 'std', 'sem', 'checked_mean', 'checked_std', 'checked_sem'], summary)


def markdown(replica_ids, properties):
    lines = ['Replica property comparison:', '',
        '[Per-replica CSV](replica_properties.csv) · [Aggregate CSV](property_summary.csv)', '',
        'Means/SD/SEM below include all available finite estimates. An asterisk marks an estimate '
        'whose checks have not passed (including smoke runs). CSV checked_* columns summarize only '
        'estimates passing checks. Missing estimates remain blank; they are never treated as zero.', '',
        'For independently constructed epoxy, the requested crosslink conversion must also be reached '
        'for a replica to count as checks-passed.', '',
        'SD is the sample standard deviation across replicas; SEM = SD / sqrt(n). '
        'Both are undefined for a single available replica. Tensor components are included in the CSV files.']
    group = None
    def fmt(value):
        return '—' if value is None else f'{value:.6g}'
    for entry in properties.values():
        if not entry['display']:
            continue
        if group != entry['group']:
            group = entry['group']
            lines += ['', f'**{group}**', '',
                '| Property | Unit | ' + ' | '.join(replica_ids) + ' | Mean | SD | SEM | n | Checks passed |',
                '|---|---|' + '---:|' * (len(replica_ids) + 5)]
        cells = []
        for replica in replica_ids:
            value = entry['replicas'][replica]
            cells.append(fmt(value['value']) + ('*' if value['value'] is not None and not value['checks_passed'] else ''))
        star = '*' if entry['n_checks_passed'] < entry['n_available'] else ''
        cells += [fmt(entry['mean']) + star, fmt(entry['std']), fmt(entry['sem']),
                  f'{entry["n_available"]}/{entry["n_expected"]}',
                  f'{entry["n_checks_passed"]}/{entry["n_expected"]}']
        lines.append(f'| {entry["label"]} | {entry["unit"]} | ' + ' | '.join(cells) + ' |')
    return lines
