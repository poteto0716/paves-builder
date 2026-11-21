"""Join physical-property estimates by replica ID, retaining missing/unchecked values."""
from __future__ import annotations

import csv
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np

from workflows.network_properties.mechanical import result_current
from workflows.network_properties.network_graph import METHOD_VERSION


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
    mode = cfg.get('input_mode', 'saved_network')
    replica_ids = [f'replica_{i:02d}' for i in range(1, cfg['replicates'] + 1)]
    properties = {}
    construction_checks = {replica: bool((_read(project / replica / 'network_origin.json') or {}).get('target_reached'))
                           for replica in replica_ids} if mode == 'crosslink' else {}
    parent = None
    if mode == 'crosslink':
        from workflows.network_properties.construction import parent_project
        parent = parent_project(project, cfg)
    from workflows.network_properties.monomer import ring_label
    ring = ring_label(cfg)
    mechanical_present = any((project / replica / 'mechanical/result.json').is_file() for replica in replica_ids)

    def add(replica, source, data, group, key, label, unit, value, checked, display=True):
        finite = isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
        value = float(value) if finite else None
        smoke = bool(cfg.get('smoke') or (data or {}).get('smoke'))
        checked = bool(value is not None and checked and not smoke and construction_checks.get(replica, True))
        status = ('pending' if data is None else 'unavailable') if value is None else (
            'smoke' if smoke else 'checks_passed' if checked else 'unchecked')
        entry = properties.setdefault(key, dict(label=label, unit=unit, group=group,
                                                 display=display, replicas={}))
        try:
            where = str(source.relative_to(project))
        except ValueError:
            where = str(source)
        entry['replicas'][replica] = dict(value=value, checks_passed=checked, status=status, source=where)

    for replica in replica_ids:
        base = project / replica
        temperature = None
        if mode == 'monomer':
            from workflows.network_properties.monomer import label as temperature_label
            source = base / 'build/done.json'
            build = _read(source)
            for key, label, unit in (('molecules', 'Monomer molecules', 'count'), ('atoms', 'Atoms', 'count'),
                                     ('reactive_sites', f'{ring} rings', 'count')):
                add(replica, source, build, 'Monomer build', f'build.{key}', label, unit, _at(build, key), True)
            temperature = _read(base / temperature_label(cfg) / 'result.json')
            temperature_source = base / temperature_label(cfg) / 'result.json'
        if mode == 'crosslink':
            source = base / 'network_origin.json'
            build = _read(source)
            rows = [('conversion', 'Crosslink conversion', '1'),
                    ('epoxide_consumption_including_activation', 'Ring consumption including activation', '1'),
                    ('activation_count', 'Start sites', 'count'),
                    ('actual_activation_fraction', 'Start sites / rings', '1'),
                    ('activation_density_nm3', 'Start sites per volume (before crosslinking)', 'nm^-3'),
                    ('activation_density_mmol_cm3', 'Start sites per volume (before crosslinking)', 'mmol/cm^3'),
                    ('kinetic_chains.number_average_epoxide_units',
                     'Kinetic chain length per active end, number average', 'ring units'),
                    ('kinetic_chains.weight_average_epoxide_units',
                     'Kinetic chain length per active end, weight average', 'ring units'),
                    ('kinetic_chains.dispersity', 'Kinetic chain length dispersity', '1'),
                    ('kinetic_chains.max_epoxide_units', 'Longest kinetic chain', 'ring units')]
            if cfg.get('crosslinking', {}).get('remove_sol_fraction'):
                rows += [('gel_epoxide_conversion', 'Ring conversion of the gel (sol removed)', '1'),
                         ('sol_mass_fraction', 'Removed sol mass fraction', '1')]
            for key, label, unit in rows:
                add(replica, source, build, 'Network construction', f'crosslinking.{key}', label, unit,
                    _at(build, *key.split('.')), _at(build, 'target_reached'))
        source = base / 'dry_result.json'
        dry = _read(source)
        for key, label, unit in COOLING:
            add(replica, source, dry, 'Cooling', f'cooling.{key}', label, unit,
                _at(dry, 'cooling', key), _at(dry, 'cooling', 'fit_checks_passed'))
        for key, label, unit in NPT:
            add(replica, source, dry, '300 K thermodynamics', f'npt.{key}', label, unit,
                _at(dry, 'npt', 'estimates', key), _at(dry, 'npt', 'sampling_checks_passed'))
        if mode == 'monomer':
            # Rings / mean 300 K NPT volume; the ring count is exact, the volume carries the sampling check.
            for unit, suffix in (('nm^-3', 'nm3'), ('mmol/cm^3', 'mmol_cm3')):
                add(replica, temperature_source, temperature, '300 K thermodynamics', f'rings.number_density_300K_{suffix}',
                    f'{ring} number density (300 K)', unit,
                    _at(temperature, 'reactive_ring_density', f'number_density_300K_{suffix}'),
                    _at(dry, 'npt', 'sampling_checks_passed'))
        if mode == 'crosslink':
            # Film retention: 300 K NPT volume of the crosslinked network (sol removed when requested,
            # after the thermal history) over the 300 K NPT volume of the monomer replica it started from.
            # Both at 300 K and 1 atm, so thermal expansion does not enter.
            start = _at(_read(parent / replica / 'dry_result.json'), 'npt', 'volume_nm3_mean')
            relaxed = _at(dry, 'npt', 'volume_nm3_mean')
            add(replica, source, dry, '300 K structure', 'structure.film_retention_volume_ratio',
                'Film retention, V(300 K, crosslinked) / V(300 K, monomer)', '1',
                relaxed / start if isinstance(relaxed, Real) and isinstance(start, Real) and start > 0 else None,
                _at(dry, 'npt', 'sampling_checks_passed'))
        add(replica, source, dry, '300 K structure', 'structure.point_probe_void_fraction',
            'Point-probe void fraction (300 K)', '1',
            _at(dry, 'structure', 'point_probe_void_fraction'),
            _at(dry, 'npt', 'sampling_checks_passed'))
        network = _at(dry, 'network')
        # A result from an older definition stays visible but never counts as checked.
        network_current = _at(network, 'method_version') == METHOD_VERSION
        for path, label, unit in (() if mode == 'monomer' else (
                (('strand_density_nu_e_mmol_cm3',), 'Effective strand density nu_e', 'mmol/cm^3'),
                (('effective_strands',), 'Effective strands', 'count'),
                (('Mc_g_mol',), 'Mc, total cell mass / strands', 'g/mol'),
                (('Mc_active_network_g_mol',), 'Mc, active network mass / strands', 'g/mol'),
                (('strand_molecular_weight_g_mol', 'mean'), 'Mean strand molecular weight', 'g/mol'),
                (('strand_heavy_atoms', 'mean'), 'Mean strand length, interior heavy atoms', '1'),
                (('strand_backbone_bonds', 'mean'), 'Mean strand backbone bonds', '1'),
                (('junctions',), 'Junctions', 'count'),
                (('junction_density_mmol_cm3',), 'Junction density', 'mmol/cm^3'),
                (('mean_junction_functionality',), 'Mean junction functionality', '1'),
                (('cycle_rank_density_mmol_cm3',), 'Cycle rank density xi/V', 'mmol/cm^3'),
                (('dangling_mass_fraction',), 'Dangling mass fraction', '1'),
                (('loop_mass_fraction',), 'Loop mass fraction', '1'),
                (('sol_mass_fraction',), 'Sol mass fraction', '1'),
                (('three_dimensional_gel_mass_fraction',), 'Three-dimensional gel mass fraction', '1'),
                (('rubbery_shear_modulus_estimates_MPa', 'affine'), 'Rubbery G estimate, affine (300 K)', 'MPa'),
                (('rubbery_shear_modulus_estimates_MPa', 'phantom'), 'Rubbery G estimate, phantom (300 K)', 'MPa'),
                (('polar_groups_mmol_cm3', 'aliphatic_OH'), 'Aliphatic OH', 'mmol/cm^3'),
                (('polar_groups_mmol_cm3', 'phenolic_OH'), 'Phenolic OH', 'mmol/cm^3'),
                (('polar_groups_mmol_cm3', 'ether_O'), 'Ether O', 'mmol/cm^3'),
                (('polar_groups_mmol_cm3', 'epoxide_rings'), 'Unreacted epoxide', 'mmol/cm^3'),
                (('polar_groups_mmol_cm3', 'oxetane_rings'), 'Unreacted oxetane', 'mmol/cm^3'),
                (('polar_groups_mmol_cm3', 'aromatic_rings'), 'Aromatic rings', 'mmol/cm^3'))):
            # Topological counts of a fixed graph; no sampling check applies.
            value = _at(network, *path)
            if path[0] == 'polar_groups_mmol_cm3' and network is not None and value is None:
                value = 0.
            add(replica, source, dry, 'Network topology (300 K volume)', 'network.' + '.'.join(path),
                label, unit, value, network_current)
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

        if mode == 'monomer':
            t = float(cfg['temperature_K'])
            group = f'{t:g} K monomer liquid (after the thermal history)'
            checked = _at(temperature, 'npt', 'sampling_checks_passed')
            for key, label, unit in NPT:
                add(replica, temperature_source, temperature, group, f'at_T.npt.{key}',
                    label.replace('(300 K)', f'({t:g} K)'), unit, _at(temperature, 'npt', 'estimates', key), checked)
            for unit, suffix in (('nm^-3', 'nm3'), ('mmol/cm^3', 'mmol_cm3')):
                add(replica, temperature_source, temperature, group, f'at_T.rings.number_density_{suffix}',
                    f'{ring} number density ({t:g} K)', unit,
                    _at(temperature, 'reactive_ring_density', f'number_density_{t:g}K_{suffix}'), checked)
            rings = _at(temperature, 'reactive_rings')
            for cut in [f'{float(c):g}' for c in cfg.get('ring_contact_cutoffs_A', [])]:
                add(replica, temperature_source, temperature, group, f'at_T.rings.coordination_{cut}A',
                    f'{ring} coordination number ({cut} A)', '1',
                    _at(rings, 'by_cutoff_A', cut, 'coordination_number', 'mean'), checked)
                add(replica, temperature_source, temperature, group, f'at_T.rings.largest_cluster_{cut}A',
                    f'Largest {ring.lower()} contact cluster fraction ({cut} A)', '1',
                    _at(rings, 'by_cutoff_A', cut, 'largest_cluster_fraction', 'mean'), checked)
            alpha = _at(rings, 'msd_exponent_alpha')
            # D needs the diffusive regime: the log-log MSD slope over the fit window near 1.
            diffusive = isinstance(alpha, Real) and .8 <= alpha <= 1.2
            add(replica, temperature_source, temperature, group, 'at_T.rings.msd_final_A2',
                f'{ring} MSD at the end of sampling', 'A^2', _at(rings, 'msd_final_A2'), checked)
            add(replica, temperature_source, temperature, group, 'at_T.rings.diffusion_coefficient_m2_s',
                f'{ring} diffusion coefficient', 'm^2/s', _at(rings, 'diffusion_coefficient_m2_s'),
                bool(checked and diffusive))
            add(replica, temperature_source, temperature, group, 'at_T.rings.msd_exponent_alpha',
                f'{ring} MSD exponent alpha', '1', alpha, checked)
            for key, label, unit in (('Cv_direct_J_kg_K', f'Cv (direct NVT, {t:g} K)', 'J/(kg K)'),
                    ('dielectric_isotropic_fixed_charge', f'Dielectric response (fixed charge, {t:g} K)', '1')):
                add(replica, temperature_source, temperature, group, f'at_T.nvt.{key}', label, unit,
                    _at(temperature, 'nvt', key), _at(temperature, 'nvt', 'sampling_checks_passed'))

        if mechanical_present:
            source = base / 'mechanical/result.json'
            mechanical = _read(source) if result_current(base) else None
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
        'For crosslinked replicas, the requested crosslink conversion must also be reached '
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
