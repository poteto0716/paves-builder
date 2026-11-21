"""Crosslink series of one monomer project: one row per start-site density (crosslinks/*)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

COLUMNS = (
    ('crosslinking.actual_activation_fraction', 'Start sites / rings'),
    ('crosslinking.activation_density_mmol_cm3', 'Start sites per volume'),
    ('crosslinking.conversion', 'Crosslink conversion'),
    ('crosslinking.gel_epoxide_conversion', 'Ring conversion of the gel'),
    ('crosslinking.sol_mass_fraction', 'Removed sol mass fraction'),
    ('cooling.tg_K', 'Tg'),
    ('npt.density_g_cm3', 'Density (300 K)'),
    ('npt.alpha_linear_1_K', 'Linear CTE (300 K)'),
    ('structure.film_retention_volume_ratio', 'Film retention'),
    ('network.strand_density_nu_e_mmol_cm3', 'Effective strand density'),
    ('network.three_dimensional_gel_mass_fraction', 'Three-dimensional gel mass fraction'),
    ('nvt.dielectric_isotropic_fixed_charge', 'Dielectric response (fixed charge)'),
    ('mechanical.young_Hill_GPa', 'Young modulus (Hill)'),
)
SHOWN = ('crosslinking.conversion', 'cooling.tg_K', 'npt.density_g_cm3', 'structure.film_retention_volume_ratio',
         'network.strand_density_nu_e_mmol_cm3')
FIGURE = (('cooling.tg_K', 'Tg / K'), ('crosslinking.conversion', 'Crosslink conversion'),
          ('npt.density_g_cm3', 'Density (300 K) / g cm$^{-3}$'),
          ('structure.film_retention_volume_ratio', 'Film retention V/V$_0$'))


def series_rows(project):
    from workflows.network_properties.network_properties import children
    from workflows.network_properties.replica_summary import collect
    from workflows.network_properties.gpu_jobs import replica_complete
    rows = []
    for child in children(project):
        cfg = json.loads((child / 'project.json').read_text())
        _, properties = collect(child, cfg)
        options = cfg['crosslinking']
        row = dict(name=child.name, activation_count=options['activation_count'],
                   activation_fraction=options['activation_fraction'],
                   replicates=cfg['replicates'],
                   completed=sum(replica_complete(child, cfg, i) for i in range(1, cfg['replicates'] + 1)),
                   properties={})
        for key, _ in COLUMNS:
            entry = properties.get(key)
            if entry is not None:
                row['properties'][key] = dict(unit=entry['unit'], mean=entry['mean'], std=entry['std'],
                    n_available=entry['n_available'], n_checks_passed=entry['n_checks_passed'])
        rows.append(row)
    rows.sort(key=lambda r: (r['activation_count'] is None, r['activation_count'] or 0, r['activation_fraction'] or 0))
    return rows


def write_series(project):
    """crosslink_series.csv, crosslink_series.png and a markdown table for the monomer report."""
    project = Path(project)
    rows = series_rows(project)
    path = project / 'crosslink_series.csv'
    fields = ['name', 'activation_count', 'activation_fraction', 'replicates', 'completed']
    units = {}
    for row in rows:
        for key, value in row['properties'].items():
            units[key] = value['unit']
    for key, _ in COLUMNS:
        if key in units:
            fields += [f'{key}.mean', f'{key}.std', f'{key}.n_available', f'{key}.n_checks_passed']
    temporary = path.with_suffix('.csv.tmp')
    with temporary.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            flat = {k: row[k] for k in fields[:5]}
            for key, value in row['properties'].items():
                flat.update({f'{key}.{stat}': value[stat] for stat in ('mean', 'std', 'n_available', 'n_checks_passed')})
            writer.writerow(flat)
    temporary.replace(path)

    def cell(row, key):
        value = row['properties'].get(key)
        if not value or value['mean'] is None:
            return '—'
        text = f'{value["mean"]:.4g}' + (f' ± {value["std"]:.2g}' if value['std'] is not None else '')
        return text + ('*' if value['n_checks_passed'] < value['n_available'] else '')

    labels = dict(COLUMNS)
    lines = ['[Series CSV](crosslink_series.csv) · mean ± replica SD; * = some estimates have not passed checks.', '',
             '| Project | Starts | Replicas | ' + ' | '.join(
                 f'{labels[k]} / {units.get(k, "")}' for k in SHOWN) + ' |',
             '|---|---:|---:|' + '---:|' * len(SHOWN)]
    for row in rows:
        starts = (f'{row["activation_count"]}' if row['activation_count'] is not None
                  else f'f = {row["activation_fraction"]:g}')
        lines.append(f'| [{row["name"]}](crosslinks/{row["name"]}/report.md) | {starts} | '
                     f'{row["completed"]}/{row["replicates"]} | ' + ' | '.join(cell(row, k) for k in SHOWN) + ' |')
    figure = _plot(project, rows)
    if figure:
        lines += ['', f'![Crosslink series]({figure})']
    return dict(csv=path.name, rows=rows, markdown=lines)


def _plot(project, rows):
    points = [r for r in rows if (r['properties'].get('crosslinking.activation_density_mmol_cm3') or {}).get('mean') is not None]
    if not points:
        return None
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5))
    x = [r['properties']['crosslinking.activation_density_mmol_cm3']['mean'] for r in points]
    for axis, (key, title) in zip(axes.flat, FIGURE):
        values = [(xi, r['properties'].get(key) or {}) for xi, r in zip(x, points)]
        values = [(xi, v) for xi, v in values if v.get('mean') is not None]
        if values:
            axis.errorbar([xi for xi, _ in values], [v['mean'] for _, v in values],
                          yerr=[v['std'] or 0 for _, v in values], fmt='o-', capsize=3)
        axis.set(xlabel='Start sites per volume / mmol cm$^{-3}$', ylabel=title)
    smoke = json.loads((project / 'project.json').read_text()).get('smoke')
    fig.suptitle(f'{project.name}: crosslink series (error bars: replica SD)'
                 + ('\nSMOKE ONLY — no physical interpretation' if smoke else ''))
    fig.tight_layout()
    path = project / 'crosslink_series.png'
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path.name
