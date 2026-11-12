"""Standalone diagnostic figures linked from report.md."""
import csv
import json
from pathlib import Path

import numpy as np


def plot_tg(base, dry=None, path=None):
    """Temperature-density and temperature-volume with the Tg fit, for one replica.

    Cooling samples, their 10 K bin means, the fitted two-asymptote curve and Tg
    are drawn together with the density/volume at the end of every heating
    phase, so an expansion far beyond thermal expansion (boiling or cavitation
    at the high temperature) is visible next to the fit. Returns the PNG path.
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    base = Path(base)
    dry = dry if dry is not None else json.loads((base/'dry_result.json').read_text())
    with (base/'cooling.csv').open() as f:
        cooling = list(csv.DictReader(f))
    t = np.array([float(r['target_temperature_K']) for r in cooling])
    rho = np.array([float(r['density_g_cm3']) for r in cooling])
    volume = np.array([float(r['volume_nm3']) for r in cooling])
    mass_g = float(np.median(rho * volume)) * 1e-21  # g/cm3 x nm3 -> g; constant topology
    to_volume = lambda density: mass_g / np.asarray(density, float) * 1e21
    fit = dry.get('cooling', {})
    history = (dry.get('thermal_history_density') or {}).get('phases', [])
    heating = [h for h in history if h['phase'][:2] in ('00', '01', '02') and h.get('final_density_g_cm3')]
    reheat = [h for h in history if h['phase'][:2] in ('04', '05') and h.get('final_density_g_cm3')]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, values, convert, label in ((axes[0], rho, lambda d: np.asarray(d, float), 'Density / g cm$^{-3}$'),
                                       (axes[1], volume, to_volume, 'Volume / nm$^3$')):
        ax.scatter(t, values, s=4, alpha=.15, c='tab:blue', label='cooling samples')
        if 'bins' in fit:
            bins = sorted(fit['bins'], key=lambda r: r['temperature_K'])
            ax.plot([b['temperature_K'] for b in bins], convert([b['density_g_cm3'] for b in bins]),
                    'o', ms=3.5, c='tab:blue', label='cooling, 10 K bin mean')
            tg, width = fit['tg_K'], fit['transition_width_K']
            slope_glass, slope_rubber = fit['density_slope_glass_g_cm3_K'], fit['density_slope_rubber_g_cm3_K']
            change = slope_rubber - slope_glass
            rho0 = fit['density_at_tg_g_cm3'] - .5 * change * width
            grid = np.linspace(t.min(), t.max(), 400)
            x = grid - tg
            curve = rho0 + slope_glass * x + .5 * change * (x + np.sqrt(x * x + width * width))
            ax.plot(grid, convert(curve), c='black', lw=1.2, label='fit')
            for slope, style, name in ((slope_glass, ':', 'glass asymptote'), (slope_rubber, '--', 'rubber asymptote')):
                ax.plot(grid, convert(fit['density_at_tg_g_cm3'] + slope * x), c='gray', ls=style, lw=.8, label=name)
            ax.axvline(tg, c='tab:red', ls='--', lw=1, label=f'Tg {tg:.0f} K')
        if heating:
            ax.plot([h['temperature_K'] for h in heating], convert([h['final_density_g_cm3'] for h in heating]),
                    's-', ms=5, c='tab:orange', lw=.8, label='heating (end of each phase)')
        if reheat:
            ax.plot([h['temperature_K'] for h in reheat], convert([h['final_density_g_cm3'] for h in reheat]),
                    '^', ms=4, c='tab:green', label='reheat to 300 K')
        ax.set(xlabel='Temperature / K', ylabel=label)
    axes[0].legend(fontsize=7, loc='lower left')
    if 'tg_K' in fit:
        title = (f'Tg {fit["tg_K"]:.0f} K, alpha_L glass {fit["alpha_linear_glass_1_K"]:.2e} / rubber '
                 f'{fit["alpha_linear_rubber_1_K"]:.2e} 1/K, R$^2$ {fit["r_squared"]:.3f}, '
                 f'checks passed: {fit.get("fit_checks_passed")}')
    else:
        title = f'Tg fit unavailable: {fit.get("reason", "")}'[:120]
    warning = fit.get('vaporization_warning')
    if warning:
        title += '\nWARNING: ' + warning
    if dry.get('smoke'):
        title = 'SMOKE ONLY - no physical interpretation\n' + title
    fig.suptitle(title, fontsize=9, color='tab:red' if warning else 'black')
    fig.tight_layout()
    path = Path(path) if path else base/'tg_density_volume.png'
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def render(project):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    project = Path(project)
    figures = []
    def save(fig, path, title):
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append((title,str(path.relative_to(project))))
    for base in sorted(project.glob('replica_*')):
        if not (base/'dry_result.json').exists():
            continue
        dry = json.loads((base/'dry_result.json').read_text())
        with (base/'cooling.csv').open() as f:
            cooling = list(csv.DictReader(f))
        fig, axes = plt.subplots(1,2,figsize=(10,3.5))
        t = np.array([float(r['target_temperature_K']) for r in cooling])
        rho = np.array([float(r['density_g_cm3']) for r in cooling])
        axes[0].scatter(t,rho,s=5,alpha=.3)
        fit = dry['cooling']
        if 'bins' in fit:
            bins = sorted(fit['bins'],key=lambda r:r['temperature_K'])
            axes[0].plot([r['temperature_K'] for r in bins], [r['fit_density_g_cm3'] for r in bins],c='black')
            axes[0].axvline(fit['tg_K'],c='tab:red',ls='--',label=f'Tg estimate {fit["tg_K"]:.1f} K')
            axes[0].legend(fontsize=8)
        axes[0].set(xlabel='Temperature / K',ylabel='Density / g cm$^{-3}$',title='Cooling: check fit and rate dependence')
        rows = json.loads((base/'06_production_300K/samples.json').read_text())
        axes[1].plot([r['time_ps']/1000 for r in rows],[r['density_g_cm3'] for r in rows],lw=.7)
        axes[1].set(xlabel='Time / ns',ylabel='Density / g cm$^{-3}$',title='300 K NPT: check stationary fluctuations')
        if dry['smoke']:
            fig.suptitle('SMOKE ONLY — no physical interpretation')
        save(fig,base/'thermal_diagnostics.png',f'{base.name}: cooling and density')
        figures.append((f'{base.name}: Tg, temperature-density and temperature-volume',
                        str(plot_tg(base, dry).relative_to(project))))
        network = dry.get('network')
        if network and network.get('effective_strands'):
            from workflows.network_properties.network_graph import plot_strand_histogram, write_strand_csv
            write_strand_csv(network, base/'network_strands.csv')
            plot_strand_histogram(network, base/'network_strands.png')
            figures.append((f'{base.name}: network strand lengths', str((base/'network_strands.png').relative_to(project))))
        for path in sorted(base.glob('solvents/*/result.json')):
            result = json.loads(path.read_text())
            fe = result['infinite_dilution_free_energy']
            d = result['fixed_loading']['diffusion']
            fig, axes = plt.subplots(1,2,figsize=(10,3.8))
            if d.get('lag_ps'):
                axes[0].loglog(d['lag_ps'],d['msd_nm2'],'o-',ms=3)
            axes[0].set(xlabel='Lag time / ps',ylabel='MSD / nm$^2$',
                        title=f'Diffusion regime verified: {d["converged"]}')
            matrix = np.asarray(fe['overlap_matrix'])
            im = axes[1].imshow(matrix,vmin=0,vmax=1,cmap='viridis')
            axes[1].set(xlabel='Lambda state',ylabel='Lambda state',
                        title=f'FE sampling verified: {fe["sampling_diagnostics_passed"]}')
            fig.colorbar(im,ax=axes[1],label='MBAR overlap')
            if result['smoke']:
                fig.suptitle('SMOKE ONLY — no physical interpretation')
            save(fig,path.parent/'diagnostics.png',f'{base.name}: {result["solvent"]["name"]}')
    return figures
