"""Standalone diagnostic figures linked from report.md."""
import csv
import json
from pathlib import Path

import numpy as np


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
