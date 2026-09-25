#!/usr/bin/env python3
"""Build the interactive results page (one self-contained HTML file) from runs/.

    python report.py --runs PROJECT/runs --out PROJECT/results.html   (or: adhesion.py report PROJECT)

Takes every stage that has finished. Hydrogens are left out of the 3D view and
coordinates are stored at 0.01 A, which keeps a 78k-atom system with a dozen
frames per trajectory under the page size limit. Numbers in the tiles and
charts are computed from all atoms.
"""
import argparse
import base64
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mdtools as T  # noqa: E402

QUANTUM = 0.01                      # A per int16 step
CONTACT_A = 3.0
HBOND_A = 2.5
COMPRESS_FRAMES = 13
AVOGADRO = 6.02214076e23


def rows(path):
    return list(csv.DictReader(open(path)))


def types_from_data(path, n):
    lines = Path(path).read_text().splitlines()
    k = lines.index('Masses')
    names = {}
    for line in lines[k + 2:]:
        if not line.strip():
            break
        w = line.split()
        names[int(w[0])] = w[-1]
    k = next(i for i, l in enumerate(lines) if l.startswith('Atoms'))
    typ = [''] * n
    for line in lines[k + 2:]:
        w = line.split('#')[0].split()
        if not w:
            break
        typ[int(w[0]) - 1] = names[int(w[2])]
    return np.array(typ)


class Interface:
    def __init__(self, runs):
        self.runs = runs
        info = json.loads((runs / '04_assemble/summary.json').read_text())
        self.nsi = info['silica_atoms']
        system = T.load_xml(runs / '04_assemble/out/interface.openmm_system.xml')
        self.mass = T.masses(system)
        self.n = len(self.mass)
        self.elem = ''.join({'H': 'H', 'C': 'C', 'O': 'O', 'Si': 'S'}.get(T.element_of(m), 'X') for m in self.mass)
        self.typ = types_from_data(runs / '04_assemble/out/interface.data', self.n)
        self.fixed = np.loadtxt(runs / '04_assemble/fixed_atoms.txt', dtype=int, ndmin=1)
        self.shown = np.array([i for i in range(self.n) if self.elem[i] != 'H'])

    def wrap(self, p, box, axes=(0, 1)):
        p = np.array(p, float)
        for a in axes:
            p[:, a] -= box[a] * np.floor(p[:, a] / box[a])
        return p

    def contacts(self, p, box):
        """Film atoms within CONTACT_A of any silica atom (periodic in x, y)."""
        q = self.wrap(p, box)
        tree = cKDTree(q[:self.nsi], boxsize=[box[0], box[1], 1e5])
        d, _ = tree.query(q[self.nsi:], k=1)
        return d, q

    def frame(self, p, box, label, subset='all'):
        q = self.wrap(p, box, (0, 1, 2) if subset == 'bulk' else (0, 1))
        if subset in ('bulk', 'polymer'):
            idx = self.shown[self.shown >= self.nsi] - self.nsi
            sel = q[idx]
            contact = []
        else:
            sel = q[self.shown]
            d, _ = self.contacts(p, box)
            contact = [int(i + self.nsi) for i in np.where(d < CONTACT_A)[0] if self.elem[i + self.nsi] != 'H']
        ints = np.round(sel / QUANTUM).astype(np.int16)
        return {'t': label, 'x': base64.b64encode(ints.tobytes()).decode(), 'c': contact}

    def metrics(self, p, box):
        d, q = self.contacts(p, box)
        A = box[0] * box[1]
        edges = np.arange(0, np.ceil(q[:, 2].max()) + 1, 1.0)
        ctr = (edges[:-1] + edges[1:]) / 2

        def prof(sel):
            h, _ = np.histogram(q[sel, 2], bins=edges, weights=self.mass[sel])
            return h / AVOGADRO / (A * 1e-24)
        film = np.arange(self.n) >= self.nsi
        pp, ps = prof(film), prof(~film)
        si_top = q[:self.nsi, 2].max()
        zp = q[self.nsi:, 2]
        lo, hi = np.percentile(zp, [20, 80])
        plateau = (ctr > lo + 5) & (ctr < hi - 5)
        rho = float(pp[plateau].mean())
        above = ctr[pp > 0.5 * rho]
        hoy = np.where(self.typ == 'hoy')[0]
        hoy = hoy[q[hoy, 2] > si_top - 6]
        po = np.where(film & np.char.startswith(self.typ.astype(str), 'o'))[0]
        tree = cKDTree(q[po], boxsize=[box[0], box[1], 1e5])
        dh, jh = tree.query(q[hoy], k=1)
        hb = dh < HBOND_A
        carbonyl = int((self.typ[po[jh[hb]]] == 'o_1').sum())
        return {
            'profile': (ctr, pp, ps), 'rho': rho, 'thickness': float(above.max() - above.min()),
            'extent': (float(above.min()), float(above.max())),
            'c3': int((d < 3).sum()), 'c5': int((d < 5).sum()), 'closest': float(d.min()),
            'hb': int(hb.sum()), 'silanols': int(len(hoy)), 'hb_carbonyl': carbonyl,
        }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--name', default='adhesion')
    a = ap.parse_args(argv)
    runs = Path(a.runs)
    I = Interface(runs)
    stages, charts, metrics = [], [], []

    st = T.load_xml(runs / '02_bulk/state.xml')
    b = T.box_A(st)
    stages.append({'id': 'bulk', 'num': '02', 'label': 'Bulk melt', 'sub': 'NPT (z only), 550→300 K', 'subset': 'polymer',
                   'box': b.tolist(), 'frames': [I.frame(T.positions_A(st), b, 'end of stage', 'bulk')]})
    st = T.load_xml(runs / '03_surface/state.xml')
    b = T.box_A(st)
    stages.append({'id': 'surface', 'num': '03', 'label': 'Free film', 'sub': 'whole molecules, vacuum in z, NVT 550→300 K',
                   'subset': 'polymer', 'box': b.tolist(), 'frames': [I.frame(T.positions_A(st), b, 'end of stage', 'polymer')]})
    st = T.load_xml(runs / '04_assemble/out/interface.openmm_state.xml')
    b = T.box_A(st)
    stages.append({'id': 'assemble', 'num': '04', 'label': 'Placed on silica', 'sub': 'film 3.0 Å above the SiO₂ surface',
                   'subset': 'all', 'box': b.tolist(), 'frames': [I.frame(T.positions_A(st), b, 'as assembled')]})
    latest = ('assemble', T.positions_A(st), b)

    if (runs / '05_compress/summary.json').exists():
        import mdtraj as md
        dcd, top = runs / '05_compress/trajectory.dcd', runs / '05_compress/topology.pdb'
        n_frames = len(md.open(str(dcd)))
        frames = []
        for k in np.unique(np.geomspace(1, n_frames, COMPRESS_FRAMES).round().astype(int) - 1):
            f = md.load_frame(str(dcd), int(k), top=str(top))
            frames.append(I.frame(f.xyz[0] * 10, b, f'frame {k + 1} of {n_frames}'))
        stages.append({'id': 'compress', 'num': '05', 'label': 'Pressed at 200 MPa', 'sub': 'wall above the film, 550 K, bottom 3 Å fixed',
                       'subset': 'all', 'box': b.tolist(), 'frames': frames})
        st = T.load_xml(runs / '05_compress/state.xml')
        latest = ('compress', T.positions_A(st), T.box_A(st))
        log = rows(runs / '05_compress/log.csv')
        tp = [float(r['time_ps']) for r in log]
        charts.append({'title': 'Wall and film top', 'sub': '05 · the wall descends at ≤ 1 Å/ps until the film pushes back',
                       'x': 'time (ps)', 'unit': 'Å', 'xd': 0, 'yd': 1,
                       'series': [{'name': 'wall position', 'short': 'wall', 'color': '--s1', 'dy': -7,
                                   'pts': [[t, float(r['wall_z_A'])] for t, r in zip(tp, log)][::4]},
                                  {'name': 'highest film atom', 'short': 'film top', 'color': '--s2', 'dy': 9,
                                   'pts': [[t, float(r['polymer_z_max_A'])] for t, r in zip(tp, log)][::4]}]})
        pr = np.array([float(r['wall_pressure_MPa']) for r in log])
        win = max(1, min(10, len(pr)))
        smooth = np.convolve(pr, np.ones(win) / win, mode='valid')
        charts.append({'title': 'Wall pressure', 'sub': f'05 · {win}-ps running mean of the instantaneous value',
                       'x': 'time (ps)', 'unit': 'MPa', 'xd': 0, 'yd': 1, 'ref': 200, 'refLabel': 'target 200 MPa',
                       'series': [{'name': 'wall pressure', 'short': 'pressure', 'color': '--s1',
                                   'pts': [[tp[k + win - 1], float(v)] for k, v in enumerate(smooth)][::4]}]})
        mean_p = float(pr[np.array(tp) > 300].mean())
        sd_p = float(pr[np.array(tp) > 300].std())

    for sid, num, label, sub in (('06_cool', '06', 'Cooled', 'wall released, 550→300 K, 0.25 fs'),
                                 ('07_relax', '07', 'Relaxed at 300 K', '300 K, 0.25 fs')):
        if (runs / sid / 'summary.json').exists():
            st = T.load_xml(runs / sid / 'state.xml')
            stages.append({'id': sid, 'num': num, 'label': label, 'sub': sub, 'subset': 'all', 'box': b.tolist(),
                           'frames': [I.frame(T.positions_A(st), b, 'end of stage')]})
            latest = (sid, T.positions_A(st), T.box_A(st))

    pulled = []
    if (runs / '09_pull/summary.json').exists():
        import mdtraj as md
        dcd, top = runs / '09_pull/trajectory.dcd', runs / '09_pull/topology.pdb'
        n_frames = len(md.open(str(dcd)))
        frames = []
        for k in np.linspace(0, n_frames - 1, 12).round().astype(int):
            f = md.load_frame(str(dcd), int(k), top=str(top))
            frames.append(I.frame(f.xyz[0] * 10, b, f'{k * 5000 * 0.25 / 1000:.0f} ps'))
        stages.append({'id': 'pull', 'num': '09', 'label': 'Pulled', 'sub': 'top 25 % on a spring, 10 m/s', 'subset': 'all',
                       'pulled': True, 'box': b.tolist(), 'frames': frames})
        pulled = np.loadtxt(runs / '09_pull/pulled_atoms.txt', dtype=int).tolist()
        pf = rows(runs / '09_pull/pull_force.csv')
        step = max(1, len(pf) // 1500)
        pf = pf[::step]
        charts.append({'title': 'Pull force on the film', 'sub': f'09 · spring 100 kcal/mol/Å², + = upward (every {step}th row)',
                       'x': 'time (ps)', 'unit': 'nN', 'xd': 1, 'yd': 3,
                       'series': [{'name': 'force', 'short': 'force', 'color': '--s1',
                                   'pts': [[float(r['time_ps']), float(r['force_nN'])] for r in pf]}]})

    name, p, box = latest
    M = I.metrics(p, box)
    if (runs / '09_pull/summary.json').exists():
        st = T.load_xml(runs / '09_pull/state.xml')
        P = I.metrics(T.positions_A(st), T.box_A(st))
        pf_all = np.loadtxt(runs / '09_pull/pull_force.csv', delimiter=',', skiprows=1)
        dt_row = float(pf_all[1, 1] - pf_all[0, 1]) if len(pf_all) > 1 else 0.00025   # ps between rows
        w = max(1, min(int(round(0.5 / dt_row)), len(pf_all)))                       # a 0.5-ps window
        smooth_f = np.convolve(pf_all[:, 6], np.ones(w) / w, mode='same')
        k = int(np.argmax(smooth_f[w // 2:len(smooth_f) - w // 2]) + w // 2) if len(pf_all) > w else int(np.argmax(pf_all[:, 6]))
        work = float(np.trapezoid(pf_all[:, 5], pf_all[:, 2]))
        area = float(box[0] * box[1])
        metrics.append({'k': 'Peak pull force', 'v': f'{smooth_f[k]:.1f}', 'u': 'nN',
                        'd': f'{w * dt_row:.1f}-ps mean, at {pf_all[k, 1]:.0f} ps (reference +{pf_all[k, 2] - pf_all[0, 2]:.1f} Å); '
                             f'spring work {work / area * 694.77:.0f} mJ/m²'})
        metrics.append({'k': 'Contact after pulling', 'v': f"{P['c3']:,}", 'u': f"of {M['c3']:,}",
                        'd': 'film atoms within 3 Å of silica, end of pull vs before; a remaining layer means failure inside the film'})
        c2, pp2, ps2 = P['profile']
        keep2 = c2 < P['extent'][1] + 15
        charts.append({'title': 'Density after pulling', 'sub': 'end of 09 · where the film thinned or separated',
                       'x': 'z (Å)', 'unit': 'g/cm³', 'xd': 1, 'yd': 3,
                       'series': [{'name': 'silica', 'short': 'SiO₂', 'color': '--s2', 'dy': 9,
                                   'pts': [[float(z), float(v)] for z, v in zip(c2[keep2], ps2[keep2])]},
                                  {'name': 'polymer film', 'short': 'polymer', 'color': '--s1', 'dy': -7,
                                   'pts': [[float(z), float(v)] for z, v in zip(c2[keep2], pp2[keep2])]}]})
    ctr, pp, ps = M['profile']
    keep = ctr < M['extent'][1] + 15
    charts.append({'title': 'Density across the interface', 'sub': f'end of {name} · 1-Å slabs over the full cell area',
                   'x': 'z (Å)', 'unit': 'g/cm³', 'xd': 1, 'yd': 3,
                   'series': [{'name': 'silica', 'short': 'SiO₂', 'color': '--s2', 'dy': 9,
                               'pts': [[float(z), float(v)] for z, v in zip(ctr[keep], ps[keep])]},
                              {'name': 'polymer film', 'short': 'polymer', 'color': '--s1', 'dy': -7,
                               'pts': [[float(z), float(v)] for z, v in zip(ctr[keep], pp[keep])]}]})
    bulk_log = rows(runs / '02_bulk/log.csv')
    charts.insert(0, {'title': 'Bulk density', 'sub': '02 · z-only NPT at 1 atm: 550 K, 550→300 K, 300 K',
                      'x': 'time (ps)', 'unit': 'g/cm³', 'xd': 0, 'yd': 4,
                      'series': [{'name': 'density', 'short': 'density', 'color': '--s1',
                                  'pts': [[float(r['time_ps']), float(r['density_g_cm3'])] for r in bulk_log][::4]}]})

    bulk_rho = float(bulk_log[-1]['density_g_cm3'])
    metrics += [
        {'k': 'Film density', 'v': f"{M['rho']:.3f}", 'u': 'g/cm³', 'd': f'interior, end of {name}; bulk at 300 K was {bulk_rho:.3f}'},
        {'k': 'Film thickness', 'v': f"{M['thickness']:.0f}", 'u': 'Å', 'd': f"z {M['extent'][0]:.0f}–{M['extent'][1]:.0f} Å (above half the interior density)"},
        {'k': 'In contact', 'v': f"{M['c3']:,}", 'u': 'atoms', 'd': f"film atoms within 3 Å of silica ({M['c5']:,} within 5 Å); closest {M['closest']:.2f} Å"},
        {'k': 'Silanol H-bonds', 'v': f"{M['hb']}", 'u': f"/ {M['silanols']}", 'd': f"top-surface silanol H within 2.5 Å of a polymer O; {M['hb_carbonyl']} to carbonyl O (PCFF o_1)"},
    ]
    if name == 'compress' or (runs / '05_compress/summary.json').exists():
        metrics.append({'k': 'Wall pressure', 'v': f'{mean_p:.0f}', 'u': f'± {sd_p:.0f} MPa',
                        'd': 'mean ± SD of the instantaneous value over 300–600 ps; target 200'})

    ie = runs / '08_interface/summary.json'
    if ie.exists():
        S = json.loads(ie.read_text())
        g, ge = S['interface_energy_mJ_m2']
        e, ee = S['E_int_kcal']
        metrics.insert(0, {'k': 'Interface energy', 'v': f'{g:.1f}', 'u': f'± {ge:.1f} mJ/m²',
                           'd': f"E_int {e:.0f} ± {ee:.0f} kcal/mol over {S['snapshots']} snapshots "
                                f"(vdW {S['E_int_vdw_kcal'][0]:.0f}, Coulomb {S['E_int_coul_kcal'][0]:.0f}); work of adhesion {-g:.1f}"})
        er = rows(runs / '08_interface/interface_energy.csv')
        charts.append({'title': 'Interface energy per snapshot', 'sub': '08 · E_int / A, E_int = E(all) − E(polymer) − E(silica), nonbonded',
                       'x': 'time (ps)', 'unit': 'mJ/m²', 'xd': 2, 'yd': 1,
                       'series': [{'name': 'interface energy', 'short': 'γ', 'color': '--s1',
                                   'pts': [[float(r['time_ps']), float(r['gamma_mJ_m2'])] for r in er]}]})
    done = [s['num'] for s in stages]
    data = {
        'meta': {'eyebrow': f'polypaves adhesion workflow · {a.name}',
                 'title': f'{a.name}: adhesion on silica',
                 'lede': (f'Polymer film (<b>{I.n - I.nsi:,} atoms</b>) pressed onto IFF amorphous '
                          f'silica (<b>{I.nsi:,} atoms</b>), PC-IFF in OpenMM. Stages finished: {", ".join(done)}. '
                          'The numbers below describe the most recent interface stage.'),
                 'hint': 'drag to rotate · scroll to zoom · z is up · hydrogens hidden'},
        'metrics': metrics, 'charts': charts,
        'notes': [
            {'h': 'How the film was made', 'p': 'The melt relaxed with only z free (1 atm), was made whole and given vacuum in z. Long chains wind through several periodic images, so the free film spreads out; the wall compacts it against the silica.'},
            {'h': 'What “in contact” means here', 'p': 'Film atoms within 3 Å of any silica atom, hydrogens included, with x and y periodic. Silanol hydrogen bonds count top-surface silanol H within 2.5 Å of a polymer oxygen.'},
            {'h': 'Limits', 'p': 'OpenMM has no slab correction, so z is periodic with a vacuum gap. Densities come from sub-nanosecond runs and a fast quench, and are typically below experimental glass densities.'},
        ],
        'n_si': int(I.nsi), 'elements': I.elem, 'shown': I.shown.tolist(), 'has_h': False,
        'fixed': I.fixed.tolist(), 'pulled': pulled, 'quantum': QUANTUM, 'stages': stages,
        'start_stage': len(stages) - 1,
    }
    blob = json.dumps(data, separators=(',', ':'))
    blob = blob.replace('</', '<\\/')           # a JSON escape; keeps </script> out of the block
    tpl = (HERE / 'report_template.html').read_text()
    script = '<script>\n' + (HERE / 'report_script.js').read_text() + '</script>\n'
    html = tpl.replace('__DATA__', blob).replace('__SCRIPT__', script)
    Path(a.out).write_text(html)
    print(f'{a.out}: {len(html) / 1e6:.2f} MB, stages {done}')
    print(json.dumps({k: v for k, v in M.items() if k != 'profile'}))


if __name__ == '__main__':
    main()
