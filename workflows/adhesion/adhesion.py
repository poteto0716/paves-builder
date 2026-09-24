#!/usr/bin/env python3
"""Polymer / silica adhesion workflow: build, relax, press, relax, measure, pull.

    python adhesion.py new NAME --monomer SMILES --dp N --chains N [options]
    python adhesion.py new NAME --polyse chemistry.polyse [options]      (copolymers etc.)
    python adhesion.py run PROJECT [--from STAGE] [--to STAGE] [--only STAGE]
    python adhesion.py status PROJECT
    python adhesion.py report PROJECT

Stages, in order:

    build      the polymer melt (polyse), cell x/y = the silica supercell
    bulk       pre-relaxation, then NPT with only z free
    surface    whole molecules, vacuum in z, NVT
    assemble   the film on the silica supercell (polyse writes the combined system)
    compress   a wall presses the film onto the silica; bottom of the silica fixed
    cool       wall released, cooled
    relax      held at the final temperature
    interface  interface energy: E(all) - E(polymer) - E(silica), averaged over snapshots
    pull       the top of the film pulled away on a spring; force every step

`run` starts at the first stage that has not finished and continues to the end.
--from STAGE reruns from STAGE, --to STAGE stops after it, --only STAGE runs one.
A finished stage is one whose runs/<NN_stage>/summary.json exists.
"""
import argparse
import copy
import json
import sys
import time
from pathlib import Path

WORKFLOW = Path(__file__).resolve().parent
BUILDER = WORKFLOW.parents[1]


# pcff-iff-mod: torsions across an angle of theta0 >= 150 deg (Si-O-Si) are switched
# off between 175 and 180 deg, where they are singular; otherwise identical.
FORCEFIELDS = {'pcff-iff': '@builder/examples/forcefields/pcff_iff_long_bulk.ff',
               'pcff-iff-mod': '@builder/examples/forcefields/pcff_iff_long_bulk_mod.ff'}


def resolve(value, project):
    """@builder/…, @workflow/… and @project/… paths, made absolute."""
    if not isinstance(value, str):
        return value
    for token, base in (('@builder/', BUILDER), ('@workflow/', WORKFLOW), ('@project/', project)):
        if value.startswith(token):
            return str(base / value[len(token):])
    return value


def load_project(path):
    project = Path(path).resolve()
    f = project / 'project.json' if project.is_dir() else project
    if not f.exists():
        sys.exit(f'no project.json at {project}; create one with: adhesion.py new NAME ...')
    cfg = json.loads(f.read_text())
    project = f.parent
    for key in ('forcefield', 'silica_mol2'):
        cfg[key] = resolve(cfg[key], project)
        if not Path(cfg[key]).exists():
            sys.exit(f'{key} not found: {cfg[key]}')
    return project, cfg


def cmd_new(a):
    project = Path(a.dir or Path.cwd() / 'projects') / a.name
    if (project / 'project.json').exists() and not a.force:
        sys.exit(f'{project}/project.json exists; use --force to overwrite it')
    cfg = json.loads((WORKFLOW / 'defaults.json').read_text())
    cfg.pop('_comment', None)
    cfg = {'name': a.name, 'system': {}, **cfg}
    y = cfg['system']
    sys.path.insert(0, str(WORKFLOW))
    import stages
    y.update({'build_density': a.density, 'seed': a.seed})
    if a.polyse:
        if a.monomer or a.dp or a.chains:
            sys.exit('give either --polyse FILE or --monomer/--dp/--chains, not both')
        try:
            kept, carried, dropped = stages.split_polyse(Path(a.polyse).read_text())
        except ValueError as e:
            sys.exit(f'{a.polyse}: {e}')
        y['polyse_lines'] = kept
        y['polyse_source'] = str(Path(a.polyse).resolve())
        for k, v in carried.items():
            y[k] = int(v) if k != 'name' else v
        if dropped:
            print(f'  set by the workflow, ignored from {Path(a.polyse).name}: {", ".join(dict.fromkeys(dropped))}')
    else:
        if not (a.monomer and a.dp and a.chains):
            sys.exit('give --monomer, --dp and --chains, or --polyse FILE (copolymers, blocks, random, ...)')
        y.update({'monomer': a.monomer, 'terminator': a.terminator, 'dp': a.dp, 'chains': a.chains})
    cfg['assemble']['supercell'] = [a.supercell[0], a.supercell[1], 1]
    cfg['forcefield'] = FORCEFIELDS[a.forcefield]
    if a.test:
        cfg['step_scale'] = 0.01
    # polyse's own parser checks the build input now, not hours into a run.
    check = dict(cfg, forcefield=resolve(cfg['forcefield'], project), silica_mol2=resolve(cfg['silica_mol2'], project))
    text, (Lx, Ly) = stages.polymer_input(check)
    project.mkdir(parents=True, exist_ok=True)
    probe = project / '.check.polyse'
    probe.write_text(text)
    try:
        from polyse.config import InputParser
        InputParser.parse(str(probe)).resolved_input()
    except Exception as e:
        sys.exit(f'polyse rejects this polymer input: {e}\n(the input it was given: {probe})')
    probe.unlink()
    (project / 'project.json').write_text(json.dumps(cfg, indent=1) + '\n')
    print(f'created {project}/project.json')
    what = (f'{Path(a.polyse).name}' if a.polyse else f'{y["monomer"]}  DP {y["dp"]} x {y["chains"]} chains')
    print(f'  polymer: {what}; substrate {a.supercell[0]} x {a.supercell[1]} silica ({Lx:.1f} x {Ly:.1f} A)'
          + ('  [test: steps x 0.01]' if a.test else ''))
    print(f'next:  python {Path(__file__).name} run {project}')


def cmd_status(a):
    project, cfg = load_project(a.project)
    sys.path.insert(0, str(WORKFLOW))
    import stages
    print(f'{cfg["name"]}  ({project})')
    for name, (d, _) in stages.STAGES.items():
        summ = project / 'runs' / d / 'summary.json'
        if summ.exists():
            s = json.loads(summ.read_text())
            extra = ''
            if 'interface_energy_mJ_m2' in s:
                g, e = s['interface_energy_mJ_m2']
                extra = f'  gamma {g:.1f} ± {e:.1f} mJ/m2'
            elif 'density_g_cm3' in s:
                extra = f'  density {s["density_g_cm3"]:.3f} g/cm3'
            print(f'  done     {name:10s} {d:13s} {s.get("wall_seconds", 0) / 60:6.1f} min{extra}')
        elif (project / 'runs' / d).exists():
            print(f'  started  {name:10s} {d:13s} (no summary: unfinished or failed; see runs/{name}.log)')
        else:
            print(f'  -        {name:10s} {d}')


def cmd_run(a):
    project, cfg = load_project(a.project)
    sys.path.insert(0, str(WORKFLOW))
    import stages
    names = list(stages.STAGES)
    runs = project / 'runs'
    runs.mkdir(exist_ok=True)
    for opt in (a.from_, a.to, a.only):
        if opt and opt not in names:
            sys.exit(f'unknown stage {opt!r}; stages are: {", ".join(names)}')
    if a.only:
        todo = [a.only]
    else:
        first = names.index(a.from_) if a.from_ else next(
            (i for i, n in enumerate(names) if not (runs / stages.STAGES[n][0] / 'summary.json').exists()), len(names))
        last = names.index(a.to) if a.to else len(names) - 1
        todo = names[first:last + 1]
    if not todo:
        print('every stage has finished; use --from STAGE to rerun')
        return
    # A stage needs the one before it.
    first_idx = names.index(todo[0])
    if first_idx > 0:
        prev = stages.STAGES[names[first_idx - 1]][0]
        if not (runs / prev / 'summary.json').exists():
            sys.exit(f'{todo[0]} needs {names[first_idx - 1]} ({prev}) to have finished first')
    print(f'{cfg["name"]}: running {", ".join(todo)}')
    for name in todo:
        d, fn = stages.STAGES[name]
        t0 = time.time()
        print(f'=== {name} ({d}) {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
        log = runs / f'{name}.log'
        with open(log, 'w') as fh:
            old_out, old_err = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = Tee(fh, old_out)
            try:
                fn(copy.deepcopy(cfg), runs)
            except SystemExit as e:
                sys.stdout, sys.stderr = old_out, old_err
                sys.exit(f'FAILED: {name}: {e} (log: {log})')
            except Exception:
                import traceback
                traceback.print_exc()
                sys.stdout, sys.stderr = old_out, old_err
                sys.exit(f'FAILED: {name} (log: {log})')
            finally:
                sys.stdout, sys.stderr = old_out, old_err
        print(f'    {name} finished in {(time.time() - t0) / 60:.1f} min', flush=True)
    print('all requested stages finished')


class Tee:
    """Write everything to the stage log; echo whole progress lines to the terminal."""
    def __init__(self, fh, echo):
        self.fh, self.echo, self.line = fh, echo, ''

    def write(self, text):
        self.fh.write(text)
        self.fh.flush()
        self.line += text
        while '\n' in self.line:
            line, self.line = self.line.split('\n', 1)
            if line.startswith('[') or 'Error' in line:
                self.echo.write(line + '\n')
                self.echo.flush()
        return len(text)

    def flush(self):
        self.fh.flush()


def cmd_report(a):
    project, cfg = load_project(a.project)
    sys.path.insert(0, str(WORKFLOW))
    import report
    report.main(['--runs', str(project / 'runs'), '--out', str(project / 'results.html'), '--name', cfg['name']])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    n = sub.add_parser('new', help='create a project for a polymer')
    n.add_argument('name')
    n.add_argument('--monomer', help="repeat unit SMILES with two * ends, e.g. '*CC(C)(C(=O)OC)*'")
    n.add_argument('--terminator', default='*C', help="chain end group (default '*C', methyl)")
    n.add_argument('--dp', type=int, help='repeat units per chain')
    n.add_argument('--chains', type=int, help='number of chains')
    n.add_argument('--polyse', metavar='FILE',
                   help='a .polyse file with the chemistry (copolymers: blocks, random, alternating, explicit '
                        'sequences, distributions, mixtures) instead of --monomer/--dp/--chains')
    n.add_argument('--supercell', type=int, nargs=2, default=[3, 2], metavar=('NX', 'NY'),
                   help='silica supercell in x and y (default 3 2 = 120.9 x 82.9 A)')
    n.add_argument('--density', type=float, default=1.0, help='build density, g/cm3 (default 1.0)')
    n.add_argument('--forcefield', choices=list(FORCEFIELDS), default='pcff-iff-mod',
                   help='pcff-iff-mod (default): torsions across a near-linear Si-O-Si switched off, '
                        'stable at 1 fs with SHAKE; pcff-iff: the published model unchanged')
    n.add_argument('--seed', type=int, default=2026)
    n.add_argument('--dir', help='parent directory for the project (default ./projects)')
    n.add_argument('--test', action='store_true', help='every stage at 1/100 of its steps, to check a setup quickly')
    n.add_argument('--force', action='store_true')
    r = sub.add_parser('run', help='run stages (resumes at the first unfinished one)')
    r.add_argument('project')
    r.add_argument('--from', dest='from_', metavar='STAGE')
    r.add_argument('--to', metavar='STAGE')
    r.add_argument('--only', metavar='STAGE')
    s = sub.add_parser('status', help='which stages have finished')
    s.add_argument('project')
    p = sub.add_parser('report', help='write results.html (interactive 3D view and charts)')
    p.add_argument('project')
    a = ap.parse_args()
    {'new': cmd_new, 'run': cmd_run, 'status': cmd_status, 'report': cmd_report}[a.cmd](a)


if __name__ == '__main__':
    main()
