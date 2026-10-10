"""Pack a SMILES solvent into a saved-charge network and run fixed-loading MD."""
import argparse
import json
from pathlib import Path
import sys


def pack(args):
    from polypaves import Solvent, pack_network
    if args.output.exists():
        raise FileExistsError(f'Output already exists: {args.output}; choose a new directory')
    result = pack_network(args.network, Solvent(args.smiles, name=args.name), args.molecules,
        seed=args.seed, cell_scale=args.cell_scale, min_distance_nm=args.min_distance_nm,
        attempts_per_molecule=args.attempts_per_molecule)
    result.save(args.output)
    print(json.dumps(dict(output=str(args.output), solvent_molecules=result.solvent_molecules,
        atoms=result.system.getNumParticles(), equilibrium_uptake_verified=False)))


def simulate(args):
    from polypaves.sorption import load_template
    from polypaves.sorption_workflow import FixedLoadingConfig, run_fixed_loading
    source = args.packed.resolve()
    data = json.loads((source / 'packed.structure.json').read_text())
    host = load_template(source, 'host')
    solvent = load_template(source, 'solvent')
    config = FixedLoadingConfig(platform=args.platform, temperature_K=args.temperature,
        pressure_atm=args.pressure_atm, timestep_fs=args.timestep_fs,
        nvt_ps=args.nvt_ps, npt_equilibration_ps=args.npt_equilibration_ps,
        production_ps=args.production_ps, sample_interval_ps=args.sample_interval_ps,
        frame_interval_ps=args.frame_interval_ps,
        minimization_iterations=args.minimization_iterations)
    # Initial loading is specified by a molecule count, never inferred as an
    # equilibrium absorption measurement. The shared engine handles checkpoints.
    result = run_fixed_loading(host, solvent, 0., args.output, config=config, seed=args.seed,
        resume=args.resume, count_override=data['solvent_molecules'], initial_state=data)
    print(json.dumps(dict(output=str(args.output), solvent_molecules=data['solvent_molecules'],
        loading_wt_percent=result['actual_loading_wt_percent'],
        fixed_loading_converged=result['fixed_loading_converged'], equilibrium_uptake_verified=False)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    preparation = commands.add_parser('pack', help='assign solvent and prepare exactly the requested count')
    preparation.add_argument('--network', type=Path, required=True)
    preparation.add_argument('--smiles', default='O', help='one neutral molecule; default: water')
    preparation.add_argument('--name', default='solvent')
    preparation.add_argument('--molecules', type=int, required=True)
    preparation.add_argument('--cell-scale', type=float, default=1.)
    preparation.add_argument('--min-distance-nm', type=float, default=.14)
    preparation.add_argument('--attempts-per-molecule', type=int, default=200)
    preparation.add_argument('--seed', type=int, default=31415)
    preparation.add_argument('--output', type=Path, required=True)
    preparation.set_defaults(action=pack)
    md = commands.add_parser('simulate', help='NVT / NPT equilibration / NPT sampling at fixed loading')
    md.add_argument('--packed', type=Path, required=True)
    md.add_argument('--output', type=Path, required=True)
    md.add_argument('--platform', choices=('CUDA', 'CPU', 'Reference'), default='CUDA')
    md.add_argument('--temperature', type=float, default=300.)
    md.add_argument('--pressure-atm', type=float, default=1.)
    md.add_argument('--timestep-fs', type=float, default=1.)
    md.add_argument('--nvt-ps', type=float, default=100.)
    md.add_argument('--npt-equilibration-ps', type=float, default=1000.)
    md.add_argument('--production-ps', type=float, default=5000.)
    md.add_argument('--sample-interval-ps', type=float, default=1.)
    md.add_argument('--frame-interval-ps', type=float, default=100.)
    md.add_argument('--minimization-iterations', type=int, default=3000)
    md.add_argument('--seed', type=int, default=31415)
    md.add_argument('--resume', action='store_true')
    md.set_defaults(action=simulate)
    args = parser.parse_args(argv)
    try:
        args.action(args)
        return 0
    except (ValueError, TypeError, RuntimeError, ImportError, OSError, KeyError) as exc:
        print(f'sorption: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
