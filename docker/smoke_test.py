"""Validate installed PAVES, solvent assignment, OVITO widgets, and headless PNG export."""
import argparse
import importlib.util
import math
from pathlib import Path
import tempfile


def check_notebook(directory, system):
    """Execute the shipped viewer with a tiny real DCD and test its export buttons."""
    import mdtraj as md
    import nbformat
    from nbclient import NotebookClient
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue('ETH', chain)
    for i, mass in enumerate(system.masses):
        element = md.element.hydrogen if mass < 2 else md.element.carbon
        topology.add_atom(f'{element.symbol}{i}', element, residue)
    xyz = np.asarray(system.positions) / 10.
    trajectory = md.Trajectory(np.array([xyz, xyz + .01]), topology,
        unitcell_lengths=np.full((2, 3), 3.), unitcell_angles=np.full((2, 3), 90.))
    trajectory[0].save_pdb(str(directory / 'topology.pdb'))
    trajectory.save_dcd(str(directory / 'trajectory.dcd'))
    (directory / 'system.data').rename(directory / 'final.data')
    notebook = nbformat.read('/workspace/workflows/adhesion/ovito_view.ipynb', as_version=4)
    code = next(cell for cell in notebook.cells if cell.cell_type == 'code')
    code.source = ('import sys\nsys.path.insert(0, "/workspace/workflows/adhesion")\n'
                   + code.source.replace("RUN_DIR = 'projects/pmma_squeeze/runs/09_pull'",
                                         f'RUN_DIR = {str(directory)!r}'))
    notebook.cells.append(nbformat.v4.new_code_cell('''
assert nframes == 2
renderer_pick.value = 'Tachyon'
size_w.value, size_h.value = 160, 120
slider.value = 1
wrap_box.value = True
render_png()
images = list(d.glob('frame_*.png'))
assert images and images[0].stat().st_size > 100, status.value
render_movie()
movies = list(d.glob('movie_*.mp4'))
assert movies and b'ftyp' in movies[0].read_bytes()[:32], status.value
'''))
    client = NotebookClient(notebook, timeout=180, kernel_name='python3',
        resources={'metadata': {'path': str(directory)}})
    client.execute()
    assert any('application/vnd.jupyter.widget-view+json' in output.get('data', {})
               for cell in notebook.cells for output in cell.get('outputs', []))
    print('PASS: existing OVITO notebook: DCD, widget, slider, wrap, PNG and MP4')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gpu', action='store_true', help='also require an actual CUDA context')
    parser.add_argument('--notebook', action='store_true', help='execute the existing viewer with a tiny DCD')
    args = parser.parse_args()
    import ipykernel
    import openmm
    from openmm import unit
    import polypaves as p
    import ovito
    from ovito.io import import_file
    from ovito.vis import TachyonRenderer, Viewport

    assert int(ipykernel.__version__.split('.')[0]) < 7
    for name in ('api', 'sorption', 'reaction.rule'):
        assert importlib.util.find_spec(f'polypaves.{name}').origin.endswith('.so')
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        system = p.pack([p.Solvent('CC')], counts=[2],
            forcefield='/workspace/examples/forcefields/opls_alkane.ff',
            cell=(30., 30., 30.), trials=1, relax=0, postrelax=0)
        assert system.n_atoms == 16 and len(system.bonds) == 14
        system.write_lammps(directory)
        pipeline = import_file(str(directory / 'system.data'), atom_style='full')
        assert pipeline.compute().particles.count == 16
        pipeline.add_to_scene()
        view = Viewport(type=Viewport.Type.Perspective, camera_dir=(1, 1, -1))
        view.zoom_all()
        widget = ovito.gui.create_ipywidget(view)
        widget.refresh()
        assert widget.get_state()
        image = directory / 'smoke.png'
        view.render_image(filename=str(image), size=(160, 120), renderer=TachyonRenderer())
        assert image.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
        pipeline.remove_from_scene()
        if args.notebook:
            check_notebook(directory, system)
    water = p.make_solvent(p.Solvent('O', name='water'))
    assert water.atom_count == 3
    platform = openmm.Platform.getPlatformByName('CUDA' if args.gpu else 'CPU')
    context = openmm.Context(water.system, openmm.VerletIntegrator(.001 * unit.picoseconds), platform)
    context.setPositions(water.positions_nm * unit.nanometer)
    state = context.getState(getEnergy=True)
    assert math.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    print(f'PASS: PAVES {p.__version__}, OVITO {ovito.version_string}, widget, PNG, '
          f'Sage water, {platform.getName()} context ({state.getPotentialEnergy()})')


if __name__ == '__main__':
    main()
