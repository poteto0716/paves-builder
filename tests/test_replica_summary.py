"""Cross-replica joins, physical units, missing values, and report-only export."""
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.network_properties.replica_summary import collect, write_tables, markdown
from workflows.network_properties.network_properties import report


class ReplicaSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.project = Path(self.temp.name)
        self.cfg = dict(replicates=3, smoke=False, mechanical=True, thermal_conductivity=True,
                        transport_axes=[0,1,2], solvents=[dict(name='water', smiles='O')])

    def put(self, path, data):
        path = self.project / path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))

    def dry(self, replica, modulus, checked):
        self.put(f'{replica}/dry_result.json', dict(smoke=False,
            npt=dict(estimates=dict(bulk_modulus_isothermal_GPa=modulus, alpha_volume_1_K=0.),
                     samples=5000, sampling_checks_passed=checked),
            cooling=dict(tg_K=400+10*modulus, fit_checks_passed=checked),
            nvt=dict(Cv_direct_J_kg_K=1000., effective_samples_energy=3000.,
                     sampling_checks_passed=checked)))

    def test_sparse_replicas_join_by_id_and_use_sample_sd(self):
        self.dry('replica_01',2.,True)
        self.dry('replica_03',4.,False)
        ids, values = collect(self.project,self.cfg)
        self.assertEqual(ids,['replica_01','replica_02','replica_03'])
        bulk = values['npt.bulk_modulus_isothermal_GPa']
        self.assertEqual(bulk['n_expected'],3)
        self.assertEqual(bulk['n_available'],2)
        self.assertEqual(bulk['mean'],3.)
        self.assertAlmostEqual(bulk['std'],math.sqrt(2))
        self.assertEqual(bulk['sem'],1.)
        self.assertEqual(bulk['replicas']['replica_02']['status'],'pending')
        self.assertEqual(bulk['replicas']['replica_03']['value'],4.)
        self.assertEqual(bulk['n_checks_passed'],1)
        self.assertEqual(bulk['checked_mean'],2.)
        self.assertIsNone(bulk['checked_std'])
        self.assertFalse(any('samples' in key or 'effective_samples' in key for key in values))
        self.assertEqual(values['npt.alpha_volume_1_K']['mean'],0.)

    def test_independent_stages_and_tensor_components(self):
        self.dry('replica_01',2.,True)
        tensor = [[float(10*i+j) for j in range(6)] for i in range(6)]
        self.put('replica_03/mechanical/result.json',dict(
            isotropic=dict(stable=True, young_Hill_GPa=8.),
            symmetric_stiffness_GPa=tensor, sampling_checks_passed=False))
        self.put('replica_02/thermal_conductivity/result.json',dict(
            conductivity_direction_mean_W_m_K=.2,thermal_diffusivity_m2_s=1e-7,
            sampling_checks_passed=False,
            axes=[dict(axis=1,conductivity_W_m_K=.3,sampling_checks_passed=False)]))
        _, values = collect(self.project,self.cfg)
        young = values['mechanical.young_Hill_GPa']
        self.assertEqual(young['replicas']['replica_03']['value'],8.)
        self.assertIsNone(young['replicas']['replica_01']['value'])
        self.assertIsNone(young['std'])
        self.assertEqual(values['mechanical.C_xx_yy_GPa']['mean'],1.)
        self.assertEqual(values['transport.conductivity_y_W_m_K']['mean'],.3)
        self.assertIsNone(values['transport.conductivity_x_W_m_K']['mean'])

    def test_partial_solvent_results_and_undefined_diffusion_are_preserved(self):
        self.put('replica_01/solvents/water/infinite_dilution_free_energy/result.json',dict(
            coupling_free_energy_kJ_mol=-10.,chemical_potential_relative_kJ_mol=-9.,
            convergence_verified=True,standard_error_kJ_mol=.2))
        self.put('replica_03/solvents/water/result.json',dict(
            infinite_dilution_free_energy=dict(coupling_free_energy_kJ_mol=-14.,
                chemical_potential_relative_kJ_mol=-13.,convergence_verified=False),
            fixed_loading=dict(fixed_loading_converged=True,diffusion=dict(
                diffusion_coefficient_m2_s=None,candidate_diffusion_coefficient_m2_s=999.,
                converged=False))))
        _, values = collect(self.project,self.cfg)
        free_energy = values['solvents.water.coupling_free_energy_kJ_mol']
        self.assertEqual(free_energy['unit'],'kJ/mol')
        self.assertEqual(free_energy['mean'],-12.)
        self.assertEqual(free_energy['checked_mean'],-10.)
        diffusion = values['solvents.water.diffusion_coefficient_m2_s']
        self.assertEqual(diffusion['n_available'],0)
        self.assertIsNone(diffusion['mean'])
        self.assertEqual(diffusion['replicas']['replica_03']['status'],'unavailable')
        self.assertEqual(diffusion['replicas']['replica_01']['status'],'pending')

    def test_nonfinite_is_excluded_and_smoke_never_passes_checks(self):
        self.dry('replica_01',float('nan'),True)
        _, values = collect(self.project,self.cfg)
        self.assertEqual(values['npt.bulk_modulus_isothermal_GPa']['n_available'],0)
        self.dry('replica_01',2.,True)
        self.cfg['smoke'] = True
        _, values = collect(self.project,self.cfg)
        bulk = values['npt.bulk_modulus_isothermal_GPa']
        self.assertEqual(bulk['n_checks_passed'],0)
        self.assertIsNone(bulk['checked_mean'])
        self.assertEqual(bulk['replicas']['replica_01']['status'],'smoke')

    def test_300K_void_fraction_and_crosslink_target_are_in_replica_summary(self):
        self.cfg['input_mode']='epoxy'
        for replica, fraction, reached in [('replica_01',.1,True),('replica_03',.3,False)]:
            self.dry(replica,2.,True)
            path=self.project/replica/'dry_result.json'
            data=json.loads(path.read_text())
            data['structure']={'point_probe_void_fraction':fraction}
            path.write_text(json.dumps(data))
            self.put(f'{replica}/network_origin.json',dict(conversion=.995 if reached else .8,
                     target_reached=reached,epoxide_consumption_including_activation=1. if reached else .9))
        _,values=collect(self.project,self.cfg)
        void=values['structure.point_probe_void_fraction']
        self.assertAlmostEqual(void['mean'],.2)
        self.assertEqual(void['unit'],'1')
        self.assertEqual(void['n_checks_passed'],1)
        self.assertEqual(void['checked_mean'],.1)
        self.assertEqual(values['npt.bulk_modulus_isothermal_GPa']['n_checks_passed'],1)
        self.assertEqual(values['crosslinking.conversion']['replicas']['replica_03']['value'],.8)

    def test_csv_exports_units_status_and_blank_missing_values(self):
        self.dry('replica_01',2.,True)
        ids, values = collect(self.project,self.cfg)
        write_tables(self.project,ids,values)
        def read(name):
            with (self.project/name).open(encoding='utf-8-sig',newline='') as fh:
                return {row['property']:row for row in csv.DictReader(fh)}
        row = read('replica_properties.csv')['npt.bulk_modulus_isothermal_GPa']
        self.assertEqual(row['unit'],'GPa')
        self.assertEqual(row['replica_01'],'2.0')
        self.assertEqual(row['replica_02'],'')
        self.assertEqual(row['replica_02_status'],'pending')
        row = read('property_summary.csv')['npt.bulk_modulus_isothermal_GPa']
        self.assertEqual(row['n_available'],'1')
        self.assertEqual(row['std'],'')
        self.assertEqual(row['checked_mean'],'2.0')
        table = '\n'.join(markdown(ids,values))
        self.assertIn('replica_01 | replica_02 | replica_03',table)
        for line in table.splitlines():
            if line.startswith('|'):
                self.assertEqual(line.count('|'),len(ids)+8)

    def test_report_refresh_requires_no_md_and_is_repeatable(self):
        self.dry('replica_01',2.,True)
        self.put('project.json',self.cfg)
        with patch('workflows.network_properties.plots.render',return_value=[]):
            report(self.project)
            before={name:(self.project/name).read_bytes() for name in (
                'summary.json','replica_properties.csv','property_summary.csv','report.md')}
            report(self.project)
        for name,content in before.items():
            self.assertEqual((self.project/name).read_bytes(),content)
        summary=json.loads(before['summary.json'])
        self.assertFalse(summary['dry_sampling_checks_passed'])
        self.assertFalse(summary['tg_fit_checks_passed'])
        self.assertIn('property_summary',summary)
        self.assertEqual(summary['property_summary']['npt.bulk_modulus_isothermal_GPa']['mean'],2.)


if __name__ == '__main__':
    unittest.main()
