"""Project-directory comparisons retain conditions and do not pool replicas across projects."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workflows.network_properties.directory_summary import collect_projects, enable_collection, MARKER
from workflows.network_properties.network_properties import DEFAULTS, report


class DirectorySummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def project(self,name,value,smoke=False):
        path=self.root/name
        path.mkdir(parents=True)
        cfg=dict(DEFAULTS,replicates=2,seed=100,mechanical=False,thermal_conductivity=False,
                 solvents=[],network='input.structure.json',input_sha256='saved-input',smoke=smoke)
        (path/'project.json').write_text(json.dumps(cfg))
        for i,offset in ((1,-1),(2,1)):
            base=path/f'replica_{i:02d}'
            base.mkdir()
            data=dict(smoke=smoke,npt=dict(estimates=dict(bulk_modulus_isothermal_GPa=value+offset),
                      sampling_checks_passed=True),cooling=dict(tg_K=400.,fit_checks_passed=True))
            (base/'dry_result.json').write_text(json.dumps(data))
        return path

    def rows(self,name):
        with (self.root/name).open(encoding='utf-8-sig',newline='') as f:
            return list(csv.DictReader(f))

    def test_nested_projects_with_same_basename_are_not_pooled(self):
        self.project('group_a/epoxy',3.)
        self.project('group_b/epoxy',9.)
        result=collect_projects(self.root)
        self.assertEqual(len(result['projects']),2)
        rows=[r for r in self.rows('projects_summary.csv') if r['property']=='npt.bulk_modulus_isothermal_GPa']
        self.assertEqual([r['project'] for r in rows],['group_a/epoxy','group_b/epoxy'])
        self.assertEqual([float(r['mean']) for r in rows],[3.,9.])
        self.assertTrue(all(r['unit']=='GPa' and r['n_available']=='2' for r in rows))
        replicas=[r for r in self.rows('replicas_summary.csv') if r['property']=='npt.bulk_modulus_isothermal_GPa']
        self.assertEqual(len(replicas),4)
        self.assertEqual([int(r['seed']) for r in replicas],[101,102,101,102])

    def test_missing_and_smoke_results_remain_distinct(self):
        p=self.project('smoke',3.,smoke=True)
        (p/'replica_02/dry_result.json').unlink()
        collect_projects(self.root)
        row=next(r for r in self.rows('projects_summary.csv') if r['property']=='npt.bulk_modulus_isothermal_GPa')
        self.assertEqual(row['n_available'],'1')
        self.assertEqual(row['n_expected'],'2')
        self.assertEqual(row['n_checks_passed'],'0')
        self.assertEqual(row['std'],'')
        replicas=[r for r in self.rows('replicas_summary.csv') if r['property']=='npt.bulk_modulus_isothermal_GPa']
        self.assertEqual([r['status'] for r in replicas],['smoke','pending'])
        self.assertEqual(replicas[1]['value'],'')

    def test_bad_directory_is_recorded_and_does_not_stop_others(self):
        self.project('good',3.)
        bad=self.root/'bad';bad.mkdir()
        (bad/'project.json').write_text('{')
        other=self.root/'other';other.mkdir()
        (other/'project.json').write_text(json.dumps({'monomer':'*CC*'}))
        result=collect_projects(self.root)
        self.assertEqual(len(result['projects']),1)
        self.assertEqual(result['errors'][0]['project'],'bad')
        self.assertIn('Directories requiring attention',(self.root/'projects_report.md').read_text())

    def test_project_report_automatically_refreshes_parent_collection(self):
        p=self.project('epoxy',3.)
        enable_collection(self.root)
        with patch('workflows.network_properties.plots.render',return_value=[]):
            report(p)
        self.assertTrue((self.root/'projects_summary.csv').exists())
        self.assertTrue((self.root/MARKER).exists())
        self.assertIn('epoxy',(self.root/'projects_report.md').read_text())
        before=(self.root/'projects_summary.csv').read_bytes()
        collect_projects(self.root)
        self.assertEqual((self.root/'projects_summary.csv').read_bytes(),before)


if __name__=='__main__':
    unittest.main()
