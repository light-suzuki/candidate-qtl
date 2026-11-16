"""Normal B continuation reuses original judgments, not another gene reviewer."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.agentic_route_b import prepare
from scripts.agentic_route_b_initial import collect_material, digest, summarize_initial
from scripts.agentic_route_b_triage import plan_light_followup
from tests.test_agentic_route_b_research import fixture, seal_prompts

ROOT = Path(__file__).resolve().parents[1]


class SharedResearchWorkflowTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.bundle, self.prompts = fixture()
        unit = self.bundle['units'][0]
        unit['reference_material_cards'] = []
        for gene in unit['genes']:
            gene.update(reference_gene='ref1', reference_species='Arabidopsis thaliana',
                        reference_function='A documented shared functional process.')
        self.source = self.root / 'route_b_semantic_bundle.json'
        self.prompts_path = self.root / 'route_b_semantic_prompts.jsonl'
        self.save()

    def save(self):
        self.prompts[0]['user'] = 'Phenotype: test\n\n' + json.dumps(self.bundle['units'][0], sort_keys=True)
        seal_prompts(self.bundle, self.prompts)
        self.source.write_text(json.dumps(self.bundle), encoding='utf8')
        self.prompts_path.write_text(''.join(json.dumps(p) + '\n' for p in self.prompts), encoding='utf8')
        self.sha = digest(self.source.read_bytes())

    def payload(self):
        material = collect_material(self.bundle)[0][0]['material_id']
        return {'case_id': 'fixture', 'bundle_sha256': self.sha, 'shared_assessments': [
            {'assessment_id': 'J1', 'material_id': material, 'member_gene_ids': ['g1', 'g2'],
             'state': 'assessed_claim', 'reason': 'Original functional judgment shared explicitly across these copies.',
             'assessed_fields': ['reference_function'],
             'research_planning': {'provisional_tier': 'T4', 'screen_ids': ['S1']}},
            {'assessment_id': 'J0', 'material_id': '', 'member_gene_ids': ['g3', 'g4'],
             'state': 'assessed_insufficient_material', 'reason': 'No functional context in the supplied input; external knowledge remains unknown.',
             'assessed_fields': []}],
            'research_screening': {'reference_screens': [{
                'screen_id': 'S1', 'reference_gene': 'ref1', 'reference_species': 'Arabidopsis thaliana',
                'member_gene_ids': ['g1', 'g2'], 'state': 'no_concrete_upper_lead',
                'queries': ['ref1 plant mutant phenotype'], 'retrieval_state': 'completed',
                'checked_species_scopes': ['focal_species', 'other_plants'], 'screened_source_ids': [],
                'reason': 'No concrete upper lead in this actual inspected scope.'}]}}

    def stage(self, payload=None, label='initial'):
        answers = self.root / (label + '.json')
        answers.write_text(json.dumps(payload or self.payload()), encoding='utf8')
        result = summarize_initial(str(self.source), self.sha, str(answers), str(self.root / label))
        summary = self.root / label / 'initial_assessment_summary.json'
        plan = plan_light_followup(str(summary), '', str(self.root / (label + '-plan')))
        return answers, result, Path(plan['plan']['path']), plan['plan']['sha256']

    def test_shared_judgments_and_embedded_screen_are_consumed_without_rewriting_genes(self):
        answers, summary, path, _ = self.stage()
        plan = json.loads(path.read_text())
        by_gene = {g['gene_id']: g for g in plan['genes']}
        self.assertEqual(2, summary['assessment_encoding']['shared_count'])
        self.assertEqual(4, summary['assessment_encoding']['expanded_count'])
        self.assertEqual(0, summary['genes_with_fields_without_read_receipt'])
        self.assertEqual({'g1', 'g2', 'g3', 'g4'}, set(by_gene))
        self.assertEqual('light_stop_eligible', by_gene['g1']['research_action'])
        self.assertEqual('not_evaluated', by_gene['g3']['research_action'])
        self.assertEqual(digest(answers.read_bytes()), plan['input']['screening']['sha256'])
        self.assertEqual(plan['input']['screening'], plan['input']['original_assessments'])

    def test_normal_prepare_consumes_plan_but_preserves_every_canonical_prompt_and_receipt(self):
        payload = self.payload()
        payload['assessments'] = [{'gene_id': 'g2', 'unit_id': 'u1', 'state': 'assessed_concrete_lead',
            'reason': 'An original copy-specific concrete developmental lead.',
            'used_material_pointers': ['/units/0/genes/1'], 'assessed_fields': ['reference_function']}]
        answers, _, plan, plan_sha = self.stage(payload)
        before = self.source.read_bytes(), self.prompts_path.read_bytes(), answers.read_bytes()
        initial, followup = self.root / 'first-view', self.root / 'followup-view'
        prepare(str(self.prompts_path), str(initial), research_bundle=str(self.source), research_bundle_sha256=self.sha)
        prepare(str(self.prompts_path), str(followup), research_bundle=str(self.source), research_bundle_sha256=self.sha,
                followup_plan=str(plan), followup_plan_sha256=plan_sha)
        manifest = json.loads((followup / 'research/manifest.json').read_text())
        self.assertEqual('followup', manifest['phase'])
        self.assertEqual(['g2'], manifest['research_action_groups']['targeted_followup'])
        self.assertEqual(['g1'], manifest['research_action_groups']['light_stop_eligible'])
        self.assertEqual([], manifest['initial_material']['pages'])
        self.assertEqual(2, manifest['initial_material']['genes_with_material'])
        self.assertEqual(['g3', 'g4'], manifest['initial_material']['genes_without_material'])
        self.assertEqual(0, manifest['initial_material']['delivered_gene_count'])
        self.assertIn('既読bundle全体', (followup / 'research/task.md').read_text())
        for file in initial.rglob('*'):
            if file.is_file() and 'research' not in file.relative_to(initial).parts and file.name != 'prepare_manifest.json':
                self.assertEqual(file.read_bytes(), (followup / file.relative_to(initial)).read_bytes())
        self.assertEqual(before, (self.source.read_bytes(), self.prompts_path.read_bytes(), answers.read_bytes()))

    def test_different_go_context_cannot_borrow_a_shared_judgment(self):
        payload = self.payload()
        self.bundle['units'][0]['genes'][1]['go_terms'] = 'different developmental context'
        self.save()
        payload['bundle_sha256'] = self.sha
        with self.assertRaisesRegex(ValueError, 'does not bind'):
            self.stage(payload)

    def test_absent_material_record_cannot_hide_available_function(self):
        payload = self.payload()
        payload['shared_assessments'][1]['member_gene_ids'] = ['g1']
        with self.assertRaisesRegex(ValueError, 'absent context'):
            self.stage(payload)

    def test_original_upper_observation_survives_a_lower_planning_override(self):
        payload = self.payload()
        payload['shared_assessments'][0]['research_planning']['provisional_tier'] = 'T1'
        payload['research_screening']['gene_assessments'] = [{
            'gene_id': 'g1', 'unit_id': 'u1', 'provisional_tier': 'T4',
            'tier_reason': 'A later provisional planning entry.', 'screen_ids': ['S1']}]
        _, _, path, _ = self.stage(payload)
        first = json.loads(path.read_text())['genes'][0]
        self.assertEqual('retain_upper_evidence', first['research_action'])
        self.assertEqual(['T1'], first['original_upper_tiers'])

    def test_tampered_original_ledger_is_rejected_before_followup_delivery(self):
        answers, _, path, sha = self.stage()
        answers.write_text(answers.read_text() + '\n')
        out = self.root / 'bad-view'
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            prepare(str(self.prompts_path), str(out), research_bundle=str(self.source), research_bundle_sha256=self.sha,
                    followup_plan=str(path), followup_plan_sha256=sha)
        self.assertFalse(out.exists())

    def test_unread_go_is_delivered_as_delta_without_repeating_read_context(self):
        self.bundle['units'][0]['genes'][1]['go_terms'] = 'specific developmental context'
        self.save()
        payload = self.payload()
        first = payload['shared_assessments'][0]
        first['member_gene_ids'] = ['g1']
        group = next(g for g in collect_material(self.bundle)[0] if 'go_terms' in g['context'])
        second = {**copy.deepcopy(first), 'assessment_id': 'J2', 'material_id': group['material_id'], 'member_gene_ids': ['g2']}
        payload['shared_assessments'].append(second)
        _, _, plan, sha = self.stage(payload)
        out = self.root / 'delta-view'
        prepare(str(self.prompts_path), str(out), research_bundle=str(self.source), research_bundle_sha256=self.sha,
                followup_plan=str(plan), followup_plan_sha256=sha)
        manifest = json.loads((out / 'research/manifest.json').read_text())
        material = manifest['initial_material']
        self.assertEqual(1, material['material_group_count'])
        self.assertEqual(2, material['full_material_group_count'])
        records = [g for p in material['pages'] for g in json.loads((out / p['file']).read_text())]
        self.assertEqual('specific developmental context', records[0]['context']['go_terms'])
        self.assertEqual(['g2'], records[0]['bindings'][0]['member_interval_genes'])

    def test_plan_hash_and_new_bundle_generation_are_not_interchangeable(self):
        _, _, plan, sha = self.stage()
        out = self.root / 'bad-hash'
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            prepare(str(self.prompts_path), str(out), research_bundle=str(self.source), research_bundle_sha256=self.sha,
                    followup_plan=str(plan), followup_plan_sha256='0' * 64)
        self.assertFalse(out.exists())
        self.bundle['units'][0]['genes'][0]['keywords'] = 'new frozen evidence'
        self.save()
        with self.assertRaises(ValueError):
            prepare(str(self.prompts_path), str(out), research_bundle=str(self.source), research_bundle_sha256=self.sha,
                    followup_plan=str(plan), followup_plan_sha256=sha)
        self.assertFalse(out.exists())

    def test_real_cli_uses_embedded_shared_ledger(self):
        answers, summary, _, _ = self.stage()
        out = self.root / 'cli-plan'
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/agentic_route_b.py'), 'plan-light-followup',
            '--summary', str(self.root / 'initial/initial_assessment_summary.json'), '--out', str(out)],
            text=True, capture_output=True)
        self.assertEqual(0, result.returncode, result.stderr)
        self.assertEqual(4, json.loads(result.stdout)['gene_count'])
        self.assertEqual(digest(answers.read_bytes()), summary['input_assessments']['sha256'])

    def test_thousand_gene_common_ledger_reduces_repeated_authored_output(self):
        unit = self.bundle['units'][0]
        genes = [f'g{i:04d}' for i in range(1000)]
        template = unit['genes'][0]
        unit['genes'] = [{**copy.deepcopy(template), 'gene_id': gid} for gid in genes]
        unit['member_interval_genes'] = genes
        self.bundle['units'] = [unit]
        self.bundle['semantic_unit_count'] = 1
        self.bundle['interval_gene_count'] = 1000
        self.save()
        mid = collect_material(self.bundle)[0][0]['material_id']
        payload = {'bundle_sha256': self.sha, 'shared_assessments': [{
            'assessment_id': 'J1000', 'material_id': mid, 'member_gene_ids': genes,
            'state': 'assessed_claim', 'reason': 'An explicitly shared original functional judgment. ' * 20,
            'assessed_fields': ['reference_function']}]}
        answers = self.root / 'compact.json'
        answers.write_text(json.dumps(payload), encoding='utf8')
        summary = summarize_initial(str(self.source), self.sha, str(answers), str(self.root / 'compact'))
        encoding = summary['assessment_encoding']
        self.assertEqual(1000, summary['gene_count'])
        self.assertEqual(1000, encoding['expanded_count'])
        self.assertLess(encoding['original_bytes'], encoding['expanded_bytes'] * 0.2)
        self.assertEqual(0, summary['genes_with_fields_without_read_receipt'])


if __name__ == '__main__':
    unittest.main()
