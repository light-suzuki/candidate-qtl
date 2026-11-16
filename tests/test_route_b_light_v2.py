"""Opt-in science records retain observations and per-member completion."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from scripts import agentic_route_b_light as light
from tests.test_route_b_light import fixture, unit, digest, checks_for, MAIN_CACHE, PAPER, PAPER_TITLE, with_database_source, run_pipeline


def observation(key='o1', **fields):
    return {'observation_id': key, 'paper_id': PAPER, 'title': PAPER_TITLE,
            'reported_gene': 'HOM1', 'reported_gene_id': 'TAIR:AT1G12345',
            'reported_species': 'Arabidopsis thaliana', 'finding': 'Homolog loss changes leaf development.', **fields}


def decision(gene='g1', reference='REF-A', **fields):
    return {'interval_gene_id': gene, 'reference_gene': reference, 'state': 'assessed',
            'basis': 'paper_evidence', 'function': 'leaf development', 'reason': 'could explain fused leaflets',
            'paper_assessments': [{'observation_id': 'o1', 'gene_correspondence': 'homolog',
                                   'checklist': checks_for(phenotype_reported='yes', equivalent_target_phenotype='yes',
                                                           focal_species='no', identity_resolved='yes'),
                                   'reason': 'This homolog copy may retain the reported developmental function.'}], **fields}


class ScienceV2Tests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def convert(self, decisions, observations=None, bundle=None):
        path = self.root/'bundle.json'
        path.write_text(json.dumps(bundle or fixture()), encoding='utf-8')
        sha = digest(path.read_bytes())
        light.light_prepare(str(path), sha, str(self.root/'prepared'), output_contract='science-v2')
        index = self.root/'prepared/light/index.json'
        results = index.parent/'b01/results.jsonl'
        results.write_text(json.dumps({'contract': light.SCIENCE_CONTRACT,
                                      'paper_observations': observations if observations is not None else [observation()],
                                      'gene_decisions': decisions}), encoding='utf-8')
        cache = self.root/'cache.json'
        cache.write_text(json.dumps(MAIN_CACHE), encoding='utf-8')
        summary = light.light_convert(str(path), sha, str(index), digest(index.read_bytes()), str(self.root/'converted'),
                                      cache_arg=str(cache), offline=True)
        responses = json.loads((self.root/'converted/agent_responses.json').read_text(encoding='utf-8'))['responses']
        coverage = json.loads((self.root/'converted/member_coverage.json').read_text(encoding='utf-8'))
        return summary, responses, coverage

    def test_no_tier_or_literature_is_needed_and_original_identifiers_survive(self):
        _, responses, coverage = self.convert([decision()])
        claim = json.loads(responses[0]['raw_text'])['claims'][0]
        self.assertEqual('T2', claim['initial_tier'])
        report = claim['paper_reports_json'][0]
        self.assertEqual('TAIR:AT1G12345', report['reported_gene_id'])
        self.assertEqual(PAPER_TITLE, report['source_title'])
        self.assertIn('copy_assessment', report['research_annotations'])
        self.assertEqual('assessed', coverage['members'][0]['state'])
        receipt = json.loads((self.root/'converted/light_convert_receipt.json').read_text(encoding='utf-8'))
        self.assertEqual(digest((self.root/'converted/member_coverage.json').read_bytes()),
                         receipt['outputs']['member_coverage']['sha256'])
        self.assertEqual(coverage['scientific_completion_gate'],
                         light.verify_light_completion(str(self.root/'converted/light_convert_receipt.json'),
                                                       digest((self.root/'converted/light_convert_receipt.json').read_bytes())))

    def test_partial_unit_retains_good_claim_and_member_failures(self):
        bundle = fixture()
        bundle['units'] = [unit('u1', [('g1', 'REF-A'), ('g2', 'REF-B'), ('g3', 'REF-C')])]
        bad = decision('g3', 'REF-C')
        bad['paper_assessments'][0]['checklist'].pop('focal_species')
        summary, responses, coverage = self.convert([decision(), decision('g2', 'REF-B', state='not_assessed'), bad], bundle=bundle)
        self.assertEqual('completed', responses[0]['status'])
        self.assertEqual(1, len(json.loads(responses[0]['raw_text'])['claims']))
        self.assertEqual({'g1': 'assessed', 'g2': 'not_assessed', 'g3': 'invalid'},
                         {row['gene_id']: row['state'] for row in coverage['members']})
        self.assertEqual(['u1'], summary['member_coverage']['partial_units'])
        self.assertFalse(summary['scientific_completion_gate']['may_claim_all_genes_assessed'])

    def test_shared_observation_keeps_copy_checklist_reasons_separate(self):
        second = decision('g2', 'REF-B')
        second['paper_assessments'][0].update(gene_correspondence='unknown', reason='copy identity unresolved',
                                             counterevidence='copy differs in expression')
        _, responses, _ = self.convert([decision(), second])
        first, second = [json.loads(row['raw_text'])['claims'][0]['paper_reports_json'][0] for row in responses[:2]]
        self.assertEqual(first['source_id'], second['source_id'])
        self.assertNotEqual(first['gene_correspondence'], second['gene_correspondence'])
        self.assertNotEqual(first['research_annotations'], second['research_annotations'])

    def test_missing_title_gene_id_and_duplicate_observation_are_local_invalid(self):
        for field in ['title', 'reported_gene_id']:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temp:
                path = Path(temp)/'results'
                obs = observation()
                obs.pop(field)
                path.write_text(json.dumps({'contract': light.SCIENCE_CONTRACT, 'paper_observations': [obs],
                                            'gene_decisions': [decision(), decision('g2', 'REF-B', state='retrieval_failed')]}), encoding='utf-8')
                rows, issues = light._science_rows(path)
                self.assertEqual('invalid', rows['g1']['_science_state'])
                self.assertEqual('retrieval_failed', rows['g2']['_science_state'])
                self.assertTrue(issues)

    def test_database_function_coexists_with_paper_and_failure_does_not_invent_basis(self):
        bundle = with_database_source(fixture(), 'u1', 'g1', 'REF-A')
        _, responses, _ = self.convert([decision(basis='source_bound_function')], bundle=bundle)
        claims = json.loads(responses[0]['raw_text'])['claims']
        self.assertEqual(2, len(claims))
        db = next(claim for claim in claims if not claim.get('paper_reports_json'))
        self.assertEqual(['function'], db['evidence_types'])
        self.assertEqual('T4', db['initial_tier'])

    def test_canonical_freeze_retains_good_member_while_sealed_gate_remains_partial(self):
        bundle = fixture()
        bundle['units'] = [unit('u1', [('g1','REF-A'), ('g2','REF-B'), ('g3','REF-C')])]
        failed = observation('missing', paper_id='PMID:99999999', title='Unretrieved source')
        failed_decision = decision('g3', 'REF-C')
        failed_decision['paper_assessments'][0]['observation_id'] = 'missing'
        payload = {'contract': light.SCIENCE_CONTRACT, 'paper_observations': [observation(), failed],
                   'gene_decisions': [decision(), decision('g2', 'REF-B', state='not_assessed'), failed_decision]}
        built = run_pipeline(self.root, bundle, payload, MAIN_CACHE, 'u1\n', output_contract='science-v2')
        self.assertTrue(built['review'].loc['g1', 'relevance_tier'].startswith('T2'))
        self.assertNotIn('g2', built['review'].index)
        self.assertEqual({'assessed': 1, 'not_assessed': 1, 'invalid': 0, 'retrieval_failed': 1},
                         built['converted']['scientific_completion_gate']['member_counts'])
        self.assertFalse(built['converted']['scientific_completion_gate']['may_claim_all_genes_assessed'])
        coverage = json.loads((self.root/'converted/member_coverage.json').read_text(encoding='utf-8'))
        receipt = json.loads((self.root/'converted/light_convert_receipt.json').read_text(encoding='utf-8'))
        self.assertEqual(coverage['scientific_completion_gate'], receipt['scientific_completion_gate'])
        self.assertEqual(digest((self.root/'converted/agent_responses.json').read_bytes()),
                         receipt['outputs']['agent_responses']['sha256'])
        gate = light.verify_light_completion(built['converted']['receipt'], built['converted']['receipt_sha256'])
        self.assertEqual('partial', gate['state'])
        (self.root/'converted/member_coverage.json').write_text('{}', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'member_coverage hash mismatch'):
            light.verify_light_completion(built['converted']['receipt'], built['converted']['receipt_sha256'])

    def test_explicit_unknown_paper_is_retained_and_never_promoted_to_db_function(self):
        item = decision()
        item['paper_assessments'][0]['checklist'] = checks_for()
        _, responses, _ = self.convert([item], bundle=with_database_source(fixture(), 'u1', 'g1', 'REF-A'))
        claims = json.loads(responses[0]['raw_text'])['claims']
        self.assertEqual(1, len(claims))
        self.assertEqual(PAPER, claims[0]['paper_reports_json'][0]['source_id'])
        self.assertNotIn('reported_function', claims[0])

    def test_missing_database_provenance_cannot_create_paperless_function_claim(self):
        _, responses, coverage = self.convert([decision(basis='source_bound_function', paper_assessments=[])], observations=[])
        self.assertEqual('not_evaluated', responses[0]['status'])
        self.assertEqual('invalid', coverage['members'][0]['state'])
        self.assertEqual('no_frozen_database_source_for_this_reference', coverage['members'][0]['reason'])

    def test_uncached_ancillary_paper_does_not_erase_explicit_identity_judgment(self):
        item = decision(basis='identity_only')
        item['paper_assessments'][0]['checklist'] = checks_for(
            identity_resolved='yes', function_supports_hypothesis='no', annotation_supports_hypothesis='no',
            equivalent_target_phenotype='no', similar_cellular_or_developmental_phenotype='no')
        source = observation(paper_id='PMID:22258224', title='An uncached ancillary paper')
        _, responses, coverage = self.convert([item], observations=[source])
        self.assertEqual('completed', responses[0]['status'])
        self.assertEqual([], json.loads(responses[0]['raw_text'])['claims'])
        member = coverage['members'][0]
        self.assertEqual('assessed', member['state'])
        self.assertEqual('identity_only', member['reason'])
        self.assertEqual([{'paper_id': 'PMID:22258224', 'state': 'not_attempted'}], member['ancillary_paper_issues'])
        checks = json.loads((self.root/'converted/light_check.json').read_text(encoding='utf-8'))
        self.assertEqual('not_attempted', checks['paper_issues'][0]['state'])

    def test_uncached_paper_basis_and_upper_checklist_still_require_retrieval(self):
        for basis, upper in [('paper_evidence', False), ('identity_only', True)]:
            with self.subTest(basis=basis, upper=upper), tempfile.TemporaryDirectory() as temp:
                self.root = Path(temp)
                item = decision(basis=basis)
                if not upper:
                    item['paper_assessments'][0]['checklist'] = checks_for()
                source = observation(paper_id='PMID:22258224', title='An uncached paper')
                _, responses, coverage = self.convert([item], observations=[source])
                self.assertEqual('not_evaluated', responses[0]['status'])
                self.assertEqual('retrieval_failed', coverage['members'][0]['state'])
                self.assertEqual('partial', coverage['scientific_completion_gate']['state'])
                self.assertNotIn('ancillary_paper_issues', coverage['members'][0])


class AffinityTests(unittest.TestCase):
    def test_exact_reference_snapshot_affinity_preserves_roster_and_splits_large_group(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            bundle = fixture()
            bundle['units'] = [unit('u1', [('g1','REF-A')]), unit('u2', [('g2','REF-B')]),
                               unit('u3', [('g3','REF-A')]), unit('u4', [('g4','REF-A')])]
            different = copy.deepcopy(bundle['units'][2])
            different['reference_material_cards'][0]['mutant_phenotype'] = 'different tissue observation'
            self.assertNotEqual(light._affinity_key(different), light._affinity_key(bundle['units'][0]))
            path = root/'bundle.json'
            path.write_text(json.dumps(bundle), encoding='utf-8')
            sha = digest(path.read_bytes())
            light.light_prepare(str(path), sha, str(root/'old'), batch_units=2)
            light.light_prepare(str(path), sha, str(root/'new'), batch_units=2, packing='reference-affinity')
            old = json.loads((root/'old/light/index.json').read_text())
            new = json.loads((root/'new/light/index.json').read_text())
            self.assertEqual(['u1','u2'], old['batches'][0]['unit_ids'])
            self.assertEqual(['u1','u3'], new['batches'][0]['unit_ids'])
            self.assertEqual(['u4','u2'], new['batches'][1]['unit_ids'])
            self.assertEqual({'g1','g2','g3','g4'}, {gene for batch in new['batches'] for gene in batch['gene_ids']})
