"""Local checks/corrections do not infer science or rewrite old generations."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts import agentic_route_b_light as light
from scripts import agentic_route_b_light_local as local
from tests.test_route_b_light import fixture, digest, gene_row, paper_row, PAPER, PAPER_TITLE, checks_for, MAIN_CACHE
from tests.test_route_b_light_v2 import decision, observation


class LocalLightTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def setup_results(self, science=True, payload=None):
        bundle = self.root/'bundle.json'
        bundle.write_text(json.dumps(fixture()), encoding='utf-8')
        sha = digest(bundle.read_bytes())
        light.light_prepare(str(bundle), sha, str(self.root/'prepared'), output_contract='science-v2' if science else 'legacy-v1')
        index = self.root/'prepared/light/index.json'
        self.args = (str(bundle), sha, str(index), digest(index.read_bytes()), 'b01')
        self.results = self.root/'results.jsonl'
        self.science = science
        if payload is None:
            payload = {'contract': light.SCIENCE_CONTRACT, 'paper_observations': [observation()],
                       'gene_decisions': [decision(), decision('g2','REF-B', state='not_assessed')]}
        self.payload = payload
        self.results.write_text(json.dumps(payload) if science else ''.join(json.dumps(row)+'\n' for row in payload), encoding='utf-8')

    def preflight(self, name='preflight.json'):
        result = local.light_preflight(*self.args, str(self.results), str(self.root/name))
        return result, json.loads(Path(result['packet']).read_text(encoding='utf-8'))

    def patch(self, operations, **overrides):
        packet = local._packet(*self.args, str(self.results))[0]
        value = {'contract': local.PATCH_CONTRACT, 'base': packet['base'], 'operations': operations, **overrides}
        path = self.root/'patch.json'
        path.write_text(json.dumps(value), encoding='utf-8')
        return path

    def apply(self, operations, **overrides):
        path = self.patch(operations, **overrides)
        return local.light_apply(*self.args, str(self.results), str(path), digest(path.read_bytes()),
                                 str(self.root/'new-results.jsonl'), str(self.root/'apply-receipt.json'))

    def operation(self, *, kind='gene_decision', key='g1', path=None, expected=None, value=None, **fields):
        return {'op': 'set', 'record': {'kind': kind, 'id': key}, 'path': path or ['reason'],
                'expected': expected if expected is not None else 'could explain fused leaflets',
                'value': value if value is not None else 'Original researcher clarified the meristem mechanism.', **fields}

    def test_local_packet_preserves_good_member_and_explicit_unassessed(self):
        self.setup_results()
        summary, packet = self.preflight()
        self.assertEqual('partial', summary['scientific_completion'])
        self.assertEqual('assessed', packet['member_states']['g1'])
        self.assertEqual('not_assessed', packet['member_states']['g2'])
        self.assertEqual(digest(self.results.read_bytes()), packet['base']['results_sha256'])
        self.assertNotIn('title', json.dumps(packet['findings']))
        self.assertIn('no metadata retrieval', packet['scope'])

    def test_missing_observation_fields_are_targeted_without_defaults(self):
        self.setup_results()
        self.payload['paper_observations'][0].pop('reported_gene_id')
        self.results.write_text(json.dumps(self.payload), encoding='utf-8')
        summary, packet = self.preflight()
        self.assertEqual('needs_correction', summary['state'])
        self.assertTrue(any(row.get('observation_id') == 'o1' and row['field'] == 'reported_gene_id' for row in packet['findings']))
        self.assertEqual('invalid', packet['member_states']['g1'])
        self.assertNotIn('reported_gene_id', json.loads(self.results.read_text())['paper_observations'][0])

    def test_nine_answers_reference_and_unknown_id_are_checked(self):
        self.setup_results()
        self.payload['gene_decisions'][0]['paper_assessments'][0]['checklist'].pop('focal_species')
        self.payload['gene_decisions'][0]['reference_gene'] = 'WRONG'
        self.payload['gene_decisions'].append(decision('unexpected','REF-A'))
        self.results.write_text(json.dumps(self.payload), encoding='utf-8')
        _, packet = self.preflight()
        self.assertTrue(any('all nine' in row['issue'] for row in packet['findings']))
        self.assertTrue(any(row['issue'] == 'reference_not_in_frozen_endpoints' for row in packet['findings']))
        self.assertTrue(any(row['issue'] == 'gene_not_in_batch' for row in packet['findings']))

    def test_missing_copy_answer_has_exact_path_and_can_be_corrected_explicitly(self):
        self.setup_results()
        answers = self.payload['gene_decisions'][0]['paper_assessments'][0]['checklist']
        answers.pop('focal_species')
        answers['identity_resolved'] = 'probably'
        self.results.write_text(json.dumps(self.payload), encoding='utf-8')
        old = self.results.read_bytes()
        _, packet = self.preflight()
        for key in ['focal_species', 'identity_resolved']:
            self.assertTrue(any(row['gene_id'] == 'g1' and row.get('observation_id') == 'o1'
                                and row['field'] == f'paper_assessments[0].checklist.{key}' for row in packet['findings']))
        result = self.apply([
            self.operation(path=['paper_assessments',0,'checklist','focal_species'], expected=local.MISSING, value='no'),
            self.operation(path=['paper_assessments',0,'checklist','identity_resolved'], expected='probably', value='yes')])
        receipt = json.loads(Path(result['receipt']).read_text())
        self.assertEqual('assessed', receipt['postpreflight']['member_states']['g1'])
        self.assertFalse(any(row['gene_id'] == 'g1' for row in receipt['postpreflight']['findings']))
        self.assertEqual(old, self.results.read_bytes())

    def test_legacy_t6_ancillary_paper_does_not_require_new_schema_completion(self):
        self.setup_results(False, [gene_row('g1','REF-A','T6',[{'paper_id':'PMID:22258224'}])])
        _, packet = self.preflight()
        self.assertEqual('assessed', packet['member_states']['g1'])
        self.assertFalse(any(row['gene_id'] == 'g1' for row in packet['findings']))

    def test_legacy_source_record_correction_adds_missing_id_and_keeps_other_gene(self):
        paper = paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes'))
        self.setup_results(False, [gene_row('g1','REF-A','T4',[paper]), gene_row('g2','REF-B','U')])
        original = self.results.read_bytes()
        result = self.apply([self.operation(kind='legacy_gene', path=['papers',0,'reported_gene_id'],
                                            expected=local.MISSING, value='TAIR:AT1G12345')])
        changed = [json.loads(line) for line in Path(result['results']).read_text().splitlines()]
        self.assertEqual('TAIR:AT1G12345', changed[0]['papers'][0]['reported_gene_id'])
        self.assertEqual(self.payload[1], changed[1])
        self.assertEqual(original, self.results.read_bytes())
        receipt = json.loads(Path(result['receipt']).read_text())
        self.assertIn('legacy_gene:g2', receipt['unchanged_record_hashes'])

    def test_explicit_scientific_leaf_correction_preserves_non_targets(self):
        self.setup_results()
        old = self.results.read_bytes()
        result = self.apply([self.operation(path=['paper_assessments',0,'checklist','shared_specific_process'],
                                            expected='unknown', value='yes')])
        payload = json.loads(Path(result['results']).read_text())
        self.assertEqual('yes', payload['gene_decisions'][0]['paper_assessments'][0]['checklist']['shared_specific_process'])
        self.assertEqual(self.payload['gene_decisions'][1], payload['gene_decisions'][1])
        self.assertEqual(self.payload['paper_observations'], payload['paper_observations'])
        self.assertEqual(old, self.results.read_bytes())
        receipt = json.loads(Path(result['receipt']).read_text())
        self.assertTrue(receipt['original_results_unchanged'])
        self.assertIn('postpreflight', receipt)

    def test_stale_base_patch_hash_and_material_tampering_are_rejected(self):
        self.setup_results()
        patch = self.patch([self.operation()], base={'results_sha256':'stale'})
        with self.assertRaisesRegex(ValueError, 'base hash'):
            local.light_apply(*self.args, str(self.results), str(patch), digest(patch.read_bytes()), str(self.root/'new'), str(self.root/'receipt'))
        with self.assertRaisesRegex(ValueError, 'patch SHA'):
            local.light_apply(*self.args, str(self.results), str(patch), 'stale', str(self.root/'new'), str(self.root/'receipt'))
        (self.root/'prepared/light/b01/material.md').write_text('tampered', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            self.preflight()

    def test_conflict_delete_unknown_ids_roster_and_whole_container_updates_rejected(self):
        cases = [self.operation(expected='stale'), self.operation(key='unknown'),
                 self.operation(path=['interval_gene_id'], expected='g1', value='g2'),
                 self.operation(path=['paper_assessments'], expected=[], value=[]),
                 self.operation(op='delete'), self.operation(path=['reference_gene'], expected='REF-A', value='WRONG')]
        for operation in cases:
            with self.subTest(operation=operation), tempfile.TemporaryDirectory() as temp:
                self.root = Path(temp)
                self.setup_results()
                with self.assertRaises(ValueError):
                    self.apply([operation])
                self.assertFalse((self.root/'new-results.jsonl').exists())

    def test_duplicate_operations_are_rejected(self):
        self.setup_results()
        with self.assertRaisesRegex(ValueError, 'conflicting'):
            self.apply([self.operation(), self.operation(expected='Original researcher clarified the meristem mechanism.')])

    def test_shared_observation_requires_explicit_all_copy_acknowledgement(self):
        self.setup_results()
        self.payload['gene_decisions'][1] = decision('g2','REF-B')
        self.results.write_text(json.dumps(self.payload), encoding='utf-8')
        operation = self.operation(kind='paper_observation', key='o1', path=['title'], expected=PAPER_TITLE, value='Corrected original title')
        with self.assertRaisesRegex(ValueError, 'affected_gene_ids'):
            self.apply([operation])
        operation['affected_gene_ids'] = ['g1','g2']
        result = self.apply([operation])
        payload = json.loads(Path(result['results']).read_text())
        self.assertEqual('Corrected original title', payload['paper_observations'][0]['title'])
        self.assertEqual(self.payload['gene_decisions'], payload['gene_decisions'])

    def test_same_paper_different_gene_observations_are_never_merged(self):
        self.setup_results()
        self.payload['paper_observations'].append(observation('o2', reported_gene='OTHER', reported_gene_id='TAIR:AT2G12345'))
        self.results.write_text(json.dumps(self.payload), encoding='utf-8')
        result = self.apply([self.operation(kind='paper_observation', key='o1', path=['finding'],
                                            expected=observation()['finding'], value='Original observation clarified.',
                                            affected_gene_ids=['g1','g2'])])
        payload = json.loads(Path(result['results']).read_text())
        self.assertEqual(2, len(payload['paper_observations']))
        self.assertEqual(self.payload['paper_observations'][1], payload['paper_observations'][1])

    def test_existing_outputs_are_never_overwritten(self):
        self.setup_results()
        self.preflight()
        with self.assertRaisesRegex(ValueError, 'must be new'):
            self.preflight()
        (self.root/'new-results.jsonl').write_text('preserve', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'new paths'):
            self.apply([self.operation()])
        self.assertEqual('preserve', (self.root/'new-results.jsonl').read_text())

    def test_absolute_existing_cli_works_outside_repository_for_preflight_and_apply(self):
        self.setup_results()
        wrapper = Path(__file__).resolve().parents[1]/'scripts/agentic_route_b.py'
        bound = ['--bundle', self.args[0], '--bundle-sha256', self.args[1], '--light-index', self.args[2],
                 '--light-index-sha256', self.args[3], '--batch-id', self.args[4], '--results', str(self.results)]
        first = subprocess.run([sys.executable, str(wrapper), 'light-preflight', *bound, '--out', str(self.root/'cli-packet.json')],
                               cwd=self.root, capture_output=True, text=True)
        self.assertEqual(0, first.returncode, first.stderr)
        self.assertIn('scientific_completion', json.loads(first.stdout))
        patch = self.patch([self.operation()])
        applied = subprocess.run([sys.executable, str(wrapper), 'light-apply', *bound, '--patch', str(patch),
                                  '--patch-sha256', digest(patch.read_bytes()), '--out-results', str(self.root/'cli-results.jsonl'),
                                  '--receipt', str(self.root/'cli-receipt.json')], cwd=self.root, capture_output=True, text=True)
        self.assertEqual(0, applied.returncode, applied.stderr)
        self.assertTrue(json.loads((self.root/'cli-receipt.json').read_text())['original_results_unchanged'])


class BatchLocalLightTests(unittest.TestCase):
    setUp = LocalLightTests.setUp
    setup_results = LocalLightTests.setup_results
    operation = LocalLightTests.operation
    def setup_batch(self):
        bundle = self.root/'bundle.json'
        value = fixture()
        value['units'] = value['units'][:2]
        bundle.write_text(json.dumps(value), encoding='utf-8')
        sha = digest(bundle.read_bytes())
        light.light_prepare(str(bundle), sha, str(self.root/'prepared'), batch_units=1, output_contract='science-v2')
        index = self.root/'prepared/light/index.json'
        self.batch_results = []
        for number, (gene, reference) in enumerate([('g1','REF-A'), ('g2','REF-B')], 1):
            path = self.root/f'results{number}.json'
            path.write_text(json.dumps({'contract':light.SCIENCE_CONTRACT, 'paper_observations':[observation()],
                                        'gene_decisions':[decision(gene, reference)]}), encoding='utf-8')
            self.batch_results.append(path)
        self.source = self.root/'original-source.json'
        self.source.write_text(json.dumps({'title':'Original unmodified title', 'unknown':'unknown',
                                           'unselected':'DO NOT DUMP THIS SOURCE'}), encoding='utf-8')
        # Deliberately noncanonical insertion order: binding keys are named.
        self.request_value = {'contract':local.REVIEW_REQUEST_CONTRACT,
                              'binding':{'light_index_sha256':digest(index.read_bytes()), 'bundle_sha256':sha,
                                         'light_index':str(index.relative_to(self.root)), 'bundle':'bundle.json'},
                              'batches':[{'batch_id':f'b{n:02d}', 'results':path.name, 'sources':[]}
                                         for n,path in enumerate(self.batch_results, 1)]}
        self.request_value['batches'][0]['sources'] = [{'source':self.source.name,
                                                      'source_sha256':digest(self.source.read_bytes()),
                                                      'pointers':['/title','/unknown']}]
        self.request = self.root/'request.json'
        self.save_request()

    def save_request(self):
        self.request.write_text(json.dumps(self.request_value), encoding='utf-8')
        self.request_sha = digest(self.request.read_bytes())

    def review_batch(self):
        self.review_result = local.light_review_batch(str(self.request), self.request_sha, str(self.root/'review'))
        return self.review_result

    def apply_batch(self, decisions):
        path = self.root/'decisions.json'
        path.write_text(json.dumps({'contract':local.DECISIONS_CONTRACT, 'batches':decisions}), encoding='utf-8')
        return local.light_apply_batch(self.review_result['review'], self.review_result['review_sha256'],
                                       str(path), digest(path.read_bytes()), str(self.root/'applied'))

    def test_batch_review_returns_no_roster_or_full_receipt_and_only_selected_source_values(self):
        self.setup_batch()
        result = self.review_batch()
        self.assertEqual('complete', result['scientific_completion'])
        self.assertNotIn('corrections', result)
        returned = json.dumps(result)
        for text in ['member_states', 'paper_assessments', 'DO NOT DUMP', 'roster_sha256']:
            self.assertNotIn(text, returned)
        self.assertEqual([{'pointer':'/title','value':'Original unmodified title'},
                          {'pointer':'/unknown','value':'unknown'}], result['sources'][0]['selections'])
        review = json.loads(Path(result['review']).read_text())
        packet = json.loads(Path(review['batches'][0]['preflight']['path']).read_text())
        self.assertEqual({'g1':'assessed'}, packet['member_states'])
        self.assertEqual(digest(Path(result['review']).read_bytes()), result['review_sha256'])

    def test_batch_review_returns_exact_missing_answers_without_filling_them(self):
        self.setup_batch()
        value = json.loads(self.batch_results[0].read_text())
        del value['gene_decisions'][0]['paper_assessments'][0]['checklist']['focal_species']
        self.batch_results[0].write_text(json.dumps(value), encoding='utf-8')
        result = self.review_batch()
        self.assertEqual(['b01'], [row['batch_id'] for row in result['corrections']])
        self.assertTrue(any(row['field'].endswith('.focal_species') for row in result['corrections'][0]['findings']))
        self.assertNotIn('focal_species', json.loads(self.batch_results[0].read_text())['gene_decisions'][0]['paper_assessments'][0]['checklist'])

    def test_batch_source_output_is_fail_closed_without_clipping(self):
        self.setup_batch()
        self.request_value['max_chars'] = 1
        self.save_request()
        result = self.review_batch()
        self.assertEqual('narrow_selection_required', result['sources']['state'])
        self.assertNotIn('Original unmodified title', json.dumps(result))

    def test_large_findings_remain_available_by_exact_hash_without_stdout_dump(self):
        self.setup_batch()
        value = json.loads(self.batch_results[0].read_text())
        del value['gene_decisions'][0]['paper_assessments'][0]['checklist']['focal_species']
        self.batch_results[0].write_text(json.dumps(value), encoding='utf-8')
        self.request_value['max_chars'] = 1
        self.save_request()
        result = self.review_batch()
        self.assertEqual('narrow_selection_required', result['corrections']['state'])
        self.assertNotIn('missing_or_invalid_original_answer', json.dumps(result))
        artifact = result['corrections_artifact']
        self.assertEqual(digest(Path(artifact['path']).read_bytes()), artifact['sha256'])
        selected = local.material.source_read(artifact['path'], artifact['sha256'], pointers=['/batches/0'])
        findings = selected['selections'][0]['value']['findings']
        self.assertEqual(result['issue_count'], len(findings))
        self.assertTrue(any(row['field'].endswith('.focal_species') for row in findings))

    def test_batch_apply_preserves_untargeted_records_and_explicit_unchanged_batch(self):
        self.setup_batch()
        originals = [path.read_bytes() for path in self.batch_results]
        self.review_batch()
        result = self.apply_batch([{'batch_id':'b01','operations':[self.operation()]},
                                   {'batch_id':'b02','unchanged':True}])
        self.assertEqual((1,1), (result['changed_batch_count'],result['unchanged_batch_count']))
        self.assertNotIn('batches', result)
        receipt = json.loads(Path(result['receipt']).read_text())
        first, second = receipt['batches']
        changed = json.loads(Path(first['results']).read_text())
        original = json.loads(originals[0])
        self.assertEqual(original['paper_observations'], changed['paper_observations'])
        expected = original['gene_decisions'][0].copy()
        expected['reason'] = self.operation()['value']
        self.assertEqual(expected, changed['gene_decisions'][0])
        self.assertEqual(str(self.batch_results[1].resolve()), second['results'])
        self.assertEqual(originals, [path.read_bytes() for path in self.batch_results])
        self.assertTrue(Path(first['receipt']).exists())

    def test_batch_rejects_tampered_request_review_results_and_source_before_output(self):
        for target in ['request','review','results','source','preflight']:
            with self.subTest(target=target), tempfile.TemporaryDirectory() as temp:
                self.root = Path(temp)
                self.setup_batch()
                self.review_batch()
                review = json.loads(Path(self.review_result['review']).read_text())
                path = {'request':self.request, 'review':Path(self.review_result['review']),
                        'results':self.batch_results[0], 'source':self.source,
                        'preflight':Path(review['batches'][0]['preflight']['path'])}[target]
                path.write_text(path.read_text()+' ', encoding='utf-8')
                with self.assertRaisesRegex(ValueError,'mismatch'):
                    self.apply_batch([{'batch_id':'b01','unchanged':True},{'batch_id':'b02','unchanged':True}])
                self.assertFalse((self.root/'applied').exists())

    def test_batch_missing_foreign_duplicate_and_implicit_noop_decisions_are_rejected(self):
        cases = [[{'batch_id':'b01','unchanged':True}],
                 [{'batch_id':'b01','unchanged':True},{'batch_id':'foreign','unchanged':True}],
                 [{'batch_id':'b01','unchanged':True},{'batch_id':'b01','unchanged':True}],
                 [{'batch_id':'b01','operations':None},{'batch_id':'b02','unchanged':True}],
                 [{'batch_id':'b01','operations':[]},{'batch_id':'b02','unchanged':True}]]
        for decisions in cases:
            with self.subTest(decisions=decisions), tempfile.TemporaryDirectory() as temp:
                self.root = Path(temp)
                self.setup_batch()
                self.review_batch()
                with self.assertRaises(ValueError):
                    self.apply_batch(decisions)
                self.assertFalse((self.root/'applied').exists())

    def test_batch_validates_later_operations_before_writing_earlier_batches(self):
        self.setup_batch()
        self.review_batch()
        with self.assertRaisesRegex(ValueError,'expected value conflict'):
            self.apply_batch([{'batch_id':'b01','operations':[self.operation()]},
                              {'batch_id':'b02','operations':[self.operation(key='g2',expected='stale')]}])
        self.assertFalse((self.root/'applied').exists())

    def test_batch_shared_observation_still_requires_all_copy_acknowledgement(self):
        self.setup_results(payload={'contract':light.SCIENCE_CONTRACT, 'paper_observations':[observation()],
                                    'gene_decisions':[decision(),decision('g2','REF-B')]})
        self.request = self.root/'request.json'
        self.request.write_text(json.dumps({'contract':local.REVIEW_REQUEST_CONTRACT,
            'binding':dict(zip(('bundle','bundle_sha256','light_index','light_index_sha256'),self.args[:4])),
            'batches':[{'batch_id':'b01','results':str(self.results)}]}), encoding='utf-8')
        self.request_sha = digest(self.request.read_bytes())
        self.review_batch()
        op = self.operation(kind='paper_observation',key='o1',path=['title'],expected=PAPER_TITLE,value='Recorded correction')
        with self.assertRaisesRegex(ValueError,'affected_gene_ids'):
            self.apply_batch([{'batch_id':'b01','operations':[op]}])
        self.assertFalse((self.root/'applied').exists())
        op['affected_gene_ids'] = ['g1','g2']
        result = self.apply_batch([{'batch_id':'b01','operations':[op]}])
        self.assertEqual(1,result['changed_batch_count'])

    def test_batch_review_rejects_tampered_request_and_foreign_batch_without_output(self):
        self.setup_batch()
        self.request.write_text(self.request.read_text()+' ', encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
            self.review_batch()
        self.request_value['batches'][1]['batch_id'] = 'foreign'
        self.save_request()
        with self.assertRaisesRegex(ValueError,'unknown light batch_id'):
            self.review_batch()
        self.assertFalse((self.root/'review').exists())


class PostprocessLightTests(unittest.TestCase):
    """One machine invocation composes the existing helpers and keeps old freezes."""
    setUp = BatchLocalLightTests.setUp
    setup_batch = BatchLocalLightTests.setup_batch
    save_request = BatchLocalLightTests.save_request
    operation = BatchLocalLightTests.operation

    def postprocess(self, **kwargs):
        cache = self.root/'cache.json'
        cache.write_text(json.dumps(MAIN_CACHE), encoding='utf-8')
        return local.light_postprocess(str(self.request), self.request_sha, str(self.root/'processed'),
                                       cache_arg=str(cache), offline=True, **kwargs)

    def test_saved_science_is_converted_without_model_calls_or_rewriting(self):
        self.setup_batch()
        original = {path: path.read_bytes() for path in [*self.batch_results, self.request, self.source]}
        result = self.postprocess()
        self.assertEqual('converted', result['state'])
        self.assertEqual(0, result['additional_model_calls'])
        self.assertEqual({'assessed':2, 'not_assessed':0, 'invalid':0, 'retrieval_failed':0}, result['member_counts'])
        receipt = json.loads(Path(result['receipt']).read_text())
        self.assertTrue(receipt['original_results_unchanged'])
        for row in receipt['result_snapshots']:
            self.assertEqual(Path(row['source']['path']).read_bytes(), Path(row['snapshot']['path']).read_bytes())
        self.assertEqual(original, {path:path.read_bytes() for path in original})
        records = json.loads(Path(result['artifacts']['agent_responses']['path']).read_text())['responses']
        claims = [json.loads(row['raw_text'])['claims'][0] for row in records]
        self.assertEqual(['g1','g2'], [row['member_interval_genes'][0] for row in claims])
        self.assertEqual(['T2','T2'], [row['initial_tier'] for row in claims])
        self.assertNotIn('paper_assessments', json.dumps(result))
        index = self.root/'prepared/light/index.json'
        for number, source in enumerate(self.batch_results, 1):
            (index.parent/f'b{number:02d}'/'results.jsonl').write_bytes(source.read_bytes())
        binding = self.request_value['binding']
        direct = light.light_convert(str(self.root/'bundle.json'), binding['bundle_sha256'], str(index),
                                     binding['light_index_sha256'], str(self.root/'direct'),
                                     cache_arg=str(self.root/'cache.json'), offline=True)
        self.assertEqual(Path(direct['agent_responses']).read_bytes(),
                         Path(result['artifacts']['agent_responses']['path']).read_bytes())
        self.assertEqual(Path(direct['member_coverage']['path']).read_bytes(),
                         Path(result['artifacts']['member_coverage']['path']).read_bytes())

    def test_partial_science_preserves_healthy_claim_and_failed_member(self):
        self.setup_batch()
        value = json.loads(self.batch_results[1].read_text())
        value['gene_decisions'][0]['state'] = 'retrieval_failed'
        self.batch_results[1].write_text(json.dumps(value), encoding='utf-8')
        result = self.postprocess()
        self.assertEqual('converted_partial', result['state'])
        self.assertFalse(result['may_claim_all_genes_assessed'])
        self.assertEqual(1, result['member_counts']['retrieval_failed'])
        responses = json.loads(Path(result['artifacts']['agent_responses']['path']).read_text())['responses']
        self.assertEqual('completed', responses[0]['status'])
        self.assertEqual('not_evaluated', responses[1]['status'])

    def test_invalid_answer_is_not_filled_and_healthy_batch_still_converts(self):
        self.setup_batch()
        value = json.loads(self.batch_results[0].read_text())
        del value['gene_decisions'][0]['paper_assessments'][0]['checklist']['focal_species']
        self.batch_results[0].write_text(json.dumps(value), encoding='utf-8')
        original = self.batch_results[0].read_bytes()
        result = self.postprocess()
        self.assertEqual('converted_partial', result['state'])
        self.assertEqual(1, result['member_counts']['invalid'])
        self.assertIn('original_findings', result)
        self.assertEqual(original, self.batch_results[0].read_bytes())
        responses = json.loads(Path(result['artifacts']['agent_responses']['path']).read_text())['responses']
        self.assertEqual('completed', responses[1]['status'])

    def test_explicit_saved_correction_is_applied_then_converted_once(self):
        self.setup_batch()
        decisions = self.root/'decisions.json'
        decisions.write_text(json.dumps({'contract':local.DECISIONS_CONTRACT, 'batches':[
            {'batch_id':'b01','operations':[self.operation()]}, {'batch_id':'b02','unchanged':True}]}), encoding='utf-8')
        original = self.batch_results[0].read_bytes()
        result = self.postprocess(decisions_arg=str(decisions), decisions_sha256=digest(decisions.read_bytes()))
        self.assertEqual('converted', result['state'])
        self.assertIn('apply_receipt', result['artifacts'])
        self.assertEqual(original, self.batch_results[0].read_bytes())
        responses = json.loads(Path(result['artifacts']['agent_responses']['path']).read_text())['responses']
        claim = json.loads(responses[0]['raw_text'])['claims'][0]
        self.assertEqual(self.operation()['value'], claim['phenotype_relation'])

    def test_missing_batch_or_stale_request_is_rejected_before_writing(self):
        self.setup_batch()
        self.request_value['batches'].pop()
        self.save_request()
        with self.assertRaisesRegex(ValueError,'every batch'):
            self.postprocess()
        self.assertFalse((self.root/'processed').exists())
        self.request.write_text(self.request.read_text()+' ', encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):
            self.postprocess()
        self.assertFalse((self.root/'processed').exists())

    def test_existing_cli_postprocess_works_outside_the_repository(self):
        self.setup_batch()
        cache = self.root/'cache.json'
        cache.write_text(json.dumps(MAIN_CACHE), encoding='utf-8')
        wrapper = Path(__file__).resolve().parents[1]/'scripts/agentic_route_b.py'
        completed = subprocess.run([sys.executable, str(wrapper), 'light-postprocess',
            '--request', str(self.request), '--request-sha256', self.request_sha,
            '--out', str(self.root/'cli-processed'), '--cache', str(cache), '--offline'],
            cwd=self.root, capture_output=True, text=True)
        self.assertEqual(0, completed.returncode, completed.stderr)
        result = json.loads(completed.stdout)
        self.assertEqual('converted', result['state'])
        self.assertEqual(0, result['additional_model_calls'])
        self.assertNotIn('paper_assessments', completed.stdout)
