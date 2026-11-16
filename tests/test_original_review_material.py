"""Original source material remains readable without becoming evidence."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.source_review_capture import content_sha256
from candidate_qtl.source_verify import source_verify


class OriginalReviewMaterialTests(unittest.TestCase):
    def run_capture(self, root, rows, records):
        table = root / 'input.tsv'
        pd.DataFrame(rows).to_csv(table, sep='\t', index=False)
        cache = root / 'cache.json'
        cache.write_text(json.dumps(records), encoding='utf-8')
        with mock.patch('candidate_qtl.source_verify._fetch_europe_pmc', side_effect=AssertionError('no new retrieval')):
            result = source_verify(str(table), str(root / 'out'), cache_arg=str(cache))
        return result, json.loads(Path(result['source_review_records']).read_text())

    def test_shared_original_paper_has_lineage_but_no_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [{'hypothesis_id': hid, 'parent_hypothesis_id': 'parent', 'source_ids': '',
                     'review_only_source_ids': 'PMID:99999; PMID:88888', 'reference_gene_name': 'NODE1'}
                    for hid in ('child_a', 'child_b')]
            record = {'id': '99999', 'pmid': '99999', 'title': 'NODE1 mutant experiment',
                      'abstractText': '  NODE1 mutants show a stem defect.\n'}
            result, capture = self.run_capture(root, rows, [record])
            source = capture['records']['PMID:99999']
            self.assertEqual(record['abstractText'], source['sections']['abstractText'])
            self.assertEqual(content_sha256(source['sections']), source['sections_sha256'])
            self.assertEqual(['child_a', 'child_b'], source['review_only_hypothesis_ids'])
            self.assertEqual(['child_a', 'child_b'], source['hypothesis_ids'])
            self.assertEqual('review_only_not_classified', source['classification_scope'])
            self.assertEqual('not_evaluated', source['semantic_binding_state'])
            self.assertEqual('none', capture['membership_effect'])
            manifest = json.loads(Path(result['manifest']).read_text())
            self.assertEqual(0, manifest['citation_count'])
            self.assertEqual(0, manifest['verified_citation_count'])
            self.assertEqual(2, manifest['review_only_material']['captured_association_count'])
            states = {item['state'] for item in manifest['review_only_material']['associations']}
            self.assertEqual({'captured_review_only', 'not_evaluated_cache_missing'}, states)
            annotations = pd.read_csv(result['mechanisms_with_verified_sources'], sep='\t', dtype=str).fillna('')
            self.assertTrue(annotations['verified_source_ids'].eq('').all())
            self.assertTrue(annotations['mechanistic_support_citation_count'].eq('0').all())
            self.assertEqual({}, json.loads(Path(result['verified_sources']).read_text()))

    def test_source_identity_mismatch_is_not_captured(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, capture = self.run_capture(Path(tmp),
                [{'hypothesis_id': 'h', 'source_ids': '', 'review_only_source_ids': 'PMID:99999'}],
                [{'cache_key': 'pmid_99999', 'id': '88888', 'pmid': '88888', 'abstractText': 'Wrong paper.'}])
            self.assertEqual({}, capture['records'])
            manifest = json.loads(Path(result['manifest']).read_text())
            self.assertEqual('not_evaluated_cached_identifier_mismatch',
                             manifest['review_only_material']['associations'][0]['state'])
            self.assertEqual(0, manifest['rejected_citation_count'])

    def test_review_only_capture_cannot_authorize_support_even_with_asserted_axes(self):
        from tests.test_agent_claim_review import ReviewFixture
        with tempfile.TemporaryDirectory() as tmp:
            fixture = ReviewFixture(Path(tmp))
            source = json.loads(fixture.sources_path.read_text())
            source['records']['PMID:12345']['review_only_hypothesis_ids'] = ['h1']
            fixture.sources_path.write_text(json.dumps(source), encoding='utf-8')
            fixture.seal()
            result = fixture.bind()
            self.assertEqual('review_required', result['binding_state'])
            self.assertIn('review-only', result['binding_reason'])

    def test_shared_paper_review_only_child_does_not_block_existing_primary_binding(self):
        from tests.test_agent_claim_review import ReviewFixture
        with tempfile.TemporaryDirectory() as tmp:
            fixture = ReviewFixture(Path(tmp))
            source = json.loads(fixture.sources_path.read_text())
            source['records']['PMID:12345']['review_only_hypothesis_ids'] = ['alternative_child']
            fixture.sources_path.write_text(json.dumps(source), encoding='utf-8')
            fixture.seal()
            self.assertEqual('verified_complete', fixture.bind()['binding_state'])

    def test_review_only_exact_accession_requires_verified_quote_and_all_binding_gates(self):
        from tests.test_agent_claim_review import ReviewFixture
        variants = [('P12345', False, 'verified_complete'),
                    ('P54321', False, 'review_required'),
                    ('P12345', True, 'not_evaluated')]
        for quoted_accession, wrong_span, expected in variants:
            with self.subTest(accession=quoted_accession, wrong_span=wrong_span), tempfile.TemporaryDirectory() as tmp:
                fixture = ReviewFixture(Path(tmp))
                text = (f'In Arabidopsis thaliana, loss of XYZ1 ({quoted_accession}) caused short roots. '
                        'Complementation restored root length.')
                sources = json.loads(fixture.sources_path.read_text())
                source = sources['records']['PMID:12345']
                source['sections']['abstractText'] = text
                source['sections_sha256'] = content_sha256(source['sections'])
                source['review_only_hypothesis_ids'] = ['h1']
                source['classification_scope'] = 'review_only_not_classified'
                source['semantic_binding_state'] = 'not_evaluated'
                fixture.sources_path.write_text(json.dumps(sources))
                fixture.record['references'] = []
                fixture.review.update(quoted_support=text, quote_start=1 if wrong_span else 0,
                    quote_end=len(text), sections_sha256=source['sections_sha256'])
                fixture.seal()
                result = fixture.checker().bind({'hypothesis_id': 'h1'}, 'P12345', fixture.record, [])
                self.assertEqual(expected, result['binding_state'])

    def test_existing_single_exact_identity_source_lane_is_unchanged(self):
        from candidate_qtl.reference_gene_resolution import reference_gene_resolve
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / 'input.tsv'
            pd.DataFrame([{'case_id': 'neutral', 'hypothesis_id': 'h', 'reference_gene_name': 'NODE1',
                           'reference_species': 'Species alpha', 'proposed_stable_id': 'Q00001',
                           'source_ids': 'PMID:99999'}]).to_csv(table, sep='\t', index=False)
            cache = root / 'record.json'
            cache.write_text(json.dumps({'primaryAccession': 'Q00001',
                'organism': {'scientificName': 'Species alpha'},
                'genes': [{'geneName': {'value': 'NODE1'}}], 'sequence': {'value': 'M' * 30}}))
            result = reference_gene_resolve(str(table), str(root / 'resolved'), record_cache_arg=str(cache))
            row = pd.read_csv(result['mechanisms_resolved'], sep='\t', dtype=str).fillna('').iloc[0]
            self.assertEqual('PMID:99999', row['source_ids'])
            self.assertEqual('', row['review_only_source_ids'])

    def test_resolver_alternatives_and_pending_keep_cached_original_for_review(self):
        from candidate_qtl.reference_gene_resolution import reference_gene_resolve
        from candidate_qtl.agent_claim_review import prepare_review_prompt
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / 'hypotheses.tsv'
            pd.DataFrame([{'case_id': 'fictional', 'hypothesis_id': 'parent',
                           'reference_gene_name': 'NODE1', 'reference_species': 'Species alpha',
                           'source_ids': 'PMID:99999', 'proposed_stable_id': 'Q00001'}]).to_csv(table, sep='\t', index=False)
            records = [{'primaryAccession': accession, 'organism': {'scientificName': 'Species alpha'},
                        'genes': [{'geneName': {'value': 'NODE1'}}], 'sequence': {'value': 'M' * 30}}
                       for accession in ('Q00002', 'Q00003')]
            for label, provider_records in (('alternatives', records), ('pending', [])):
                with self.subTest(label=label):
                    cache = root / (label + '.json')
                    cache.write_text(json.dumps({'results': provider_records}), encoding='utf-8')
                    resolved = reference_gene_resolve(str(table), str(root / label), record_cache_arg=str(cache))
                    rows = pd.read_csv(resolved['mechanisms_resolved'], sep='\t', dtype=str).fillna('')
                    self.assertEqual(2 if label == 'alternatives' else 1, len(rows))
                    self.assertTrue(rows['review_only_source_ids'].eq('PMID:99999').all())
                    self.assertTrue(rows['source_ids'].eq('').all())
                    capture_root = root / (label + '-capture')
                    capture_root.mkdir()
                    captured, sources = self.run_capture(capture_root, rows.to_dict('records'),
                        [{'id': '99999', 'pmid': '99999', 'title': 'NODE1 experiment',
                          'abstractText': 'Original NODE1 perturbation material.'}])
                    self.assertEqual(len(rows), len(sources['records']['PMID:99999']['review_only_hypothesis_ids']))
                    prompt = prepare_review_prompt(captured['mechanisms_with_verified_sources'],
                                                   captured['source_review_records'], str(root / (label + '-prompt')))
                    self.assertIn('Original NODE1 perturbation material.', prompt.read_text())
                    self.assertEqual({}, json.loads(Path(captured['verified_sources']).read_text()))


if __name__ == '__main__':
    unittest.main()
