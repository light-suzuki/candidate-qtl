"""Large frozen provenance must survive scope freeze and selection losslessly."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.boundary_tier_review import _read, consume_boundary_selection, freeze_boundary_scope
from candidate_qtl.relevance_tiers import T3


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


class LargeBoundaryFieldsTests(unittest.TestCase):
    def setUp(self):
        self.original_limit = csv.field_size_limit()
        csv.field_size_limit(131072)
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        csv.field_size_limit(self.original_limit)
        self.temp.cleanup()

    def test_large_unicode_quoted_field_is_read_exactly(self):
        payload = ('evidence\tquoted "value"\n根拠\n' * 10000)
        source = self.root / 'large.tsv'
        rows = [{'candidate_gene_id': 'G1', 'provenance': payload}]
        write(source, rows)
        expected_sha = sha(source)
        self.assertEqual(_read(source), (list(rows[0]), rows))
        self.assertEqual(sha(source), expected_sha)

    def test_scope_and_consumer_preserve_large_provenance_and_all_copies(self):
        evidence = json.dumps({'observations': ['literal evidence "line"\n' * 10000]}, ensure_ascii=False)
        roster = [dict(candidate_gene_id=f'G{i:02}', relevance_tier=T3, support_state='review',
                       copy_state='paralog_ambiguous', relevance_evidence_json=evidence if i == 0 else '{}')
                  for i in range(26)]
        source = self.root / 'inventory.tsv'
        write(source, roster)
        scope = freeze_boundary_scope(str(source), str(self.root / 'scope'), inventory_sha256=sha(source), case_id='neutral')
        boundary = Path(scope['artifacts']['boundary']['path'])
        self.assertEqual(_read(boundary)[1], roster)
        proposal = self.root / 'proposals.tsv'
        write(proposal, [dict(candidate_gene_id=r['candidate_gene_id'], decision='adopt', reason='neutral retain all copies',
                             source_wet_sha256=sha(boundary), executor_id='proposer', executor_model='model-a') for r in roster])
        independent = self.root / 'independent.tsv'
        write(independent, [dict(candidate_gene_id=r['candidate_gene_id'], decision='adopt', reason='independent neutral retain',
                                source_wet_sha256=sha(boundary), proposal_sha256=sha(proposal),
                                reviewer_id='reviewer', reviewer_model='model-b') for r in roster])
        result = consume_boundary_selection(scope['manifest'], str(proposal), str(self.root / 'selected'),
                                            scope_sha256=scope['manifest_sha256'], case_id='neutral',
                                            executor_id='proposer', executor_model='model-a',
                                            reviewer_id='reviewer', reviewer_model='model-b',
                                            independent_review_arg=str(independent))
        selected = _read(Path(result['wet_candidates']['path']))[1]
        self.assertEqual(len(selected), 26)
        self.assertEqual(selected[0]['relevance_evidence_json'], evidence)
        self.assertTrue(all(r['copy_state'] == 'paralog_ambiguous' and r['support_state'] == 'review' for r in selected))

    def test_malformed_roster_guards_still_reject_missing_or_extra_fields(self):
        for suffix in ('G1\n', 'G1\tx\textra\n'):
            path = self.root / 'bad.tsv'
            path.write_text('candidate_gene_id\tprovenance\n' + suffix, encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'malformed candidate table row'):
                _read(path)
