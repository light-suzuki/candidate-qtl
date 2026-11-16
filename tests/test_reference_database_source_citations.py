"""Database identities must not suppress retrieval of linked publications."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.source_verify import source_verify
from test_reference_gene_resolution import _mechanism_row, _uniprot_record, _write_cache, _write_mechanisms


class DatabaseSourceCitationTests(unittest.TestCase):
    def test_database_only_source_reaches_source_text_consumer(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            mechanisms = root/'mechanisms.tsv'
            _write_mechanisms(mechanisms, [_mechanism_row(
                'synthetic', source_ids='UniProt:Q9ZZZ9', proposed_stable_id='Q9ZZZ9',
                reference_gene_name='')])
            cache = root/'records.json'
            _write_cache(cache, [_uniprot_record('Q9ZZZ9', 'Arabidopsis thaliana', 'MAAA', pmid='12345')])
            result = reference_gene_resolve(str(mechanisms), str(root/'resolved'), record_cache_arg=str(cache))
            rows = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('')
            self.assertEqual('UniProt:Q9ZZZ9;PMID:12345', rows.iloc[0]['source_ids'])
            ledger = pd.read_csv(result['ledger'], sep='\t').fillna('')
            self.assertEqual('UniProt:Q9ZZZ9', ledger.iloc[0]['input_source_ids'])
            literature = root/'literature.json'
            literature.write_text(json.dumps([{'id':'12345', 'pmid':'12345', 'source':'MED',
                'title':'Synthetic primary experiment', 'abstractText':'Synthetic captured abstract.'}]))
            source_verify(result['mechanisms_resolved'], str(root/'source'), cache_arg=str(literature))
            captured = json.loads((root/'source/source_review_records.json').read_text())
            self.assertEqual(1, len(captured['records']))

    def test_real_literature_and_unknown_sources_remain_unchanged(self):
        for source in ['PMID:12345;UniProt:Q9ZZZ9', 'unrecognised source']:
            with self.subTest(source=source), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                path = root/'mechanisms.tsv'
                _write_mechanisms(path, [_mechanism_row('synthetic', source_ids=source, proposed_stable_id='Q9ZZZ9')])
                result = reference_gene_resolve(str(path), str(root/'out'))
                row = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('').iloc[0]
                self.assertEqual(source, row['source_ids'])

    def test_offline_unavailable_database_source_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root/'mechanisms.tsv'
            _write_mechanisms(path, [_mechanism_row('synthetic', source_ids='UniProt:Q9ZZZ9', proposed_stable_id='Q9ZZZ9')])
            result = reference_gene_resolve(str(path), str(root/'out'))
            row = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('').iloc[0]
            ledger = pd.read_csv(result['ledger'], sep='\t').fillna('').iloc[0]
            self.assertEqual('UniProt:Q9ZZZ9', row['source_ids'])
            self.assertEqual('not_attempted', ledger['evaluation_status'])
            self.assertEqual('False', str(ledger['function_source_backed']))

    def test_source_only_identity_can_resolve_without_gene_name(self):
        for source, tair in [('UniProt:Q9ZZZ9', None), ('TAIR:AT1G00010', 'AT1G00010')]:
            with self.subTest(source=source), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                path, cache = root/'mechanisms.tsv', root/'records.json'
                _write_mechanisms(path, [_mechanism_row('synthetic', source_ids=source)])
                _write_cache(cache, [_uniprot_record('Q9ZZZ9', 'Arabidopsis thaliana', 'MAAA', tair=tair)])
                result = reference_gene_resolve(str(path), str(root/'out'), record_cache_arg=str(cache))
                row = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('').iloc[0]
                self.assertEqual('Q9ZZZ9', row['proposed_stable_id'])
                self.assertIn(source, row['source_ids'])
                self.assertIn('PMID:12345', row['source_ids'])

    def test_conflicting_identities_are_never_chosen_or_queried(self):
        for source, supplied in [('UniProt:Q9ZZZ9;UniProt:Q9ZZZ8', ''), ('UniProt:Q9ZZZ9', 'Q9ZZZ8')]:
            with self.subTest(source=source), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                path = root/'mechanisms.tsv'
                _write_mechanisms(path, [_mechanism_row('synthetic', source_ids=source, proposed_stable_id=supplied)])
                with mock.patch('candidate_qtl.reference_gene_resolution._fetch_uniprot') as fetch:
                    result = reference_gene_resolve(str(path), str(root/'out'), allow_external_api=True)
                fetch.assert_not_called()
                row = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('').iloc[0]
                ledger = pd.read_csv(result['ledger'], sep='\t').fillna('').iloc[0]
                self.assertEqual(source, row['source_ids'])
                self.assertEqual(supplied, row['proposed_stable_id'])
                self.assertEqual('database_source_identity_unresolved', ledger['review_reason'])

    def test_missing_publications_and_provider_failure_do_not_create_support(self):
        for unavailable in [False, True]:
            with self.subTest(unavailable=unavailable), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                path = root/'mechanisms.tsv'
                _write_mechanisms(path, [_mechanism_row('synthetic', source_ids='UniProt:Q9ZZZ9', proposed_stable_id='Q9ZZZ9')])
                with mock.patch('candidate_qtl.reference_gene_resolution._fetch_uniprot',
                                side_effect=TimeoutError('synthetic timeout') if unavailable else None,
                                return_value=_uniprot_record('Q9ZZZ9', 'Arabidopsis thaliana', 'MAAA', pmid='')):
                    result = reference_gene_resolve(str(path), str(root/'out'), allow_external_api=True)
                row = pd.read_csv(result['mechanisms_resolved'], sep='\t').fillna('').iloc[0]
                ledger = pd.read_csv(result['ledger'], sep='\t').fillna('').iloc[0]
                self.assertEqual('UniProt:Q9ZZZ9', row['source_ids'])
                self.assertEqual('not_resolved', ledger['resolution_origin'])
                if unavailable:
                    self.assertEqual('not_evaluated', ledger['evaluation_status'])
