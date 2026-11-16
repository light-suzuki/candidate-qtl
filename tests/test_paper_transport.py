"""Bibliography compression must not merge scientific observations in A or B."""
import copy
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.paper_reports import imported_paper_fields, normalize_paper_reports
from candidate_qtl.relevance_tiers import annotate_relevance
from tests import test_initial_paper_prompt_examples as a_fixtures
from tests.test_initial_paper_prompt_examples import _filled_hypothesis, _valid_report
from tests.test_route_b_light import ConvertHarness, MAIN_CACHE, PAPER, PAPER_TITLE, checks_for, fixture, gene_row, paper_row, run_pipeline, with_database_source


def envelope(reports):
    return {"bibliography": {"p1": {"source_id": reports[0]["source_id"], "source_title": reports[0]["source_title"]}},
            "observations": [{**{k: v for k, v in r.items() if k not in {"source_id", "source_title"}}, "paper_ref": "p1"} for r in reports]}


class PaperTransportTests(unittest.TestCase):
    def test_review_gene_trait_context_differences_survive_lossless_normalization(self):
        first = _valid_report()
        second = {**first, "reported_gene": "OTHER", "reported_gene_id": "NCBI Gene:222222",
                  "reported_phenotype": "root growth", "reported_stage": "seedling",
                  "paper_checklist": {"function_supports_hypothesis": "yes"}}
        third = {**first, "reported_phenotype": "flower development", "reported_stage": "adult"}
        reports = [first, second, third]
        self.assertEqual(normalize_paper_reports(reports), normalize_paper_reports(envelope(reports)))
        self.assertEqual(3, len(json.loads(normalize_paper_reports(envelope(reports)))))
        self.assertEqual("NCBI Gene:222222", second["reported_gene_id"])

    def test_catalog_cannot_supply_biology_or_override_ids(self):
        good = envelope([_valid_report()])
        for change in [lambda v: v['bibliography']['p1'].update(reported_gene='WRONG'),
                       lambda v: v['observations'][0].update(paper_ref='missing'),
                       lambda v: v['observations'][0].update(source_id='PMID:999'),
                       lambda v: v['observations'][0].update(source_title='Another title')]:
            value = copy.deepcopy(good)
            change(value)
            with self.assertRaises(ValueError):
                normalize_paper_reports(value)
        self.assertEqual(good, envelope([_valid_report()]))

    def test_database_function_and_annotation_remain_separate_without_paper(self):
        base = {"candidate_gene_id": "G1", "reference_gene": "REF1", "source_id": "UniProt:P12345",
                "phenotype_relation": "could alter meristem size", "paper_reports_json": []}
        for kind, field, text, expected in [('function', 'reported_function', 'meristem maintenance', 'T4'),
                                           ('annotation', 'reported_annotation', 'specific receptor domain', 'T5')]:
            row = {**base, 'evidence_types': [kind], field: text}
            row.update(imported_paper_fields(row, base['source_id']))
            actual = annotate_relevance([{'candidate_gene_id': 'G1'}], [row])[0]
            self.assertTrue(actual['relevance_tier'].startswith(expected), actual)
        for fields in [{'reported_function': ''}, {'source_id': ''}, {'source_id': 'invented'},
                       {'phenotype_relation': ''}, {'evidence_types': ['annotation'], 'reported_annotation': ''}]:
            row = {**base, 'evidence_types': ['function'], 'reported_function': 'meristem maintenance', **fields}
            actual = annotate_relevance([{'candidate_gene_id': 'G1'}], [row])[0]
            self.assertFalse(actual['relevance_tier'].startswith('T4'), actual)

    def test_route_a_actual_import_accepts_catalog_and_preserves_distinct_review_observations(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            example = a_fixtures.InitialPaperPromptExampleTests()._route_a_example(root)
            reports = [_valid_report(), {**_valid_report(), 'reported_gene': 'OTHER',
                       'reported_gene_id': 'NCBI Gene:222222', 'reported_phenotype': 'root growth'}]
            hypothesis = _filled_hypothesis(example['hypotheses'][0], 'FAMILY_A', envelope(reports))
            payload = root / 'payload.json'
            payload.write_text(json.dumps({'hypotheses': [hypothesis]}), encoding='utf-8')
            with redirect_stdout(StringIO()):
                code = main(['mechanism-import', '--payload', str(payload), '--out-dir', str(root/'import'),
                             '--investigation', str(root/'inv/phenotype_investigation.tsv')])
            self.assertEqual(0, code)
            rows = pd.read_csv(root/'import/mechanism_hypotheses.tsv', sep='\t', dtype=str).fillna('')
            stored = json.loads(rows.iloc[0]['paper_reports_json'])
            self.assertEqual({'NCBI Gene:111111', 'NCBI Gene:222222'}, {r['reported_gene_id'] for r in stored})
            self.assertEqual({'Synthetic paper PMID:123'}, {r['source_title'] for r in stored})


class LightPreservationTests(ConvertHarness):
    def test_inapplicable_paper_does_not_erase_independent_database_function(self):
        for correspondence in ['unknown', 'different_gene']:
            for database in [True, False]:
                with self.subTest(correspondence=correspondence, database=database):
                    root = self.root / f'{correspondence}-{database}'
                    root.mkdir()
                    value = fixture()
                    if database:
                        with_database_source(value, 'u1', 'g1', 'REF-A')
                    paper = paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes'),
                                      gene_correspondence=correspondence, reported_gene_id='TAIR:AT1G12345')
                    built = run_pipeline(root, value, [gene_row('g1', 'REF-A', 'T4', [paper])], MAIN_CACHE, 'u1\n')
                    tiers = built['review'].loc[['g1'], 'relevance_tier'].tolist()
                    self.assertEqual(database, any(tier.startswith('T4') for tier in tiers), tiers)
                    proposals = pd.read_csv(built['imported']['proposals'], sep='\t', dtype=str).fillna('')
                    reports = [report for encoded in proposals['paper_reports_json'] for report in json.loads(encoded or '[]')]
                    self.assertEqual(1, len(reports))
                    self.assertEqual(correspondence, reports[0]['gene_correspondence'])
                    self.assertEqual('TAIR:AT1G12345', reports[0]['reported_gene_id'])
                    self.assertEqual(PAPER_TITLE, reports[0]['source_title'])
                    self.assertEqual(checks_for(function_supports_hypothesis='yes'), json.loads(reports[0]['paper_checklist']))
                    if database:
                        db = proposals[proposals['source_id'] == 'UniProtKB:Q9XXXX']
                        self.assertEqual(1, len(db))
                        self.assertEqual('', db.iloc[0]['paper_reports_json'])
                        self.assertEqual({'unverified_route_b_review'}, set(built['review'].loc[['g1'], 'support_status']))

    def test_absolute_paths_full_cards_and_context_specific_sharing(self):
        value = fixture()
        unit = value['units'][0]
        card = unit['reference_material_cards'][0]
        tail = 'specific meristem process after 400 characters'
        card['function_description'] = 'x' * 500 + tail
        unit['reference_material_cards'].append(copy.deepcopy(card))
        different = copy.deepcopy(card)
        different['go_terms'] = 'distinct root process'
        unit['reference_material_cards'].append(different)
        self.convert([], mutate=lambda bundle: bundle.update(value), cache=MAIN_CACHE, unit_ids='u1\n')
        material = (self.root/'light/light/b01/material.md').read_text(encoding='utf-8')
        task = (self.root/'light/light/b01/task.md').read_text(encoding='utf-8')
        self.assertEqual(2, material.count(tail))
        self.assertIn('distinct root process', material)
        self.assertIn(str((self.root/'light/light/b01/material.md').resolve()), task)
        self.assertIn(str((self.root/'light/light/b01/results.jsonl').resolve()), task)
        self.assertIn('reported_gene_id', task)
        full = json.loads((self.root/'light/light/b01/material.full.json').read_text(encoding='utf-8'))
        self.assertEqual([unit], full['units'])

    def test_full_material_tampering_is_rejected(self):
        from scripts import agentic_route_b_light as light
        from tests.test_route_b_light import digest
        self.convert([], cache=MAIN_CACHE, unit_ids='u1\n')
        full = self.root/'light/light/b01/material.full.json'
        full.write_text('{}', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            light.light_convert(str(self.bundle), self.sha, str(self.index), digest(self.index.read_bytes()),
                                str(self.root/'tampered'), offline=True)

    def test_light_bibliography_expansion_and_missing_ref_is_local(self):
        paper = paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes'), reported_gene_id='TAIR:AT1G12345')
        observation = {k:v for k,v in paper.items() if k not in {'paper_id','title'}}
        rows = [{'bibliography': {'p1': {'source_id': PAPER, 'source_title': PAPER_TITLE}}},
                gene_row('g1','REF-A','T4',[{**observation,'paper_ref':'p1'}]),
                gene_row('g2','REF-B','T4',[{**observation,'paper_ref':'missing'}])]
        self.convert(rows, cache=MAIN_CACHE, unit_ids='u1\nu2\n')
        stored = self.claims('u1')[0]['paper_reports_json'][0]
        self.assertEqual('TAIR:AT1G12345', stored['reported_gene_id'])
        self.assertEqual(PAPER_TITLE, stored['source_title'])
        self.assertEqual('not_evaluated', self.responses()['route-b-semantic-u2']['status'])

    def test_six_distinct_observations_of_one_review_survive_to_candidate_freeze(self):
        reports = [paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes', identity_resolved='yes'),
                            reported_gene=f'GENE{i}', reported_gene_id=f'NCBI Gene:{100+i}',
                            reported_phenotype=f'phenotype {i}', finding=f'Independent observation {i}') for i in range(6)]
        built = run_pipeline(self.root, fixture(), [gene_row('g1','REF-A','T4',reports)], MAIN_CACHE, 'u1\n')
        proposals = pd.read_csv(built['imported']['proposals'], sep='\t', dtype=str).fillna('')
        stored = json.loads(proposals.iloc[0]['paper_reports_json'])
        self.assertEqual(6, len(stored))
        self.assertEqual({f'NCBI Gene:{100+i}' for i in range(6)}, {r['reported_gene_id'] for r in stored})
        self.assertTrue(built['review'].loc['g1','relevance_tier'].startswith('T4'))
        discoveries = json.loads((self.root/'converted/discoveries.json').read_text(encoding='utf-8'))
        self.assertEqual(6, len(discoveries['discovered_sources']))

    def test_light_preserves_optional_observation_context_and_followup(self):
        fields = {'phenotype_comparison': 'different', 'underlying_source_ids': 'PMID:456',
                  'corroboration_note': 'independent root observation', 'followup_question': 'inspect figure 2',
                  'followup_state': 'completed', 'reported_function': 'root maintenance',
                  'reported_annotation': 'specific receptor domain'}
        paper = paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes'), **fields)
        self.convert([gene_row('g1', 'REF-A', 'T4', [paper])], cache=MAIN_CACHE, unit_ids='u1\n')
        stored = self.claims('u1')[0]['paper_reports_json'][0]
        self.assertEqual(fields, {key: stored[key] for key in fields})

    def test_distinct_counterevidence_survives_light_and_candidate_freeze(self):
        reports = [paper_row(PAPER, PAPER_TITLE, checks_for(function_supports_hypothesis='yes', identity_resolved='yes'),
                            counterevidence=reason, read_scope='Results', paper_checklist_reasons={'function_supports_hypothesis': reason})
                   for reason in ['root effect limited to seedlings', 'adult shoots were unaffected']]
        built = run_pipeline(self.root, fixture(), [gene_row('g1', 'REF-A', 'T4', reports)], MAIN_CACHE, 'u1\n')
        proposals = pd.read_csv(built['imported']['proposals'], sep='\t', dtype=str).fillna('')
        stored = json.loads(proposals.iloc[0]['paper_reports_json'])
        self.assertEqual(2, len(stored))
        self.assertEqual({report['counterevidence'] for report in reports},
                         {report['research_annotations']['counterevidence'] for report in stored})
        self.assertEqual({'Results'}, {report['research_annotations']['read_scope'] for report in stored})
        self.assertTrue(built['review'].loc['g1', 'relevance_tier'].startswith('T4'))
