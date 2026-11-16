"""Synthetic review conditions and native-context tests; no scientific ranking."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from candidate_qtl import second_stage_checks as checks
from candidate_qtl import second_stage_recheck as recheck
from tests import test_second_stage_recheck as recheck_tests
from tests import test_second_stage_guided as guided_tests


def table(path, fields, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    return str(path)


def hit(**changes):
    return dict(gene_id='gA', protein_id='pA', subject_id='rA',
                query_residue_identity=60, subject_residue_identity=60,
                query_coverage=.9, subject_coverage=.9, evalue=1e-20) | changes


class ReviewProfileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.profiles = checks.load_profiles(None)
    def states(self, row):
        return {r['profile_id']: r for r in checks.alignment_checks(row, self.profiles, 10)}
    def test_stricter_view_does_not_change_broad_view_or_hit(self):
        row = hit(query_coverage=.7); before = dict(row)
        result = self.states(row)
        self.assertEqual('meets_review_conditions', result['broad_alignment_review']['state'])
        self.assertEqual('below_review_conditions', result['strict_full_length_review']['state'])
        self.assertEqual(['query_coverage'], result['strict_full_length_review']['unmet_conditions'])
        self.assertEqual(row, before)
    def test_both_coverages_and_both_identities_are_checked(self):
        for key, value in [('query_coverage', .2), ('subject_coverage', .2),
                           ('query_residue_identity', 5), ('subject_residue_identity', 5), ('evalue', 1)]:
            with self.subTest(key=key):
                self.assertIn(key, self.states(hit(**{key: value}))['strict_full_length_review']['unmet_conditions'])
    def test_high_local_identity_does_not_override_low_subject_coverage(self):
        result = self.states(hit(query_residue_identity=100, subject_residue_identity=100,
                                 query_coverage=1, subject_coverage=.1))
        self.assertEqual('below_review_conditions', result['strict_full_length_review']['state'])
        self.assertIn('subject_coverage', result['strict_full_length_review']['unmet_conditions'])
    def test_boundary_is_a_measurement_comparison_not_gene_deletion(self):
        row = hit(query_residue_identity=50, subject_residue_identity=50,
                  query_coverage=.8, subject_coverage=.8, evalue=1e-10)
        result = self.states(row)['strict_full_length_review']
        self.assertEqual('meets_review_conditions', result['state'])
        self.assertEqual('none', result['candidate_membership_effect'])
    def test_search_limits_are_not_complete_negative_search_coverage(self):
        result = checks.alignment_checks(hit(), self.profiles, 1e-40)
        self.assertTrue(all(r['search_more_restrictive_than_review_evalue'] for r in result))
    def test_custom_profiles_are_explicit_and_order_stable(self):
        p = self.root/'profiles.json'; p.write_text(json.dumps(list(reversed(self.profiles))))
        self.assertEqual(self.profiles, checks.load_profiles(p))
    def test_invalid_profiles_fail_closed_without_guessing(self):
        invalid = [[], {}, [True], [self.profiles[0], self.profiles[0]]]
        for key, value in [('min_identity', True), ('min_identity', 101), ('min_query_coverage', -1),
                           ('max_evalue', float('nan')), ('max_evalue', float('inf')), ('basis', '')]:
            invalid.append([dict(self.profiles[0]) | {key: value}])
        invalid.append([dict(self.profiles[0]) | {'score': 100}])
        p = self.root/'profiles.json'
        for payload in invalid:
            with self.subTest(payload=str(payload)):
                p.write_text(json.dumps(payload))
                with self.assertRaises(ValueError): checks.load_profiles(p)


class FamilyChecksTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.states = {'interval_atlas': 'available', 'family_copy': 'available'}
    def test_coarse_bucket_and_domain_are_not_functional_family(self):
        atlas = [{'family_id': 'RLK', 'discovery_bucket_method': 'symbol_stem', 'domain_ids': 'PF00069'}]
        result = checks.family_checks('gA', atlas, [], self.states)
        self.assertIn('coarse_bucket_is_not_evolutionary_family_verification', result['review_flags'])
        self.assertIn('domain_only_is_not_functional_family_identity', result['review_flags'])
        self.assertFalse(result['new_family_assignment_performed'])
    def test_target_mismatch_is_reported_not_a_refinement_verdict(self):
        rows = [{'family_id': 'F', 'target_id': 'T', 'target_verification_status': 'curated',
                 'subfamily_vs_target': 'target_mismatch', 'domain_vs_target': 'target_match'}]
        r = checks.family_checks('gA', [], rows, self.states)
        self.assertEqual('reported_target_mismatch', r['target_relative_checks'][0]['subfamily_vs_target'])
        self.assertEqual('none', r['candidate_membership_effect'])
        self.assertFalse(r['independent_source_verification_performed'])
        self.assertNotIn('support_falsified', json.dumps(r))
    def test_unverified_target_cannot_turn_match_or_mismatch_into_verified_relation(self):
        for status in ('pending_review', 'unknown', 'false', ''):
            row = {'target_id': 'T', 'target_verification_status': status,
                   'subfamily_vs_target': 'target_mismatch', 'domain_vs_target': 'target_match'}
            r = checks.family_checks('gA', [], [row], self.states)
            self.assertEqual('not_evaluated_target_not_bound', r['target_relative_checks'][0]['subfamily_vs_target'])
            self.assertNotIn('target_relative_mismatch_requires_review_not_automatic_exclusion', r['review_flags'])
    def test_target_status_without_identifier_is_not_bound(self):
        r = checks.family_checks('gA', [], [{'target_verification_status': 'verified',
                                           'domain_vs_target': 'target_mismatch'}], self.states)
        self.assertEqual('not_evaluated_target_not_bound', r['target_relative_checks'][0]['domain_vs_target'])
    def test_modal_domain_difference_is_separate_from_target_comparison(self):
        row = {'domain_architecture_match': 'architecture_mismatch', 'target_id': 'T',
               'target_verification_status': 'verified', 'domain_vs_target': 'target_match'}
        r = checks.family_checks('gA', [], [row], self.states)
        self.assertIn('family_modal_domain_difference_is_qc_not_selection', r['review_flags'])
        self.assertEqual('reported_target_match', r['target_relative_checks'][0]['domain_vs_target'])
    def test_ambiguity_and_partial_model_are_not_negative_evidence(self):
        a = [{'gene_model_completeness': 'partial_model', 'subfamily_source': 'orthogroup_proxy'}]
        c = [{'orthology_class': 'one_to_many'}]
        r = checks.family_checks('gA', a, c, self.states)
        self.assertIn('copy_or_orthology_ambiguity_must_remain', r['review_flags'])
        self.assertIn('model_incompleteness_may_explain_alignment_or_domain_difference', r['review_flags'])
        self.assertEqual('none', r['candidate_membership_effect'])
    def test_no_family_data_is_not_evaluated_not_rejected(self):
        r = checks.family_checks('gA', [], [], self.states)
        self.assertEqual('not_evaluated', r['state'])
        self.assertEqual('none', r['candidate_membership_effect'])
    def test_all_gene_family_target_rows_survive_without_source_mixing(self):
        p = self.root/'copy.tsv'
        rows = [{'candidate_gene_id': 'gA', 'family_id': 'F1', 'target_id': 'T1',
                 'target_verification_status': 'curated', 'subfamily_vs_target': 'target_match'},
                {'candidate_gene_id': 'gA', 'family_id': 'F2', 'target_id': 'T2',
                 'target_verification_status': 'pending_review', 'subfamily_vs_target': 'target_mismatch'}]
        table(p, list(rows[0]), rows)
        first, state = checks.load_context(p, 'family_copy', {'gA'})
        table(p, list(rows[0]), list(reversed(rows)) + [rows[0]])
        second, _ = checks.load_context(p, 'family_copy', {'gA'})
        self.assertEqual(first, second); self.assertEqual(2, len(first['gA']))
        r = checks.family_checks('gA', [], first['gA'], self.states)
        self.assertEqual('reported_target_match', r['target_relative_checks'][0]['subfamily_vs_target'])
        self.assertEqual('not_evaluated_target_not_bound', r['target_relative_checks'][1]['subfamily_vs_target'])
    def test_native_ids_are_exact_and_rank_fields_are_not_consumed(self):
        p = self.root/'atlas.tsv'
        table(p, ['gene_id', 'family_id', 'rank', 'dry_score'],
              [{'gene_id': 'gA.t1', 'family_id': 'F', 'rank': '1', 'dry_score': '99'},
               {'gene_id': 'gA', 'family_id': 'G', 'rank': '2', 'dry_score': '98'}])
        result, _ = checks.load_context(p, 'interval_atlas', {'gA'})
        self.assertEqual(['gA'], list(result))
        self.assertEqual('G', result['gA'][0]['family_id'])
        self.assertNotIn('rank', result['gA'][0]); self.assertNotIn('dry_score', result['gA'][0])
    def test_missing_and_empty_and_invalid_context_distinguished(self):
        self.assertEqual('not_evaluated_not_supplied', checks.load_context(None, 'interval_atlas', set())[1])
        p = self.root/'atlas.tsv'
        self.assertEqual('not_evaluated_missing_file', checks.load_context(p, 'interval_atlas', set())[1])
        p.write_text('gene_id\n')
        self.assertEqual(({}, 'available'), checks.load_context(p, 'interval_atlas', set()))
        p.write_text('symbol\ngA\n')
        with self.assertRaises(ValueError): checks.load_context(p, 'interval_atlas', {'gA'})
    def test_conflicting_atlas_annotations_are_preserved(self):
        p = self.root/'atlas.tsv'
        table(p, ['gene_id', 'evolutionary_family_id'],
              [{'gene_id': 'gA', 'evolutionary_family_id': 'F1'}, {'gene_id': 'gA', 'evolutionary_family_id': 'F2'}])
        a, _ = checks.load_context(p, 'interval_atlas', {'gA'})
        r = checks.family_checks('gA', a['gA'], [], self.states)
        self.assertEqual(2, len(r['atlas_context']))
        self.assertIn('multiple_atlas_assignments_require_review', r['review_flags'])


class ChecksIntegrationTests(unittest.TestCase):
    """Reuse the existing subprocess fixture without duplicating its test suite."""
    def setUp(self):
        self.fixture = recheck_tests.FreshRecheckTests('test_fresh_searches_execute_and_first_stage_is_unchanged')
        self.fixture.setUp(); self.addCleanup(self.fixture.doCleanups)
        self.f = self.fixture
    def test_real_recheck_delivers_checks_and_family_to_every_prompt(self):
        f = self.f
        atlas = f.root/'atlas.tsv'
        table(atlas, ['gene_id', 'family_id', 'evolutionary_family_id', 'subfamily_id'],
              [{'gene_id': 'gA', 'family_id': 'F', 'evolutionary_family_id': '', 'subfamily_id': ''}])
        copy_path = f.root/'copy.tsv'
        table(copy_path, ['candidate_gene_id', 'family_id', 'target_id', 'target_verification_status', 'subfamily_vs_target'],
              [{'candidate_gene_id': 'gA', 'family_id': 'F', 'target_id': 'T',
                'target_verification_status': 'curated', 'subfamily_vs_target': 'target_mismatch'}])
        before = f.request.read_bytes()
        r = f.run_recheck(interval_atlas_arg=str(atlas), family_copy_arg=str(copy_path))
        self.assertEqual(3, r['prompt_count']); self.assertEqual(before, f.request.read_bytes())
        payload = json.loads(f.prompts()[0]['user'])
        self.assertEqual(2, len(payload['reference_hits']))
        self.assertTrue(all(len(h['review_conditions']) == 2 for h in payload['reference_hits']))
        self.assertEqual('reported_target_mismatch', payload['family_review']['target_relative_checks'][0]['subfamily_vs_target'])
        self.assertIn('family_copy', r['inputs']); self.assertTrue(r['inputs']['family_copy']['sha256'])
        self.assertTrue(Path(r['review_checks']).is_file())
        self.assertTrue((f.out/'alignment_review_checks.tsv').is_file())
        self.assertFalse((f.out/'wet_candidates.tsv').exists())
    def test_profile_failure_keeps_every_returned_hit_and_all_genes(self):
        r = self.f.run_recheck()
        for p in self.f.prompts()[:2]:
            data = json.loads(p['user'])
            self.assertEqual(2, len(data['reference_hits']))
            strict = [v for h in data['reference_hits'] for v in h['review_conditions'] if v['profile_id'] == 'strict_full_length_review']
            self.assertTrue(all(x['state'] == 'below_review_conditions' for x in strict))
        self.assertEqual(['gA', 'gB', 'missing'], r['candidate_gene_ids'])
    def test_custom_profile_changes_reassessment_without_research_reranking(self):
        p = self.f.root/'profiles.json'; data = [dict(checks.DEFAULT_PROFILES[0])]
        data[0]['max_evalue'] = 1; data[0]['min_identity'] = 0; data[0]['basis'] = 'custom-inspection-only'
        p.write_text(json.dumps(data))
        r = self.f.run_recheck(review_profiles_arg=str(p))
        self.assertEqual(data, r['settings']['review_profiles'])
        self.assertEqual(3, r['candidate_count'])
    def test_reuse_refuses_changed_family_input_or_review_profiles(self):
        p = self.f.root/'profiles.json'; p.write_text(json.dumps(list(checks.DEFAULT_PROFILES)))
        atlas = self.f.root/'atlas.tsv'; table(atlas, ['gene_id','family_id'], [{'gene_id':'gA','family_id':'F'}])
        self.f.run_recheck(review_profiles_arg=str(p), interval_atlas_arg=str(atlas))
        atlas.write_text('gene_id\tfamily_id\ngA\tG\n')
        with self.assertRaisesRegex(ValueError, 'changed'):
            self.f.run_recheck(review_profiles_arg=str(p), interval_atlas_arg=str(atlas), reuse_existing=True)
    def test_missing_family_file_keeps_prompts_and_states(self):
        r = self.f.run_recheck(family_copy_arg=str(self.f.root/'absent.tsv'))
        self.assertEqual('not_evaluated_missing_file', r['family_context_states']['family_copy'])
        self.assertEqual(3, r['prompt_count'])
    def test_invalid_profile_fails_before_blast(self):
        p = self.f.root/'profiles.json'; p.write_text('[]')
        with mock.patch.object(recheck, '_blast') as run:
            with self.assertRaises(ValueError): self.f.run_recheck(review_profiles_arg=str(p))
        run.assert_not_called(); self.assertFalse(self.f.out.exists())
    def test_standalone_cli_new_flags_reach_real_outputs(self):
        import contextlib
        import io
        f = self.f
        profile = f.root/'profiles.json'; profile.write_text(json.dumps(list(checks.DEFAULT_PROFILES)))
        atlas = f.root/'atlas.tsv'; table(atlas, ['gene_id','family_id'], [{'gene_id':'gA','family_id':'F'}])
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            code = recheck.main(['--request', str(f.request), '--out-dir', str(f.out),
                                 '--phenotype', 'observed morphology', '--proteins', str(f.proteins),
                                 '--reference-proteins', str(f.references), '--blastp', f.binary,
                                 '--review-profiles', str(profile), '--interval-atlas', str(atlas)])
        self.assertEqual(0, code)
        result = json.loads(stream.getvalue())
        self.assertEqual('available', result['family_context_states']['interval_atlas'])
        self.assertEqual(checks.load_profiles(profile), result['settings']['review_profiles'])
    def test_changed_custom_profile_requires_new_revision(self):
        f = self.f
        profile = f.root/'profiles.json'; data = [dict(p) for p in checks.DEFAULT_PROFILES]
        profile.write_text(json.dumps(data))
        f.run_recheck(review_profiles_arg=str(profile))
        data[1]['min_identity'] = 70
        profile.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'changed'):
            f.run_recheck(review_profiles_arg=str(profile), reuse_existing=True)

    def test_old_report_contract_cannot_masquerade_as_new_checks(self):
        self.f.run_recheck()
        p = self.f.out/'second_stage_recheck_manifest.json'; d = json.loads(p.read_text())
        d['contract'] = 'second-stage-fresh-homology-1.0'; p.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError, 'changed'):
            self.f.run_recheck(reuse_existing=True)


class NativeGuidedChecksTests(unittest.TestCase):
    def test_actual_guided_handler_supplies_existing_native_family_artifacts(self):
        f = guided_tests.GuidedIntegrationTests('test_actual_cli_existing_arguments_reach_recheck_and_persist_link')
        f.setUp(); self.addCleanup(f.doCleanups)
        atlas = f.output/'b01_route_b/interval_gene_atlas.tsv'; atlas.parent.mkdir()
        table(atlas, ['gene_id', 'family_id'], [{'gene_id':'g00','family_id':'F'}])
        cp = f.output/'family_copy_discrimination.tsv'
        table(cp, ['candidate_gene_id','family_id','orthology_class'],
              [{'candidate_gene_id':'g00','family_id':'F','orthology_class':'one_to_many'}])
        f.canonical.return_value['family_copy'] = str(cp)
        with mock.patch.object(recheck.shutil, 'which', return_value=None):
            r = f.cli.HANDLERS['guided-profile'](f.args())
        path = Path(r['second_stage_recheck']['review_checks'])
        data = json.loads(path.read_text())
        self.assertEqual('available', data['context_states']['interval_atlas'])
        self.assertEqual('available', data['context_states']['family_copy'])
        self.assertIn('copy_or_orthology_ambiguity_must_remain', data['families']['g00']['review_flags'])
        self.assertEqual(26, len(data['families']))

