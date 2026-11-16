import json
import tempfile
from pathlib import Path
import pandas as pd
import test_route_b_semantic_import as import_fixtures
import test_route_b_semantic_transfer_verify as transfer_fixtures
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates


import unittest


class RouteBDeclaredMemberSubsetTests(unittest.TestCase):
    def test_import_preserves_declared_subset_and_legacy_omission(self):
        helper = import_fixtures.RouteBSemanticImportTests()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = helper._bundle(root)
            bundle = json.loads(Path(built['bundle']).read_text())
            for unit in bundle['units']:
                unit['member_interval_genes'].append('G3')
                unit['genes'].append({**unit['genes'][0], 'gene_id': 'G3'})
            Path(built['bundle']).write_text(json.dumps(bundle))
            for name, extra, expected in [('subset', {'member_interval_genes': ['G2', 'G1']}, 'G1;G2'),
                                           ('string', {'member_interval_genes': 'G2;G1'}, 'G1;G2'),
                                           ('legacy', {}, 'G1;G2;G3')]:
                response = root / (name + '.json')
                helper._response(Path(built['prompts']), [{'claim_type': 'mutant_phenotype', 'reference_gene': 'REF1',
                    'source_ids': ['PMID:123'], 'evidence_scope': 'reference_gene', **extra}], response)
                result = import_route_b_semantic_proposals(built['bundle'], built['prompts'], str(response), str(root / name))
                assert set(pd.read_csv(result['proposals'], sep='\t')['member_interval_genes']) == {expected}
            for i, invalid in enumerate([[], '', None, ['OUTSIDE'], ['G1', 1]]):
                response = root / f'bad{i}.json'
                helper._response(Path(built['prompts']), [{'claim_type': 'mutant_phenotype', 'member_interval_genes': invalid}], response)
                result = import_route_b_semantic_proposals(built['bundle'], built['prompts'], str(response), str(root / f'bad{i}'))
                assert result['claim_row_count'] == 0
                assert result['rejected_claim_row_count'] == 1


    def test_transfer_and_freeze_keep_only_declared_members_for_support_and_review(self):
        helper = transfer_fixtures.RouteBSemanticTransferVerifyTests()
        for claim_type in ['mutant_phenotype', 'mechanistic_review_seed']:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result = helper._run(root, transfer_fixtures._claim(claim_type=claim_type, member_interval_genes='G1;G2'),
                                     [('G1', '1:many'), ('G2', '1:many'), ('G3', '1:many')], 'supported_mechanistic')
                frozen = freeze_route_b_semantic_candidates(result['claims'], str(root / 'frozen'), case_id='case403')
                support = pd.read_csv(frozen['candidates'], sep='\t')
                review = pd.read_csv(frozen['review_candidates'], sep='\t')
                target = support if claim_type == 'mutant_phenotype' else review
                assert set(target['candidate_gene_id']) == {'G1', 'G2'}
                other = review if claim_type == 'mutant_phenotype' else support
                assert other.empty


    def test_transfer_rejects_outside_member_without_support_or_review(self):
        helper = transfer_fixtures.RouteBSemanticTransferVerifyTests()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = helper._run(root, transfer_fixtures._claim(member_interval_genes='G1;OUTSIDE'), [('G1', '1:many'), ('G2', '1:many')], 'supported_mechanistic')
            frozen = freeze_route_b_semantic_candidates(result['claims'], str(root / 'frozen'), case_id='case403')
            assert pd.read_csv(frozen['candidates'], sep='\t').empty
            assert pd.read_csv(frozen['review_candidates'], sep='\t').empty

    def test_support_freeze_honors_hash_bound_transfer_subset_and_rejects_escape(self):
        import hashlib
        helper = transfer_fixtures.RouteBSemanticTransferVerifyTests()
        for declared, valid in [('G1', True), ('G3', False)]:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result = helper._run(root, transfer_fixtures._claim(), [('G1', '1:many'), ('G2', '1:many')], 'supported_mechanistic')
                transfer_path = Path(result['transfer_verification'])
                table = pd.read_csv(transfer_path, sep='\t').fillna('')
                table['supported_member_genes'] = declared
                table.to_csv(transfer_path, sep='\t', index=False)
                manifest_path = Path(result['claims']).with_name('semantic_claim_transfer_verification_manifest.json')
                manifest = json.loads(manifest_path.read_text())
                manifest['outputs']['transfer_verification_sha256'] = hashlib.sha256(transfer_path.read_bytes()).hexdigest()
                manifest_path.write_text(json.dumps(manifest))
                if valid:
                    frozen = freeze_route_b_semantic_candidates(result['claims'], str(root / 'frozen'), case_id='case403')
                    assert set(pd.read_csv(frozen['candidates'], sep='\t')['candidate_gene_id']) == {'G1'}
                else:
                    with self.assertRaisesRegex(ValueError, 'non-empty claim subset'):
                        freeze_route_b_semantic_candidates(result['claims'], str(root / 'frozen'), case_id='case403')



    def test_supported_relevance_does_not_leak_between_claim_subsets(self):
        import hashlib
        helper = transfer_fixtures.RouteBSemanticTransferVerifyTests()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = helper._run(root, transfer_fixtures._claim(), [('G1', '1:many'), ('G2', '1:many')], 'supported_mechanistic')
            path = Path(result['claims'])
            claims = pd.read_csv(path, sep='\t').fillna('')
            claims['initial_tier'] = 'T2'
            claims['initial_tier_reason'] = 'reported phenotype'
            second = claims.iloc[0].copy()
            second['claim_id'] = 'claim-2'
            second['initial_tier'] = 'T5'
            second['initial_tier_reason'] = 'annotation hypothesis'
            second['source_id'] = 'PMID:7654321'
            claims = pd.concat([claims, second.to_frame().T], ignore_index=True)
            claims.to_csv(path, sep='\t', index=False)
            transfer_path = Path(result['transfer_verification'])
            transfer = pd.read_csv(transfer_path, sep='\t').fillna('')
            second_transfer = transfer.iloc[0].copy()
            second_transfer['claim_id'] = 'claim-2'
            second_transfer['supported_member_genes'] = 'G2'
            transfer['supported_member_genes'] = 'G1'
            pd.concat([transfer, second_transfer.to_frame().T], ignore_index=True).to_csv(transfer_path, sep='\t', index=False)
            manifest_path = path.with_name('semantic_claim_transfer_verification_manifest.json')
            manifest = json.loads(manifest_path.read_text())
            manifest['artifact_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest['outputs']['transfer_verification_sha256'] = hashlib.sha256(transfer_path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            frozen = freeze_route_b_semantic_candidates(str(path), str(root / 'frozen'), case_id='case403')
            frame = pd.read_csv(frozen['candidates'], sep='\t').set_index('candidate_gene_id')
            assert frame.loc['G1', 'initial_tier'] == 'T2'
            assert frame.loc['G2', 'initial_tier'] == 'T5'
            assert frame.loc['G1', 'claim_ids'] == 'claim-1'
            assert frame.loc['G2', 'claim_ids'] == 'claim-2'
            assert 'claim-1' not in frame.loc['G2', 'relevance_evidence_json']
            assert 'claim-2' in frame.loc['G2', 'relevance_evidence_json']

    def test_current_prompt_documents_subset_without_changing_legacy_policy(self):
        from candidate_qtl.route_b_semantic_bundle import SYSTEM_PROMPT, CONTENT_SYSTEM_PROMPT
        from candidate_qtl.paper_reports import LITERATURE_SEARCH_GUIDANCE
        assert SYSTEM_PROMPT.count(LITERATURE_SEARCH_GUIDANCE) == 1
        assert 'Optional claim.member_interval_genes' in SYSTEM_PROMPT
        assert 'omit it only when the claim applies to the whole frozen unit' in SYSTEM_PROMPT
        assert 'Optional claim.member_interval_genes' not in CONTENT_SYSTEM_PROMPT
        assert LITERATURE_SEARCH_GUIDANCE not in CONTENT_SYSTEM_PROMPT

    def test_explicit_subset_does_not_erase_competing_strong_copy(self):
        helper = transfer_fixtures.RouteBSemanticTransferVerifyTests()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = helper._run(root, transfer_fixtures._claim(evidence_scope='interval_gene', member_interval_genes='G1'),
                                 [('G1', '1:1'), ('G2', '1:1')], 'supported_mechanistic')
            final = pd.read_csv(result['claims'], sep='\t')
            transfer = pd.read_csv(result['transfer_verification'], sep='\t')
            assert final.iloc[0]['verification_status'] == 'not_evaluated'
            assert transfer.iloc[0]['transfer_status'] == 'family_supported_copy_unresolved'
            assert transfer.iloc[0]['supported_member_genes'] == 'G1'
