"""Prompt examples must obey the real producer/validator contracts."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl import phenotype_investigation as pi
from candidate_qtl.mechanism_payload_adapter import validate_payload_file
from candidate_qtl.route_c_direct_pipeline import DIRECT_PROPOSAL_COLUMNS, direct_rescue_freeze, render_direct_proposal_prompt, render_harness_native_direct_prompt
from candidate_qtl.route_c_direct_rescue import SPECIES_DISTANCE, EVIDENCE_DIRECTNESS
from tests.test_initial_paper_prompt_examples import _example_after, _filled_hypothesis, _valid_report


class PromptResearchContractTests(unittest.TestCase):
    def test_c_prompt_is_schema_generated_and_examples_pass_real_freeze(self):
        text = render_direct_proposal_prompt()
        template = Path(__file__).resolve().parents[1] / 'workflow/route_c_practical_prompt.txt'
        self.assertEqual(text, template.read_text(encoding='utf-8'))
        self.assertIn('\t'.join(DIRECT_PROPOSAL_COLUMNS), text)
        for value in SPECIES_DISTANCE | EVIDENCE_DIRECTNESS:
            self.assertIn(value, text)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # All canonical enum pairs go through the actual freeze consumer,
            # offline, with no scientific claim acceptance required.
            rows = []
            for distance in sorted(SPECIES_DISTANCE):
                for directness in sorted(EVIDENCE_DIRECTNESS):
                    row = dict.fromkeys(DIRECT_PROPOSAL_COLUMNS, '')
                    row.update(record_id=f'r{len(rows)}', query_phenotype='synthetic phenotype',
                               focal_species='Plant alpha', reported_phenotype='synthetic phenotype',
                               reference_species='Plant alpha' if distance.startswith(('S0_', 'S1_')) else 'Plant beta',
                               reported_gene_symbol='SYNTHETIC', species_distance=distance,
                               evidence_directness=directness, evidence_type='synthetic observation', source_id='PMID:12345678')
                    rows.append(row)
            path = root / 'proposals.tsv'
            pd.DataFrame(rows, columns=DIRECT_PROPOSAL_COLUMNS).to_csv(path, sep='\t', index=False)
            result = direct_rescue_freeze(str(path), str(root / 'freeze'))
            frozen = pd.read_csv(result['frozen'], sep='\t').fillna('')
            self.assertEqual(len(rows), len(frozen))
            self.assertNotIn('verified', set(frozen.source_verification_state))
            for field, value in [('species_distance', 'S1'), ('evidence_directness', 'D0')]:
                invalid = deepcopy(rows[:1]); invalid[0][field] = value
                pd.DataFrame(invalid, columns=DIRECT_PROPOSAL_COLUMNS).to_csv(path, sep='\t', index=False)
                with self.subTest(field=field), self.assertRaises(ValueError):
                    direct_rescue_freeze(str(path), str(root / ('bad_' + field)))

    def test_c_native_prompt_stays_blind_broad_and_uses_existing_schema(self):
        native = render_harness_native_direct_prompt()
        historical = render_direct_proposal_prompt()
        # Do not edit the historical benchmark contract to change normal research.
        self.assertTrue(native.endswith(historical))
        for fragment in (
            "Which genes, classical loci or variants",
            "Do not read Route A/B results, QTL intervals",
            "search/fetch tools actually available",
            "are optional, not prerequisites",
            "not_evaluated rather than completed_empty",
            "Do not pre-filter by gene identity",
            "Preserve multiple source records",
            "Do not suppress overlap",
            "project to the QTL only after freeze",
        ):
            self.assertIn(fragment, native)
        self.assertIn("\t".join(DIRECT_PROPOSAL_COLUMNS), native)

    def test_a_generated_prompt_obeys_mode_and_paper_annotation_location(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for strict in (False, True):
                result = pi.phenotype_investigate(str(root / str(strict)), 'synthetic', 'Plant alpha', 'synthetic phenotype', benchmark_mode=strict)
                text = Path(result['prompt']).read_text(encoding='utf-8')
                self.assertNotIn('- model-species evidence: a named mutant or functional study in a species other', text)
                self.assertNotIn('do not make known focal causal genes the main retrieval target', text)
                self.assertIn('direct focal-species', text)
                self.assertIn('record queries', text)
                self.assertEqual(strict, 'STRICT BENCHMARK POLICY:' in text)
                manifest = json.loads(Path(result['manifest']).read_text(encoding='utf-8'))
                self.assertEqual(strict, manifest['payload_contract']['reference_target_contract']['target_must_be_non_focal'])
                template = _example_after(text, 'Return JSON only:')['hypotheses'][0]
                report = _valid_report()
                report['research_annotations'] = {'read_scope': 'abstract', 'queries': ['synthetic phenotype Plant alpha']}
                hypothesis = _filled_hypothesis(template, 'SYNTHETIC_FAMILY', [report])
                payload = root / 'payload.json'
                payload.write_text(json.dumps({'hypotheses': [hypothesis]}), encoding='utf-8')
                validate_payload_file(payload, None, '')
                hypothesis['research_annotations'] = {'read_scope': 'abstract'}
                payload.write_text(json.dumps({'hypotheses': [hypothesis]}), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'unknown JSON key'):
                    validate_payload_file(payload, None, '')


if __name__ == '__main__':
    unittest.main()
