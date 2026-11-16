"""Original-judgment Tier policy and legacy replay; fixed synthetic evidence."""
import unittest

from candidate_qtl.relevance_tiers import NEW_TIERS, T2, annotate_relevance, classify_relevance
from candidate_qtl.boundary_tier_review import EVIDENCE_SCOPED_POLICY, partition_relevance_tiers


def observation(tier="T2", **changes):
    row = {"candidate_gene_id": "G1", "candidate_admission_policy": "initial_bibliography_v1",
           "initial_tier": tier, "initial_tier_reason": "original reference observation and conditional bridge",
           "source_id": "PMID:1234567", "reference_gene": "REF1", "biological_function": "meristem maintenance",
           "mechanism": "receptor signaling", "phenotype_relation": "could alter meristem size",
           "reported_phenotype": "enlarged meristem", "reported_process": "meristem cell proliferation",
           "phenotype_comparison": "same_cellular_process", "reference_species": "Plant alpha",
           "focal_species": "Plant alpha"}
    row.update(changes)
    return row


class InitialTierTests(unittest.TestCase):
    def test_report_function_annotation_identity_and_unknown_have_distinct_meanings(self):
        examples = [observation("T1", phenotype_comparison="same_equivalent"), observation("T2"),
                    observation("T3", phenotype_comparison="shared_process", reported_organ="root"),
                    observation("T4", reported_phenotype=""), observation("T5", reported_phenotype=""),
                    observation("T6", biological_function="kinase", phenotype_relation=""),
                    {"candidate_admission_policy": "initial_bibliography_v1"}]
        self.assertEqual(list(NEW_TIERS), [classify_relevance(row)[0] for row in examples])

    def test_real_source_does_not_turn_function_or_unresolved_report_into_t2(self):
        for row in [observation("T2", reported_phenotype=""), observation("T2", phenotype_comparison="different"),
                    observation("T1", reference_species="Plant beta", phenotype_comparison="same_equivalent"),
                    observation("", reported_phenotype="")]:
            self.assertEqual(NEW_TIERS[3], classify_relevance(row)[0])
        self.assertEqual(NEW_TIERS[5], classify_relevance(observation(biological_function="transcription factor"))[0])

    def test_copies_and_initial_reason_survive_annotation_and_boundary_consumer(self):
        claims = [observation("T2", candidate_gene_id=f"G{i}") for i in range(30)]
        rows = annotate_relevance([{"candidate_gene_id": f"G{i}", "family_id": "F"} for i in range(1000)], claims)
        self.assertEqual(30, sum(row["relevance_tier"] == NEW_TIERS[1] for row in rows))
        self.assertEqual(970, sum(row["relevance_tier"] == NEW_TIERS[5] for row in rows))
        plan = partition_relevance_tiers(rows, policy_version=EVIDENCE_SCOPED_POLICY,
                                        eligible_gene_ids={f"G{i}" for i in range(30)})
        self.assertEqual(30, len(plan["boundary"]))
        self.assertEqual(970, len(plan["lower"]))
        self.assertTrue(all(row["initial_tier_reason"] for row in plan["boundary"]))

    def test_original_full_labels_keep_their_reported_function_meanings(self):
        examples = [(NEW_TIERS[0], {"phenotype_comparison": "same_equivalent"}),
                    (NEW_TIERS[1], {}),
                    (NEW_TIERS[2], {"phenotype_comparison": "shared_process"}),
                    (NEW_TIERS[4], {"reported_phenotype": ""})]
        for label, fields in examples:
            with self.subTest(label=label):
                row = annotate_relevance([{"candidate_gene_id": "G1"}], [observation(label, **fields)])[0]
                self.assertEqual(label, row["relevance_tier"])
                self.assertEqual(label, row["initial_tier"])
        self.assertEqual(NEW_TIERS[3], classify_relevance(observation("T2_trait_or_mechanism_linked"))[0])

    def test_old_unversioned_evidence_keeps_old_tier_meaning(self):
        row = observation()
        del row["candidate_admission_policy"]
        self.assertEqual(T2, classify_relevance(row)[0])
        annotated = annotate_relevance([{"candidate_gene_id": "G1"}], [row])[0]
        self.assertEqual("legacy_relevance_v1", annotated["relevance_policy"])
        self.assertEqual(T2, annotated["relevance_tier"])

    def test_rehydrated_t1_application_does_not_rescue_other_reference_copy(self):
        import json
        claim = observation("T1", phenotype_comparison="same_equivalent", candidate_gene_id="", member_interval_genes="G1;G2")
        frozen = {"candidate_gene_id": "G1", "candidate_admission_policy": "initial_bibliography_v1",
                  "relevance_evidence_json": json.dumps([claim])}
        first = annotate_relevance([{"candidate_gene_id": "G1"}], [frozen])
        rows = annotate_relevance([{"candidate_gene_id": "G1"}, {"candidate_gene_id": "G2", "family_id": "F"}], first)
        self.assertEqual(NEW_TIERS[0], rows[0]["relevance_tier"])
        self.assertEqual(NEW_TIERS[5], rows[1]["relevance_tier"])


if __name__ == "__main__":
    unittest.main()
