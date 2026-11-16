"""Source-free synthetic checks of relevance, not biological accuracy."""
import ast
import copy
import hashlib
import itertools
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from candidate_qtl.relevance_tiers import (
    T1, T2, T3, T4, UNRESOLVED, CONTRACT, RELEVANCE_COLUMNS,
    annotate_relevance, classify_relevance,
)


def known():
    return dict(candidate_gene_id="g", focal_species="Species one", reference_species="Species one",
        species_distance="S1_same_species", trait_relation_state="direct_target_trait",
        query_phenotype="seed color", reported_phenotype="seed color", source_id="PMID:1",
        source_reference_gene_id="REF1", source_verification_state="verified",
        gene_identity_state="verified_gene", direct_rescue_eligibility="eligible_direct_record",
        projection_state="inside_current_qtl", projection_basis="exact_focal_id",
        claim_binding_projection_authorized="True", claim_binding_status="verified_complete",
        evidence_directness="D0_causal_perturbation_or_rescue", route="c")


def bridge(**kwargs):
    return dict(candidate_gene_id="g", family_id="F", reference_gene="R", source_id="PMID:2",
        claim_type="mechanistic_process_bridge", biological_function="pigment enzyme",
        mechanism="the reaction supplies the seed pigment", phenotype_relation="seed color changes",
        route="b", **kwargs)


class RelevanceTests(unittest.TestCase):
    def test_same_species_bound_target_trait_only(self):
        self.assertEqual(classify_relevance(known())[0],T1)
        for field in ("claim_binding_status","claim_binding_projection_authorized","source_id",
                      "source_reference_gene_id","focal_species","gene_identity_state"):
            row=known();row[field]=""
            self.assertNotEqual(classify_relevance(row)[0],T1,field)

    def test_related_trait_is_not_equivalent_trait(self):
        row=known();row["trait_relation_state"]="direct_related_trait"
        self.assertEqual(classify_relevance(row)[0],T2)

    def test_other_species_direct_report_is_t2_not_t1(self):
        row=known();row["reference_species"]="Species two";row["species_distance"]="S3_close_clade"
        self.assertEqual(classify_relevance(row)[0],T2)

    def test_d2_existing_association_vocabulary(self):
        row=known();row.update(reference_species="Species two",evidence_directness="D2_direct_trait_gene_association")
        self.assertEqual(classify_relevance(row)[0],T2)

    def test_verified_homology_not_t1_or_t3(self):
        row=dict(family_id="kinase",support_status="verified_semantic_support")
        self.assertEqual(classify_relevance(row)[0],T4)

    def test_concrete_bridge_and_review_seed(self):
        self.assertEqual(classify_relevance(bridge())[0],T2)
        row=bridge();row["claim_type"]="mechanistic_review_seed"
        self.assertEqual(classify_relevance(row)[0],T3)
        self.assertIn("seed color changes",classify_relevance(row)[2])

    def test_function_alone_is_not_trait_relevance(self):
        self.assertEqual(classify_relevance({"biological_function":"kinase"})[0],T4)
        self.assertEqual(classify_relevance({"family_id":"F","phenotype_relation":"seed color"})[0],T4)
        self.assertEqual(classify_relevance({})[0],UNRESOLVED)

    def test_a_mechanism_uses_same_tier(self):
        row=dict(candidate_gene_id="g",reference_gene="R",route="a",
                 mechanism_statement="enzyme makes seed pigment",expected_role="seed color changes")
        self.assertEqual(classify_relevance(row)[0],T2)

    def test_mechanism_field_aliases_have_equal_untyped_bridge_tier(self):
        base=dict(candidate_gene_id="g",family_id="F",reference_gene="R",
                  phenotype_relation="seed color changes")
        a={**base,"mechanism_statement":"enzyme makes seed pigment"}
        b={**base,"mechanistic_link":"enzyme makes seed pigment"}
        c={**base,"mechanism":"enzyme makes seed pigment"}
        self.assertEqual(classify_relevance(a)[0],T2)
        self.assertEqual(classify_relevance(b)[0],T2)
        self.assertEqual(classify_relevance(c)[0],T2)

    def test_weak_or_unknown_claim_type_does_not_gain_t2_from_mechanism_alias(self):
        for claim_type in ("mechanistic_review_seed", "unknown"):
            row=bridge();row["claim_type"]=claim_type
            self.assertEqual(classify_relevance(row)[0],T3,claim_type)

    def test_mechanism_bridge_remains_bound_to_its_member_and_source(self):
        claim=bridge();claim["candidate_gene_id"]="g1"
        rows=annotate_relevance([{"candidate_gene_id":"g1"},{"candidate_gene_id":"g2"}],[claim])
        by_id={row["candidate_gene_id"]:row for row in rows}
        self.assertEqual(by_id["g1"]["relevance_tier"],T2)
        self.assertEqual(by_id["g1"]["source_ids"],"PMID:2")
        self.assertEqual(by_id["g2"]["relevance_tier"],UNRESOLVED)
        self.assertEqual(json.loads(by_id["g1"]["relevance_evidence_json"])[0]["candidate_gene_id"],"g1")

    def test_copy_ambiguity_does_not_demote(self):
        a=bridge(mapping_tier="tier1_source_backed_unique",copy_state="preferred")
        b=bridge(mapping_tier="tier3_paralog_ambiguous",copy_state="co_preferred")
        self.assertEqual(classify_relevance(a)[0],classify_relevance(b)[0])
        result=annotate_relevance([dict(candidate_gene_id="g",support_status="unverified_route_b_review")],[b])[0]
        self.assertEqual(result["copy_state"],"co_preferred")
        self.assertEqual(result["mapping_tier"],"tier3_paralog_ambiguous")
        self.assertEqual(result["support_state"],"review")
        self.assertEqual(result["relevance_tier"],T2)

    def test_failed_verification_not_negative(self):
        a=bridge(verification_status="verified"); b=bridge(verification_status="provider_failed")
        self.assertEqual(classify_relevance(a)[0],classify_relevance(b)[0])
        b["verification_status"]="contradicted"
        self.assertEqual(classify_relevance(b)[0],T4)

    def test_no_cross_claim_stitching_into_t1(self):
        a=known();b=known();a["source_id"]="";b["gene_identity_state"]=""
        row=annotate_relevance([{"candidate_gene_id":"g"}],[a,b])[0]
        self.assertNotEqual(row["relevance_tier"],T1)

    def test_legacy_tier_is_retained_not_renamed(self):
        a=dict(candidate_gene_id="g",family_id="F",support_tier="tier1_verified_support")
        row=annotate_relevance([a])[0]
        self.assertEqual(row["legacy_support_tier"],"tier1_verified_support")
        self.assertEqual(row["support_tier"],"tier1_verified_support")
        self.assertEqual(row["support_state"],"verified")
        self.assertEqual(row["relevance_tier"],T4)
        self.assertEqual(row["relevance_migration_state"],"legacy_identity_only_no_trait_proof")

    def test_all_paralog_members_retained(self):
        claim=bridge();del claim["candidate_gene_id"];claim["member_interval_genes"]="g1;g2;g3"
        rows=annotate_relevance([dict(candidate_gene_id=f"g{i}") for i in (1,2,3)],[claim])
        self.assertEqual(len(rows),3)
        self.assertTrue(all(r["relevance_tier"]==T2 for r in rows))

    def test_unbound_claim_rejected(self):
        with self.assertRaisesRegex(ValueError,"bound interval"):
            annotate_relevance([{"candidate_gene_id":"g"}],[{"family_id":"F"}])

    def test_invalid_or_duplicate_roster_rejected(self):
        for rows in ([{}],[{"gene_id":"g"},{"gene_id":"g"}]):
            with self.assertRaisesRegex(ValueError,"unique"):
                annotate_relevance(rows)

    def test_membership_and_input_immutable(self):
        rows=[dict(candidate_gene_id="g"),dict(candidate_gene_id="z")]
        evidence=[bridge(),{**bridge(),"candidate_gene_id":"outside"}]
        before=copy.deepcopy((rows,evidence))
        out=annotate_relevance(rows,evidence)
        self.assertEqual((rows,evidence),before)
        self.assertEqual({r["candidate_gene_id"] for r in out},{"g","z"})

    def test_permutation_and_duplicate_invariance(self):
        claims=[bridge(),known(),{**bridge(),"source_id":"PMID:3"}]
        expected=annotate_relevance([{"candidate_gene_id":"g"}],claims)
        for order in itertools.permutations(claims):
            self.assertEqual(annotate_relevance([{"candidate_gene_id":"g"}],list(order)+[claims[0]]),expected)
        self.assertEqual(len(json.loads(expected[0]["relevance_evidence_json"])),3)

    def test_actual_b_freeze_retains_relevance_and_claim_pairing(self):
        # Actual producer function, isolated from unrelated optional dependencies.
        source=Path(__file__).resolve().parents[1]/"src/candidate_qtl/route_b_semantic_candidates.py"
        tree=ast.parse(source.read_text())
        tree.body=[node for node in tree.body if not (isinstance(node,ast.ImportFrom) and node.level)]
        ns={"runtime_path":Path,"normalize_semantic_claims":lambda frame:frame,
            "RELEVANCE_CONTRACT":CONTRACT,"RELEVANCE_COLUMNS":RELEVANCE_COLUMNS,
            "annotate_relevance":annotate_relevance}
        exec(compile(tree,str(source),"exec"),ns)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);bundle=root/"bundle.json"
            bundle.write_text(json.dumps({"route":"interval_first","units":[{
                "unit_id":"u","member_interval_genes":["g"],"identity":{"family_ids":["F"]}}]}))
            claims=[]
            for i,status in enumerate(("verified","provider_failed")):
                row=bridge(verification_status=status)
                row.pop("candidate_gene_id");row.update(member_interval_genes="g",route="interval_first",
                    claim_id=f"c{i}",semantic_unit_id="u",evidence_scope="family",verification_reason="synthetic")
                if i: row["claim_type"]="mechanistic_review_seed"
                claims.append(row)
            path=root/"claims.tsv";pd.DataFrame(claims).to_csv(path,sep="\t",index=False)
            transfer=root/"semantic_claim_transfer_verification.tsv"
            pd.DataFrame([
                {"claim_id":"c0","supported_member_genes":"g"},
                {"claim_id":"c1","supported_member_genes":""},
            ]).to_csv(transfer,sep="\t",index=False)
            (root/"semantic_claim_transfer_verification_manifest.json").write_text(json.dumps({
                "stage":"route_b_semantic_transfer_verification","contract":"route-b-semantic-transfer-verification-1.0.0",
                "case_id":"synthetic","artifact_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                "outputs":{"transfer_verification":str(transfer),
                           "transfer_verification_sha256":hashlib.sha256(transfer.read_bytes()).hexdigest()},
                "ancestors":{"bundle":{"path":str(bundle),"sha256":hashlib.sha256(bundle.read_bytes()).hexdigest()}}}))
            result=ns["freeze_route_b_semantic_candidates"](str(path),str(root/"out"),case_id="synthetic")
            supported=pd.read_csv(result["candidates"],sep="\t").iloc[0]
            review=pd.read_csv(result["review_candidates"],sep="\t").iloc[0]
            self.assertEqual(supported["relevance_tier"],T2)
            self.assertEqual(review["relevance_tier"],T3)
            self.assertEqual(supported["support_state"],"verified")
            self.assertEqual(review["support_state"],"review")
            self.assertEqual({item["claim_id"] for item in json.loads(supported["relevance_evidence_json"])},{"c0"})
            self.assertEqual({item["claim_id"] for item in json.loads(review["relevance_evidence_json"])},{"c1"})
            manifest=json.loads(Path(result["manifest"]).read_text())
            self.assertEqual(manifest["relevance_contract"],CONTRACT)
            self.assertEqual(manifest["artifact_sha256"],hashlib.sha256(Path(result["candidates"]).read_bytes()).hexdigest())

if __name__=="__main__":
    unittest.main()
