from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.supported_member_union import build_supported_member_union
from candidate_qtl.family_candidate_bridge import family_candidate_bridge
from candidate_qtl.relevance_tiers import T2, annotate_relevance
from tests._issue664_helpers import manifest_for


class SupportedRelevanceEvidenceTests(unittest.TestCase):
    def test_real_bridge_schema_recovers_original_hypothesis_at_t2(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HA", "gene_family": "FA", "reference_gene_name": "REF-A",
                "proposed_stable_id": "ACC-A", "member_interval_genes": "",
                "mechanism_statement": "controls trait A", "expected_role": "trait1",
            }]).to_csv(hypotheses, sep="\t", index=False, lineterminator="\n")
            status = root / "status.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HA", "family_status": "concrete_interval_gene_matched",
                "interval_gene_ids": "GA", "source_ids": "PMID:A",
            }]).to_csv(status, sep="\t", index=False, lineterminator="\n")
            bait = root / "bait.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HA", "bait_id": "ACC-A", "function_source_ids": "PMID:A",
                "source_provider": "fixture",
            }]).to_csv(bait, sep="\t", index=False, lineterminator="\n")
            bridge = family_candidate_bridge(
                str(root / "bridge"), route_a_hypotheses_arg=str(hypotheses),
                route_a_family_status_arg=str(status), route_a_bait_panel_arg=str(bait), case_id="case-real",
            )
            bridge_frame = pd.read_csv(bridge["bridge_tsv"], sep="\t", dtype=str).fillna("")
            self.assertNotIn("hypothesis_id", bridge_frame.columns)
            result = build_supported_member_union(
                str(root / "union"), route_a_bridge_arg=bridge["bridge_tsv"],
                route_a_manifest_arg=bridge["manifest"], case_id="case-real",
            )
            row = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).iloc[0]
            evidence = json.loads(row.relevance_evidence_json)
            self.assertIn("controls trait A", {item.get("mechanism_statement") for item in evidence})
            self.assertEqual(T2, annotate_relevance([{"candidate_gene_id": "GA"}], evidence)[0]["relevance_tier"])

    def test_real_bridge_symbol_hypothesis_binds_verified_accession_bait(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HS", "gene_family": "FS", "reference_gene_name": "SYMBOL",
                "proposed_stable_id": "", "member_interval_genes": "",
                "mechanism_statement": "controls symbol trait", "expected_role": "trait1",
            }]).to_csv(hypotheses, sep="\t", index=False, lineterminator="\n")
            status = root / "status.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HS", "family_status": "concrete_interval_gene_matched",
                "interval_gene_ids": "GS", "source_ids": "PMID:S",
            }]).to_csv(status, sep="\t", index=False, lineterminator="\n")
            bait = root / "bait.tsv"
            pd.DataFrame([{
                "hypothesis_id": "HS", "bait_id": "ACC-S", "reference_gene_name": "SYMBOL",
                "function_source_ids": "PMID:S", "source_provider": "fixture",
            }]).to_csv(bait, sep="\t", index=False, lineterminator="\n")
            bridge = family_candidate_bridge(
                str(root / "bridge"), route_a_hypotheses_arg=str(hypotheses),
                route_a_family_status_arg=str(status), route_a_bait_panel_arg=str(bait), case_id="case-symbol",
            )
            result = build_supported_member_union(
                str(root / "union"), route_a_bridge_arg=bridge["bridge_tsv"],
                route_a_manifest_arg=bridge["manifest"], case_id="case-symbol",
            )
            row = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).iloc[0]
            evidence = json.loads(row.relevance_evidence_json)
            self.assertIn("controls symbol trait", {item.get("mechanism_statement") for item in evidence})
            self.assertEqual(T2, annotate_relevance([{"candidate_gene_id": "GS"}], evidence)[0]["relevance_tier"])

    def _claims(self, path: Path, order=(0, 1, 2)) -> Path:
        rows = [
            {
                "claim_id": "C1", "verification_status": "verified",
                "member_interval_genes": "G1;G2", "family_id": "F1",
                "reference_gene": "REF1", "source_id": "PMID:1",
                "mechanism": "cell proliferation", "phenotype_relation": "trait1",
            },
            {
                "claim_id": "C2", "verification_status": "verified",
                "member_interval_genes": "G1", "family_id": "F1",
                "reference_gene": "REF2", "source_id": "PMID:2",
                "mechanism": "organ patterning", "phenotype_relation": "trait2",
            },
            {
                "claim_id": "C3", "verification_status": "not_evaluated",
                "member_interval_genes": "G3", "family_id": "F2",
                "reference_gene": "REF3", "source_id": "PMID:3",
            },
        ]
        pd.DataFrame([rows[index] for index in order]).to_csv(
            path, sep="\t", index=False, lineterminator="\n")
        manifest_for(path, route="b_claims", case_id="case-evidence")
        return path

    def _bridge(self, root: Path) -> Path:
        path = root / "bridge.tsv"
        pd.DataFrame([{
            "gene_id": "GA", "route": "a", "mapping_tier": "tier1_source_backed_unique",
            "provider_status": "complete", "review_status": "verified",
            "family_id": "FA", "reference_accession": "ACC-A",
        }]).to_csv(path, sep="\t", index=False, lineterminator="\n")
        hypothesis = root / "hypotheses.tsv"
        pd.DataFrame([{
            "hypothesis_id": "HA", "gene_family": "FA", "reference_gene_name": "REF-A",
            "proposed_stable_id": "ACC-A", "member_interval_genes": "GA",
            "mechanism_statement": "controls trait1 development",
            "expected_role": "trait1",
        }]).to_csv(hypothesis, sep="\t", index=False, lineterminator="\n")
        bait = root / "bait.tsv"
        pd.DataFrame([{
            "hypothesis_id": "HA", "bait_id": "ACC-A", "function_source_ids": "PMID:A",
            "source_provider": "fixture",
        }]).to_csv(bait, sep="\t", index=False, lineterminator="\n")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest = {
            "schema_version": "1.0", "command": "family-candidate-bridge",
            "contract": "family-candidate-bridge-1.0.0", "case_id": "case-evidence",
            "inputs": {
                "route_a_hypotheses": {"present": True, "path": str(hypothesis),
                    "sha256": hashlib.sha256(hypothesis.read_bytes()).hexdigest(), "strength": "strong"},
                "route_a_bait_panel": {"present": True, "path": str(bait),
                    "sha256": hashlib.sha256(bait.read_bytes()).hexdigest(), "strength": "strong"},
            },
            "outputs": {"family_candidate_bridge": str(path),
                        "family_candidate_bridge_tsv_sha256": digest},
        }
        manifest_path = root / "bridge_manifest.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        return path

    def test_verified_multiclaims_are_retained_and_member_bound(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = self._claims(root / "claims.tsv")
            result = build_supported_member_union(
                str(root / "out"), route_b_claims_arg=str(claims),
                route_b_claims_manifest_arg=str(root / "claims_manifest.json"),
                case_id="case-evidence",
            )
            frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
            evidence = {row.candidate_gene_id: json.loads(row.relevance_evidence_json)
                        for row in frame.itertuples()}
            self.assertEqual({"C1", "C2"}, {row["claim_id"] for row in evidence["G1"]})
            self.assertEqual({"C1"}, {row["claim_id"] for row in evidence["G2"]})
            self.assertNotIn("G3", set(frame.candidate_gene_id))

    def test_route_a_original_hypothesis_and_panel_are_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = self._bridge(root)
            result = build_supported_member_union(
                str(root / "out"), route_a_bridge_arg=str(bridge),
                route_a_manifest_arg=str(root / "bridge_manifest.json"),
                case_id="case-evidence",
            )
            row = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).iloc[0]
            origins = {item["evidence_origin"] for item in json.loads(row.relevance_evidence_json)}
            self.assertEqual({"route_a_hypothesis", "route_a_bait_panel"}, origins)
            annotated = annotate_relevance(
                [{"candidate_gene_id": "GA"}],
                json.loads(row.relevance_evidence_json),
            )
            self.assertEqual(T2, annotated[0]["relevance_tier"])

    def test_route_a_reference_binding_excludes_other_reference_same_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = self._bridge(root)
            hypothesis = root / "hypotheses.tsv"
            pd.DataFrame([
                {"hypothesis_id": "HA", "gene_family": "FA", "proposed_stable_id": "ACC-A",
                 "mechanism_statement": "controls trait A", "expected_role": "trait1"},
                {"hypothesis_id": "HA", "gene_family": "FB", "proposed_stable_id": "ACC-B",
                 "mechanism_statement": "controls unrelated trait B", "expected_role": "trait2"},
            ]).to_csv(hypothesis, sep="\t", index=False, lineterminator="\n")
            bait = root / "bait.tsv"
            pd.DataFrame([
                {"hypothesis_id": "HA", "bait_id": "ACC-A", "function_source_ids": "PMID:A"},
                {"hypothesis_id": "HA", "bait_id": "ACC-B", "function_source_ids": "PMID:B"},
            ]).to_csv(bait, sep="\t", index=False, lineterminator="\n")
            manifest_path = root / "bridge_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            for name, path in (("route_a_hypotheses", hypothesis), ("route_a_bait_panel", bait)):
                manifest["inputs"][name]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            result = build_supported_member_union(
                str(root / "out"), route_a_bridge_arg=str(bridge),
                route_a_manifest_arg=str(manifest_path), case_id="case-evidence",
            )
            row = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).iloc[0]
            evidence = json.loads(row.relevance_evidence_json)
            self.assertEqual({"ACC-A"}, {item.get("proposed_stable_id") for item in evidence if item.get("proposed_stable_id")})
            self.assertNotIn("controls unrelated trait B", {item.get("mechanism_statement") for item in evidence})

    def test_claim_tamper_fails_before_union(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = self._claims(root / "claims.tsv")
            claims.write_text(claims.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash"):
                build_supported_member_union(
                    str(root / "out"), route_b_claims_arg=str(claims),
                    route_b_claims_manifest_arg=str(root / "claims_manifest.json"),
                    case_id="case-evidence",
                )

    def test_claim_row_order_does_not_change_evidence_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self._claims(root / "first.tsv", (0, 1, 2))
            second = self._claims(root / "second.tsv", (1, 0, 2))
            first_result = build_supported_member_union(
                str(root / "out1"), route_b_claims_arg=str(first),
                route_b_claims_manifest_arg=str(root / "first_manifest.json"),
                case_id="case-evidence",
            )
            second_result = build_supported_member_union(
                str(root / "out2"), route_b_claims_arg=str(second),
                route_b_claims_manifest_arg=str(root / "second_manifest.json"),
                case_id="case-evidence",
            )
            self.assertEqual(
                Path(first_result["candidate_set"]).read_bytes(),
                Path(second_result["candidate_set"]).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
