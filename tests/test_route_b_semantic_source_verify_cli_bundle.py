"""CLI reachability for Route-B frozen curated-function source verification (#645)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


def _write_proposal(path: Path) -> None:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": "claim-cli-bundle",
        "route": "interval_first",
        "semantic_unit_id": "unit-cli-bundle",
        "family_id": "EV1",
        "orthogroup_id": "OG1",
        "subfamily_or_clade": "SUB1",
        "member_interval_genes": "G1;G2",
        "claim_type": "direct_functional_bridge",
        "biological_function": "circadian clock regulation",
        "mechanism": "source-backed reference function",
        "phenotype_relation": "flowering-time phenotype",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "reference_gene",
        "source_id": "UniProt:REF1",
        "alternative_hypothesis": "related family copy",
        "uncertainty_note": "focal transfer remains unresolved",
        "provider": "fixture",
        "model": "fixture-model",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": "proposed",
        "verification_reason": "source_id_present_in_frozen_bundle; semantic_relation_verification_pending",
        "rank_contribution": "none",
    })
    pd.DataFrame([row], columns=CLAIM_COLUMNS).to_csv(
        path, sep="\t", index=False, lineterminator="\n"
    )


def _write_bundle(path: Path) -> None:
    path.write_text(json.dumps({
        "schema_version": "1.0",
        "case_id": "case645",
        "route": "interval_first",
        "units": [{
            "unit_id": "unit-cli-bundle",
            "unit_type": "family",
            "unit_anchor": "EV1",
            "member_interval_genes": ["G1", "G2"],
            "identity": {
                "family_ids": ["EV1"],
                "orthogroup_ids": ["OG1"],
                "subfamily_ids": ["SUB1"],
            },
            "genes": [{
                "gene_id": "G1",
                "reference_gene": "REF1",
                "orthology_reference_gene": "REF1",
                "orthology_relationship": "best_hit_homology",
                "lookup_status": "source_backed",
                "reference_source_id": "UniProt:REF1",
                "reference_function": "circadian clock regulation",
                "protein_family": "EV1 family",
                "go_terms": "GO:0042752",
                "source_database": "UniProtKB",
                "source_url": "https://rest.uniprot.org/uniprotkb/REF1",
            }],
            "source_provenance": [{
                "gene_id": "G1",
                "reference_gene": "REF1",
                "source_ids": ["UniProt:REF1"],
                "source_database": "UniProtKB",
                "source_url": "https://rest.uniprot.org/uniprotkb/REF1",
                "provenance_detail": "frozen_uniprot.tsv",
                "retrieved_at": "2026-08-27T00:00:00Z",
            }],
            "stable_source_ids": ["UniProt:REF1"],
            "rank_contribution": "none",
        }],
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


class RouteBSemanticSourceVerifyCliBundleTests(unittest.TestCase):
    def test_cli_bundle_enables_curated_function_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            bundle = root / "bundle.json"
            out = root / "verify"
            _write_proposal(proposals)
            _write_bundle(bundle)

            self.assertEqual(0, main([
                "route-b-semantic-source-verify",
                "--proposals", str(proposals),
                "--bundle", str(bundle),
                "--out-dir", str(out),
                "--case-id", "case645",
            ]))

            frame = pd.read_csv(
                out / "semantic_claim_source_verification.tsv",
                sep="\t", dtype=str,
            ).fillna("")
            self.assertEqual(1, len(frame))
            row = frame.iloc[0]
            self.assertEqual("supported_curated_function", row["source_relation_status"])
            self.assertEqual("source_backed_curated_function", row["claim_support_status"])
            self.assertEqual("frozen_route_b_bundle", row["retrieval_scope"])
            self.assertEqual("True", row["transfer_verification_required"])
            self.assertEqual("proposed", row["final_semantic_verification_status"])
            self.assertEqual("none", row["rank_contribution"])


if __name__ == "__main__":
    unittest.main()
