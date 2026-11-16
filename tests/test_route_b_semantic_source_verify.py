from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import source_verify as source_verify_module
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


def _proposal(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": "claim-1",
        "route": "interval_first",
        "semantic_unit_id": "bsem-source-fixture",
        "family_id": "EV1",
        "orthogroup_id": "OG1",
        "subfamily_or_clade": "SUB1",
        "member_interval_genes": "G1;G2",
        "claim_type": "mutant_phenotype",
        "biological_function": "meristem maintenance",
        "mechanism": "receptor signaling",
        "phenotype_relation": "fasciated shoot meristem",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "reference_gene",
        "source_id": "PMID:1234567",
        "source_url": "",
        "alternative_hypothesis": "related receptor copy",
        "uncertainty_note": "transfer to focal copy remains unresolved",
        "provider": "fixture",
        "model": "fixture-model",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": "proposed",
        "verification_reason": "source_id_present_in_frozen_bundle; semantic_relation_verification_pending",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


def _write_proposals(path: Path, rows: list[dict[str, str]]) -> None:
    pd.DataFrame(rows, columns=CLAIM_COLUMNS).to_csv(path, sep="\t", index=False, lineterminator="\n")


def _cache(root: Path, *, relation: str = "mechanistic") -> Path:
    cache = root / "pmc_cache"
    cache.mkdir()
    if relation == "mechanistic":
        title = "REF1 controls shoot meristem development"
        abstract = "Loss-of-function knockout mutants of REF1 display a fasciated shoot meristem phenotype."
    elif relation == "associative":
        title = "REF1 expression in shoot meristems"
        abstract = "RNA-seq expression analysis found REF1 transcript abundance associated with meristem development."
    else:
        title = "REF1 in plant development"
        abstract = "REF1 is discussed as part of background knowledge about plant development."
    (cache / "med_1234567.json").write_text(json.dumps({
        "cache_key": "pmid_1234567",
        "id": "1234567",
        "pmid": "1234567",
        "title": title,
        "abstractText": abstract,
        "journalTitle": "Plant Journal",
        "pubYear": "2026",
        "doi": "10.1000/ref1",
        "authorString": "Example A, Example B",
    }), encoding="utf-8")
    return cache


def _bundle(root: Path, *, function: str = "circadian clock regulation", source_id: str = "UniProt:REF1") -> Path:
    path = root / "bundle.json"
    path.write_text(json.dumps({
        "schema_version": "1.0",
        "case_id": "case403",
        "route": "interval_first",
        "units": [{
            "unit_id": "bsem-source-fixture",
            "unit_type": "family",
            "unit_anchor": "EV1",
            "member_interval_genes": ["G1", "G2"],
            "identity": {"family_ids": ["EV1"], "orthogroup_ids": ["OG1"], "subfamily_ids": ["SUB1"]},
            "genes": [{
                "gene_id": "G1",
                "reference_gene": "REF1",
                "orthology_reference_gene": "REF1",
                "orthology_relationship": "best_hit_homology",
                "lookup_status": "source_backed",
                "reference_source_id": source_id,
                "reference_function": function,
                "protein_family": "EV1 family" if function else "",
                "go_terms": "GO:0042752" if function else "",
                "source_database": "UniProtKB",
                "source_url": "https://rest.uniprot.org/uniprotkb/REF1",
            }],
            "source_provenance": [{
                "gene_id": "G1",
                "reference_gene": "REF1",
                "source_ids": [source_id],
                "source_database": "UniProtKB",
                "source_url": "https://rest.uniprot.org/uniprotkb/REF1",
                "provenance_detail": "frozen_uniprot.tsv",
                "retrieved_at": "2026-08-21T00:00:00Z",
            }],
            "stable_source_ids": [source_id],
            "rank_contribution": "none",
        }],
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


class RouteBSemanticSourceVerifyTests(unittest.TestCase):
    def test_frozen_curated_database_function_is_a_weaker_explicit_support_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(
                source_id="UniProt:REF1", claim_type="direct_functional_bridge"
            )])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), bundle_arg=str(_bundle(root)), case_id="case403"
            )
            row = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("supported_curated_function", row["source_relation_status"])
            self.assertEqual("source_backed_curated_function", row["claim_support_status"])
            self.assertEqual("curated_function_annotation", row["evidence_relation"])
            self.assertEqual("frozen_route_b_bundle", row["retrieval_scope"])
            self.assertEqual("not_literature_source", row["citation_verification_status"])
            self.assertEqual("False", row["citation_verified"])
            self.assertEqual("True", row["transfer_verification_required"])
            self.assertEqual("proposed", row["final_semantic_verification_status"])

    def test_database_accession_without_non_generic_frozen_function_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(
                source_id="UniProt:REF1", claim_type="direct_functional_bridge"
            )])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), bundle_arg=str(_bundle(root, function=""))
            )
            row = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("not_evaluated", row["source_relation_status"])
            self.assertEqual("not_literature_source", row["citation_verification_status"])
            self.assertEqual("False", row["transfer_verification_required"])

    def test_database_accession_with_mismatched_frozen_source_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(
                source_id="UniProt:OTHER", claim_type="direct_functional_bridge"
            )])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"),
                bundle_arg=str(_bundle(root, source_id="UniProt:REF1")),
            )
            row = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("not_evaluated", row["source_relation_status"])
            self.assertEqual("False", row["transfer_verification_required"])

    def test_generic_functional_relevance_proposal_cannot_use_curated_membership_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(
                source_id="UniProt:REF1", claim_type="functional_relevance_proposal"
            )])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), bundle_arg=str(_bundle(root))
            )
            row = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("not_evaluated", row["source_relation_status"])
            self.assertEqual("False", row["transfer_verification_required"])

    def test_direct_genetic_source_support_does_not_auto_verify_cross_species_transfer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal()])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), cache_arg=str(_cache(root)), case_id="case403"
            )
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            row = frame.iloc[0]
            self.assertEqual("supported_mechanistic", row["source_relation_status"])
            self.assertEqual("supported_direct_genetic_test", row["claim_support_status"])
            self.assertEqual("True", row["citation_verified"])
            self.assertEqual("title_abstract", row["retrieval_scope"])
            self.assertEqual("proposed", row["final_semantic_verification_status"])
            self.assertEqual("True", row["transfer_verification_required"])
            self.assertEqual("none", row["rank_contribution"])
            self.assertEqual(0, result["final_verified_claim_count"])
            self.assertEqual("cache_hit", row["cache_investigation_state"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(
                {"cache_hit": 1, "new_investigation": 0, "provider_failure": 0, "not_applicable": 0},
                manifest["cache_investigation_state_counts"],
            )
            self.assertEqual(1, manifest["cache_accounting"]["cache_hit_count"])

    def test_expression_association_stays_explicitly_weaker_than_mechanistic_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal()])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), cache_arg=str(_cache(root, relation="associative"))
            )
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("supported_associative", frame.iloc[0]["source_relation_status"])
            self.assertEqual("associative_support_only", frame.iloc[0]["claim_support_status"])
            self.assertEqual("proposed", frame.iloc[0]["final_semantic_verification_status"])

    def test_offline_without_cache_is_not_evaluated_not_unsupported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal()])
            result = verify_route_b_semantic_sources(str(proposals), str(root / "verify"))
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", frame.iloc[0]["source_relation_status"])
            self.assertEqual("review_required", frame.iloc[0]["citation_verification_status"])
            self.assertEqual("proposed", frame.iloc[0]["final_semantic_verification_status"])
            # No cache and no provider allowed: the claim was never investigated.
            # That must be recorded as provider_failure, and provider_failure
            # must never be converted into unsupported/negative evidence.
            self.assertEqual("provider_failure", frame.iloc[0]["cache_investigation_state"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["cache_investigation_state_counts"]["provider_failure"])
            self.assertEqual(1, manifest["cache_accounting"]["provider_failure_count"])
            self.assertNotEqual("unsupported", frame.iloc[0]["source_relation_status"])

    def test_proposal_already_rejected_by_frozen_bundle_import_remains_unsupported_without_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(
                verification_status="unsupported",
                verification_reason="source_id_not_in_frozen_bundle:PMID:9999999",
                source_id="PMID:9999999",
            )])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), cache_arg=str(_cache(root))
            )
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("unsupported", frame.iloc[0]["source_relation_status"])
            self.assertEqual("unsupported", frame.iloc[0]["final_semantic_verification_status"])
            self.assertEqual("", frame.iloc[0]["source_provider"])
            # Rejected before any lookup ran; this is not a provider outage.
            self.assertEqual("not_applicable", frame.iloc[0]["cache_investigation_state"])

    def test_empty_proposal_table_is_valid_abstention(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            pd.DataFrame(columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False, lineterminator="\n")
            result = verify_route_b_semantic_sources(str(proposals), str(root / "verify"))
            self.assertEqual(0, result["claim_count"])
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(frame))
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            # Schema-checkable even for the empty-abstention path: keys never
            # disappear just because a bucket is zero this run.
            self.assertEqual(
                {"cache_hit": 0, "new_investigation": 0, "provider_failure": 0, "not_applicable": 0},
                manifest["cache_investigation_state_counts"],
            )

    def test_cache_hit_and_provider_failure_coexist_in_one_canonical_run(self) -> None:
        """Two claims in the same call: one resolves from cache, one has no
        provider available. Both cache-investigation states must be visible
        together from the real entrypoint, not asserted separately from
        synthetic rows."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [
                _proposal(claim_id="claim-cached", source_id="PMID:1234567"),
                _proposal(claim_id="claim-uncached", source_id="PMID:9999998"),
            ])
            result = verify_route_b_semantic_sources(
                str(proposals), str(root / "verify"), cache_arg=str(_cache(root))
            )
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            states = dict(zip(frame["claim_id"], frame["cache_investigation_state"]))
            self.assertEqual("cache_hit", states["claim-cached"])
            self.assertEqual("provider_failure", states["claim-uncached"])
            # provider_failure must never leak into a negative relation status.
            uncached_row = frame[frame["claim_id"] == "claim-uncached"].iloc[0]
            self.assertEqual("not_evaluated", uncached_row["source_relation_status"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(
                {"cache_hit": 1, "new_investigation": 0, "provider_failure": 1, "not_applicable": 0},
                manifest["cache_investigation_state_counts"],
            )

    def test_cache_hit_and_new_investigation_coexist_in_one_canonical_run(self) -> None:
        """Two claims in the same call, live provider allowed: one is served
        from cache, one triggers a fresh live lookup this run. The live call is
        monkeypatched at the module boundary used by ``source_verify`` so the
        test never touches the network, matching the pattern used for other
        allow-external-api tests in this repository (see test_reference_sources.py)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [
                _proposal(claim_id="claim-cached", source_id="PMID:1234567"),
                _proposal(claim_id="claim-live", source_id="PMID:5555555"),
            ])

            def _fake_fetch(query: str, timeout: int = 30):
                return {
                    "id": "5555555", "pmid": "5555555",
                    "title": "REF1 controls shoot meristem development",
                    "abstractText": "Loss-of-function knockout mutants of REF1 display a fasciated shoot meristem.",
                    "journalTitle": "Plant Journal", "pubYear": "2026", "doi": "10.1000/ref1-live",
                    "authorString": "Example C",
                }

            original_fetch = source_verify_module._fetch_europe_pmc
            source_verify_module._fetch_europe_pmc = _fake_fetch
            try:
                result = verify_route_b_semantic_sources(
                    str(proposals), str(root / "verify"),
                    cache_arg=str(_cache(root)), allow_external_api=True,
                )
            finally:
                source_verify_module._fetch_europe_pmc = original_fetch

            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            states = dict(zip(frame["claim_id"], frame["cache_investigation_state"]))
            self.assertEqual("cache_hit", states["claim-cached"])
            self.assertEqual("new_investigation", states["claim-live"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(
                {"cache_hit": 1, "new_investigation": 1, "provider_failure": 0, "not_applicable": 0},
                manifest["cache_investigation_state_counts"],
            )
            self.assertEqual(1, manifest["cache_accounting"]["new_investigation_count"])

    def test_ncbi_fallback_is_recorded_as_new_investigation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal(source_id="PMID:5555555")])

            original_europe = source_verify_module._fetch_europe_pmc
            original_ncbi = source_verify_module._fetch_ncbi_eutils
            source_verify_module._fetch_europe_pmc = lambda query, timeout=30: None
            source_verify_module._fetch_ncbi_eutils = lambda citation_type, identifier, timeout=30: {
                "id": "5555555", "pmid": "5555555",
                "title": "REF1 controls shoot meristem development",
                "abstractText": "Loss-of-function knockout mutants of REF1 display a fasciated shoot meristem.",
                "journalTitle": "Plant Journal", "pubYear": "2026", "doi": "10.1000/ref1-ncbi",
                "authorString": "Example N",
                "_source_provider": "ncbi_eutils_pubmed",
                "_source_url": "https://eutils.ncbi.nlm.nih.gov/fixture",
            }
            try:
                result = verify_route_b_semantic_sources(
                    str(proposals), str(root / "verify"), allow_external_api=True,
                )
            finally:
                source_verify_module._fetch_europe_pmc = original_europe
                source_verify_module._fetch_ncbi_eutils = original_ncbi

            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("ncbi_eutils_pubmed", frame.iloc[0]["source_provider"])
            self.assertEqual("new_investigation", frame.iloc[0]["cache_investigation_state"])

    def test_live_provider_failure_is_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proposals = root / "proposals.tsv"
            _write_proposals(proposals, [_proposal()])
            original_europe = source_verify_module._fetch_europe_pmc
            original_ncbi = source_verify_module._fetch_ncbi_eutils
            source_verify_module._fetch_europe_pmc = lambda query, timeout=30: (_ for _ in ()).throw(RuntimeError("down"))
            source_verify_module._fetch_ncbi_eutils = lambda citation_type, identifier, timeout=30: (_ for _ in ()).throw(RuntimeError("down"))
            try:
                result = verify_route_b_semantic_sources(
                    str(proposals), str(root / "verify"), allow_external_api=True,
                )
            finally:
                source_verify_module._fetch_europe_pmc = original_europe
                source_verify_module._fetch_ncbi_eutils = original_ncbi
            frame = pd.read_csv(result["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", frame.iloc[0]["source_relation_status"])
            self.assertEqual("provider_failure", frame.iloc[0]["cache_investigation_state"])


if __name__ == "__main__":
    unittest.main()
