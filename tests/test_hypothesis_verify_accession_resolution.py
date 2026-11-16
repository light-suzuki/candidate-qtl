"""A proposal that obeyed "do not invent an accession" must still be verifiable (issue #492).

The investigation prompt tells the responder:

    Do not invent accessions or source ids. Leave the field empty instead.
    An empty field is recoverable; a fabricated one is not.

Nothing recovered it. A proposal naming a real gene with an empty
`proposed_stable_id` was rejected `no_stable_accession_proposed`, while the
accession sat one query away on the UniProt API this stage may already call. All
37 agent proposals across ten loci were rejected that way.

Resolution is not an exemption: an accession found from a gene name goes through
the same existence / species / sequence / source-backed checks as one the
proposal supplied. These tests pin both halves -- that a name resolves, and that
resolving does not skip a gate. No network: the API is stubbed.
"""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
import urllib.error
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from candidate_qtl.llm_pathway import (
    ACCESSION_RESOLUTION_STATES, _explicit_agi_loci, _searchable_symbol, bait_panel_freeze, hypothesis_verify,
)
from tests.bound_claim_fixtures import bound_claim

HYPOTHESIS_COLUMNS = ["hypothesis_id", "case_id", "gene_family", "reference_gene_name",
                      "reference_species", "proposed_stable_id"]


def _entry(accession: str, species: str, tair_locus: str = "", *, include_references: bool = True) -> dict:
    entry = {
        "primaryAccession": accession,
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "entryAudit": {"entryVersion": 7},
        "organism": {"scientificName": species},
        "sequence": {"value": "MAAGSVDPLA"},
        "references": ([{"citation": {"citationCrossReferences": [
            {"database": "PubMed", "id": "12345"}]}}] if include_references else []),
    }
    if tair_locus:
        entry["uniProtKBCrossReferences"] = [{"database": "TAIR", "id": tair_locus}]
    return entry


class _Response(BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _api(search_hits: list[str], entry_species: str = "Arabidopsis thaliana", seen: list[str] | None = None):
    """Stub both UniProt calls: the symbol search and the entry fetch."""
    def fake(request, timeout=60):  # noqa: ARG001
        url = request.full_url
        if seen is not None:
            seen.append(url)
        if "/search" in url:
            payload = {"results": [{"primaryAccession": a} for a in search_hits]}
        else:
            accession = url.rsplit("/", 1)[-1].removesuffix(".json")
            payload = _entry(accession, entry_species)
        return _Response(json.dumps(payload).encode("utf-8"))
    return fake


class AccessionResolutionTests(unittest.TestCase):
    @staticmethod
    def _source_records(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
        records: dict[str, dict[str, object]] = {}
        for row in rows:
            source_ids = sorted({
                source_id.strip()
                for source_id in str(row.get("verified_source_ids", "")).replace(",", ";").split(";")
                if source_id.strip()
            })
            expected_rows: list[dict[str, str]] = []
            for source_id in source_ids:
                provenance = {
                    "citation": source_id,
                    "resolved_id": source_id.removeprefix("PMID:"),
                    "provider": "local_record_cache",
                    "url": "frozen-literature-cache.json",
                    "retrieved_at": "2026-08-22T00:00:00+00:00",
                    "claim_support_status": (
                        "supported_direct_genetic_test"
                        if row.get("claim_support_summary") == "mechanistic_support"
                        else "associative_support_only"
                    ),
                }
                record = records.setdefault(source_id, {
                    "citation": source_id,
                    "candidate_qtl_hypothesis_verifications": {},
                })
                verifications = record["candidate_qtl_hypothesis_verifications"]
                assert isinstance(verifications, dict)
                verifications[str(row["hypothesis_id"])] = provenance
                expected_rows.append(provenance)
            current = str(row.get("publication_source_provenance", ""))
            try:
                current_rows = json.loads(current) if current else []
            except json.JSONDecodeError:
                current_rows = []
            if source_ids and current_rows and all(
                isinstance(item, dict) and set(item) == {"citation"} for item in current_rows
            ):
                row["publication_source_provenance"] = json.dumps(
                    expected_rows, sort_keys=True, separators=(",", ":"))
        return records

    @staticmethod
    def _attach_source_hashes(rows: list[dict[str, str]],
                               source_records: dict[str, dict[str, object]]) -> None:
        for row in rows:
            source_ids = sorted({
                value.strip() for value in str(row.get("verified_source_ids", "")).split(";")
                if value.strip()
            })
            if row.get("publication_records_sha256") == "AUTO":
                row["publication_records_sha256"] = hashlib.sha256(json.dumps(
                    [source_records[source_id] for source_id in source_ids],
                    ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                ).encode("utf-8")).hexdigest()

    def _verify(self, rows: list[dict[str, str]], api, *, allow_external_api: bool = True,
                cache_records: list[dict] | None = None) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_records = self._source_records(rows)
            self._attach_source_hashes(rows, source_records)
            hypotheses = root / "hypotheses.tsv"
            extra_columns = sorted(set().union(*(row.keys() for row in rows)) - set(HYPOTHESIS_COLUMNS))
            pd.DataFrame(rows, columns=HYPOTHESIS_COLUMNS + extra_columns).to_csv(
                hypotheses, sep="\t", index=False)
            cache = root / "cache"
            cache.mkdir()
            for index, record in enumerate(cache_records or []):
                (cache / f"record-{index}.json").write_text(json.dumps(record), encoding="utf-8")
            source_records_path = root / "source_records.json"
            source_records_path.write_text(json.dumps(source_records), encoding="utf-8")
            with patch("candidate_qtl.llm_pathway.urllib.request.urlopen", api):
                hypothesis_verify(str(hypotheses), str(root / "out"), "Pisum sativum",
                                  record_cache_arg=str(cache) if cache_records is not None else None,
                                  allow_external_api=allow_external_api,
                                  source_records_arg=(str(source_records_path)
                                                      if source_records else None))
            return pd.read_csv(root / "out" / "proposal_validation.tsv", sep="\t", dtype=str).fillna("")

    def _row(self, **overrides: str) -> dict[str, str]:
        row = {"hypothesis_id": "mech_001", "case_id": "c", "gene_family": "PEBP",
               "reference_gene_name": "TFL1 (TERMINAL FLOWER 1)",
               "reference_species": "Arabidopsis thaliana", "proposed_stable_id": ""}
        row.update(overrides)
        return row

    def test_a_gene_name_resolves_and_then_verifies(self) -> None:
        result = self._verify([self._row()], _api(["P93003"])).iloc[0]
        self.assertEqual("resolved_from_gene_name", result["accession_resolution_status"])
        self.assertEqual("P93003", result["proposed_stable_id"])
        self.assertEqual("verified", result["verification_status"])
        # The resolved accession was actually fetched and checked, not assumed.
        self.assertEqual("True", result["sequence_retrieved"])
        self.assertEqual("True", result["species_match"])
        self.assertEqual("Arabidopsis thaliana", result["retrieved_species"])

    def test_resolution_does_not_skip_the_species_gate(self) -> None:
        # The search returns an accession, but the entry belongs to another species.
        result = self._verify([self._row()], _api(["P93003"], entry_species="Zea mays")).iloc[0]
        self.assertEqual("resolved_from_gene_name", result["accession_resolution_status"])
        self.assertEqual("rejected", result["verification_status"])
        self.assertEqual("declared_species_does_not_match_accession", result["rejection_reason"])

    def test_later_explicit_species_alternative_satisfies_the_species_gate(self) -> None:
        declared = "Citrullus lanatus / Arabidopsis thaliana"
        result = self._verify([
            self._row(reference_species=declared, proposed_stable_id="P93003")
        ], _api([])).iloc[0]
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("True", result["species_match"])
        self.assertEqual(declared, result["declared_species"])

    def test_explicit_species_alternatives_still_reject_a_true_mismatch(self) -> None:
        result = self._verify([
            self._row(
                reference_species="Citrullus lanatus; Arabidopsis thaliana | Lotus japonicus",
                proposed_stable_id="P93003",
            )
        ], _api([], entry_species="Zea mays")).iloc[0]
        self.assertEqual("rejected", result["verification_status"])
        self.assertEqual("False", result["species_match"])
        self.assertEqual("declared_species_does_not_match_accession", result["rejection_reason"])

    def test_an_ambiguous_symbol_is_not_tie_broken(self) -> None:
        result = self._verify([self._row()], _api(["P93003", "Q9ZNV5"])).iloc[0]
        self.assertEqual("ambiguous_gene_name", result["accession_resolution_status"])
        self.assertEqual("2", result["accession_resolution_hit_count"])
        self.assertEqual("", result["proposed_stable_id"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("gene_name_resolves_to_multiple_accessions", result["review_reason"])

    def test_an_unresolvable_name_keeps_the_existing_rejection(self) -> None:
        result = self._verify([self._row(reference_gene_name="IPT family")], _api([])).iloc[0]
        self.assertEqual("unresolved_gene_name", result["accession_resolution_status"])
        self.assertEqual("rejected", result["verification_status"])
        self.assertEqual("no_stable_accession_proposed", result["rejection_reason"])

    def test_non_object_uniprot_search_is_provider_unavailable(self) -> None:
        def json_array(request, timeout=60):  # noqa: ARG001
            return _Response(b"[]")

        result = self._verify([self._row()], json_array).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])

    def test_missing_uniprot_search_results_are_provider_unavailable(self) -> None:
        def empty_object(request, timeout=60):  # noqa: ARG001
            return _Response(b"{}")

        result = self._verify([self._row()], empty_object).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])

    def test_null_uniprot_search_results_are_provider_unavailable(self) -> None:
        def null_results(request, timeout=60):  # noqa: ARG001
            return _Response(b'{"results":null}')

        result = self._verify(
            [self._row(proposed_stable_id="TAIR:AT5G03840")], null_results,
        ).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])

    def test_offline_does_not_attempt_resolution(self) -> None:
        seen: list[str] = []
        result = self._verify([self._row()], _api(["P93003"], seen=seen), allow_external_api=False).iloc[0]
        self.assertEqual("not_attempted_offline", result["accession_resolution_status"])
        self.assertEqual([], seen, "offline runs must not call the API")
        self.assertEqual("rejected", result["verification_status"])

    def test_a_supplied_accession_is_not_re_resolved(self) -> None:
        seen: list[str] = []
        result = self._verify([self._row(proposed_stable_id="P93003")], _api([], seen=seen)).iloc[0]
        self.assertEqual("supplied_by_proposal", result["accession_resolution_status"])
        self.assertEqual("verified", result["verification_status"])
        self.assertFalse([url for url in seen if "/search" in url], "no symbol search for a supplied accession")

    def test_one_declared_source_accession_is_resolved_after_literature_verification(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="",
            source_ids="PMID:12345;UniProtKB:Q6UD73",
            verified_source_ids="PMID:12345",
            verified_citation_count="1",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], _api([]), cache_records=[
            _entry("Q6UD73", "Arabidopsis thaliana", include_references=False)
        ]).iloc[0]
        self.assertEqual("resolved_from_declared_source_accession", result["accession_resolution_status"])
        self.assertEqual("Q6UD73", result["proposed_stable_id"])
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("True", result["proposal_literature_verified"])
        self.assertEqual("PMID:12345", result["proposal_verified_source_ids"])
        self.assertEqual("PMID:12345;UniProt:Q6UD73", result["function_source_ids"])
        self.assertEqual("7", result["source_record_version"])

    def test_verified_paper_and_resolved_accession_freeze_both_provenance_axes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = root / "hypotheses.tsv"
            rows = [self._row(
                reference_gene_name="",
                source_ids="PMID:12345;UniProtKB:Q6UD73",
                verified_source_ids="PMID:12345",
                verified_citation_count="1",
                claim_support_summary="mechanistic_support",
                mechanistic_support_citation_count="1",
                publication_source_provenance='[{"citation":"PMID:12345"}]',
                publication_records_sha256="AUTO",
            )]
            source_records = self._source_records(rows)
            self._attach_source_hashes(rows, source_records)
            pd.DataFrame(rows).to_csv(hypotheses, sep="\t", index=False)
            source_records_path = root / "source_records.json"
            source_records_path.write_text(json.dumps(source_records), encoding="utf-8")
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(json.dumps(
                _entry("Q6UD73", "Arabidopsis thaliana", include_references=False)
            ), encoding="utf-8")
            verified_dir = root / "verified"
            bound_path = root / "bound.tsv"
            pd.DataFrame([bound_claim(
                route="a", claim_id="mech_001", citation="PMID:12345",
                subject_gene="Q6UD73", organism="Arabidopsis thaliana",
                phenotype="flowering architecture",
            )], columns=BOUND_CLAIM_COLUMNS).to_csv(
                bound_path, sep="\t", index=False, lineterminator="\n"
            )
            hypothesis_verify(str(hypotheses), str(verified_dir), "Pisum sativum",
                              record_cache_arg=str(cache), allow_external_api=False,
                              source_records_arg=str(source_records_path),
                              bound_claims_arg=str(bound_path))
            panel_dir = root / "panel"
            bait_panel_freeze(
                str(hypotheses), str(verified_dir / "proposal_validation.tsv"),
                str(verified_dir / "verified_sequences.json"), str(panel_dir), case_id="c",
            )
            panel = pd.read_csv(panel_dir / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["Q6UD73"], panel["bait_id"].tolist())
            self.assertEqual("PMID:12345", panel.iloc[0]["proposal_verified_source_ids"])
            self.assertEqual("resolved_from_declared_source_accession",
                             panel.iloc[0]["accession_resolution_status"])
            self.assertEqual("7", panel.iloc[0]["source_record_version"])
            manifest = json.loads((panel_dir / "bait_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("PMID:12345", manifest["provenance"]["Q6UD73"]["proposal_verified_source_ids"])

    def test_publication_hash_must_match_the_frozen_source_records(self) -> None:
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="b" * 64,
        )], _api([])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_publication_provenance_must_match_its_hypothesis_binding(self) -> None:
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance=json.dumps([{
                "citation": "PMID:12345", "resolved_id": "12345",
                "provider": "forged-provider", "url": "https://example.invalid/forged",
                "retrieved_at": "2099-01-01T00:00:00Z",
                "claim_support_status": "supported_direct_genetic_test",
            }]),
            publication_records_sha256="AUTO",
        )], _api([])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_publication_provenance_rejects_non_object_extra_elements(self) -> None:
        provenance = {
            "citation": "PMID:12345", "resolved_id": "12345",
            "provider": "local_record_cache", "url": "frozen-literature-cache.json",
            "retrieved_at": "2026-08-22T00:00:00+00:00",
            "claim_support_status": "supported_direct_genetic_test",
        }
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance=json.dumps([provenance, "garbage"]),
            publication_records_sha256="AUTO",
        )], _api([])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_multiple_declared_source_accessions_are_not_tie_broken(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="",
            source_ids="UniProtKB:Q6UD73;UniProt:P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], _api([]), cache_records=[
            _entry("Q6UD73", "Arabidopsis thaliana"),
            _entry("P93003", "Arabidopsis thaliana"),
        ]).iloc[0]
        self.assertEqual("ambiguous_declared_source_accessions", result["accession_resolution_status"])
        self.assertEqual("", result["proposed_stable_id"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("declared_source_ids_contain_multiple_accessions", result["review_reason"])

    def test_accession_alone_does_not_replace_verified_proposal_literature(self) -> None:
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            source_ids="UniProt:P93003",
            verified_source_ids="",
            verified_citation_count="0",
        )], _api([]), cache_records=[_entry("P93003", "Arabidopsis thaliana")]).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_associative_paper_plus_curated_function_passes_the_dual_path(self) -> None:
        """Associative support alone never passed; with an authoritative curated
        function it satisfies the two-point combination (issue #611 sibling)."""
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="associative_support_only",
            mechanistic_support_citation_count="0",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], _api([]), cache_records=[_entry("P93003", "Arabidopsis thaliana")]).iloc[0]
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("", result["review_reason"])

    def test_associative_paper_without_curated_function_still_fails(self) -> None:
        entry = _entry("P93003", "Arabidopsis thaliana", include_references=False)
        entry["entryType"] = "UniProtKB unreviewed (TrEMBL)"
        result = self._verify([self._row(
            proposed_stable_id="P93003",
            verified_source_ids="PMID:12345",
            claim_support_summary="associative_support_only",
            mechanistic_support_citation_count="0",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], _api([]), cache_records=[entry]).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_uniprot_provider_failure_is_review_required_not_an_exception(self) -> None:
        def unavailable(request, timeout=60):  # noqa: ARG001
            raise OSError("provider offline")

        result = self._verify([self._row(proposed_stable_id="P93003")], unavailable).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])
        self.assertIn("provider offline", result["accession_provider_error"])

    def test_invalid_uniprot_json_is_review_required_not_an_exception(self) -> None:
        def invalid_json(request, timeout=60):  # noqa: ARG001
            return _Response(b"{not-json")

        result = self._verify([self._row(proposed_stable_id="P93003")], invalid_json).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])
        self.assertTrue(result["accession_provider_error"])

    def test_non_object_uniprot_json_is_review_required_not_an_exception(self) -> None:
        def json_array(request, timeout=60):  # noqa: ARG001
            return _Response(b'[{"unexpected":"shape"}]')

        result = self._verify([self._row(proposed_stable_id="P93003")], json_array).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])
        self.assertIn("invalid record payload", result["accession_provider_error"])

    def test_missing_uniprot_entry_fields_are_review_required(self) -> None:
        def empty_object(request, timeout=60):  # noqa: ARG001
            return _Response(b"{}")

        result = self._verify([self._row(proposed_stable_id="P93003")], empty_object).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])

    def test_duplicate_accessions_keep_hypothesis_level_manifest_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [self._row(
                hypothesis_id=hypothesis_id,
                proposed_stable_id="P93003",
                verified_source_ids="PMID:12345",
                claim_support_summary="mechanistic_support",
                mechanistic_support_citation_count="1",
                publication_source_provenance='[{"citation":"PMID:12345"}]',
                publication_records_sha256="AUTO",
            ) for hypothesis_id in ("mech_001", "mech_002")]
            source_records = self._source_records(rows)
            self._attach_source_hashes(rows, source_records)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame(rows).to_csv(hypotheses, sep="\t", index=False)
            source_records_path = root / "source_records.json"
            source_records_path.write_text(json.dumps(source_records), encoding="utf-8")
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(json.dumps(
                _entry("P93003", "Arabidopsis thaliana")
            ), encoding="utf-8")
            verified_dir = root / "verified"
            bound_path = root / "bound.tsv"
            pd.DataFrame([
                bound_claim(
                    route="a", claim_id=str(row["hypothesis_id"]), citation="PMID:12345",
                    subject_gene="P93003", organism="Arabidopsis thaliana",
                    phenotype="flowering architecture",
                )
                for row in rows
            ], columns=BOUND_CLAIM_COLUMNS).to_csv(
                bound_path, sep="\t", index=False, lineterminator="\n"
            )
            hypothesis_verify(
                str(hypotheses), str(verified_dir), "Pisum sativum",
                record_cache_arg=str(cache), allow_external_api=False,
                source_records_arg=str(source_records_path),
                bound_claims_arg=str(bound_path),
            )
            panel_dir = root / "panel"
            bait_panel_freeze(
                str(hypotheses), str(verified_dir / "proposal_validation.tsv"),
                str(verified_dir / "verified_sequences.json"), str(panel_dir), case_id="c",
            )
            manifest = json.loads((panel_dir / "bait_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["bait_count"])
            self.assertEqual(
                {"P93003::mech_001", "P93003::mech_002"},
                set(manifest["provenance"]),
            )

    def test_a_supplied_tair_locus_resolves_then_verifies(self) -> None:
        seen: list[str] = []
        result = self._verify(
            [self._row(proposed_stable_id="TAIR:AT5G03840")],
            _api(["P93003"], seen=seen),
        ).iloc[0]
        self.assertEqual("resolved_from_tair_locus", result["accession_resolution_status"])
        self.assertEqual("P93003", result["proposed_stable_id"])
        self.assertEqual("verified", result["verification_status"])
        search_urls = [url for url in seen if "/search" in url]
        self.assertEqual(1, len(search_urls))
        self.assertIn("xref%3Atair-AT5G03840", search_urls[0])

    def test_an_ambiguous_tair_locus_is_not_tie_broken(self) -> None:
        result = self._verify(
            [self._row(proposed_stable_id="AT5G03840")],
            _api(["P93003", "Q9ZNV5"]),
        ).iloc[0]
        self.assertEqual("ambiguous_tair_locus", result["accession_resolution_status"])
        self.assertEqual("", result["proposed_stable_id"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("tair_locus_resolves_to_multiple_accessions", result["review_reason"])

    def test_a_tair_locus_stays_offline_without_network_resolution(self) -> None:
        seen: list[str] = []
        result = self._verify(
            [self._row(proposed_stable_id="AT5G03840")],
            _api(["P93003"], seen=seen),
            allow_external_api=False,
        ).iloc[0]
        self.assertEqual("not_attempted_offline", result["accession_resolution_status"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual([], seen)

    def test_a_tair_locus_resolves_from_a_frozen_record_cache_offline(self) -> None:
        seen: list[str] = []
        result = self._verify(
            [self._row(proposed_stable_id="TAIR:AT5G03840")],
            _api([], seen=seen),
            allow_external_api=False,
            cache_records=[_entry("P93003", "Arabidopsis thaliana", "AT5G03840")],
        ).iloc[0]
        self.assertEqual("resolved_from_tair_locus_cache", result["accession_resolution_status"])
        self.assertIn(result["accession_resolution_status"], ACCESSION_RESOLUTION_STATES)
        self.assertEqual("record_cache:xref:tair-AT5G03840", result["accession_resolution_query"])
        self.assertEqual("1", result["accession_resolution_hit_count"])
        self.assertEqual("P93003", result["proposed_stable_id"])
        self.assertEqual("local_record_cache", result["source_provider"])
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual([], seen, "offline cache resolution must not call the API")

    def test_an_ambiguous_tair_locus_in_cache_is_not_tie_broken(self) -> None:
        seen: list[str] = []
        records = [
            _entry("P93003", "Arabidopsis thaliana", "AT5G03840"),
            _entry("Q9ZNV5", "Arabidopsis thaliana", "AT5G03840"),
        ]
        result = self._verify(
            [self._row(proposed_stable_id="AT5G03840")],
            _api([], seen=seen),
            allow_external_api=False,
            cache_records=records,
        ).iloc[0]
        self.assertEqual("ambiguous_tair_locus", result["accession_resolution_status"])
        self.assertEqual("2", result["accession_resolution_hit_count"])
        self.assertEqual("", result["proposed_stable_id"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("tair_locus_resolves_to_multiple_accessions", result["review_reason"])
        self.assertEqual([], seen)

    def test_an_unmatched_tair_locus_remains_offline_review_required(self) -> None:
        seen: list[str] = []
        result = self._verify(
            [self._row(proposed_stable_id="AT5G03840")],
            _api([], seen=seen),
            allow_external_api=False,
            cache_records=[_entry("P93003", "Arabidopsis thaliana", "AT1G01010")],
        ).iloc[0]
        self.assertEqual("not_attempted_offline", result["accession_resolution_status"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual([], seen)

    def test_the_query_used_is_recorded(self) -> None:
        result = self._verify([self._row()], _api(["P93003"])).iloc[0]
        self.assertIn("gene:TFL1", result["accession_resolution_query"])
        self.assertIn("Arabidopsis thaliana", result["accession_resolution_query"])
        self.assertIn("reviewed:true", result["accession_resolution_query"])

    def test_symbol_extraction_from_free_text(self) -> None:
        self.assertEqual("TFL1", _searchable_symbol("TFL1 (TERMINAL FLOWER 1)"))
        self.assertEqual("IPT", _searchable_symbol("IPT (isopentenyltransferase) family"))
        self.assertEqual("", _searchable_symbol(""))
        # An accession is not a symbol: searching gene:P93003 would match nothing useful.
        self.assertEqual("", _searchable_symbol("P93003"))

    def test_explicit_agi_in_description_reaches_tair_resolver(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="LONG DESCRIPTIVE NAME (XYZ1; at1g12345)")],
            _api(["P93003"])).iloc[0]
        self.assertEqual("resolved_from_tair_locus", result["accession_resolution_status"])
        self.assertIn("xref:tair-AT1G12345", result["accession_resolution_query"])
        self.assertEqual("P93003", result["proposed_stable_id"])

    def test_explicit_agi_is_found_when_symbol_comes_first(self) -> None:
        self.assertEqual(["AT1G12345"], _explicit_agi_loci("XYZ1 (description; at1g12345)"))

    def test_multiple_or_conflicting_explicit_agi_fail_closed(self) -> None:
        multiple = self._verify([self._row(
            reference_gene_name="description (AT1G12345; AT2G23456)")], _api(["P93003"])).iloc[0]
        self.assertEqual("review_required", multiple["verification_status"])
        self.assertEqual("multiple_explicit_agi_in_reference_gene_name", multiple["review_reason"])
        conflict = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)", proposed_stable_id="AT2G23456")],
            _api(["P93003"])).iloc[0]
        self.assertEqual("review_required", conflict["verification_status"])
        self.assertEqual("explicit_agi_conflicts_with_proposed_stable_id", conflict["review_reason"])
        unresolved = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)", proposed_stable_id="P93003")],
            _api(["Q6UD73"])).iloc[0]
        self.assertEqual("review_required", unresolved["verification_status"])
        self.assertEqual("explicit_agi_requires_identifier_correspondence_review", unresolved["review_reason"])

    def test_live_agi_correspondence_matching_proposed_accession_can_pass(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)", proposed_stable_id="P93003")],
            _api(["P93003"])).iloc[0]
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("resolved_from_agi_correspondence", result["accession_resolution_status"])

    def test_live_agi_correspondence_provider_failure_is_review(self) -> None:
        def unavailable(request, timeout=60):  # noqa: ARG001
            raise urllib.error.URLError("offline")
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)", proposed_stable_id="P93003")], unavailable).iloc[0]
        self.assertEqual("provider_unavailable", result["accession_resolution_status"])
        self.assertEqual("review_required", result["verification_status"])

    def test_live_agi_correspondence_cardinality_is_not_tie_broken(self) -> None:
        seen: list[str] = []
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)", proposed_stable_id="P93003")],
            _api(["P93003", "Q6UD73"], seen=seen)).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("unresolved_reference_gene_correspondence", result["accession_resolution_status"])
        self.assertEqual("2", result["accession_resolution_hit_count"])
        self.assertEqual(1, len(seen))

    def test_matching_mixed_case_proposed_agi_is_normalized(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (at1g12345)", proposed_stable_id="At1g12345")],
            _api(["P93003"])).iloc[0]
        self.assertEqual("resolved_from_tair_locus", result["accession_resolution_status"])
        self.assertEqual("P93003", result["proposed_stable_id"])

    def test_single_source_id_requires_agi_correspondence(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)",
            source_ids="UniProt:Q6UD73")], _api(["P93003"])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("explicit_agi_requires_identifier_correspondence_review", result["review_reason"])

    def test_cached_source_id_with_matching_agi_xref_can_pass(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)",
            source_ids="UniProt:P93003")], _api([]),
            cache_records=[_entry("P93003", "Arabidopsis thaliana", "AT1G12345")]).iloc[0]
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("resolved_from_agi_correspondence", result["accession_resolution_status"])

    def test_cached_source_id_with_different_agi_xref_stays_review(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)",
            source_ids="UniProt:P93003")], _api([]),
            cache_records=[_entry("P93003", "Arabidopsis thaliana", "AT2G23456")]).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("explicit_agi_requires_identifier_correspondence_review", result["review_reason"])

    def test_cached_agi_correspondence_cardinality_is_not_tie_broken(self) -> None:
        seen: list[str] = []
        records = [
            _entry("P93003", "Arabidopsis thaliana", "AT1G12345"),
            _entry("Q6UD73", "Arabidopsis thaliana", "AT1G12345"),
        ]
        result = self._verify([self._row(
            reference_gene_name="XYZ1 (AT1G12345)",
            source_ids="UniProt:P93003")], _api([], seen=seen), cache_records=records).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("2", result["accession_resolution_hit_count"])
        self.assertEqual([], seen)

    def test_explicit_agi_boundaries_reject_embedded_or_malformed_tokens(self) -> None:
        self.assertEqual([], _explicit_agi_loci(
            "XAT1G12345Y AT1G123456 AT1G1234 foo_AT1G12345 AT1G12345.1x AT1G12345.1.2"))
        self.assertEqual(["AT1G12345", "AT1G12345.1"],
                         _explicit_agi_loci("AT1G12345, AT1G12345.1"))

    def test_explicit_agi_does_not_bypass_multiple_source_id_gate(self) -> None:
        result = self._verify([self._row(
            reference_gene_name="description (AT1G12345)",
            source_ids="UniProt:P93003;UniProt:Q6UD73")], _api(["P93003"])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("declared_source_ids_contain_multiple_accessions", result["review_reason"])

    def test_multiple_explicit_agi_stays_ambiguous_with_proposed_agi(self) -> None:
        result = self._verify([self._row(
            proposed_stable_id="AT1G01010",
            reference_gene_name="description (AT1G12345; AT2G23456)")], _api([])).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("multiple_explicit_agi_in_reference_gene_name", result["review_reason"])

    def test_explicit_agi_provider_failure_is_review_not_negative(self) -> None:
        def unavailable(request, timeout=60):  # noqa: ARG001
            raise urllib.error.URLError("offline")

        result = self._verify([self._row(
            reference_gene_name="description (at1g12345)")], unavailable).iloc[0]
        self.assertEqual("provider_unavailable", result["accession_resolution_status"])
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("accession_provider_unavailable", result["review_reason"])


if __name__ == "__main__":
    unittest.main()
