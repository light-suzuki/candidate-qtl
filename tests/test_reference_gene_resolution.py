"""Reference-gene citation resolution in front of source-verify (issue #649).

All tests run fully offline: any test that needs --allow-external-api mocks
candidate_qtl.reference_gene_resolution's imported UniProt helpers directly
(_fetch_uniprot / _resolve_accession) instead of touching the network. Offline
resolution (proposed_stable_id already present + a --records cache hit, or a
TAIR cross-reference cache lookup) needs no mocking at all -- it exercises the
same offline record-cache code path hypothesis-verify already uses.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS
from candidate_qtl.reference_gene_resolution import reference_gene_resolve


def _uniprot_record(accession: str, species: str, sequence: str, *, reviewed: bool = True,
                    pmid: str = "12345", tair: str | None = None,
                    gene_name: str | None = None, synonyms: list[str] | None = None) -> dict:
    record: dict = {
        "primaryAccession": accession,
        "uniProtkbId": f"{accession}_TEST",
        "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
        "organism": {"scientificName": species},
        "sequence": {"value": sequence, "length": len(sequence)},
    }
    if pmid:
        record["references"] = [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": pmid}]}}]
    if tair:
        record["uniProtKBCrossReferences"] = [{"database": "TAIR", "id": tair}]
    if gene_name or synonyms:
        gene: dict = {}
        if gene_name:
            gene["geneName"] = {"value": gene_name}
        if synonyms:
            gene["synonyms"] = [{"value": synonym} for synonym in synonyms]
        record["genes"] = [gene]
    return record


def _mechanism_row(hypothesis_id: str, **overrides: object) -> dict[str, object]:
    row = {column: "" for column in MECHANISM_COLUMNS}
    row.update({
        "hypothesis_id": hypothesis_id, "case_id": "fasciation", "pathway_or_process": "CLAVATA-WUSCHEL",
        "gene_family": "CLE receptor kinase", "mechanism_statement": "a lesion enlarges the meristem",
        "expected_role": "receptor kinase restricting stem cell pool", "direction_consistency": "consistent",
        "model_species_evidence": "Arabidopsis clv1 mutant enlarges the meristem",
        "model_species_mutant_phenotype": "fasciated shoot apex",
        "reference_species": "Arabidopsis thaliana", "selection_rationale": "reviewed reference in a non-focal species",
        "counter_hypothesis": "meristem size can also change through hormone transport",
        "discriminating_observation": "measure meristem width over time",
        "confidence": "medium", "uncertainty": "family member identity in the focal species",
        "payload_source": "llm", "model": "test-model", "prompt_version": "phenotype-investigation-1.2.0",
        "parser_version": "mechanism-hypothesis-parser-1.2.0",
        "proposal_round": "1", "interval_informed": "False", "review_status": "pending_verification",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


def _write_mechanisms(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows)[MECHANISM_COLUMNS].to_csv(path, sep="\t", index=False, lineterminator="\n")


def _write_cache(path: Path, records: list[dict]) -> None:
    path.write_text(json.dumps({"results": records}), encoding="utf-8")


class ModelSuppliedTests(unittest.TestCase):
    def test_legacy_28_column_table_is_upgraded_with_blank_causal_layer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "legacy_mechanisms.tsv"
            row = _mechanism_row("legacy_001", source_ids="PMID:11111111")
            pd.DataFrame([row]).drop(columns=["causal_layer"]).to_csv(
                mechanisms, sep="\t", index=False, lineterminator="\n")
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("causal_layer", resolved.columns)
            self.assertEqual("", resolved.loc[0, "causal_layer"])

    def test_model_supplied_source_ids_are_preserved_untouched(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            original = _mechanism_row("mech_001", source_ids="PMID:11111111", proposed_stable_id="")
            _write_mechanisms(mechanisms, [original])
            out = root / "out"
            result = reference_gene_resolve(str(mechanisms), str(out))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            for column in MECHANISM_COLUMNS:
                self.assertEqual(str(original[column]), str(row[column]), column)
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("model_supplied", ledger_row["resolution_origin"])
            self.assertEqual("PMID:11111111", ledger_row["resolved_source_ids"])
            self.assertEqual(0, result["not_resolved_count"])
            self.assertEqual(1, result["model_supplied_count"])
            manifest = json.loads((out / "reference_gene_resolution_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["resolution_origin_counts"]["model_supplied"])
            self.assertEqual(0, manifest["review_required_count"])


class ProviderResolvedTests(unittest.TestCase):
    def test_supplied_accession_gets_citations_filled_from_offline_cache(self) -> None:
        """Rule 2: an already-supplied proposed_stable_id is used as-is; only its citations are fetched."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="P0CLV1", reference_gene_name="",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("P0CLV1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, pmid="99887766")])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("P0CLV1", row["proposed_stable_id"])
            self.assertIn("PMID:99887766", row["source_ids"])
            provenance = json.loads(row["retrieval_provenance"])
            self.assertEqual("local_record_cache", provenance["provider"])
            self.assertEqual("supplied_by_proposal", provenance["resolution_status"])
            self.assertTrue(provenance["retrieved_at"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("provider_resolved", ledger_row["resolution_origin"])
            self.assertEqual("local_record_cache", ledger_row["source_provider"])
            self.assertEqual("supplied_by_proposal", ledger_row["accession_resolution_status"])
            self.assertEqual("True", str(ledger_row["record_found"]))
            self.assertEqual("True", str(ledger_row["function_source_backed"]))
            ledger_provenance = json.loads(ledger_row["retrieval_provenance"])
            self.assertEqual(provenance, ledger_provenance)

    def test_gene_name_resolves_to_exactly_one_accession(self) -> None:
        """Rule 3: a reference gene name resolving to one reviewed record fills proposed_stable_id too."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="", reference_gene_name="CLV1",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("Q9XYZ1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, pmid="55443322", gene_name="CLV1")])
            out = root / "out"
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                            return_value=(["Q9XYZ1"], "resolved_from_gene_name",
                                         'gene:CLV1 AND organism_name:"Arabidopsis thaliana"', 1)) as patched:
                reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True, record_cache_arg=str(cache))
            patched.assert_called_once()
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("Q9XYZ1", row["proposed_stable_id"])
            self.assertIn("PMID:55443322", row["source_ids"])
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("provider_resolved", ledger_row["resolution_origin"])
            self.assertEqual("resolved_from_gene_name", ledger_row["accession_resolution_status"])
            self.assertEqual("1", str(ledger_row["accession_resolution_hit_count"]))


class AmbiguousResolutionTests(unittest.TestCase):
    def test_tair_locus_resolving_to_multiple_accessions_is_never_tie_broken(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="AT1G69120", reference_gene_name="",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [
                _uniprot_record("Q1111A", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, tair="AT1G69120"),
                _uniprot_record("Q2222B", "Arabidopsis thaliana", "MGTAAALPD" * 12, tair="AT1G69120"),
            ])
            out = root / "out"
            result = reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            # Left completely empty: never guessed between Q1111A and Q2222B.
            self.assertEqual("AT1G69120", row["proposed_stable_id"])
            self.assertEqual("", row["source_ids"])
            self.assertNotIn(row["proposed_stable_id"], {"Q1111A", "Q2222B"})

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("not_resolved", ledger_row["resolution_origin"])
            self.assertEqual("ambiguous_tair_locus", ledger_row["accession_resolution_status"])
            self.assertEqual("2", str(ledger_row["accession_resolution_hit_count"]))
            self.assertEqual("", ledger_row["resolved_stable_id"])

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[reviews["target"] == "mech_001"].iloc[0]
            self.assertEqual("reference_gene_resolution_required:declared_tair_locus_resolves_to_multiple_accessions",
                            review_row["reason"])
            self.assertIn("2", review_row["required_next_action"])
            self.assertEqual(1, result["review_required_count"])


class ProviderFailureTests(unittest.TestCase):
    def test_provider_failure_is_not_evaluated_never_a_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="P0FAI1", reference_species="Arabidopsis thaliana")])
            out = root / "out"
            with mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                            side_effect=RuntimeError("UniProt API HTTP 503 for P0FAI1")):
                result = reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True)
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("not_evaluated", ledger_row["evaluation_status"])
            self.assertEqual("not_resolved", ledger_row["resolution_origin"])
            self.assertIn("503", ledger_row["provider_error"])
            self.assertNotEqual("rejected", ledger_row["resolution_origin"])

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[reviews["target"] == "mech_001"].iloc[0]
            self.assertEqual("reference_gene_resolution_required:accession_provider_unavailable", review_row["reason"])

            manifest = json.loads((out / "reference_gene_resolution_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["mech_001"], manifest["not_evaluated_hypothesis_ids"])
            self.assertEqual(1, manifest["evaluation_status_counts"]["not_evaluated"])
            self.assertEqual(1, result["not_evaluated_count"])
            # A provider failure never shows up as a negative/rejected result anywhere.
            self.assertNotIn("rejected", json.dumps(manifest))


class OfflineWithoutApiTests(unittest.TestCase):
    def test_accession_missing_from_offline_cache_is_not_attempted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="P0OFF1", reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("P0OTH1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12)])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out), allow_external_api=False, record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("P0OFF1", row["proposed_stable_id"])
            self.assertEqual("", row["source_ids"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("not_attempted", ledger_row["evaluation_status"])
            self.assertEqual("not_attempted_offline", ledger_row["accession_resolution_status"])

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[reviews["target"] == "mech_001"].iloc[0]
            self.assertEqual("reference_gene_resolution_required:not_attempted_offline", review_row["reason"])

    def test_cli_reference_gene_resolve_matches_direct_call(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row("mech_001", source_ids="PMID:22223333")])
            out = root / "out"
            self.assertEqual(0, main(["reference-gene-resolve", "--mechanisms", str(mechanisms), "--out", str(out)]))
            self.assertTrue((out / "mechanism_hypotheses_resolved.tsv").exists())
            self.assertTrue((out / "reference_gene_resolution.tsv").exists())
            self.assertTrue((out / "review_required.tsv").exists())
            manifest = json.loads((out / "reference_gene_resolution_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("reference-gene-resolve", manifest["command"])


class NoSourceBackedCitationTests(unittest.TestCase):
    def test_record_without_source_backed_citation_stays_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="P0NOC1", reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("P0NOC1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12,
                                                 reviewed=False, pmid="")])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("", row["source_ids"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("not_resolved", ledger_row["resolution_origin"])
            self.assertEqual("evaluated", ledger_row["evaluation_status"])
            self.assertEqual("True", str(ledger_row["record_found"]))
            self.assertEqual("False", str(ledger_row["function_source_backed"]))

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[reviews["target"] == "mech_001"].iloc[0]
            self.assertEqual("reference_gene_resolution_required:reference_gene_function_is_not_source_backed",
                            review_row["reason"])


class MalformedProposedStableIdTests(unittest.TestCase):
    def test_malformed_proposed_stable_id_is_never_queried_and_falls_back_to_gene_name(self) -> None:
        """Amendment 1, criterion 1: a fabricated pseudo-accession is ignored, not queried, and
        the hypothesis is instead resolved from reference_gene_name + reference_species."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="Medicago_truncatula__NIN",
                reference_gene_name="NIN", reference_species="Medicago truncatula")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "Q9XYZ2", "Medicago truncatula", "MKLLTAVWQ" * 12, pmid="44556677",
                gene_name="NIN")])
            out = root / "out"
            with mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot") as fetch:
                result = reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            fetch.assert_not_called()  # the malformed string is never queried against any provider

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            # The clean table never carries the malformed string forward.
            self.assertEqual("Q9XYZ2", row["proposed_stable_id"])
            self.assertNotIn("Medicago_truncatula__NIN", row["proposed_stable_id"])
            self.assertIn("PMID:44556677", row["source_ids"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("Medicago_truncatula__NIN", ledger_row["input_proposed_stable_id"])
            self.assertEqual("malformed_proposed_stable_id_ignored", ledger_row["malformed_proposed_stable_id_status"])
            self.assertEqual("resolved_from_gene_name_cache", ledger_row["accession_resolution_status"])
            self.assertEqual("provider_resolved", ledger_row["resolution_origin"])
            self.assertEqual("Q9XYZ2", ledger_row["resolved_stable_id"])

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[(reviews["target"] == "mech_001")
                                 & (reviews["reason"] == "reference_gene_resolution_required:malformed_proposed_stable_id_ignored")].iloc[0]
            self.assertIn("Medicago_truncatula__NIN", review_row["required_next_action"])
            self.assertEqual(1, result["malformed_proposed_stable_id_count"])
            manifest = json.loads((out / "reference_gene_resolution_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["malformed_proposed_stable_id_counts"]["malformed_proposed_stable_id_ignored"])

    def test_malformed_value_naming_focal_species_gets_distinct_status(self) -> None:
        """Amendment 1, criterion 3: a fabricated pseudo-accession that names the focal species
        (e.g. ``Pisum_sativum__NIN``) is flagged with its own status, separate from ordinary
        malformation, and the focal species itself is threaded in as a caller-supplied parameter
        rather than a hardcoded string."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="Pisum_sativum__NIN",
                reference_gene_name="NIN", reference_species="Medicago truncatula")])
            out = root / "out"
            result = reference_gene_resolve(str(mechanisms), str(out), focal_species="Pisum sativum")
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("Pisum_sativum__NIN", ledger_row["input_proposed_stable_id"])
            self.assertEqual("malformed_proposed_stable_id_names_focal_species",
                            ledger_row["malformed_proposed_stable_id_status"])
            self.assertNotEqual("malformed_proposed_stable_id_ignored", ledger_row["malformed_proposed_stable_id_status"])

            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review_row = reviews[(reviews["target"] == "mech_001")
                                 & (reviews["reason"] == "reference_gene_resolution_required:malformed_proposed_stable_id_names_focal_species")].iloc[0]
            self.assertIn("focal species", review_row["required_next_action"])
            self.assertEqual(1, result["malformed_proposed_stable_id_count"])
            manifest = json.loads((out / "reference_gene_resolution_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["mech_001"], manifest["malformed_proposed_stable_id_names_focal_species_hypothesis_ids"])
            self.assertEqual("Pisum sativum", manifest["focal_species"])

    def test_malformed_value_not_naming_focal_species_stays_ordinary(self) -> None:
        """Without a focal_species match, the value is ordinary malformation, not the
        focal-species-specific status -- and with no focal_species argument at all, the
        distinction can never fire (conservative default: never guess at blinding)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="Pisum_sativum__NIN",
                reference_gene_name="", reference_species="")])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out))  # no focal_species supplied
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("malformed_proposed_stable_id_ignored", ledger_row["malformed_proposed_stable_id_status"])

    def test_well_formed_supplied_accession_is_unaffected_by_malformation_check(self) -> None:
        """Amendment 1, criterion 2: a well-formed model-supplied accession keeps today's
        behaviour exactly, even when a focal_species argument is supplied."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="P0CLV1", reference_gene_name="",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("P0CLV1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, pmid="99887766")])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache), focal_species="Pisum sativum")
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("", ledger_row["malformed_proposed_stable_id_status"])
            self.assertEqual("supplied_by_proposal", ledger_row["accession_resolution_status"])
            self.assertEqual("provider_resolved", ledger_row["resolution_origin"])
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("P0CLV1", row["proposed_stable_id"])


class NamedReferenceIdentityTests(unittest.TestCase):
    def test_exact_name_species_resolution_overrides_unrelated_supplied_accession(self) -> None:
        """A valid-looking model accession is still only an input claim when a
        complete reference name/species target is present."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "named_001", proposed_stable_id="P0WRONG", reference_gene_name="GENERIC1",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "Q9RIGHT", "Arabidopsis thaliana", "M" * 120, pmid="77889900", gene_name="GENERIC1")])
            out = root / "out"

            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved.iloc[0]
            self.assertEqual("Q9RIGHT", row["proposed_stable_id"])
            self.assertIn("PMID:77889900", row["source_ids"])
            self.assertNotIn("P0WRONG", row["source_ids"])
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            audit = ledger.iloc[0]
            self.assertEqual("P0WRONG", audit["input_proposed_stable_id"])
            self.assertEqual("Q9RIGHT", audit["resolved_stable_id"])
            self.assertEqual("resolved_from_gene_name_cache", audit["accession_resolution_status"])

    def test_named_identity_mismatch_does_not_verify_supplied_accession(self) -> None:
        """When the exact name/species lookup has no cache hit, the supplied
        accession remains audit input and is not allowed to carry publications.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "named_002", proposed_stable_id="P0WRONG", reference_gene_name="GENERIC2",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "P0WRONG", "Arabidopsis thaliana", "M" * 120, pmid="11223344", gene_name="OTHER1")])
            out = root / "out"

            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved.iloc[0]
            self.assertEqual("", row["proposed_stable_id"])
            self.assertEqual("", row["source_ids"])
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("P0WRONG", ledger.iloc[0]["input_proposed_stable_id"])
            self.assertEqual("not_resolved", ledger.iloc[0]["resolution_origin"])

    def test_missing_provider_gene_metadata_does_not_authenticate_supplied_accession(self) -> None:
        """A cache record without genes is insufficient to validate a named
        hypothesis's supplied accession; metadata absence is not identity proof."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "named_003", proposed_stable_id="P0WRONG", reference_gene_name="GENERIC3",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "P0WRONG", "Arabidopsis thaliana", "M" * 120, pmid="55667788")])
            out = root / "out"

            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("", resolved.iloc[0]["proposed_stable_id"])
            self.assertEqual("", resolved.iloc[0]["source_ids"])
            self.assertEqual("name_species_unresolved", resolved.iloc[0]["reference_identity_status"])

    def test_named_target_with_model_literature_source_still_resolves_identity(self) -> None:
        """A model PMID must not make a guessed accession authoritative; the
        PMID remains input audit while replacement sources come from the match."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "named_004", source_ids="PMID:10101010", proposed_stable_id="P0WRONG",
                reference_gene_name="GENERIC4", reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "Q9RIGHT", "Arabidopsis thaliana", "M" * 120, pmid="20202020", gene_name="GENERIC4")])
            out = root / "out"

            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved.iloc[0]
            self.assertEqual("Q9RIGHT", row["proposed_stable_id"])
            self.assertIn("PMID:20202020", row["source_ids"])
            self.assertNotIn("10101010", row["source_ids"])
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            audit = ledger.iloc[0]
            self.assertEqual("P0WRONG", audit["input_proposed_stable_id"])
            self.assertEqual("PMID:10101010", audit["input_source_ids"])


class GeneSymbolOfflineIndexTests(unittest.TestCase):
    def test_bare_gene_name_resolves_offline_from_cache_gene_symbol_index(self) -> None:
        """Amendment 2, criterion 4: a hypothesis with only reference_gene_name +
        reference_species resolves fully offline against a --records cache through the new
        gene-symbol index, and the ledger names the cache as provider."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="", reference_gene_name="YUC1",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "Q9XYZ3", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, pmid="33445566",
                gene_name="YUC1", synonyms=["YUCCA1"])])
            out = root / "out"
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession") as live_call:
                result = reference_gene_resolve(str(mechanisms), str(out), allow_external_api=False,
                                                 record_cache_arg=str(cache))
            live_call.assert_not_called()  # fully offline: no network path touched

            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("Q9XYZ3", row["proposed_stable_id"])
            self.assertIn("PMID:33445566", row["source_ids"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("provider_resolved", ledger_row["resolution_origin"])
            self.assertEqual("resolved_from_gene_name_cache", ledger_row["accession_resolution_status"])
            self.assertEqual("local_record_cache", ledger_row["source_provider"])
            self.assertEqual(1, result["provider_resolved_count"])

    def test_ambiguous_gene_symbol_cache_is_never_tie_broken(self) -> None:
        """Amendment 2, criterion 5: two cached records claiming the same symbol for the same
        species yield the ambiguity status, never a pick."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="", reference_gene_name="DUP1",
                reference_species="Arabidopsis thaliana")])
            cache = root / "records.json"
            _write_cache(cache, [
                _uniprot_record("Q9XYZ1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, gene_name="DUP1"),
                _uniprot_record("Q9XYZ2", "Arabidopsis thaliana", "MGTAAALPD" * 12, gene_name="DUP1"),
            ])
            out = root / "out"
            result = reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            children = resolved[resolved["parent_hypothesis_id"] == "mech_001"]
            self.assertEqual({"Q9XYZ1", "Q9XYZ2"}, set(children["proposed_stable_id"]))

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"mech_001__ref_q9xyz1", "mech_001__ref_q9xyz2"}, set(ledger["hypothesis_id"]))
            self.assertEqual(2, result["provider_resolved_count"])

    def test_gene_symbol_only_cached_for_a_different_species_does_not_match(self) -> None:
        """Amendment 2, criterion 6: a symbol present in the cache only for a different species
        never matches -- offline resolution stays not_attempted_offline, not a cross-species pick."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            _write_mechanisms(mechanisms, [_mechanism_row(
                "mech_001", source_ids="", proposed_stable_id="", reference_gene_name="NIN",
                reference_species="Medicago truncatula")])
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record(
                "Q3333E", "Lotus japonicus", "MKLLTAVWQ" * 12, gene_name="NIN")])
            out = root / "out"
            reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
            row = resolved[resolved["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("", row["proposed_stable_id"])
            self.assertEqual("", row["source_ids"])

            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            ledger_row = ledger[ledger["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("not_attempted_offline", ledger_row["accession_resolution_status"])
            self.assertEqual("not_attempted", ledger_row["evaluation_status"])
            self.assertEqual("not_resolved", ledger_row["resolution_origin"])


class MultiRowMixTests(unittest.TestCase):
    def test_mixed_batch_keeps_each_row_independent(self) -> None:
        """One mechanism table exercising every origin at once, matching a real a02 batch."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanism_hypotheses.tsv"
            rows = [
                _mechanism_row("mech_001", source_ids="PMID:11111111"),
                _mechanism_row("mech_002", source_ids="", proposed_stable_id="P0CLV1", reference_gene_name=""),
                _mechanism_row("mech_003", source_ids="", proposed_stable_id="", reference_gene_name=""),
            ]
            _write_mechanisms(mechanisms, rows)
            cache = root / "records.json"
            _write_cache(cache, [_uniprot_record("P0CLV1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12, pmid="99887766")])
            out = root / "out"
            result = reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
            self.assertEqual(3, result["hypothesis_count"])
            self.assertEqual(1, result["model_supplied_count"])
            self.assertEqual(1, result["provider_resolved_count"])
            self.assertEqual(1, result["not_resolved_count"])
            ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("no_gene_name_to_resolve", ledger.set_index("hypothesis_id").loc["mech_003", "accession_resolution_status"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("reference_gene_resolution_required:no_reference_gene_to_resolve",
                            reviews.set_index("target").loc["mech_003", "reason"])


if __name__ == "__main__":
    unittest.main()
