"""TEST_ONLY packet-to-canonical adapter; synthetic fixtures, no LLM/API/truth."""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path[:0] = [str(HERE), str(REPO), str(REPO / "src")]

import packet_to_canonical as adapter  # noqa: E402
from tests.bound_claim_fixtures import bound_claim  # noqa: E402
from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS  # noqa: E402


ATLAS_ROWS = [
    {"gene_id": "G1", "family_id": "FAM1", "evolutionary_family_id": "EV1",
     "orthogroup": "OG1", "subfamily_id": "SUB1",
     "orthology_reference_gene": "REF1", "domain_architecture": "kinase",
     "broad_function_bins": "signaling"},
    {"gene_id": "G2", "family_id": "FAM1", "evolutionary_family_id": "EV1",
     "orthogroup": "OG1", "subfamily_id": "SUB1",
     "orthology_reference_gene": "REF1", "domain_architecture": "kinase",
     "broad_function_bins": "signaling"},
]

LITERATURE_ROWS = [{
    "interval_gene_ids": "G1;G2", "reference_gene": "REF1",
    "reference_function": "meristem receptor", "query_status": "retrieved",
    "article_rank": "1", "source_id": "PMID:123",
    "title": "Reference receptor controls meristem development",
    "doi": "10.1/example", "journal": "Plant", "year": "2020",
    "matched_phenotype_terms": "fasciation",
    "claim_type": "literature_link_candidate",
}]


def _entry(accession: str, species: str) -> dict:
    return {
        "primaryAccession": accession,
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "entryAudit": {"entryVersion": 7},
        "organism": {"scientificName": species},
        "sequence": {"value": "MAAGSVDPLA"},
        "references": [{"citation": {"citationCrossReferences": [
            {"database": "PubMed", "id": "12345"}]}}],
    }


def _payload_row(**overrides: str) -> dict[str, str]:
    row = {
        "hypothesis_id": "H1",
        "causal_layer": "regulation",
        "pathway_or_process": "meristem maintenance",
        "gene_family": "LRR-RLK",
        "mechanism_statement": "receptor lesion blocks meristem signaling",
        "expected_role": "positive regulator",
        "direction_consistency": "loss reduces meristem activity",
        "model_species_evidence": "mutant study",
        "model_species_mutant_phenotype": "fasciated meristem",
        "source_ids": "PMID:12345;UniProtKB:Q6UD73",
        "proposed_stable_id": "",
        "reference_species": "Arabidopsis thaliana",
        "reference_gene_name": "",
        "selection_rationale": "family bridge",
        "counter_hypothesis": "paralog redundancy",
        "discriminating_observation": "reporter assay",
        "confidence": "low",
        "uncertainty": "fixture",
    }
    row.update(overrides)
    return row


def _source_records(rows: list[dict[str, str]]) -> dict[str, dict]:
    records: dict[str, dict] = {}
    for row in rows:
        source_ids = sorted({
            part.strip() for part in str(row.get("verified_source_ids", "")).split(";")
            if part.strip()
        })
        expected_rows: list[dict[str, str]] = []
        for source_id in source_ids:
            record = records.setdefault(source_id, {
                "citation": source_id,
                "candidate_qtl_hypothesis_verifications": {},
            })
            verifications = record["candidate_qtl_hypothesis_verifications"]
            assert isinstance(verifications, dict)
            verifications[str(row["hypothesis_id"])] = {
                "citation": source_id,
                "resolved_id": source_id.removeprefix("PMID:"),
                "provider": "local_record_cache",
                "url": "frozen-literature-cache.json",
                "retrieved_at": "2026-08-22T00:00:00+00:00",
                "claim_support_status": "supported_direct_genetic_test"
                if row.get("claim_support_summary") == "mechanistic_support"
                else "associative_support_only",
            }
            expected_rows.append(verifications[str(row["hypothesis_id"])])
        current = str(row.get("publication_source_provenance", ""))
        try:
            current_rows = json.loads(current) if current else []
        except json.JSONDecodeError:
            current_rows = []
        if source_ids and current_rows and all(
            isinstance(item, dict) and set(item) == {"citation"}
            for item in current_rows
        ):
            row["publication_source_provenance"] = json.dumps(
                expected_rows, sort_keys=True, separators=(",", ":"))
    return records


def _stamp_hashes(rows: list[dict[str, str]],
                  source_records: dict[str, dict]) -> None:
    for row in rows:
        if row.get("publication_records_sha256") == "AUTO":
            ids = sorted({
                part.strip() for part in str(row.get("verified_source_ids", "")).split(";")
                if part.strip()
            })
            row["publication_records_sha256"] = hashlib.sha256(json.dumps(
                [source_records[i] for i in ids],
                ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            ).encode("utf-8")).hexdigest()


class PacketToCanonicalTests(unittest.TestCase):
    def test_a_import_normalizes_envelope_and_preserves_raw_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            raw = json.dumps({
                "hypotheses": [_payload_row()],
                "research_log": [{"query": "family function",
                                  "direct_answer_exposure": "no"}],
                "direct_answer_exposure": "no",
                "limitations": ["synthetic fixture"],
            }, ensure_ascii=False, indent=2).encode("utf-8") + b"\n"
            payload.write_bytes(raw)
            result = adapter.import_a_payload(
                payload, root / "mech", case_id="case900",
                model="TEST_ONLY/model")
            self.assertEqual(raw, Path(result["raw_a_payload"]).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(),
                             result["raw_a_payload_sha256"])
            provenance = json.loads(Path(
                result["a_payload_envelope_provenance"]).read_text(encoding="utf-8"))
            self.assertEqual(result["raw_a_payload_sha256"],
                             provenance["original_payload_sha256"])
            self.assertEqual("no", provenance["retained_metadata"][
                "direct_answer_exposure"])
            frame = pd.read_csv(result["hypotheses"], sep="\t",
                                dtype=str).fillna("")
            self.assertEqual(1, len(frame))

    def test_a_import_rejects_unexpected_envelope_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({
                "hypotheses": [_payload_row()],
                "unexpected_metadata": "must not be silently dropped",
            }), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unknown top-level key"):
                adapter.import_a_payload(
                    payload, root / "mech", case_id="case900",
                    model="TEST_ONLY/model")

    def test_a_import_accepts_single_limitations_string(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({
                "hypotheses": [_payload_row()],
                "limitations": "Only a conditional hypothesis.",
            }), encoding="utf-8")
            result = adapter.import_a_payload(
                payload, root / "mech", case_id="case900",
                model="TEST_ONLY/model")
            provenance = json.loads(Path(
                result["a_payload_envelope_provenance"]).read_text(encoding="utf-8"))
            self.assertEqual("Only a conditional hypothesis.",
                             provenance["retained_metadata"]["limitations"])

    def test_a_import_preserves_structured_limitations_without_changing_hypotheses(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = [_payload_row()]
            limitations = {"reading_scope": "abstract only", "uncertainty": {
                "mapping": ["multiple copies", "function transfer unresolved"], "new_search": False}}
            results = []
            for label, metadata in (("text", "abstract only"), ("structured", limitations)):
                payload = root / f"{label}.json"
                raw = (json.dumps({"hypotheses": hypotheses, "limitations": metadata}, ensure_ascii=False) + "\n").encode()
                payload.write_bytes(raw)
                result = adapter.import_a_payload(payload, root / label, case_id="case900", model="TEST_ONLY/model")
                provenance = json.loads(Path(result["a_payload_envelope_provenance"]).read_text(encoding="utf-8"))
                self.assertEqual(metadata, provenance["retained_metadata"]["limitations"])
                self.assertEqual(raw, payload.read_bytes())
                self.assertEqual(raw, Path(result["raw_a_payload"]).read_bytes())
                self.assertEqual(hashlib.sha256(raw).hexdigest(), provenance["original_payload_sha256"])
                results.append((provenance["normalized_hypotheses_sha256"], Path(result["hypotheses"]).read_bytes()))
            self.assertEqual(results[0], results[1])

    def test_a_import_structured_limitations_cannot_hide_forbidden_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": [_payload_row()],
                "limitations": {"details": [{"candidate_answer": "leak"}]}}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "forbidden nested key"):
                adapter.import_a_payload(payload, root / "imported", case_id="case900", model="TEST_ONLY/model")
            self.assertFalse((root / "imported").exists())

    def test_a_import_rejects_nested_leakage_keys_in_envelope_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({
                "hypotheses": [_payload_row()],
                "research_log": [{"query": "family function",
                                  "candidate_answer": "leak"}],
            }), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "forbidden nested key"):
                adapter.import_a_payload(
                    payload, root / "mech", case_id="case900",
                    model="TEST_ONLY/model")

    def test_a_import_preserves_proposal_without_invention(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": [_payload_row()]}),
                               encoding="utf-8")
            result = adapter.import_a_payload(
                str(payload), str(root / "mech"), case_id="case900",
                model="TEST_ONLY/model")
            frame = pd.read_csv(result["hypotheses"], sep="\t",
                                dtype=str).fillna("")
            self.assertEqual(1, len(frame))
            self.assertEqual("LRR-RLK", frame.iloc[0]["gene_family"])
            self.assertEqual("Arabidopsis thaliana",
                             frame.iloc[0]["reference_species"])

    def test_a_offline_verify_freezes_review_only_without_bound_claims(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": [_payload_row(
                proposed_stable_id="Q6UD73",
            )]}), encoding="utf-8")
            imported = adapter.import_a_payload(
                str(payload), str(root / "mech"), case_id="case900",
                model="TEST_ONLY/model")
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(
                json.dumps(_entry("Q6UD73", "Arabidopsis thaliana")),
                encoding="utf-8")
            result = adapter.verify_and_freeze_a(
                str(imported["hypotheses"]), str(root / "a"),
                focal_species="Pisum sativum", record_cache=str(cache),
                case_id="case900")
            panel = pd.read_csv(result["verified_panel"], sep="\t",
                                dtype=str).fillna("")
            self.assertEqual(0, len(panel))  # no bound authorization: no support
            projection = pd.read_csv(result["projection_panel"], sep="\t",
                                     dtype=str).fillna("")
            self.assertEqual(1, len(projection))
            blob = " ".join(projection.iloc[0].astype(str).tolist())
            self.assertIn("review_only", blob)

    def test_a_verified_panel_requires_bound_claim_authorization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [_payload_row(
                reference_gene_name="",
                proposed_stable_id="Q6UD73",
            )]
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": rows}),
                               encoding="utf-8")
            imported = adapter.import_a_payload(
                str(payload), str(root / "mech"), case_id="case900",
                model="TEST_ONLY/model")
            # Verify-gate columns join after import from frozen source
            # records; the import schema itself carries no gate fields.
            gate = {
                "verified_source_ids": "PMID:12345",
                "verified_citation_count": "1",
                "claim_support_summary": "mechanistic_support",
                "mechanistic_support_citation_count": "1",
                "publication_source_provenance": '[{"citation":"PMID:12345"}]',
                "publication_records_sha256": "AUTO",
            }
            frame = pd.read_csv(str(imported["hypotheses"]), sep="\t",
                                dtype=str).fillna("")
            for key, value in gate.items():
                frame[key] = value
            gated = frame.to_dict("records")
            source_records = _source_records(gated)
            _stamp_hashes(gated, source_records)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame(gated).to_csv(hypotheses, sep="\t", index=False)
            (root / "source_records.json").write_text(
                json.dumps(source_records), encoding="utf-8")
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(
                json.dumps(_entry("Q6UD73", "Arabidopsis thaliana")),
                encoding="utf-8")
            bound = root / "bound.tsv"
            pd.DataFrame([bound_claim(
                route="a", claim_id="H1", citation="PMID:12345",
                subject_gene="Q6UD73", organism="Arabidopsis thaliana",
                phenotype="meristem maintenance",
            )], columns=BOUND_CLAIM_COLUMNS).to_csv(
                bound, sep="\t", index=False, lineterminator="\n")
            result = adapter.verify_and_freeze_a(
                str(hypotheses), str(root / "a"),
                focal_species="Pisum sativum", record_cache=str(cache),
                source_records=str(root / "source_records.json"),
                bound_claims=str(bound), case_id="case900")
            panel = pd.read_csv(result["verified_panel"], sep="\t",
                                dtype=str).fillna("")
            self.assertEqual(["Q6UD73"], panel["bait_id"].tolist())

    def _packet_b(self, root: Path,
                  gene_ids: list[str] | None = None) -> dict[str, object]:
        packet = root / "packet-b"
        packet.mkdir()
        ids = gene_ids or ["G1", "G2"]
        with (packet / "candidates.tsv").open("w", encoding="utf-8",
                                              newline="") as handle:
            handle.write("gene_id\tchrom\tstart\tend\n")
            for i, gene_id in enumerate(ids):
                handle.write(f"{gene_id}\tchr1\t{1000 + i}\t{2000 + i}\n")
        (packet / "case.json").write_text(json.dumps(
            {"case_id": "case900", "interval": "chr1:1-9999", "candidate_count": len(ids)}),
            encoding="utf-8")
        (packet / "phenotype.json").write_text(json.dumps(
            {"case_id": "case900", "focal_species": "Pisum sativum",
             "phenotype": "fasciated meristem"}), encoding="utf-8")
        return adapter.read_packet_b(str(packet))

    def _bundle(self, root: Path) -> dict[str, object]:
        src = root / "src"
        src.mkdir()
        pd.DataFrame(ATLAS_ROWS).to_csv(src / "atlas.tsv", sep="\t",
                                        index=False, lineterminator="\n")
        pd.DataFrame(LITERATURE_ROWS).to_csv(src / "literature.tsv", sep="\t",
                                             index=False, lineterminator="\n")
        return src

    def _responses(self, prompts_path: str, claims: list[dict],
                   out_path: Path) -> None:
        entries = []
        for line in Path(prompts_path).read_text(encoding="utf-8").splitlines():
            prompt = json.loads(line)
            entries.append({
                "prompt_id": prompt["prompt_id"],
                "gene_id": prompt["gene_id"],
                "claims": claims,
                "model_json": {"claims": claims},
                "raw_text": json.dumps({"claims": claims}),
                "parse_status": "json",
                "review_status": "not_reviewed",
            })
        out_path.write_text(json.dumps({"candidates": entries},
                                       ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")

    def test_b_packet_requires_case_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            packet = self._packet_b(root)
            case_path = Path(str(packet["dir"])) / "case.json"
            metadata = json.loads(case_path.read_text(encoding="utf-8"))
            metadata.pop("case_id")
            case_path.write_text(json.dumps(metadata), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "case_id"):
                adapter.read_packet_b(str(packet["dir"]))

    def test_b_packet_builds_phenotype_free_atlas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            packet = self._packet_b(root)
            self.assertEqual(["G1", "G2"], packet["gene_ids"])
            result = adapter.build_atlas(packet, str(root / "atlas"),
                                         case_id="case900")
            frame = pd.read_csv(str(result["atlas"]), sep="\t",
                                dtype=str).fillna("")
            # Every packet gene lands in the atlas; atlas takes no phenotype.
            self.assertEqual({"G1", "G2"}, set(frame["gene_id"].tolist()))

    def test_b_offline_chain_freezes_review_without_support_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = self._bundle(root)
            bundle = adapter.build_bundle(
                str(src / "atlas.tsv"), str(root / "bundle"),
                phenotype="fasciated meristem",
                literature=str(src / "literature.tsv"), case_id="case900")
            # Source-free retention is the unversioned replay contract, not
            # initial_bibliography_v1 admission. Seal this choice before import.
            from tests._route_b_policy_fixtures import unversioned_bundle
            unversioned_bundle(bundle)
            responses = root / "responses.json"
            self._responses(str(bundle["prompts"]), [{
                "claim_type": "mutant_phenotype_bridge",
                "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling",
                "phenotype_relation": "could explain fasciation",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "family",
                "source_ids": ["PMID:123"],
                "alternative_hypothesis": "related receptor copy",
                "uncertainty_note": "transfer unresolved",
            }], responses)
            result = adapter.import_verify_freeze_b(
                str(bundle["bundle"]), str(bundle["prompts"]), str(responses),
                str(root / "b"), case_id="case900")
            # Offline with no source cache: nothing verified, nothing invented.
            # Both copies land in the review lane as an unresolved set.
            verified = pd.read_csv(result["candidates"], sep="\t",
                                   dtype=str).fillna("")
            self.assertEqual(0, len(verified))
            review = pd.read_csv(result["review_candidates"], sep="\t",
                                 dtype=str).fillna("")
            self.assertEqual({"G1", "G2"},
                             set(review["candidate_gene_id"].tolist()))
            self.assertEqual({"unverified_route_b_review"},
                             set(review["support_status"].tolist()))
            self.assertEqual({"none"},
                             set(review["rank_contribution"].tolist()))

    def test_b_cached_citation_is_checked_without_becoming_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = self._bundle(root)
            bundle = adapter.build_bundle(
                src / "atlas.tsv", root / "bundle", phenotype="fasciation",
                literature=str(src / "literature.tsv"), case_id="case900")
            responses = root / "responses.json"
            self._responses(str(bundle["prompts"]), [{
                "claim_type": "mutant_phenotype_bridge",
                "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling",
                "phenotype_relation": "could explain fasciation",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "family", "source_ids": ["PMID:123"],
                "alternative_hypothesis": "related receptor copy",
                "uncertainty_note": "transfer unresolved",
            }], responses)
            cache = root / "cache"
            cache.mkdir()
            (cache / "med_123.json").write_text(json.dumps({
                "id": "123", "pmid": "123",
                "title": "REF1 controls meristem development",
                "abstractText": "Loss-of-function knockout mutants of REF1 display fasciation.",
                "journalTitle": "Fixture", "pubYear": "2020",
            }), encoding="utf-8")
            result = adapter.import_verify_freeze_b(
                bundle["bundle"], bundle["prompts"], responses, root / "b",
                case_id="case900", source_cache=str(cache))
            source = pd.read_csv(result["source_verification"], sep="\t", dtype=str)
            self.assertEqual({"true"}, set(source["citation_verified"].str.lower()))
            self.assertEqual(0, len(pd.read_csv(result["candidates"], sep="\t")))
            self.assertEqual({"G1", "G2"}, set(pd.read_csv(
                result["review_candidates"], sep="\t")["candidate_gene_id"]))

    def test_b_review_requires_real_corresponding_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = self._bundle(root)
            bundle = adapter.build_bundle(
                src / "atlas.tsv", root / "bundle", phenotype="fasciation",
                literature=str(src / "literature.tsv"), case_id="case900")
            responses = root / "responses.json"
            self._responses(str(bundle["prompts"]), [], responses)
            with self.assertRaisesRegex(ValueError, "requires --source-review-records"):
                adapter.import_verify_freeze_b(
                    bundle["bundle"], bundle["prompts"], responses, root / "b",
                    agent_claim_review=str(root / "invented-review.tsv"))

    def test_handoff_keeps_support_review_fill_separate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            packet = self._packet_b(root, gene_ids=["G1", "G2", "G3", "G4"])
            atlas = adapter.build_atlas(packet, str(root / "atlas"),
                                        case_id="case900")
            # B chain on the family-annotated fixture atlas; its G1/G2 members
            # sit inside the packet atlas universe, so the handoff can verify
            # interval membership end to end.
            src = self._bundle(root)
            bundle = adapter.build_bundle(
                str(src / "atlas.tsv"), str(root / "bundle"),
                phenotype="fasciated meristem",
                literature=str(src / "literature.tsv"), case_id="case900")
            from tests._route_b_policy_fixtures import unversioned_bundle
            unversioned_bundle(bundle)
            responses = root / "responses.json"
            self._responses(str(bundle["prompts"]), [{
                "claim_type": "mutant_phenotype_bridge",
                "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling",
                "phenotype_relation": "could explain fasciation",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "family",
                "source_ids": ["PMID:123"],
                "alternative_hypothesis": "related receptor copy",
                "uncertainty_note": "transfer unresolved",
            }], responses)
            frozen = adapter.import_verify_freeze_b(
                str(bundle["bundle"]), str(bundle["prompts"]), str(responses),
                str(root / "b"), case_id="case900")
            result = adapter.handoff_to_wet(
                str(root / "wet"), case_id="case900", wet_budget=4,
                route_b_claims=str(frozen["transfer_claims"]),
                route_b_claims_manifest=str(frozen["transfer_manifest"]),
                interval_atlas=str(atlas["atlas"]),
                review_lanes=[str(frozen["review_candidates"])],
                review_manifests=[str(frozen["review_manifest"])])
            wet = pd.read_csv(result["wet_candidates"], sep="\t",
                              dtype=str).fillna("")
            id_column = ("gene_id" if "gene_id" in wet.columns
                         else "candidate_gene_id")
            self.assertLessEqual(len(wet), 4)  # budget is a target, never padded
            self.assertIn("G1", wet[id_column].tolist())
            self.assertIn("G2", wet[id_column].tolist())
            for _, row in wet.iterrows():
                basis = str(row.get("membership_basis", ""))
                self.assertTrue(basis)  # every row labeled
                if basis == "interval_review_coverage":
                    support = " ".join(str(value) for key, value in
                                       row.items() if "support" in key.lower())
                    self.assertNotIn("verified", support.lower())


if __name__ == "__main__":
    unittest.main()
