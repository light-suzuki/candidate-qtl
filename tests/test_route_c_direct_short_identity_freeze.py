"""Exact short identity tokens through real, offline Direct-C source freeze."""
from __future__ import annotations

import json
import csv
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_direct_pipeline import DIRECT_PROPOSAL_COLUMNS, direct_rescue_freeze


class ShortIdentityFreezeTests(unittest.TestCase):
    def test_literal_missing_value_like_identifier_survives_freeze_projection_and_members(self):
        from candidate_qtl.route_c_direct_pipeline import direct_rescue_project
        from candidate_qtl.route_c_direct_members import build_route_c_direct_members

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            quote = "The NA mutation caused altered growth in Synthetic species."
            result, _, _ = self._freeze(root, "NA", [self._record(quote)], extraction=[{
                "citation": "PMID:12345678", "evidence_relation": "direct_genetic_test", "quoted_support": quote}])
            with Path(result["frozen"]).open(encoding="utf-8", newline="") as stream:
                frozen = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual("NA", frozen["reported_gene_symbol"])
            self.assertEqual("verified_gene", frozen["gene_identity_state"])
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "NA", "chrom": "chr1", "start": "1", "end": "100"}]).to_csv(atlas, sep="\t", index=False)
            projected = direct_rescue_project(result["frozen"], result["manifest"], str(atlas), "neutral", str(root / "projection"))
            members = build_route_c_direct_members(result["frozen"], result["manifest"], projected["projection"],
                                                  projected["manifest"], str(atlas), str(root / "members"), case_id="neutral")
            with Path(members["members"]).open(encoding="utf-8", newline="") as stream:
                member = next(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual("NA", member["candidate_gene_id"])
            self.assertEqual(0, members["supported_count"])

    def _freeze(self, root: Path, symbol: str, records: list[dict], *, citations: str = "PMID:12345678", extraction: list[dict] | None = None):
        proposal = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
        proposal.update({
            "record_id": "neutral-c1", "query_phenotype": "altered growth", "reported_phenotype": "altered growth",
            "focal_species": "Synthetic species", "reference_species": "Synthetic species",
            "reported_gene_symbol": symbol, "species_distance": "S1_same_species",
            "evidence_directness": "D0_causal_perturbation_or_rescue", "evidence_type": "neutral experiment", "source_id": citations,
        })
        proposals, cache = root / "proposals.tsv", root / "cache.json"
        pd.DataFrame([proposal], columns=DIRECT_PROPOSAL_COLUMNS).to_csv(proposals, sep="\t", index=False)
        cache.write_text(json.dumps(records), encoding="utf-8")
        extraction_path = None
        if extraction is not None:
            extraction_path = root / "extraction.json"
            extraction_path.write_text(json.dumps(extraction), encoding="utf-8")
        result = direct_rescue_freeze(str(proposals), str(root / "freeze"), str(cache),
                                      claim_extraction_arg=str(extraction_path) if extraction_path else None)
        frozen = pd.read_csv(result["frozen"], sep="\t", dtype=str).fillna("").iloc[0]
        validation = pd.read_csv(root / "freeze/source_verification/source_validation.tsv", sep="\t", dtype=str).fillna("")
        return result, frozen, validation

    def _record(self, text: str, *, title: str = "Neutral source", pmid: str = "12345678", **extra) -> dict:
        return {"id": pmid, "pmid": pmid, "title": title, "abstractText": text, "pubYear": "2020", **extra}

    def test_explicit_two_letter_symbols_match_literal_frozen_source_after_accepted_extraction(self):
        for symbol in ("XY", "ZZ"):
            with self.subTest(symbol=symbol), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                quote = f"The {symbol} allele caused altered growth in Synthetic species."
                result, frozen, validation = self._freeze(root, symbol, [self._record(quote)], extraction=[{
                    "citation": "PMID:12345678", "evidence_relation": "direct_genetic_test", "quoted_support": quote, "locator": "abstract",
                }])
                self.assertEqual("accepted", validation.iloc[0]["extraction_status"])
                self.assertNotIn(symbol.casefold(), validation.iloc[0]["claim_terms_found"].split(";"))
                self.assertEqual("verified_gene", frozen["gene_identity_state"])
                self.assertEqual("eligible_direct_record", frozen["direct_rescue_eligibility"])
                self.assertEqual("none", frozen["rank_contribution"])
                manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
                self.assertEqual(sha256(root / "freeze/source_verification/source_review_records.json"), manifest["source_review_records_sha256"])
                self.assertNotIn("claim_binding", manifest)

    def test_absent_wrong_case_and_embedded_tokens_do_not_verify_short_identity(self):
        for token in ("ZZ", "xy", "AXYB", "XY_suffix", "prefixXY"):
            with self.subTest(token=token), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                quote = f"The {token} allele caused altered growth in Synthetic species."
                _, frozen, validation = self._freeze(root, "XY", [self._record(quote)], extraction=[{
                    "citation": "PMID:12345678", "evidence_relation": "direct_genetic_test", "quoted_support": quote,
                }])
                self.assertEqual("accepted", validation.iloc[0]["extraction_status"])
                self.assertEqual("verified", frozen["source_verification_state"])
                self.assertEqual("unresolved", frozen["gene_identity_state"])
                self.assertEqual("pending_verification", frozen["direct_rescue_eligibility"])

    def test_two_letter_symbol_in_title_can_match_but_journal_or_keywords_cannot(self):
        for in_title in (True, False):
            with self.subTest(in_title=in_title), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                quote = "The mutant caused altered growth; complementation rescued the phenotype."
                record = self._record(quote, title="XY mutant analysis" if in_title else "Neutral mutant analysis", journalTitle="XY", keywordList="XY")
                _, frozen, _ = self._freeze(root, "XY", [record])
                self.assertEqual("verified_gene" if in_title else "unresolved", frozen["gene_identity_state"])

    def test_identity_from_an_unrelated_citation_cannot_complete_an_accepted_relation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            quote = "The ZZ allele caused altered growth in Synthetic species."
            records = [self._record(quote), self._record("XY is a reported identifier.", pmid="87654321")]
            _, frozen, validation = self._freeze(root, "XY", records, citations="PMID:12345678;PMID:87654321", extraction=[{
                "citation": "PMID:12345678", "evidence_relation": "direct_genetic_test", "quoted_support": quote,
            }])
            self.assertEqual("accepted", validation.iloc[0]["extraction_status"])
            self.assertEqual("unresolved", frozen["gene_identity_state"])
            self.assertEqual("pending_verification", frozen["direct_rescue_eligibility"])

    def test_matched_identity_source_supplies_frozen_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [self._record("ZZ knockout caused altered growth.", title="Unrelated accepted source"),
                       self._record("XY knockout caused altered growth.", title="Matching accepted source", pmid="87654321")]
            _, frozen, _ = self._freeze(root, "XY", records, citations="PMID:12345678;PMID:87654321")
            self.assertEqual("verified_gene", frozen["gene_identity_state"])
            self.assertEqual("87654321", frozen["pmid"])
            self.assertEqual("Matching accepted source", frozen["source_title"])

    def test_unresolved_citation_or_short_symbol_without_relation_stays_pending(self):
        for missing in (True, False):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _, frozen, validation = self._freeze(root, "XY", [] if missing else [self._record("XY is a reported identifier with altered growth context.")])
                self.assertNotEqual("eligible_direct_record", frozen["direct_rescue_eligibility"])
                self.assertEqual("unresolved", frozen["gene_identity_state"])
                self.assertNotIn("supported_direct_genetic_test", set(validation["claim_support_status"]))


if __name__ == "__main__":
    unittest.main()
