"""What a source supports, separated from whether it exists (issue #226)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.claim_evidence import classify_claim, gate_extraction
from candidate_qtl.cli import main
from candidate_qtl import source_verify as source_verify_module

MECHANISMS = "\n".join([
    "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids",
    "mech_direct\tmeristem maintenance\tCLV\tCLV1\tPMID:1000001",
    "mech_expression\tmeristem maintenance\tWUS\tWUS\tPMID:1000002",
    "mech_mention\tmeristem maintenance\tCLE\tCLE40\tPMID:1000003",
    "mech_review\tmeristem maintenance\tCRN\tCRN\tPMID:1000004",
]) + "\n"

RECORDS = {
    # a knockout plus a complementation: the strongest relation
    "1000001": {
        "title": "CLV1 controls meristem size in Arabidopsis",
        "abstractText": "A clv1 knockout enlarged the shoot apical meristem, and complementation with the "
                        "wild-type gene rescued the phenotype.",
    },
    # the gene is on-topic, but only its transcript level was measured
    "1000002": {
        "title": "WUS expression during meristem development",
        "abstractText": "Transcript level of WUS was up-regulated in the central zone as measured by qRT-PCR.",
    },
    # the gene is named and nothing is done to it
    "1000003": {
        "title": "A catalogue of Arabidopsis peptide genes including CLE40",
        "abstractText": "We list peptide-encoding loci annotated in the reference genome.",
    },
    "1000004": {
        "title": "Receptor kinases in plant development: a review of recent advances",
        "abstractText": "This review summarises signalling by CRN and related receptor kinases.",
    },
}


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _cache(root: Path, pmids: tuple[str, ...] = tuple(RECORDS)) -> Path:
    cache = root / "pmc_cache"
    cache.mkdir(exist_ok=True)
    for pmid in pmids:
        payload = {"cache_key": f"pmid_{pmid}", "id": pmid, "pmid": pmid,
                   "journalTitle": "Plant Cell", "pubYear": "2020", "authorString": "A B, C D",
                   **RECORDS[pmid]}
        (cache / f"med_{pmid}.json").write_text(json.dumps(payload), encoding="utf-8")
    return cache


class ClassifierTests(unittest.TestCase):
    def test_a_complementation_is_a_direct_genetic_test(self) -> None:
        relation, status, cues = classify_claim(
            "A clv1 knockout enlarged the meristem and complementation rescued the phenotype.",
            ["CLV1"], ["CLV1"], "title_abstract")
        self.assertEqual("direct_genetic_test", relation)
        self.assertEqual("supported_direct_genetic_test", status)
        self.assertIn("complementation", cues)

    def test_an_expression_measurement_is_associative_only(self) -> None:
        relation, status, _ = classify_claim(
            "Transcript level of WUS was up-regulated as measured by qRT-PCR.", ["WUS"], ["WUS"], "title_abstract")
        self.assertEqual("expression_association", relation)
        self.assertEqual("associative_support_only", status)

    def test_naming_a_gene_without_an_experiment_is_review_required(self) -> None:
        relation, status, _ = classify_claim(
            "We list peptide-encoding loci including CLE40.", ["CLE40"], ["CLE40"], "title_abstract")
        self.assertEqual("mention_without_mechanism", relation)
        self.assertEqual("review_required", status)

    def test_a_source_that_never_names_the_subject_is_an_unsupported_mention(self) -> None:
        relation, status, _ = classify_claim(
            "A knockout of an unrelated locus changed root length.", [], ["CLV1"], "title_abstract")
        self.assertEqual("no_claim_term_in_source", relation)
        self.assertEqual("unsupported_mention", status)

    def test_nothing_retrievable_is_not_evaluated_never_a_verdict(self) -> None:
        relation, status, _ = classify_claim("", [], ["CLV1"], "none")
        self.assertEqual("not_evaluated_no_text", relation)
        self.assertEqual("not_evaluated", status)

    def test_cached_cues_preserve_each_claim_gate_and_do_not_share_mutable_results(self) -> None:
        cache = {}
        for text in (
            "CLV1 knockout and complementation rescued the phenotype. RNAi and GWAS were also reported.",
            "CLV1 transcript level was up-regulated. This review gives an overview of the work.",
            "CLV1 is listed in a catalogue.",
            " \n ",
        ):
            for found, terms, scope in (
                (["CLV1"], ["CLV1"], "title_abstract"),
                ([], ["UNRELATED"], "title_abstract"),
                ([], [], "title_only"),
                (["CLV1"], ["CLV1"], "none"),
                (["CLV1"], ["CLV1"], "full_text_cached"),
            ):
                with self.subTest(text=text, terms=terms, scope=scope):
                    expected = classify_claim(text, found, terms, scope)
                    actual = classify_claim(text, found, terms, scope, cue_cache=cache)
                    self.assertEqual(expected, actual)
                    actual[2].append("caller mutation")
                    self.assertEqual(expected, classify_claim(text, found, terms, scope, cue_cache=cache))


class SourceVerifyClaimTests(unittest.TestCase):
    def _run(self, root: Path, *extra: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        out = root / "sources"
        self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(_write(root, "mech.tsv", MECHANISMS)),
                                  "--out-dir", str(out), "--cache", str(_cache(root)), *extra]))
        validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("")
        status = pd.read_csv(out / "hypothesis_source_status.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "source_verification_manifest.json").read_text(encoding="utf-8"))
        return validation, status, manifest

    def test_the_two_axes_are_reported_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, status, manifest = self._run(root)
            by_citation = validation.set_index("citation")
            self.assertEqual("True", by_citation.loc["PMID:1000001", "citation_verified"])
            self.assertEqual("supported_direct_genetic_test",
                             by_citation.loc["PMID:1000001", "claim_support_status"])
            # a real paper that only names the gene: the citation is fine, the claim is not
            self.assertEqual("True", by_citation.loc["PMID:1000003", "citation_verified"])
            self.assertEqual("review_required", by_citation.loc["PMID:1000003", "claim_support_status"])
            self.assertEqual("title_abstract", by_citation.loc["PMID:1000003", "retrieval_scope"])
            self.assertIn("claim_support_status", manifest["claim_axes"])

    def test_a_review_and_an_expression_paper_are_not_mechanistic_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, status, _ = self._run(root)
            by_citation = validation.set_index("citation")
            self.assertEqual("associative_support_only", by_citation.loc["PMID:1000002", "claim_support_status"])
            self.assertEqual("background_or_review_only", by_citation.loc["PMID:1000004", "claim_support_status"])
            by_hypothesis = status.set_index("hypothesis_id")
            self.assertEqual("mechanistic_support", by_hypothesis.loc["mech_direct", "claim_support_summary"])
            self.assertEqual("associative_support_only", by_hypothesis.loc["mech_expression", "claim_support_summary"])
            self.assertEqual("1", by_hypothesis.loc["mech_direct", "mechanistic_support_citation_count"])
            self.assertEqual("0", by_hypothesis.loc["mech_review", "mechanistic_support_citation_count"])

    def test_a_weaker_relation_is_surfaced_for_review_without_being_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, _status, _ = self._run(root)
            reviews = pd.read_csv(root / "sources" / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("claim_support_is_associative_support_only", set(reviews["reason"]))
            # and the citation itself is still verified, not rejected
            self.assertEqual("verified",
                             validation.set_index("citation").loc["PMID:1000002", "verification_status"])


class ExtractionGateTests(unittest.TestCase):
    def test_repeated_quotes_keep_citation_identity_and_refresh_on_a_new_call(self) -> None:
        rows = [
            {"citation": "PMID:1", "evidence_relation": "direct_genetic_test", "quoted_support": quote}
            for quote in ("STRASSE knockout", "mutant phenotype", "never reported")
        ] + [
            {"citation": "PMID:2", "evidence_relation": "direct_genetic_test",
             "quoted_support": "STRASSE knockout"},
        ]
        sources = {"PMID:1": "Straße  knockout\n changed the mutant\tphenotype.",
                   "PMID:2": "A different paper."}
        self.assertEqual(
            ["accepted", "accepted", "rejected_quote_not_in_retrieved_text", "rejected_quote_not_in_retrieved_text"],
            [row["extraction_status"] for row in gate_extraction(rows, sources)])
        sources["PMID:1"] = "A revised record without the quoted material."
        self.assertEqual(
            ["rejected_quote_not_in_retrieved_text"] * 4,
            [row["extraction_status"] for row in gate_extraction(rows, sources)])

    def test_a_quote_that_is_not_in_the_record_is_rejected(self) -> None:
        gated = gate_extraction(
            [{"citation": "PMID:1", "evidence_relation": "direct_genetic_test",
              "quoted_support": "complementation rescued the phenotype"},
             {"citation": "PMID:1", "evidence_relation": "direct_genetic_test",
              "quoted_support": "a sentence the paper never contained"},
             {"citation": "PMID:9", "evidence_relation": "direct_genetic_test", "quoted_support": "x"},
             {"citation": "PMID:1", "evidence_relation": "made_up_relation", "quoted_support": "complementation"},
             {"citation": "", "evidence_relation": "", "quoted_support": ""}],
            {"PMID:1": "A knockout was made and complementation rescued the phenotype."})
        self.assertEqual(["accepted", "rejected_quote_not_in_retrieved_text", "rejected_unknown_citation",
                          "rejected_unknown_relation", "rejected_missing_field"],
                         [row["extraction_status"] for row in gated])

    def test_an_accepted_extraction_annotates_but_an_unlocatable_one_is_reviewed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            extraction = _write(root, "extraction.tsv", "\n".join([
                "citation\tclaim_type\tevidence_relation\tquoted_support\tlocator",
                "PMID:1000003\tcausal\tdirect_genetic_test\tthe knockout reproduced the fasciated phenotype\tabstract",
            ]) + "\n")
            out = root / "sources"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(_write(root, "mech.tsv", MECHANISMS)),
                                      "--out-dir", str(out), "--cache", str(_cache(root)),
                                      "--claim-extraction", str(extraction)]))
            gate = pd.read_csv(out / "claim_extraction_gate.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("rejected_quote_not_in_retrieved_text", gate.iloc[0]["extraction_status"])
            validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").set_index("citation")
            # the fabricated quote did not promote the mention
            self.assertEqual("review_required", validation.loc["PMID:1000003", "claim_support_status"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("claim_extraction_rejected_quote_not_in_retrieved_text", set(reviews["reason"]))

    def test_an_accepted_extraction_promotes_the_hypothesis_rollup_consistently(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            extraction = _write(root, "extraction.tsv", "\n".join([
                "citation\tclaim_type\tevidence_relation\tquoted_support\tlocator",
                "PMID:1000003\tcausal\tdirect_genetic_test\t"
                "We list peptide-encoding loci annotated in the reference genome.\tabstract",
            ]) + "\n")
            out = root / "sources"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(_write(root, "mech.tsv", MECHANISMS)),
                                      "--out-dir", str(out), "--cache", str(_cache(root)),
                                      "--claim-extraction", str(extraction)]))
            validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").set_index("citation")
            self.assertEqual("supported_direct_genetic_test", validation.loc["PMID:1000003", "claim_support_status"])
            status = pd.read_csv(out / "hypothesis_source_status.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            # the accepted extraction re-opens the frozen roll-up
            self.assertEqual("1", status.loc["mech_mention", "mechanistic_support_citation_count"])
            self.assertEqual("mechanistic_support", status.loc["mech_mention", "claim_support_summary"])
            provenance = json.loads(status.loc["mech_mention", "publication_source_provenance"])
            self.assertEqual(["supported_direct_genetic_test"],
                             [row["claim_support_status"] for row in provenance])
            # the frozen record graph and its published digest agree with the promotion
            records = json.loads((out / "verified_sources.json").read_text(encoding="utf-8"))
            verification = records["PMID:1000003"]["candidate_qtl_hypothesis_verifications"]["mech_mention"]
            self.assertEqual("supported_direct_genetic_test", verification["claim_support_status"])
            digest = hashlib.sha256(json.dumps(
                [records["PMID:1000003"]], ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            ).encode("utf-8")).hexdigest()
            self.assertEqual(digest, status.loc["mech_mention", "publication_records_sha256"])

    def test_an_accepted_extraction_replaces_the_stale_weak_review_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            extraction = _write(root, "extraction.tsv", "\n".join([
                "citation\tclaim_type\tevidence_relation\tquoted_support\tlocator",
                "PMID:1000002\tcausal\tdirect_genetic_test\t"
                "Transcript level of WUS was up-regulated in the central zone as measured by qRT-PCR.\tabstract",
            ]) + "\n")
            out = root / "sources"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(_write(root, "mech.tsv", MECHANISMS)),
                                      "--out-dir", str(out), "--cache", str(_cache(root)),
                                      "--claim-extraction", str(extraction)]))
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("claim_support_is_associative_support_only", set(reviews["reason"]))
            status = pd.read_csv(out / "hypothesis_source_status.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            self.assertEqual("1", status.loc["mech_expression", "mechanistic_support_citation_count"])
            self.assertEqual("0", status.loc["mech_expression", "associative_support_citation_count"])
            self.assertEqual("mechanistic_support", status.loc["mech_expression", "claim_support_summary"])


class ClaimBenchmarkTests(unittest.TestCase):
    CASES = "\n".join([
        "case_id\tcitation\treference_gene_name\tgene_family\tpathway_or_process\texpected_claim_support_status\texpected_citation_verified",
        "true_support\tPMID:1000001\tCLV1\tCLV\tmeristem maintenance\tsupported_direct_genetic_test\tTrue",
        "weak_support\tPMID:1000002\tWUS\tWUS\tmeristem maintenance\tassociative_support_only\tTrue",
        "irrelevant_mention\tPMID:1000003\tCLE40\tCLE\tmeristem maintenance\treview_required\tTrue",
        "fabricated\tPMID:9999999\tFAKE1\tFAKE\tmeristem maintenance\tnot_evaluated\tFalse",
    ]) + "\n"

    def test_all_four_label_classes_are_scored_with_their_own_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "bench"
            self.assertEqual(0, main(["claim-verification-benchmark",
                                      "--cases", str(_write(root, "cases.tsv", self.CASES)),
                                      "--out-dir", str(out), "--cache", str(_cache(root)),
                                      "--cache-is-complete"]))
            cases = pd.read_csv(out / "claim_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("correct", cases.loc["true_support", "outcome"])
            self.assertEqual("correct", cases.loc["weak_support", "outcome"])
            self.assertEqual("correct", cases.loc["irrelevant_mention", "outcome"])
            self.assertEqual("False", cases.loc["fabricated", "citation_verified"])
            manifest = json.loads((out / "claim_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(4, manifest["case_count"])
            self.assertEqual(4, manifest["evaluated"])
            self.assertEqual(1.0, manifest["accuracy"])

    def test_without_a_provider_a_case_is_not_evaluated_rather_than_wrong(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "bench"
            self.assertEqual(0, main(["claim-verification-benchmark",
                                      "--cases", str(_write(root, "cases.tsv", self.CASES)),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "claim_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(0, manifest["evaluated"])
            self.assertEqual(4, manifest["not_evaluated"])
            self.assertIsNone(manifest["accuracy"])

    def test_ncbi_fallback_provider_is_recorded_in_benchmark_manifest(self) -> None:
        one_case = self.CASES.splitlines()[0] + "\n" + self.CASES.splitlines()[1] + "\n"
        record = dict(RECORDS["1000001"], _source_provider="ncbi_eutils_pubmed")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "bench"
            with patch.object(source_verify_module, "_fetch_europe_pmc", return_value=None), \
                 patch.object(source_verify_module, "_fetch_ncbi_eutils", return_value=record):
                self.assertEqual(0, main(["claim-verification-benchmark",
                                          "--cases", str(_write(root, "cases.tsv", one_case)),
                                          "--out-dir", str(out), "--allow-external-api"]))
            manifest = json.loads((out / "claim_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("ncbi_eutils_pubmed", manifest["provider"])
            self.assertEqual(["ncbi_eutils_pubmed"], manifest["providers_used"])

    def test_repeated_doi_with_different_case_is_fetched_once(self) -> None:
        header = self.CASES.splitlines()[0]
        cases = "\n".join([
            header,
            "doi-upper\tDOI:10.1000/EXAMPLE\tCLV1\tCLV\tmeristem maintenance\tsupported_direct_genetic_test\tTrue",
            "doi-lower\tDOI:10.1000/example\tCLV1\tCLV\tmeristem maintenance\tsupported_direct_genetic_test\tTrue",
        ]) + "\n"
        record = {
            "doi": "10.1000/example", "title": "CLV1 meristem study",
            "abstractText": "Loss-of-function of CLV1 caused fasciation in the meristem.",
            "journalTitle": "Plant Journal", "pubYear": "2024", "authorString": "Example",
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "bench"
            with patch.object(source_verify_module, "_fetch_europe_pmc", return_value=record) as fetch:
                self.assertEqual(0, main(["claim-verification-benchmark",
                                          "--cases", str(_write(root, "cases.tsv", cases)),
                                          "--out-dir", str(out), "--allow-external-api"]))
            self.assertEqual(1, fetch.call_count)


if __name__ == "__main__":
    unittest.main()
