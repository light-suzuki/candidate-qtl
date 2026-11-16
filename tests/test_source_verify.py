"""Literature citations are verified, not just counted (issue #215 follow-up)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.source_verify import classify_citation
from candidate_qtl import source_verify as source_verify_module

MECHANISMS = "\n".join([
    "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids",
    "mech_001\tstarch biosynthesis\tSBE\tSBEI\tPMID:1978859",
    "mech_002\tsucrose transport\tSUT\tSUC2\tPMID:29999999;UniProt:Q9SXK1",
    "mech_003\tcell wall biosynthesis\tCESA\tCESA9\tnot-a-citation",
    "mech_004\tmeristem maintenance\tCLV\tCLV1\t",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _cache(root: Path) -> Path:
    """A Europe PMC-shaped cache: only PMID:1978859 exists, and its abstract
    mentions the branching enzyme it was cited for."""
    cache = root / "pmc_cache"
    cache.mkdir()
    (cache / "med_1978859.json").write_text(json.dumps({
        "cache_key": "pmid_1978859", "id": "1978859", "pmid": "1978859",
        "title": "Molecular analysis of the SBE locus and starch branching enzyme I in pea",
        "abstractText": "Loss of starch branching enzyme (SBE) activity gives the wrinkled seed phenotype.",
        "journalTitle": "Cell", "pubYear": "1990", "doi": "10.1016/0092-8674(90)90441-h",
        "authorString": "Bhattacharyya MK, Smith AM, Ellis TH",
    }), encoding="utf-8")
    return cache


class CitationClassificationTests(unittest.TestCase):
    def test_citation_forms(self) -> None:
        self.assertEqual(("pmid", "1978859"), classify_citation("PMID:1978859"))
        self.assertEqual(("pmid", "1978859"), classify_citation("PubMed:1978859"))
        self.assertEqual(("pmid", "1978859"), classify_citation(" pubmed 1978859 "))
        self.assertEqual(("pmid", "1978859"), classify_citation("1978859"))
        self.assertEqual(("pmcid", "PMC123456"), classify_citation("pmc123456"))
        self.assertEqual(("pmcid", "PMC123456"), classify_citation("PMC:PMC123456"))
        self.assertEqual(("pmcid", "PMC123456"), classify_citation(" pmc: pmc123456 "))
        self.assertEqual(("pmcid", "PMC123456"), classify_citation("PMCID:PMC123456"))
        self.assertEqual(("doi", "10.1016/0092-8674(90)90441-h"),
                         classify_citation("https://doi.org/10.1016/0092-8674(90)90441-h"))
        self.assertEqual(("doi", "10.1000/example"), classify_citation("DOI:10.1000/EXAMPLE"))
        self.assertEqual("database_accession", classify_citation("UniProt:Q9SXK1")[0])
        self.assertEqual("database_accession", classify_citation("UniProtKB:Q6UD73")[0])
        self.assertEqual("unrecognised", classify_citation("Smith et al. 1990")[0])
        self.assertEqual("unrecognised", classify_citation("PubMed:19x8859")[0])
        self.assertEqual("unrecognised", classify_citation("PubMed:1978859 extra")[0])
        self.assertEqual("unrecognised", classify_citation("PMC:PMC123456 extra")[0])
        self.assertEqual("unrecognised", classify_citation("PMC:123456")[0])



class SourceVerifyTests(unittest.TestCase):
    def test_pmc_namespace_alias_still_requires_a_publication_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = "\n".join([
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids",
                "mech_alias\tstarch biosynthesis\tSBE\tSBEI\tPMC:PMC123456",
            ]) + "\n"
            validation, _ = self._run_with_content(root, mechanisms)
            row = validation.iloc[0]
            self.assertEqual("PMC:PMC123456", row["citation"])
            self.assertEqual("pmcid", row["citation_type"])
            self.assertEqual("PMCID:PMC123456", row["citation_query"])
            self.assertEqual("review_required", row["verification_status"])
            self.assertEqual("verification_not_performed_offline", row["review_reason"])

    def _run(self, root: Path, *extra: str) -> tuple[pd.DataFrame, dict]:
        out = root / "sources"
        self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(_write(root, "mech.tsv", MECHANISMS)),
                                  "--out-dir", str(out), *extra]))
        validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "source_verification_manifest.json").read_text(encoding="utf-8"))
        return validation, manifest

    def test_a_fabricated_citation_is_rejected_and_a_real_one_verified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, manifest = self._run(root, "--cache", str(_cache(root)))
            by_citation = validation.set_index("citation")
            self.assertEqual("verified", by_citation.loc["PMID:1978859", "verification_status"])
            self.assertEqual("supported", by_citation.loc["PMID:1978859", "claim_support"])
            self.assertIn("SBE", by_citation.loc["PMID:1978859", "claim_terms_found"])
            self.assertEqual("1990", by_citation.loc["PMID:1978859", "year"])
            # not in the cache and no API allowed -> review_required, never a rejection
            self.assertEqual("review_required", by_citation.loc["PMID:29999999", "verification_status"])
            self.assertEqual("verification_not_performed_offline",
                             by_citation.loc["PMID:29999999", "review_reason"])
            # an accession is not a literature citation
            self.assertEqual("not_literature_source", by_citation.loc["UniProt:Q9SXK1", "verification_status"])
            # a free-text "citation" cannot be checked at all
            self.assertEqual("rejected", by_citation.loc["not-a-citation", "verification_status"])
            self.assertEqual("citation_is_not_a_recognisable_pmid_pmcid_or_doi",
                             by_citation.loc["not-a-citation", "rejection_reason"])
            self.assertEqual("local_record_cache", manifest["provider"])
            self.assertEqual("none", manifest["external_transmission"])

    def test_pubmed_alias_in_multi_source_row_resolves_as_pmid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = "\n".join([
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids",
                "mech_alias\tstarch biosynthesis\tSBE\tSBEI\tPubMed:1978859; UniProt:Q9SXK1",
            ]) + "\n"
            validation, _ = self._run_with_content(root, mechanisms, "--cache", str(_cache(root)))
            by_citation = validation.set_index("citation")
            self.assertEqual("pmid", by_citation.loc["PubMed:1978859", "citation_type"])
            self.assertEqual("verified", by_citation.loc["PubMed:1978859", "verification_status"])
            self.assertEqual("not_literature_source", by_citation.loc["UniProt:Q9SXK1", "verification_status"])

    def _run_with_content(self, root: Path, mechanisms: str, *extra: str) -> tuple[pd.DataFrame, dict]:
        out = root / "sources"
        self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms",
                                  str(_write(root, "mech.tsv", mechanisms)),
                                  "--out-dir", str(out), *extra]))
        validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "source_verification_manifest.json").read_text(encoding="utf-8"))
        return validation, manifest

    def test_europe_pmc_invalid_json_shape_is_provider_failure(self) -> None:
        class Response(BytesIO):
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.close()
                return False

        with patch.object(source_verify_module.urllib.request, "urlopen",
                          return_value=Response(b"[]")):
            with self.assertRaisesRegex(RuntimeError, "invalid payload shape"):
                source_verify_module._fetch_europe_pmc("EXT_ID:1")
        with patch.object(source_verify_module.urllib.request, "urlopen",
                          return_value=Response(b"{}")):
            with self.assertRaisesRegex(RuntimeError, "invalid result list"):
                source_verify_module._fetch_europe_pmc("EXT_ID:1")

    def test_ncbi_esearch_invalid_json_shape_is_provider_failure(self) -> None:
        class Response(BytesIO):
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.close()
                return False

        with patch.object(source_verify_module.urllib.request, "urlopen",
                          return_value=Response(b'{"esearchresult":{"idlist":null}}')):
            with self.assertRaisesRegex(RuntimeError, "invalid result list"):
                source_verify_module._fetch_ncbi_eutils("pmcid", "PMC123456")
        with patch.object(source_verify_module.urllib.request, "urlopen",
                          return_value=Response(b"{}")):
            with self.assertRaisesRegex(RuntimeError, "invalid result list"):
                source_verify_module._fetch_ncbi_eutils("pmcid", "PMC123456")

    def test_hypothesis_level_status_and_downstream_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, _ = self._run(root, "--cache", str(_cache(root)))
            status = pd.read_csv(root / "sources" / "hypothesis_source_status.tsv", sep="\t", dtype=str) \
                .fillna("").set_index("hypothesis_id")
            self.assertEqual("all_verified", status.loc["mech_001", "source_verification_status"])
            self.assertEqual("none_verified", status.loc["mech_002", "source_verification_status"])
            self.assertEqual("none_verified", status.loc["mech_003", "source_verification_status"])
            self.assertEqual("no_literature_citation", status.loc["mech_004", "source_verification_status"])
            self.assertEqual("PMID:1978859", status.loc["mech_001", "verified_source_ids"])
            verified_sources = json.loads(
                (root / "sources" / "verified_sources.json").read_text(encoding="utf-8"))
            source_record = verified_sources["PMID:1978859"]
            binding = source_record["candidate_qtl_hypothesis_verifications"]["mech_001"]
            self.assertEqual("local_record_cache", binding["provider"])
            self.assertEqual(
                validation.set_index("citation").loc["PMID:1978859", "claim_support_status"],
                binding["claim_support_status"],
            )
            expected_hash = hashlib.sha256(json.dumps(
                [source_record], ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            ).encode("utf-8")).hexdigest()
            self.assertEqual(expected_hash, status.loc["mech_001", "publication_records_sha256"])
            annotated = pd.read_csv(root / "sources" / "mechanisms_with_verified_sources.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("source_verification_status", annotated.columns)
            self.assertEqual(4, len(annotated))

    def test_shared_sources_keep_claim_terms_and_hash_the_complete_extracted_graph(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = _write(root, "cache.json", json.dumps({
                "id": "12345", "pmid": "12345", "title": "ABC1 and DEF2 transcript level",
                "abstractText": "ABC1 and DEF2 were up-regulated.",
                "fullText": "ONLYBODY was also discussed.",
            }))
            mechanisms = "\n".join([
                "hypothesis_id\treference_gene_name\tsource_ids",
                "h1\tABC1\tPMID:12345",
                "h2\tDEF2\tPMID:12345",
                "h3\tONLYBODY\tPMID:12345",
                "h4\t\tPMID:12345",
                "h5\tABC1\tPMID:12345;PMID:12345",
            ]) + "\n"
            extraction = _write(root, "extraction.json", json.dumps([{
                "citation": "PMID:12345", "evidence_relation": "direct_genetic_test",
                "quoted_support": "ABC1 and DEF2 were up-regulated.",
            }]))
            validation, _ = self._run_with_content(
                root, mechanisms, "--cache", str(cache), "--claim-extraction", str(extraction))
            by_hypothesis = validation.set_index("hypothesis_id")
            self.assertEqual("ABC1", by_hypothesis.loc["h1", "claim_terms_found"])
            self.assertEqual("DEF2", by_hypothesis.loc["h2", "claim_terms_found"])
            # The claim-term gate deliberately sees title/abstract/metadata,
            # even when the captured full text contains another subject.
            self.assertEqual("unsupported", by_hypothesis.loc["h3", "claim_support"])
            self.assertEqual("review_required", by_hypothesis.loc["h3", "verification_status"])
            self.assertEqual("not_checked_no_claim_term", by_hypothesis.loc["h4", "claim_support"])
            out = root / "sources"
            status = pd.read_csv(out / "hypothesis_source_status.tsv", sep="\t", dtype=str).fillna("")
            records = json.loads((out / "verified_sources.json").read_text(encoding="utf-8"))
            verifications = records["PMID:12345"]["candidate_qtl_hypothesis_verifications"]
            self.assertEqual({"h1", "h2", "h4", "h5"}, set(verifications))
            self.assertEqual({"supported_direct_genetic_test"},
                             {row["claim_support_status"] for row in verifications.values()})
            for row in status.to_dict("records"):
                ids = sorted(value for value in row["verified_source_ids"].split(";") if value)
                expected = hashlib.sha256(json.dumps(
                    [records[value] for value in ids], ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                ).encode("utf-8")).hexdigest() if ids else ""
                self.assertEqual(expected, row["publication_records_sha256"])
            digests = status.set_index("hypothesis_id")["publication_records_sha256"]
            self.assertEqual(digests["h1"], digests["h2"])
            self.assertNotEqual(digests["h1"], digests["h5"])

    def test_source_preparation_is_not_reused_across_invocations(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = "hypothesis_id\treference_gene_name\tsource_ids\nh1\tABC1\tPMID:12345\n"
            cache = root / "cache.json"
            for title, expected in (
                ("ABC1 knockout", "supported_direct_genetic_test"),
                ("ABC1 transcript level", "associative_support_only"),
                ("Different subject knockout", "unsupported_mention"),
            ):
                with self.subTest(title=title):
                    cache.write_text(json.dumps({"id": "12345", "pmid": "12345", "title": title}), encoding="utf-8")
                    validation, _ = self._run_with_content(root, mechanisms, "--cache", str(cache))
                    self.assertEqual(expected, validation.iloc[0]["claim_support_status"])

    def test_unsourced_and_unresolved_hypotheses_are_reviewed_not_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._run(root, "--cache", str(_cache(root)))
            reasons = set(pd.read_csv(root / "sources" / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("hypothesis_without_literature_citation", reasons)
            self.assertIn("hypothesis_citations_all_unresolved", reasons)
            self.assertIn("citation_unparsable", reasons)
            # the citation table has a row per citation (mech_004 has none), while
            # every hypothesis keeps a status row: nothing is dropped for being unsourced
            validation = pd.read_csv(root / "sources" / "source_validation.tsv", sep="\t", dtype=str).fillna("")
            status = pd.read_csv(root / "sources" / "hypothesis_source_status.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(set(validation["hypothesis_id"])))
            self.assertEqual(4, len(status))

    def test_a_citation_that_does_not_mention_the_claim_is_reviewed_not_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = _cache(root)
            (cache / "med_2222222.json").write_text(json.dumps({
                "cache_key": "pmid_2222222", "id": "2222222", "pmid": "2222222",
                "title": "Nitrogen fixation in root nodules", "abstractText": "Unrelated to the cited claim.",
                "journalTitle": "Plant Cell", "pubYear": "2001", "doi": "10.1/x", "authorString": "A B",
            }), encoding="utf-8")
            mechanisms = _write(root, "mech2.tsv", "\n".join([
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids",
                "mech_001\tstarch biosynthesis\tSBE\tSBEI\tPMID:2222222",
            ]) + "\n")
            out = root / "sources2"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanisms), "--out-dir", str(out),
                                      "--cache", str(cache)]))
            row = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("review_required", row["verification_status"])
            self.assertEqual("unsupported", row["claim_support"])
            self.assertEqual("citation_resolves_but_does_not_mention_the_claimed_subject", row["review_reason"])
            self.assertEqual("", row["rejection_reason"])

    def test_offline_without_a_cache_never_rejects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            validation, manifest = self._run(root)
            statuses = set(validation[validation["citation_type"].isin(["pmid", "doi", "pmcid"])]["verification_status"])
            self.assertEqual({"review_required"}, statuses)
            self.assertEqual("none", manifest["provider"])
            self.assertEqual(0, manifest["verified_citation_count"])

    def test_cache_doi_matching_is_case_insensitive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            cache.mkdir()
            (cache / "doi.json").write_text(json.dumps({
                "doi": "10.1000/EXAMPLE", "title": "CLV1 meristem study",
                "abstractText": "Loss-of-function of CLV1 caused fasciation in the meristem.",
                "journalTitle": "Plant Journal", "pubYear": "2024", "authorString": "Example",
            }), encoding="utf-8")
            mechanisms = _write(
                root, "doi.tsv",
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids\n"
                "h1\tmeristem\tCLV\tCLV1\tDOI:10.1000/example\n",
            )
            out = root / "out"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanisms),
                                      "--out-dir", str(out), "--cache", str(cache)]))
            row = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("verified", row["verification_status"])
            self.assertEqual("local_record_cache", row["source_provider"])

    def test_transport_timeout_is_converted_to_provider_failure(self) -> None:
        with patch.object(source_verify_module.urllib.request, "urlopen", side_effect=TimeoutError("timed out")):
            with self.assertRaisesRegex(RuntimeError, "Europe PMC connection failed"):
                source_verify_module._fetch_europe_pmc("EXT_ID:1 AND SRC:MED")
            with self.assertRaisesRegex(RuntimeError, "NCBI EFetch connection failed"):
                source_verify_module._fetch_ncbi_eutils("pmid", "1")

    def test_ncbi_fallback_supplies_abstract_and_provider_provenance(self) -> None:
        class Response(BytesIO):
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.close()
                return False

        xml = b"""<PubmedArticleSet><PubmedArticle><MedlineCitation><Article><ArticleTitle>CLV1 meristem study</ArticleTitle><Abstract><AbstractText>Loss-of-function of CLV1 caused fasciation.</AbstractText></Abstract><Journal><Title>Plant Journal</Title></Journal><PubDate><Year>2024</Year></PubDate><AuthorList><Author><LastName>Example</LastName></Author></AuthorList></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType='pubmed'>29999999</ArticleId><ArticleId IdType='doi'>10.1000/example</ArticleId></ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>"""
        mechanisms = "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids\n" \
            "h1\tmeristem\tCLV\tCLV1\tPMID:29999999\n"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanism_path = _write(root, "ncbi.tsv", mechanisms)
            out = root / "out"
            with patch.object(source_verify_module, "_fetch_europe_pmc", return_value=None), \
                 patch.object(source_verify_module.urllib.request, "urlopen", return_value=Response(xml)):
                self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanism_path),
                                          "--out-dir", str(out), "--allow-external-api"]))
            row = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("ncbi_eutils_pubmed", row["source_provider"])
            self.assertEqual("verified", row["verification_status"])
            self.assertEqual("2024", row["year"])
            self.assertEqual("title_abstract", row["retrieval_scope"])
            manifest = json.loads((out / "source_verification_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["ncbi_eutils_pubmed"], manifest["providers_used"])

    def test_ncbi_pmc_fallback_uses_pmc_provider_provenance(self) -> None:
        class Response(BytesIO):
            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.close()
                return False

        search = b'{"esearchresult":{"idlist":["123456"]}}'
        article = b"""<pmc-articleset><article><front><article-meta><article-id pub-id-type='pmc'>PMC123456</article-id><article-id pub-id-type='pmid'>29999999</article-id><title-group><article-title>CLV1 meristem study</article-title></title-group><abstract>Loss-of-function of CLV1 caused fasciation.</abstract></article-meta></front></article></pmc-articleset>"""
        with patch.object(
            source_verify_module.urllib.request, "urlopen",
            side_effect=[Response(search), Response(article)],
        ):
            record = source_verify_module._fetch_ncbi_eutils("pmcid", "PMC123456")
        self.assertIsNotNone(record)
        self.assertEqual("ncbi_eutils_pmc", record["_source_provider"])
        self.assertEqual("PMC123456", record["pmcid"])

    def test_repeated_citation_is_fetched_once_per_run(self) -> None:
        mechanisms = (
            "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids\n"
            "h1\tmeristem\tCLV\tCLV1\tPMID:29999999\n"
            "h2\tmeristem maintenance\tCLV\tCLV1\tPMID:29999999\n"
        )
        record = {
            "id": "29999999", "pmid": "29999999", "title": "CLV1 meristem study",
            "abstractText": "Loss-of-function of CLV1 caused fasciation in the meristem.",
            "journalTitle": "Plant Journal", "pubYear": "2024", "authorString": "Example",
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanism_path = _write(root, "repeated.tsv", mechanisms)
            out = root / "out"
            with patch.object(source_verify_module, "_fetch_europe_pmc", return_value=record) as fetch:
                self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanism_path),
                                          "--out-dir", str(out), "--allow-external-api"]))
            self.assertEqual(1, fetch.call_count)
            validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(validation))
            self.assertEqual({"europe_pmc_rest"}, set(validation["source_provider"]))

    def test_clean_ncbi_no_hit_after_europe_failure_is_not_provider_failure(self) -> None:
        mechanisms = (
            "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tsource_ids\n"
            "h1\tmeristem\tCLV\tCLV1\tPMID:29999999\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanism_path = _write(root, "nohit.tsv", mechanisms)
            out = root / "out"
            with patch.object(source_verify_module, "_fetch_europe_pmc", side_effect=RuntimeError("down")), \
                 patch.object(source_verify_module, "_fetch_ncbi_eutils", return_value=None):
                self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanism_path),
                                          "--out-dir", str(out), "--allow-external-api"]))
            row = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("rejected", row["verification_status"])
            self.assertEqual("citation_not_found_in_literature_database", row["rejection_reason"])
            self.assertNotEqual("literature_provider_unavailable", row["review_reason"])
            self.assertEqual("europe_pmc_rest+ncbi_eutils_pubmed", row["source_provider"])
            manifest = json.loads((out / "source_verification_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["europe_pmc_rest", "ncbi_eutils_pubmed"], manifest["providers_used"])


class AgentPayloadSourceTests(unittest.TestCase):
    def test_agent_is_a_first_class_payload_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload.tsv", "\n".join([
                "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
                "discriminating_observation\tmodel_species_evidence\tproposed_stable_id\treference_species",
                "starch biosynthesis\tSBE\tbranching enzyme loss gives wrinkled seed\tPMID:1978859\t"
                "sucrose transport could do the same\tamylopectin chain length\tArabidopsis sbe mutant\tP30924\tArabidopsis thaliana",
            ]) + "\n")
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out),
                                      "--payload-source", "agent"]))
            mechanisms = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"agent"}, set(mechanisms["payload_source"]))
            self.assertEqual({"none"}, set(mechanisms["rank_contribution"]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("agent", manifest["payload_source"])


if __name__ == "__main__":
    unittest.main()
