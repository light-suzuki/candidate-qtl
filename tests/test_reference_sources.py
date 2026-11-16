"""Phenotype-rich reference records, and what they change downstream (issue #232).

The failure these tests pin down is the real one from the monitored run: a
description-only corpus answered "fasciated stem with enlarged shoot apical
meristem" with NAC proteins -- because "No Apical Meristem" is a domain name --
while the CLAVATA machinery that actually produces the phenotype scored nothing.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.reference_sources import (
    clean_comment, evidence_codes, evidence_type, expand_preset, load_refcorpus_presets,
    pubmed_ids, read_uniprot_table, uniprot_fetch,
)

# A UniProt TSV export in the exact shape reference-corpus-fetch writes.
UNIPROT = "\n".join([
    "\t".join(["Entry", "Entry Name", "Gene Names (primary)", "Gene Names", "Protein names", "Organism",
               "Organism (ID)", "Protein families", "Function [CC]", "Disruption phenotype",
               "Developmental stage", "Tissue specificity", "Pathway", "Gene Ontology IDs",
               "Gene Ontology (biological process)", "PubMed ID", "TAIR"]),
    "\t".join([
        "P30924", "CLV1_ARATH", "CLV1", "CLV1 AT1G75820", "Receptor protein kinase CLAVATA1",
        "Arabidopsis thaliana", "3702", "Protein kinase superfamily",
        "FUNCTION: Receptor kinase restricting the shoot apical meristem stem cell pool. "
        "{ECO:0000269|PubMed:10394359}.",
        "DISRUPTION PHENOTYPE: Enlarged shoot apical meristem and fasciated stems with extra floral organs. "
        "{ECO:0000269|PubMed:8402901}.",
        "", "Expressed in the shoot apical meristem.", "PATHWAY: CLAVATA-WUSCHEL meristem maintenance signalling.",
        "GO:0004672; GO:0010080", "meristem maintenance", "8402901; 10394359; 27862469", "AT1G75820",
    ]),
    "\t".join([
        "Q9SXK2", "NAC01_ARATH", "NAC001", "NAC001", "NAC domain-containing protein 1 (No Apical Meristem protein)",
        "Arabidopsis thaliana", "3702", "NAC transcription factor family",
        "FUNCTION: Transcription factor. {ECO:0000250|UniProtKB:Q9ZQ21}.",
        "", "", "", "", "GO:0003700", "regulation of transcription", "27862469", "AT1G01010",
    ]),
    # No gene symbol: kept in the ledger rather than silently dropped.
    "\t".join([
        "A0A123", "UNK_ARATH", "", "", "Uncharacterized protein", "Arabidopsis thaliana", "3702", "",
        "", "DISRUPTION PHENOTYPE: No visible phenotype. {ECO:0000269|PubMed:11111111}.",
        "", "", "", "", "", "11111111", "",
    ]),
]) + "\n"

# The same two accessions as a description-only corpus: the "before" state.
THIN_DUMP = "\n".join([
    "reference_gene_id\treference_species\ttaxon_id\tgene_symbol\tdescription\tsource_id",
    "P30924\tArabidopsis thaliana\t3702\tCLV1\tReceptor protein kinase CLAVATA1\tUniProt:P30924",
    "Q9SXK2\tArabidopsis thaliana\t3702\tNAC001\tNAC domain-containing protein 1 (No Apical Meristem protein)\tUniProt:Q9SXK2",
]) + "\n"

# Multi-taxon fixtures for the versioned preset: Medicago rows
# carry the preset's required terms, Lotus rows carry the optional ones (with
# "nin" lowercased on purpose, so the casefold match rule is exercised).
# Accessions are illustrative, not verified UniProt mappings.
UNIPROT_MEDICAGO = "\n".join([
    "\t".join(["Entry", "Entry Name", "Gene Names (primary)", "Gene Names", "Protein names", "Organism",
               "Organism (ID)", "Protein families", "Function [CC]", "Disruption phenotype",
               "Developmental stage", "Tissue specificity", "Pathway", "Gene Ontology IDs",
               "Gene Ontology (biological process)", "PubMed ID", "TAIR"]),
    "\t".join([
        "Q8L3K5", "NFP_MEDTR", "NFP", "NFP", "Nodulation factor receptor kinase NFP",
        "Medicago truncatula", "3880", "LysM domain-containing receptor-like kinase family",
        "FUNCTION: Perception of rhizobial Nod factors. {ECO:0000269|PubMed:12566564}.",
        "DISRUPTION PHENOTYPE: Nodulation abolished; no root nodule primordia. {ECO:0000269|PubMed:12566564}.",
        "", "", "PATHWAY: Nod factor signalling.", "GO:0004888", "nodulation", "12566564", "",
    ]),
    "\t".join([
        "Q6RIC4", "LYK3_MEDTR", "LYK3", "LYK3", "LysM domain receptor-like kinase 3",
        "Medicago truncatula", "3880", "LysM domain-containing receptor-like kinase family",
        "FUNCTION: Nod factor receptor signalling. {ECO:0000269|PubMed:17293529}.",
        "DISRUPTION PHENOTYPE: Reduced nodulation. {ECO:0000269|PubMed:17293529}.",
        "", "", "PATHWAY: Nod factor signalling.", "GO:0004888", "nodulation", "17293529", "",
    ]),
    "\t".join([
        "Q5S3S4", "NSP2_MEDTR", "NSP2", "NSP2", "Nodulation signaling pathway 2 protein",
        "Medicago truncatula", "3880", "GRAS transcription factor family",
        "FUNCTION: Transcription factor for nodulation signaling. {ECO:0000269|PubMed:16954775}.",
        "DISRUPTION PHENOTYPE: No nodule formation. {ECO:0000269|PubMed:16954775}.",
        "", "", "PATHWAY: Nodulation signaling.", "GO:0003700", "nodulation", "16954775", "",
    ]),
    "\t".join([
        "A0A072V3W5", "CLE13_MEDTR", "CLE13", "CLE13", "CLAVATA3/ESR-related protein 13",
        "Medicago truncatula", "3880", "CLE peptide family",
        "FUNCTION: CLE peptide signalling. {ECO:0000250}.",
        "", "", "", "", "GO:0005576", "cell communication", "", "",
    ]),
]) + "\n"

UNIPROT_LOTUS = "\n".join([
    "\t".join(["Entry", "Entry Name", "Gene Names (primary)", "Gene Names", "Protein names", "Organism",
               "Organism (ID)", "Protein families", "Function [CC]", "Disruption phenotype",
               "Developmental stage", "Tissue specificity", "Pathway", "Gene Ontology IDs",
               "Gene Ontology (biological process)", "PubMed ID", "TAIR"]),
    "\t".join([
        "Q9S7B1", "NIN_LOTJA", "nin", "nin", "Nodule inception protein",
        "Lotus japonicus", "34305", "RWP-RK transcription factor family",
        "FUNCTION: Transcription factor required for nodule organogenesis. {ECO:0000269|PubMed:10398678}.",
        "DISRUPTION PHENOTYPE: No nodules. {ECO:0000269|PubMed:10398678}.",
        "", "", "PATHWAY: Nodule organogenesis.", "GO:0003700", "nodulation", "10398678", "",
    ]),
    "\t".join([
        "Q8LKS5", "SYMRK_LOTJA", "SYMRK", "SYMRK", "Symbiosis receptor-like kinase",
        "Lotus japonicus", "34305", "Receptor-like kinase family",
        "FUNCTION: Perception of symbiotic signals. {ECO:0000269|PubMed:11923494}.",
        "DISRUPTION PHENOTYPE: Nodulation and mycorrhization abolished. {ECO:0000269|PubMed:11923494}.",
        "", "", "PATHWAY: Symbiosis signalling.", "GO:0005524", "signal transduction", "11923494", "",
    ]),
]) + "\n"

FASCIATION = "fasciated stem with enlarged shoot apical meristem"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


class CommentParsingTests(unittest.TestCase):
    def test_comment_text_loses_the_block_prefix_and_the_eco_braces(self) -> None:
        parsed = clean_comment("DISRUPTION PHENOTYPE: Fasciated stems. {ECO:0000269|PubMed:8402901}.")
        self.assertEqual("Fasciated stems", parsed)

    def test_experimental_evidence_outranks_similarity(self) -> None:
        codes = evidence_codes("{ECO:0000250|UniProtKB:Q1, ECO:0000269|PubMed:8402901}")
        self.assertEqual(["ECO:0000250", "ECO:0000269"], codes)
        self.assertEqual("experimental", evidence_type(codes))
        self.assertEqual("sequence_similarity", evidence_type(["ECO:0000250"]))
        self.assertEqual("unspecified", evidence_type([]))

    def test_only_the_pubmed_ids_the_phenotype_block_cites_are_phenotype_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(Path(tmp), "uniprot_3702.tsv", UNIPROT)
            rows, excluded, meta = read_uniprot_table(path)
            clv = next(row for row in rows if row["reference_gene_id"] == "P30924")
            # the entry cites three papers; only 8402901 supports the phenotype
            self.assertEqual("PMID:8402901", clv["phenotype_source_ids"])
            self.assertEqual("PMID:8402901;PMID:10394359;PMID:27862469", clv["literature_source_ids"])
            self.assertEqual("experimental", clv["evidence_type"])
            self.assertEqual(2, meta["rows"])
            self.assertEqual(1, meta["excluded"])
            self.assertEqual("no_gene_symbol_in_uniprot_record", excluded[0]["exclusion_reason"])
            self.assertEqual(["8402901"], pubmed_ids("{ECO:0000269|PubMed:8402901}"))


class CorpusBuildTests(unittest.TestCase):
    def _build(self, root: Path, version: str, *extra: str) -> tuple[pd.DataFrame, dict]:
        out = root / "corpus"
        self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(out),
                                  "--corpus-version", version, *extra]))
        corpus = pd.read_csv(out / version / "reference_corpus.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / version / "reference_corpus_manifest.json").read_text(encoding="utf-8"))
        return corpus, manifest

    def test_curated_phenotype_and_pathway_reach_the_corpus(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = _write(root, "uniprot_3702.tsv", UNIPROT)
            corpus, manifest = self._build(root, "rich", "--from-uniprot", str(table))
            clv = corpus.set_index("reference_gene_id").loc["P30924"]
            self.assertIn("fasciated stems", str(clv["mutant_phenotype"]).casefold())
            self.assertIn("CLAVATA-WUSCHEL", str(clv["pathway"]))
            self.assertEqual("Protein kinase superfamily", clv["protein_family"])
            self.assertEqual("UniProtKB", clv["source_database"])
            self.assertEqual(1, manifest["phenotype_row_count"])
            self.assertEqual(1, manifest["species_count"])
            self.assertEqual(0.5, manifest["field_coverage"]["mutant_phenotype"]["coverage"])
            coverage = pd.read_csv(root / "corpus" / "rich" / "reference_corpus_coverage.tsv",
                                   sep="\t", dtype=str).fillna("")
            self.assertIn("mutant_phenotype", set(coverage["field"]))
            self.assertIn("Arabidopsis thaliana", set(coverage["reference_species"]))

    def test_an_absent_phenotype_is_knowledge_unavailable_not_negative_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = _write(root, "uniprot_3702.tsv", UNIPROT)
            corpus, manifest = self._build(root, "rich", "--from-uniprot", str(table))
            nac = corpus.set_index("reference_gene_id").loc["Q9SXK2"]
            self.assertEqual("", nac["mutant_phenotype"])
            self.assertIn("knowledge_unavailable", manifest["missing_phenotype_meaning"])

    def test_the_curated_row_wins_over_the_bare_description_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus, _ = self._build(
                root, "merged",
                "--from-dump", str(_write(root, "dump.tsv", THIN_DUMP)),
                "--from-uniprot", str(_write(root, "uniprot_3702.tsv", UNIPROT)))
            self.assertEqual(2, len(corpus))
            clv = corpus.set_index("reference_gene_id").loc["P30924"]
            self.assertIn("fasciated", str(clv["mutant_phenotype"]).casefold())

    def test_an_orthogroup_map_fills_the_orthogroup_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mapping = _write(root, "og.tsv", "reference_gene_id\torthogroup_id\nP30924\tOG0000123\n")
            corpus, _ = self._build(root, "mapped",
                                    "--from-uniprot", str(_write(root, "uniprot_3702.tsv", UNIPROT)),
                                    "--orthogroup-map", str(mapping))
            corpus = corpus.set_index("reference_gene_id")
            self.assertEqual("OG0000123", corpus.loc["P30924", "orthogroup_id"])
            self.assertEqual("", corpus.loc["Q9SXK2", "orthogroup_id"])


class RetrievalImpactTests(unittest.TestCase):
    """The point of the corpus work: what the phenotype actually retrieves."""

    def _retrieve(self, root: Path, label: str, *build_args: str) -> pd.DataFrame:
        out = root / label
        self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(out / "corpus"),
                                  "--corpus-version", label, *build_args]))
        corpus_path = out / "corpus" / label / "reference_corpus.tsv"
        self.assertEqual(0, main(["retrieval-index-build", "--corpus", str(corpus_path),
                                  "--database-version", label, "--out-dir", str(out / "index")]))
        self.assertEqual(0, main(["retrieval-run", "--index", str(out / "index" / label),
                                  "--phenotype", FASCIATION, "--out-dir", str(out / "run"), "--top-k", "10"]))
        return pd.read_csv(out / "run" / "route_c_local_retrieval.tsv", sep="\t", dtype=str).fillna("")

    def test_the_description_only_corpus_ranks_the_domain_name_and_misses_clavata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            thin = self._retrieve(root, "thin", "--from-dump", str(_write(root, "dump.tsv", THIN_DUMP)))
            # "No Apical Meristem" is a domain name, and it is all the thin corpus has
            self.assertEqual("Q9SXK2", thin.iloc[0]["reference_gene_id"])

    def test_the_phenotype_rich_corpus_puts_clavata_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rich = self._retrieve(root, "rich", "--from-uniprot", str(_write(root, "uniprot_3702.tsv", UNIPROT)))
            self.assertEqual("P30924", rich.iloc[0]["reference_gene_id"])


class FetchManifestTests(unittest.TestCase):
    """Species are usually added one run at a time; the earlier provenance must survive."""

    def _fake_request(self, taxon: str):
        header = UNIPROT.splitlines()[0]
        body = UNIPROT.splitlines()[1]

        def request(url: str, timeout: int = 120):
            return "\n".join([header, body]) + "\n", {"x-uniprot-release": "2026_02",
                                                      "x-total-results": "1", "link": ""}
        return request

    def test_a_second_taxon_does_not_erase_the_first_ones_release(self) -> None:
        from candidate_qtl import reference_sources

        original = reference_sources._request
        reference_sources._request = self._fake_request("3702")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                uniprot_fetch(["3702"], str(root), allow_external_api=True)
                result = uniprot_fetch(["3880"], str(root), allow_external_api=True)
                self.assertEqual(["3702", "3880"], result["taxa"])
                manifest = json.loads((root / "reference_fetch_manifest.json").read_text(encoding="utf-8"))
                releases = {str(item["taxon_id"]): str(item["source_release"]) for item in manifest["sources"]}
                self.assertEqual({"3702": "2026_02", "3880": "2026_02"}, releases)
                # and the release reaches the corpus rows parsed from that table
                rows, _excluded, meta = read_uniprot_table(root / "uniprot_3702.tsv")
                self.assertEqual("2026_02", meta["source_release"])
                self.assertEqual("2026_02", rows[0]["source_release"])
        finally:
            reference_sources._request = original


class PresetTests(unittest.TestCase):
    """Versioned corpus presets (issue #299): expansion, validation, fetch integration."""

    def _fetch(self, root: Path, taxon_args: list[str], capture, preset: str | None = None) -> dict:
        from candidate_qtl import reference_sources

        original = reference_sources._request
        reference_sources._request = capture
        try:
            return uniprot_fetch(taxon_args, str(root), allow_external_api=True, preset_arg=preset)
        finally:
            reference_sources._request = original

    def test_expand_preset_resolves_the_configured_taxa_and_reviewed_flags(self) -> None:
        taxa = expand_preset("legume_v1")
        self.assertEqual(["3702", "3880", "34305"], [str(taxon["taxon_id"]) for taxon in taxa])
        self.assertEqual("Arabidopsis thaliana", taxa[0]["species"])
        self.assertEqual([True, False, False], [taxon["reviewed"] for taxon in taxa])

    def test_an_unknown_preset_names_the_available_ids(self) -> None:
        with self.assertRaises(ValueError) as error:
            expand_preset("nope")
        self.assertIn("unknown reference corpus preset", str(error.exception))
        self.assertIn("legume_v1", str(error.exception))

    def test_preset_loading_validates_the_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = _write(root, "bad.json", json.dumps({
                "schema_version": "1.0",
                "presets": [{"preset_id": "dup", "provider": "x", "taxa": [{"taxon_id": 1, "species": "s", "reviewed": True}]},
                            {"preset_id": "dup", "provider": "x", "taxa": [{"taxon_id": 2, "species": "s", "reviewed": True}]}],
            }))
            with self.assertRaises(ValueError) as error:
                load_refcorpus_presets(str(bad))
            self.assertIn("duplicate preset_id 'dup'", str(error.exception))
            broken = _write(root, "broken.json", json.dumps({
                "schema_version": "1.0",
                "presets": [{"preset_id": "p", "provider": "x",
                             "taxa": [{"taxon_id": "not-a-number", "species": "s", "reviewed": True}]}],
            }))
            with self.assertRaises(ValueError) as error:
                load_refcorpus_presets(str(broken))
            self.assertIn("numeric 'taxon_id'", str(error.exception))

    def test_preset_alone_fetches_all_taxa_with_their_reviewed_flags(self) -> None:
        class Capture:
            def __init__(self) -> None:
                self.urls: list[str] = []

            def __call__(self, url: str, timeout: int = 120):
                self.urls.append(url)
                return "\n".join([UNIPROT.splitlines()[0], UNIPROT.splitlines()[1]]) + "\n", {
                    "x-uniprot-release": "2026_02", "x-total-results": "1", "link": ""}

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            capture = Capture()
            result = self._fetch(root, [], capture=capture, preset="legume_v1")
            self.assertEqual(["34305", "3702", "3880"], result["taxa"])
            manifest = json.loads((root / "reference_fetch_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("legume_v1", manifest["preset"])
            by_taxon = {str(item["taxon_id"]): item for item in manifest["sources"]}
            # per-taxon reviewed flags reach the manifest and the query
            self.assertTrue(by_taxon["3702"]["reviewed"])
            self.assertFalse(by_taxon["3880"]["reviewed"])
            self.assertFalse(by_taxon["34305"]["reviewed"])
            self.assertEqual("legume_v1", by_taxon["34305"]["preset_id"])
            self.assertEqual("organism_id:3702 AND reviewed:true", by_taxon["3702"]["query"])
            self.assertEqual("organism_id:3880", by_taxon["3880"]["query"])
            self.assertTrue(by_taxon["3702"]["reviewed_only"])
            self.assertFalse(by_taxon["3880"]["reviewed_only"])
            urls_3702 = [url for url in capture.urls if "organism_id%3A3702" in url]
            self.assertTrue(urls_3702 and "reviewed%3Atrue" in urls_3702[0])
            urls_3880 = [url for url in capture.urls if "organism_id%3A3880" in url]
            self.assertTrue(urls_3880 and "reviewed%3Atrue" not in urls_3880[0])

    def test_explicit_taxon_unions_with_the_preset_and_wins_on_overlap(self) -> None:
        from candidate_qtl import reference_sources

        def fake_request(url: str, timeout: int = 120):
            return "\n".join([UNIPROT.splitlines()[0], UNIPROT.splitlines()[1]]) + "\n", {
                "x-uniprot-release": "2026_02", "x-total-results": "1", "link": ""}

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = reference_sources._request
            reference_sources._request = fake_request
            try:
                result = uniprot_fetch(["3880"], str(root), allow_external_api=True, preset_arg="legume_v1")
            finally:
                reference_sources._request = original
            self.assertEqual(["34305", "3702", "3880"], result["taxa"])
            manifest = json.loads((root / "reference_fetch_manifest.json").read_text(encoding="utf-8"))
            by_taxon = {str(item["taxon_id"]): item for item in manifest["sources"]}
            # the explicit --taxon 3880 wins completely: no preset provenance and
            # the global reviewed-only semantics; the preset fills in the rest
            self.assertIsNone(by_taxon["3880"]["preset_id"])
            self.assertTrue(by_taxon["3880"]["reviewed"])
            self.assertTrue(by_taxon["3880"]["reviewed_only"])
            self.assertEqual("organism_id:3880 AND reviewed:true", by_taxon["3880"]["query"])
            self.assertEqual("organism_id:34305", by_taxon["34305"]["query"])
            self.assertEqual("legume_v1", by_taxon["34305"]["preset_id"])

    def test_fetch_still_requires_a_taxon_or_a_preset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError) as error:
                uniprot_fetch([], tmp, allow_external_api=True)
            self.assertIn("--taxon or --preset", str(error.exception))

    def test_preset_fetch_manifest_merge_preserves_earlier_releases(self) -> None:
        from candidate_qtl import reference_sources

        def fake_request(url: str, timeout: int = 120):
            return "\n".join([UNIPROT.splitlines()[0], UNIPROT.splitlines()[1]]) + "\n", {
                "x-uniprot-release": "2026_02", "x-total-results": "1", "link": ""}

        original = reference_sources._request
        reference_sources._request = fake_request
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                uniprot_fetch(["3702"], str(root), allow_external_api=True)
                result = uniprot_fetch([], str(root), allow_external_api=True, preset_arg="legume_v1")
                self.assertEqual(["34305", "3702", "3880"], result["taxa"])
                manifest = json.loads((root / "reference_fetch_manifest.json").read_text(encoding="utf-8"))
                releases = {str(item["taxon_id"]): str(item["source_release"]) for item in manifest["sources"]}
                self.assertEqual({"3702": "2026_02", "3880": "2026_02", "34305": "2026_02"}, releases)
        finally:
            reference_sources._request = original


class PresetAcceptanceTests(unittest.TestCase):
    """The fail-closed vocabulary gate on reference-corpus-build (issue #299)."""

    def _build(self, root: Path, *extra: str) -> int:
        return main(["reference-corpus-build", "--out-dir", str(root / "corpus"),
                     "--corpus-version", "legume_test", *extra])

    def _manifest(self, root: Path) -> dict:
        return json.loads((root / "corpus" / "legume_test" / "reference_corpus_manifest.json")
                          .read_text(encoding="utf-8"))

    def test_all_required_terms_present_passes_and_records_every_term(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = [("uniprot_3702.tsv", UNIPROT),
                      ("uniprot_3880.tsv", UNIPROT_MEDICAGO),
                      ("uniprot_34305.tsv", UNIPROT_LOTUS)]
            extra = []
            for name, content in tables:
                extra += ["--from-uniprot", str(_write(root, name, content))]
            self.assertEqual(0, self._build(root, *extra, "--preset-acceptance", "legume_v1"))
            acceptance = self._manifest(root)["preset_acceptance"]
            self.assertEqual("legume_v1", acceptance["preset_id"])
            self.assertTrue(acceptance["fail_closed_on_required"])
            for term in ("NFP", "LYK3", "NSP2", "CLE13"):
                record = acceptance["required_terms"][term]
                self.assertTrue(record["found"], term)
                self.assertGreaterEqual(record["count"], 1)
            # optional terms are recorded, not gated: "nin" matches casefolded,
            # CLE-RS1/LHK1 are absent from the fixture tables and that is fine
            self.assertTrue(acceptance["optional_terms"]["NIN"]["found"])
            self.assertEqual(1, acceptance["optional_terms"]["NIN"]["count"])
            self.assertTrue(acceptance["optional_terms"]["SYMRK"]["found"])
            self.assertFalse(acceptance["optional_terms"]["CLE-RS1"]["found"])
            self.assertFalse(acceptance["optional_terms"]["LHK1"]["found"])

    def test_a_missing_required_term_fails_closed_before_anything_is_written(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError) as error:
                self._build(root, "--from-uniprot", str(_write(root, "uniprot_3702.tsv", UNIPROT)),
                            "--preset-acceptance", "legume_v1")
            message = str(error.exception)
            self.assertIn("preset acceptance failed", message)
            self.assertIn("NFP", message)
            self.assertIn("LYK3", message)
            # the previous frozen corpus cannot be replaced by a failed build
            self.assertFalse((root / "corpus" / "legume_test").exists())

    def test_acceptance_matches_usable_symbols_casefolded_and_counts(self) -> None:
        from candidate_qtl.reference_corpus import preset_acceptance

        frame = pd.DataFrame({"gene_symbol": ["nfp", "NFP", "NAC001"]})
        result = preset_acceptance(frame, {"preset_id": "x", "required_terms": ["NFP"], "optional_terms": []})
        self.assertEqual({"found": True, "count": 2}, result["required_terms"]["NFP"])
        # an accession masquerading as a gene symbol is never a symbol
        accessioned = preset_acceptance(pd.DataFrame({"gene_symbol": ["P30924"]}),
                                        {"preset_id": "x", "required_terms": ["P30924"], "optional_terms": []})
        self.assertEqual({"found": False, "count": 0}, accessioned["required_terms"]["P30924"])


class FetchGateTests(unittest.TestCase):
    def test_no_network_call_without_an_explicit_allow_flag(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError) as error:
                uniprot_fetch(["3702"], tmp, allow_external_api=False)
            self.assertIn("--allow-external-api", str(error.exception))

    def test_the_fetch_command_has_no_phenotype_or_query_argument(self) -> None:
        from candidate_qtl.cli import build_parser
        parser = build_parser()
        action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
        options: set[str] = set()
        for sub_action in action.choices["reference-corpus-fetch"]._actions:
            options.update(sub_action.option_strings)
        for forbidden in ("--phenotype", "--query", "--interval", "--candidates", "--case-id"):
            self.assertNotIn(forbidden, options)


class CorpusCompareTests(unittest.TestCase):
    def _corpora(self, root: Path) -> tuple[Path, Path]:
        self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(root / "thin"),
                                  "--corpus-version", "thin",
                                  "--from-dump", str(_write(root, "dump.tsv", THIN_DUMP))]))
        self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(root / "rich"),
                                  "--corpus-version", "rich",
                                  "--from-uniprot", str(_write(root, "uniprot_3702.tsv", UNIPROT))]))
        return (root / "thin" / "thin" / "reference_corpus.tsv",
                root / "rich" / "rich" / "reference_corpus.tsv")

    def test_coverage_delta_shows_where_the_enriched_corpus_adds_biology(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            thin, rich = self._corpora(root)
            out = root / "compare"
            self.assertEqual(0, main(["reference-corpus-compare", "--corpus", f"thin={thin}",
                                      "--corpus", f"rich={rich}", "--baseline-label", "thin",
                                      "--out-dir", str(out)]))
            coverage = pd.read_csv(out / "corpus_coverage_comparison.tsv", sep="\t", dtype=str).fillna("")
            row = coverage[(coverage["corpus_label"] == "rich") & (coverage["scope"] == "all")
                           & (coverage["field"] == "mutant_phenotype")].iloc[0]
            self.assertEqual("0.0", str(row["baseline_coverage"]))
            self.assertEqual(0.5, float(row["coverage_delta"]))
            summary = pd.read_csv(out / "corpus_summary.tsv", sep="\t", dtype=str).fillna("").set_index("corpus_label")
            self.assertEqual("0", summary.loc["thin", "phenotype_rows"])
            self.assertEqual("1", summary.loc["rich", "phenotype_rows"])
            self.assertEqual("2", summary.loc["rich", "shared_accessions_with_baseline"])

    def test_retrieval_recall_is_compared_per_corpus_with_its_own_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            thin, rich = self._corpora(root)
            paths = {}
            for label, corpus in (("thin", thin), ("rich", rich)):
                index_out = root / f"index_{label}"
                self.assertEqual(0, main(["retrieval-index-build", "--corpus", str(corpus),
                                          "--database-version", label, "--out-dir", str(index_out)]))
                run_out = root / f"run_{label}"
                self.assertEqual(0, main(["retrieval-run", "--index", str(index_out / label),
                                          "--phenotype", FASCIATION, "--out-dir", str(run_out), "--top-k", "10"]))
                paths[label] = run_out / "route_c_local_retrieval.tsv"
            manifest = _write(root, "retrieval_manifest.tsv", "\n".join([
                "corpus_label\tcase_id\tretrieval_tsv\ttruth_family_id\ttruth_gene_id",
                f"thin\tfasciation\t{paths['thin']}\tCLV\tP30924",
                f"rich\tfasciation\t{paths['rich']}\tCLV\tP30924",
            ]) + "\n")
            out = root / "compare"
            self.assertEqual(0, main(["reference-corpus-compare", "--corpus", f"thin={thin}",
                                      "--corpus", f"rich={rich}", "--baseline-label", "thin",
                                      "--retrieval-manifest", str(manifest), "--out-dir", str(out)]))
            summary = pd.read_csv(out / "corpus_retrieval_summary.tsv", sep="\t", dtype=str).fillna("")
            row = summary[(summary["metric"] == "family_recall_5") & (summary["corpus_label"] == "rich")].iloc[0]
            self.assertEqual("1", row["evaluated"])
            self.assertEqual("1", row["hits"])
            self.assertEqual(1.0, float(row["recall"]))
            gene = summary[(summary["metric"] == "gene_recall_20") & (summary["corpus_label"] == "rich")].iloc[0]
            self.assertEqual("1", gene["evaluated"])
            manifest_json = json.loads((out / "reference_corpus_comparison_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("evaluated", manifest_json["retrieval_status"])

    def test_two_corpora_are_required_and_the_baseline_must_exist(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            thin, rich = self._corpora(root)
            from candidate_qtl.reference_corpus_compare import reference_corpus_compare
            with self.assertRaises(ValueError):
                reference_corpus_compare([f"thin={thin}"], str(root / "x"))
            with self.assertRaises(ValueError):
                reference_corpus_compare([f"thin={thin}", f"rich={rich}"], str(root / "x"), baseline_label="nope")


if __name__ == "__main__":
    unittest.main()
