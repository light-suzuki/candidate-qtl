"""The reference resource Routes B and C share: retrieve, then verify proposals."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.reference_corpus import _parse_defline

DUMP = "\n".join([
    "reference_gene_id\treference_species\ttaxon_id\tgene_symbol\tdescription\tsource_id",
    "P30924\tArabidopsis thaliana\t3702\tCLV1\tCLAVATA1 receptor kinase\tUniProt:P30924",
    "Q9SXK1\tArabidopsis thaliana\t3702\tSUC2\tsucrose transporter 2\tUniProt:Q9SXK1",
    "Q9SWW6\tArabidopsis thaliana\t3702\tCESA9\tcellulose synthase A9\tUniProt:Q9SWW6",
]) + "\n"

# Two of these are legitimate additions; two must not survive.
PAYLOAD = [
    {"reference_gene_id": "P30924", "gene_symbol": "CLV1",
     "mutant_phenotype": "fasciated stem and enlarged meristem",
     "pathway": "meristem_maintenance", "source_id": "PMID:8402901"},
    {"reference_gene_id": "Q9SXK1", "gene_symbol": "SUC2",
     "mutant_phenotype": "reduced seed filling", "source_id": "PMID:8402901"},
    {"reference_gene_id": "P99999", "gene_symbol": "FAKE1",
     "mutant_phenotype": "invented phenotype", "source_id": "PMID:8402901"},
    {"reference_gene_id": "Q9SWW6", "gene_symbol": "CESA9",
     "mutant_phenotype": "collapsed xylem", "source_id": "PMID:7777777"},
]


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _citation_cache(root: Path) -> Path:
    cache = root / "pmc"
    cache.mkdir()
    (cache / "med_8402901.json").write_text(json.dumps({
        "pmid": "8402901", "id": "8402901", "title": "CLAVATA1 controls meristem size",
    }), encoding="utf-8")
    return cache


class DeflineTests(unittest.TestCase):
    def test_uniprot_defline_is_parsed_into_fields(self) -> None:
        parsed = _parse_defline("E4 SUMO-protein ligase PIAL1 OS=Arabidopsis thaliana OX=3702 GN=PIAL1 PE=1 SV=1")
        self.assertEqual("E4 SUMO-protein ligase PIAL1", parsed["description"])
        self.assertEqual("Arabidopsis thaliana", parsed["species"])
        self.assertEqual("3702", parsed["taxon"])
        self.assertEqual("PIAL1", parsed["symbol"])

    def test_defline_without_gene_symbol_still_parses(self) -> None:
        parsed = _parse_defline("Uncharacterized protein OS=Arabidopsis thaliana OX=3702 PE=4 SV=1")
        self.assertEqual("", parsed["symbol"])
        self.assertEqual("Arabidopsis thaliana", parsed["species"])


class CorpusBuildTests(unittest.TestCase):
    def _build(self, root: Path, *extra: str) -> tuple[pd.DataFrame, dict]:
        out = root / "corpus"
        payload = _write(root, "payload.json", json.dumps({"records": PAYLOAD}))
        self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(out),
                                  "--corpus-version", "plant_ref_test",
                                  "--from-dump", str(_write(root, "dump.tsv", DUMP)),
                                  "--from-payload", str(payload),
                                  "--citation-cache", str(_citation_cache(root)), *extra]))
        corpus = pd.read_csv(out / "plant_ref_test" / "reference_corpus.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "plant_ref_test" / "reference_corpus_manifest.json").read_text(encoding="utf-8"))
        return corpus, manifest

    def test_a_fabricated_accession_never_enters_the_corpus(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus, manifest = self._build(root)
            self.assertNotIn("P99999", set(corpus["reference_gene_id"]))
            excluded = pd.read_csv(root / "corpus" / "plant_ref_test" / "reference_corpus_excluded.tsv",
                                   sep="\t", dtype=str).fillna("").set_index("reference_gene_id")
            self.assertEqual("accession_not_in_any_retrieved_source", excluded.loc["P99999", "exclusion_reason"])
            self.assertEqual("llm_unverified", excluded.loc["P99999", "provenance_tier"])

    def test_an_unresolvable_citation_blocks_promotion_but_keeps_the_retrieved_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus, _ = self._build(root)
            cesa = corpus.set_index("reference_gene_id").loc["Q9SWW6"]
            # the proposal cited a PMID that does not resolve locally, so the added
            # phenotype text is dropped -- but the retrieved skeleton row survives
            self.assertEqual("retrieved", cesa["provenance_tier"])
            self.assertEqual("", cesa["mutant_phenotype"])
            excluded = pd.read_csv(root / "corpus" / "plant_ref_test" / "reference_corpus_excluded.tsv",
                                   sep="\t", dtype=str).fillna("").set_index("reference_gene_id")
            self.assertEqual("no_citation_resolved_against_the_local_cache",
                             excluded.loc["Q9SWW6", "exclusion_reason"])

    def test_a_verified_proposal_is_promoted_and_marked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus, manifest = self._build(root)
            clv = corpus.set_index("reference_gene_id").loc["P30924"]
            self.assertEqual("llm_verified", clv["provenance_tier"])
            self.assertEqual("fasciated stem and enlarged meristem", clv["mutant_phenotype"])
            self.assertEqual("PMID:8402901", clv["verified_source_ids"])
            # species is taken from the authoritative record, not from the proposal
            self.assertEqual("Arabidopsis thaliana", clv["reference_species"])
            self.assertEqual("CLV", clv["family_id"])
            self.assertTrue(manifest["contains_model_proposed_rows"])
            self.assertEqual(2, manifest["model_proposed_rows"])
            self.assertEqual(3, manifest["row_count"])

    def test_the_corpus_is_case_independent_by_construction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _corpus, manifest = self._build(root)
            self.assertTrue(manifest["case_independent"])
            self.assertFalse(manifest["phenotype_used"])
            self.assertEqual([], manifest["interval_arguments"])
            self.assertTrue(manifest["corpus_sha256"])
        parser = build_parser()
        action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
        options = set()
        for sub_action in action.choices["reference-corpus-build"]._actions:
            options.update(sub_action.option_strings)
        for forbidden in ("--phenotype", "--interval", "--candidates", "--case-id"):
            self.assertNotIn(forbidden, options)

    def test_benchmark_mode_bars_focal_lineage_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "corpus"
            payload = _write(root, "payload.json", json.dumps({"records": [
                {"reference_gene_id": "P30924", "gene_symbol": "CLV1",
                 "mutant_phenotype": "matches the Pisum sativum fasciata locus",
                 "source_id": "PMID:8402901"},
            ]}))
            self.assertEqual(0, main(["reference-corpus-build", "--out-dir", str(out),
                                      "--corpus-version", "v1", "--from-dump", str(_write(root, "dump.tsv", DUMP)),
                                      "--from-payload", str(payload),
                                      "--citation-cache", str(_citation_cache(root)), "--benchmark-mode"]))
            excluded = pd.read_csv(out / "v1" / "reference_corpus_excluded.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("focal_lineage_or_answer_shaped_marker_in_benchmark_mode",
                          set(excluded["exclusion_reason"]))

    def test_the_corpus_feeds_both_routes_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus, _ = self._build(root)
            corpus_path = root / "corpus" / "plant_ref_test" / "reference_corpus.tsv"
            # Route C: index build accepts it as-is
            index_out = root / "index"
            self.assertEqual(0, main(["retrieval-index-build", "--corpus", str(corpus_path),
                                      "--database-version", "plant_ref_test", "--out-dir", str(index_out)]))
            index = json.loads((index_out / "plant_ref_test" / "index_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, index["document_count"])
            # Route B: the same file is a reference-function resource
            candidates = _write(root, "candidates.tsv",
                                "gene_id\tchrom\tstart\tend\tdescription\nPsat_G1\tchr5\t1\t100\treceptor kinase\n")
            orthology = _write(root, "orthology.tsv",
                               "gene_id\treference_gene_id\torthogroup\trelationship\tprovider\n"
                               "Psat_G1\tP30924\tOG1\tbest_hit_homology\ttest\n")
            route_b = root / "route_b"
            self.assertEqual(0, main(["route-b-prepare", "--candidates", str(candidates),
                                      "--orthology", str(orthology), "--reference-function", str(corpus_path),
                                      "--out-dir", str(route_b), "--no-auto-local-db"]))
            manifest = json.loads((route_b / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("source_backed_reverse_lookup", manifest["route_b_readiness"])


class PresetGateTests(unittest.TestCase):
    """reference-corpus-build --preset-acceptance against the packaged presets (issue #299)."""

    def _build(self, root: Path, *extra: str) -> int:
        return main(["reference-corpus-build", "--out-dir", str(root / "corpus"),
                     "--corpus-version", "v1",
                     "--from-dump", str(_write(root, "dump.tsv", DUMP)), *extra])

    def test_an_unknown_preset_id_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError) as error:
                self._build(root, "--preset-acceptance", "nope")
            self.assertIn("unknown reference corpus preset", str(error.exception))
            self.assertIn("legume_v1", str(error.exception))

    def test_a_corpus_without_the_preset_vocabulary_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError) as error:
                self._build(root, "--preset-acceptance", "legume_v1")
            message = str(error.exception)
            self.assertIn("preset acceptance failed", message)
            for term in ("NFP", "LYK3", "NSP2", "CLE13"):
                self.assertIn(term, message)
            self.assertFalse((root / "corpus" / "v1").exists())


if __name__ == "__main__":
    unittest.main()
