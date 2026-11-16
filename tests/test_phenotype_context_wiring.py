"""Issue #361: phenotype_context.yaml wired into the real Route A/C flows.

The context file is an audit input, never a query input: the free text is
always preserved, no structured field is inserted into any prompt or retrieval
query, and rank/score/weight outputs stay byte-identical.  The manifest block
records the honest ``query_usage`` (``stored_only`` in the current wiring) and
the deterministic context fingerprint.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import yaml

from candidate_qtl.cli import main

CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_CLV1\tchr5\t1000\t3000\tCLAVATA1 like receptor kinase acting in the meristem",
    "Psat_UNK\tchr5\t9000\t11000\thypothetical protein",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\torthogroup\trelationship\tprovider",
    "Psat_CLV1\tAtCLV1\tOG0002\tbest_hit_homology\tlocal_blastp:test",
]) + "\n"

REFERENCE_CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase restricting the stem cell pool\t"
    "fasciated stem and enlarged meristem\tmeristem_maintenance\tPMID:8402901",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed\tstarch_biosynthesis\tPMID:1978859",
]) + "\n"

PAYLOAD = "\n".join([
    "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
    "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
    "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
    "meristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
    "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
    "P30924\tArabidopsis thaliana\tCLV1\treviewed receptor kinase with a fasciation phenotype",
    "hormone transport\tPIN\tauxin redistribution changes stem geometry\tPMID:1978859\t"
    "meristem lesions give the same shape\tassay auxin flux\tArabidopsis pin mutant\tfused organs\t"
    "Q9SXK1\tArabidopsis thaliana\tPIN1\tcanonical auxin efflux carrier",
    "cell wall biosynthesis\tCESA\tweakened wall changes stem shape\tPMID:12172025\t"
    "wall phenotypes rarely fasciate\tsection the stem\tArabidopsis irx mutant\tcollapsed xylem\t"
    "Q9SWW6\tArabidopsis thaliana\tCESA9\tsecondary wall cellulose synthase",
]) + "\n"

PHENOTYPE_TEXT = "fasciated stem with enlarged meristem"

CONTEXT_STRUCTURED = (
    "case_id: case_fas\n"
    "phenotype_context:\n"
    "  trait_ontology_id: TO:0002616\n"
    "  trait_ontology_term: shoot development trait\n"
    "  developmental_stage: seedling\n"
    "  tissue: shoot apical meristem\n"
    "  environment_or_treatment: long day\n"
)

CONTEXT_SEEDLING = (
    "phenotype_context:\n"
    "  developmental_stage: seedling\n"
    "  environment_or_treatment: control\n"
)

CONTEXT_ADULT = (
    "phenotype_context:\n"
    "  developmental_stage: adult\n"
    "  environment_or_treatment: drought\n"
)

CONTEXT_UNMAPPED = (
    "phenotype_context:\n"
    "  unmapped_term: 頂端花束\n"
    "  developmental_stage: adult reproductive\n"
)


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _record_cache(root: Path) -> Path:
    cache = root / "records"
    cache.mkdir()
    for accession in ("P30924", "Q9SXK1", "Q9SWW6"):
        (cache / f"{accession}.json").write_text(json.dumps({
            "primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Arabidopsis thaliana"},
            "sequence": {"value": "MASTPRQTVSEK" * 8},
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "8402901"}]}}],
        }), encoding="utf-8")
    return cache


def _source_cache(root: Path) -> Path:
    cache = root / "pmc"
    cache.mkdir()
    for pmid, title in (("8402901", "CLAVATA1 receptor kinase controls meristem size in Arabidopsis"),
                        ("1978859", "Starch branching enzyme SBE and the wrinkled seed character"),
                        ("12172025", "CESA cellulose synthase and secondary wall formation")):
        (cache / f"med_{pmid}.json").write_text(json.dumps({
            "cache_key": f"pmid_{pmid}", "id": pmid, "pmid": pmid, "title": title,
            "abstractText": title, "journalTitle": "Plant Cell", "pubYear": "1993",
            "doi": f"10.1/{pmid}", "authorString": "A B, C D",
        }), encoding="utf-8")
    return cache


class PhenotypeContextWiringTests(unittest.TestCase):
    def _build_index(self, root: Path) -> Path:
        index_out = root / "index"
        self.assertEqual(0, main(["retrieval-index-build",
                                  "--corpus", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                                  "--database-version", "plant_ref_test", "--out-dir", str(index_out)]))
        return index_out / "plant_ref_test" / "index_manifest.json"

    def _case_and_resources(self, root: Path) -> list[str]:
        case = root / "case.yaml"
        case.write_text(yaml.safe_dump({
            "case_id": "case_fas",
            "trait": {"name": "fasciated stem"},
            "source": {"source_id": "fixture", "target_interval": "chr5:1-20000"},
            "target": {"assembly": "fixture"},
            "inputs": {"candidates_resource": "candidates"},
            "ranking_policy": "deterministic_v1",
        }, sort_keys=False), encoding="utf-8")
        resources = root / "resources.yaml"
        resources.write_text(yaml.safe_dump({
            "resources": {"candidates": {"type": "candidates_tsv", "path": "candidates.tsv"}}
        }, sort_keys=False), encoding="utf-8")
        return ["--case", str(case), "--resources", str(resources)]

    def _profile_args(self, root: Path, *, out: Path, bridge_case_run: bool = False,
                      context_path: Path | None = None) -> list[str]:
        args = [
            "guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
            "--case-id", "case_fas", "--focal-species", "Pisum sativum",
            "--phenotype", PHENOTYPE_TEXT,
            "--mechanism-payload", str(_write(root, "payload.tsv", PAYLOAD)), "--payload-source", "agent",
            "--record-cache", str(_record_cache(root)), "--source-cache", str(_source_cache(root)),
            "--retrieval-index", str(self._build_index(root)),
            "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
            "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
            "--reference-function", str(root / "corpus.tsv"),
        ]
        if context_path is not None:
            args += ["--phenotype-context", str(context_path)]
        if bridge_case_run:
            args += self._case_and_resources(root)
        return args

    def _run_profile(self, root: Path, *, bridge_case_run: bool = False,
                     context_path: Path | None = None) -> Path:
        out = root / "triple"
        self.assertEqual(0, main(self._profile_args(
            root, out=out, bridge_case_run=bridge_case_run, context_path=context_path)))
        return out

    def _status_block(self, path: Path) -> dict:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data["structured_context"]

    def test_legacy_free_text_profile_records_missing_context_and_preserves_rank(self) -> None:
        """No context file: every Route A/C manifest says structured_context_missing,
        query_usage stays honest (stored_only), and ranked_candidates.tsv is
        byte-identical to a plain case-run baseline."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root, bridge_case_run=True)

            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("structured_context_missing", profile["structured_context"]["status"])
            self.assertEqual("raw_text", profile["structured_context"]["term_usage"])
            self.assertEqual("stored_only", profile["structured_context"]["query_usage"]["usage"])

            for status_path in (
                out / "a01_phenotype_investigation" / "phenotype_investigation_manifest.json",
                out / "a06_route_a_status" / "phenotype_first_status.json",
                out / "c01_retrieval_run" / "route_c_local_retrieval_status.json",
                out / "c02_retrieval_freeze" / "route_c_retrieval_status.json",
            ):
                block = self._status_block(status_path)
                self.assertEqual("structured_context_missing", block["status"], status_path)
                self.assertEqual("raw_text", block["term_usage"], status_path)
                self.assertEqual("none", block["present_fields"], status_path)
                self.assertEqual("none", block["rank_contribution"], status_path)
                self.assertEqual("stored_only", block["query_usage"]["usage"], status_path)
                self.assertEqual(64, len(block["context_fingerprint"]), status_path)
                self.assertEqual(PHENOTYPE_TEXT, block["phenotype_text"], status_path)

            baseline = root / "baseline_case"
            self.assertEqual(0, main(["case-run", "--case", str(root / "case.yaml"),
                                      "--resources", str(root / "resources.yaml"),
                                      "--out-dir", str(baseline)]))
            self.assertEqual(
                (baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                (out / "z03_case_run" / "rank" / "ranked_candidates.tsv").read_bytes(),
            )

    def test_context_file_wires_into_route_a_and_c_manifests_with_honest_usage(self) -> None:
        """A context file with stage/tissue/environment + ontology term is recorded
        in the Route A/C manifests (used status, normalized term, fingerprint,
        present fields), query_usage stays stored_only, and rank bytes are unchanged."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            context_path = _write(root, "phenotype_context.yaml", CONTEXT_STRUCTURED)
            out = self._run_profile(root, bridge_case_run=True, context_path=context_path)

            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            profile_block = profile["structured_context"]
            self.assertEqual("structured_context_used", profile_block["status"])
            self.assertEqual("normalized_term", profile_block["term_usage"])
            self.assertIn("developmental_stage", profile_block["present_fields"])
            self.assertIn("environment_or_treatment", profile_block["present_fields"])
            self.assertEqual("stored_only", profile_block["query_usage"]["usage"])
            self.assertEqual("none", profile_block["rank_contribution"])

            fingerprints: set[str] = set()
            for status_path in (
                out / "a01_phenotype_investigation" / "phenotype_investigation_manifest.json",
                out / "a06_route_a_status" / "phenotype_first_status.json",
                out / "c01_retrieval_run" / "route_c_local_retrieval_status.json",
                out / "c02_retrieval_freeze" / "route_c_retrieval_status.json",
            ):
                block = self._status_block(status_path)
                self.assertEqual("structured_context_used", block["status"], status_path)
                self.assertEqual("normalized_term", block["term_usage"], status_path)
                self.assertIn("trait_ontology_id", block["present_fields"], status_path)
                self.assertIn("tissue", block["present_fields"], status_path)
                self.assertIn("environment_or_treatment", block["present_fields"], status_path)
                self.assertEqual("stored_only", block["query_usage"]["usage"], status_path)
                self.assertEqual("none", block["rank_contribution"], status_path)
                self.assertEqual(PHENOTYPE_TEXT, block["phenotype_text"], status_path)
                fingerprints.add(block["context_fingerprint"])
            self.assertEqual(1, len(fingerprints), "one context, one fingerprint across stages")

            c02 = json.loads((out / "c02_retrieval_freeze" / "route_c_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertTrue(c02["phenotype_used"], "the free-text flag stays truthful")
            self.assertEqual("none", c02["rank_contribution"])

            baseline = root / "baseline_case"
            self.assertEqual(0, main(["case-run", "--case", str(root / "case.yaml"),
                                      "--resources", str(root / "resources.yaml"),
                                      "--out-dir", str(baseline)]))
            self.assertEqual(
                (baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                (out / "z03_case_run" / "rank" / "ranked_candidates.tsv").read_bytes(),
            )

    def test_same_free_text_different_stage_environment_distinct_fingerprints(self) -> None:
        """Same --phenotype, different structured context: distinct fingerprints in
        the Route C status, and the retrieval output itself is unchanged."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._build_index(root)
            first = root / "run_seedling"
            second = root / "run_adult"
            for out, content in ((first, CONTEXT_SEEDLING), (second, CONTEXT_ADULT)):
                context_path = _write(root, f"context_{out.name}.yaml", content)
                self.assertEqual(0, main(["retrieval-run", "--index", str(index),
                                          "--phenotype", PHENOTYPE_TEXT, "--case-id", "case_fas",
                                          "--out-dir", str(out), "--phenotype-context", str(context_path)]))
            first_block = self._status_block(first / "route_c_local_retrieval_status.json")
            second_block = self._status_block(second / "route_c_local_retrieval_status.json")
            self.assertNotEqual(first_block["context_fingerprint"], second_block["context_fingerprint"])
            self.assertEqual("structured_context_used", first_block["status"])
            self.assertEqual("structured_context_used", second_block["status"])
            self.assertEqual("raw_text", first_block["term_usage"])
            self.assertEqual(
                (first / "route_c_local_retrieval.tsv").read_bytes(),
                (second / "route_c_local_retrieval.tsv").read_bytes(),
                "the retrieval query is the free text; context must not change it",
            )

    def test_unmapped_term_fixture_preserves_original_text(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._build_index(root)
            context_path = _write(root, "phenotype_context.yaml", CONTEXT_UNMAPPED)
            out = root / "run_unmapped"
            self.assertEqual(0, main(["retrieval-run", "--index", str(index),
                                      "--phenotype", PHENOTYPE_TEXT, "--case-id", "case_fas",
                                      "--out-dir", str(out), "--phenotype-context", str(context_path)]))
            block = self._status_block(out / "route_c_local_retrieval_status.json")
            self.assertEqual("unmapped_term", block["term_usage"])
            self.assertEqual("頂端花束", block["unmapped_term"])
            self.assertEqual(PHENOTYPE_TEXT, block["phenotype_text"])
            self.assertEqual("structured_context_used", block["status"])
            self.assertIn("unmapped_term", block["present_fields"])

    def test_conflicting_context_case_id_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            context_path = _write(root, "phenotype_context.yaml",
                                  "case_id: another_case\nphenotype_context:\n  tissue: leaf\n")
            with self.assertRaisesRegex(ValueError, "does not match"):
                main(["guided-profile", "--profile", "triple-route", "--out-dir", str(root / "plan"),
                      "--case-id", "case_fas", "--phenotype", PHENOTYPE_TEXT,
                      "--phenotype-context", str(context_path)])
            index = self._build_index(root)
            with self.assertRaisesRegex(ValueError, "does not match"):
                main(["retrieval-run", "--index", str(index), "--phenotype", PHENOTYPE_TEXT,
                      "--case-id", "case_fas", "--out-dir", str(root / "c1"),
                      "--phenotype-context", str(context_path)])

    def test_missing_context_file_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(ValueError, "not found"):
                main(["guided-profile", "--profile", "triple-route", "--out-dir", str(root / "plan"),
                      "--case-id", "case_fas", "--phenotype", PHENOTYPE_TEXT,
                      "--phenotype-context", str(root / "missing.yaml")])

    def test_deterministic_rerun_records_the_same_fingerprint_and_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._build_index(root)
            context_path = _write(root, "phenotype_context.yaml", CONTEXT_STRUCTURED)
            first = root / "run_a"
            second = root / "run_b"
            for out in (first, second):
                self.assertEqual(0, main(["retrieval-run", "--index", str(index),
                                          "--phenotype", PHENOTYPE_TEXT, "--case-id", "case_fas",
                                          "--out-dir", str(out), "--phenotype-context", str(context_path)]))
            first_block = self._status_block(first / "route_c_local_retrieval_status.json")
            second_block = self._status_block(second / "route_c_local_retrieval_status.json")
            self.assertEqual(first_block, second_block)
            self.assertEqual("normalized_term", first_block["term_usage"])
            self.assertEqual("stored_only", first_block["query_usage"]["usage"])


if __name__ == "__main__":
    unittest.main()
