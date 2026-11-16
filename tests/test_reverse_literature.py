"""Route B stage 2: interval gene -> reference homolog -> literature."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.core import safe_id
from candidate_qtl.reverse_literature import _query_for, _phenotype_terms

LOOKUP = "\n".join([
    "gene_id\treference_gene\treference_function",
    "Psat_A\tCLV1\tReceptor protein kinase CLAVATA1",
    "Psat_B\tCLV1\tReceptor protein kinase CLAVATA1",
    "Psat_C\tWUS\tProtein WUSCHEL",
    "Psat_D\t\t",
]) + "\n"

PHENOTYPE = "fasciated stem with enlarged meristem"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _cache(root: Path) -> Path:
    cache = root / "queries"
    cache.mkdir()
    terms = _phenotype_terms(PHENOTYPE)
    hit_query = _query_for("CLV1", terms)
    (cache / "clv1.json").write_text(json.dumps({"query": hit_query, "results": [
        {"pmid": "8402901", "title": "CLAVATA1 controls meristem size and fasciated stem development",
         "abstractText": "clv1 mutants show a fasciated stem and an enlarged meristem.",
         "journalTitle": "Cell", "pubYear": "1993", "doi": "10.1/x", "authorString": "Clark SE, Running MP"},
    ]}), encoding="utf-8")
    # WUS was queried and genuinely returned nothing
    (cache / "wus.json").write_text(json.dumps({"query": _query_for("WUS", terms), "results": []}), encoding="utf-8")
    return cache


class QueryTests(unittest.TestCase):
    def test_query_focuses_the_reference_gene_by_phenotype_terms(self) -> None:
        query = _query_for("CLV1", _phenotype_terms(PHENOTYPE))
        self.assertTrue(query.startswith('"CLV1" AND ('))
        self.assertIn('"fasciated"', query)
        self.assertIn('"meristem"', query)
        # generic words carry no signal and must not widen the query
        self.assertNotIn('"with"', query)

    def test_cache_key_is_stable(self) -> None:
        query = _query_for("CLV1", _phenotype_terms(PHENOTYPE))
        self.assertEqual(safe_id(query), safe_id(_query_for("CLV1", _phenotype_terms(PHENOTYPE))))


class ReverseLiteratureTests(unittest.TestCase):
    def _run(self, root: Path, *extra: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        out = root / "revlit"
        self.assertEqual(0, main(["reverse-literature", "--lookup", str(_write(root, "lookup.tsv", LOOKUP)),
                                  "--out-dir", str(out), "--phenotype", PHENOTYPE, *extra]))
        literature = pd.read_csv(out / "reverse_literature.tsv", sep="\t", dtype=str).fillna("")
        status = pd.read_csv(out / "reverse_literature_status.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "reverse_literature_manifest.json").read_text(encoding="utf-8"))
        return literature, status, manifest

    def test_articles_are_linked_to_the_interval_genes_that_share_the_homolog(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature, status, manifest = self._run(root, "--cache", str(_cache(root)))
            self.assertEqual(1, len(literature))
            row = literature.iloc[0]
            self.assertEqual("CLV1", row["reference_gene"])
            # both interval genes hitting the same reference gene share the article
            self.assertEqual("Psat_A;Psat_B", row["interval_gene_ids"])
            self.assertEqual("PMID:8402901", row["source_id"])
            self.assertEqual("1993", row["year"])
            self.assertIn("fasciated", row["matched_phenotype_terms"])
            self.assertIn("meristem", row["matched_phenotype_terms"])
            # retrieval is not relevance
            self.assertEqual("literature_link_candidate", row["claim_type"])
            self.assertEqual("none", row["rank_contribution"])
            self.assertEqual("local_query_cache", manifest["provider"])
            self.assertEqual("none", manifest["external_transmission"])

    def test_a_genuine_no_hit_is_distinct_from_a_query_that_never_ran(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _literature, status, manifest = self._run(root, "--cache", str(_cache(root)))
            by_gene = status.set_index("reference_gene")
            self.assertEqual("retrieved", by_gene.loc["CLV1", "query_status"])
            self.assertEqual("no_hit", by_gene.loc["WUS", "query_status"])
            self.assertEqual(1, manifest["no_hit"])
            self.assertEqual(0, manifest["not_evaluated_offline"])

    def test_offline_without_a_cache_is_not_evaluated_never_no_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature, status, manifest = self._run(root)
            self.assertEqual(0, len(literature))
            self.assertEqual({"not_evaluated_offline", "skipped_no_reference_gene"}, set(status["query_status"]))
            self.assertEqual(2, manifest["not_evaluated_offline"])
            self.assertEqual(0, manifest["no_hit"])
            note = status.set_index("reference_gene").loc["CLV1", "note"]
            self.assertIn("absence of literature is not established", note)

    def test_genes_without_a_characterised_relative_stay_an_explicit_bucket(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _literature, status, manifest = self._run(root, "--cache", str(_cache(root)))
            unresolved = status[status["reference_gene"] == ""].iloc[0]
            self.assertEqual("skipped_no_reference_gene", unresolved["query_status"])
            self.assertEqual("1", unresolved["interval_gene_count"])
            self.assertEqual(1, manifest["interval_genes_without_reference"])


if __name__ == "__main__":
    unittest.main()
