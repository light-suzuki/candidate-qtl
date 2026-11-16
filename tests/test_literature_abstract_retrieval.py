"""The claim classifier reads abstracts, so the fetchers must ask for them (issue #484).

Europe PMC returns `abstractText` only for `resultType=core`. Both fetchers asked
for `lite`, so the classifier -- whose cue lexicon is built from phrases that live
in abstracts ("loss-of-function", "complementation", "T-DNA insertion") -- was
handed titles alone. Every Route B claim then fell through to `review_required`,
and the transfer gate, which requires `supported_mechanistic`, could never pass.

No network here: the request URL is captured, and the classifier is exercised on
text directly.
"""
from __future__ import annotations

import json
import tempfile
import unittest
import urllib.parse
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import reverse_literature, source_verify
from candidate_qtl.claim_evidence import classifiable_text, classify_claim
from candidate_qtl.cli import main

ABSTRACT = ("Loss-of-function alleles of the gene abolished the phenotype, and complementation "
            "restored it in transgenic lines.")
TITLE = "A gene with a role in inflorescence architecture"


class _Response(BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _captured_urlopen(seen: list[str]):
    def fake(request, timeout=30):  # noqa: ARG001
        seen.append(request.full_url)
        return _Response(json.dumps({"resultList": {"result": []}}).encode("utf-8"))
    return fake


class LiteratureAbstractRetrievalTests(unittest.TestCase):
    def test_reverse_literature_requests_abstracts(self) -> None:
        seen: list[str] = []
        with patch("candidate_qtl.reverse_literature.urllib.request.urlopen", _captured_urlopen(seen)):
            reverse_literature._fetch("TFL1", 3)
        self.assertEqual(1, len(seen))
        query = urllib.parse.parse_qs(urllib.parse.urlparse(seen[0]).query)
        self.assertEqual(["core"], query["resultType"])

    def test_source_verify_requests_abstracts(self) -> None:
        seen: list[str] = []
        with patch("candidate_qtl.source_verify.urllib.request.urlopen", _captured_urlopen(seen)):
            source_verify._fetch_europe_pmc("PMID:1")
        self.assertEqual(1, len(seen))
        query = urllib.parse.parse_qs(urllib.parse.urlparse(seen[0]).query)
        self.assertEqual(["core"], query["resultType"])

    def test_mechanistic_cues_are_only_reachable_with_an_abstract(self) -> None:
        # Why the result type matters: the same record classifies differently
        # depending on whether the abstract came back.
        with_abstract = {"title": TITLE, "abstractText": ABSTRACT}
        title_only = {"title": TITLE}
        terms = ["inflorescence"]

        _, support_with, cues = classify_claim(classifiable_text(with_abstract), ["inflorescence"], terms, "title_abstract")
        self.assertEqual("supported_direct_genetic_test", support_with)
        self.assertIn("loss-of-function", cues)

        _, support_without, _ = classify_claim(classifiable_text(title_only), ["inflorescence"], terms, "title_only")
        self.assertNotEqual("supported_direct_genetic_test", support_without)

    def test_literature_table_records_whether_an_abstract_arrived(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lookup = root / "lookup.tsv"
            lookup.write_text("gene_id\treference_gene\treference_function\n"
                              "G1\tTFL1\tPEBP protein\n", encoding="utf-8")
            cache = root / "cache"
            cache.mkdir()
            terms = reverse_literature._phenotype_terms("determinate inflorescence")
            (cache / "hit.json").write_text(json.dumps({
                "query": reverse_literature._query_for("TFL1", terms),
                "results": [{"pmid": "1", "title": TITLE, "abstractText": ABSTRACT,
                             "journalTitle": "J", "pubYear": "2026", "authorString": "A"},
                            {"pmid": "2", "title": TITLE, "journalTitle": "J",
                             "pubYear": "2026", "authorString": "B"}],
            }), encoding="utf-8")
            out = root / "revlit"
            self.assertEqual(0, main(["reverse-literature", "--lookup", str(lookup),
                                      "--phenotype", "determinate inflorescence",
                                      "--cache", str(cache), "--limit", "2",
                                      "--out-dir", str(out), "--case-id", "abs"]))
            literature = pd.read_csv(out / "reverse_literature.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["True", "False"], list(literature["abstract_available"]))
            manifest = json.loads((out / "reverse_literature_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("core", manifest["result_type"])
            self.assertEqual(1, manifest["abstract_available_article_count"])


if __name__ == "__main__":
    unittest.main()
