"""reverse-literature must search on something the literature index knows (issue #480).

The family index writes the BLAST subject verbatim, so a reference gene reaches
this stage as a UniProt FASTA defline (``sp|P93003|TFL1_ARATH``). Europe PMC
returns nothing for that string -- 0 hits, where ``TFL1`` returns thousands --
so every reference came back ``no_hit``, the Route B semantic bundle was handed
no citation to verify, and the lane could never freeze a verified candidate.

These tests run entirely from a query-keyed cache: no network, and the assertion
is about which identifier is searched, not about what Europe PMC happens to hold.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.core import safe_id
from candidate_qtl.reverse_literature import _phenotype_terms, _query_for, _searchable_reference
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle

PHENOTYPE = "determinate inflorescence"

LOOKUP = "\n".join([
    "gene_id\treference_gene\treference_gene_symbol\treference_function",
    # The B0 lookup may carry a bare accession and its symbol in separate
    # columns; the accession itself has no symbol to recover from a defline.
    "G_TFL\tP93003\tTFL1\tphosphatidylethanolamine-binding protein",
    "G_TR\ttr|A0A1I9LLD6|A0A1I9LLD6_ARATH\tMTR_3g010000\tuncharacterised",
    "G_BAD\tweird|x\t\tunparsable identifier",
    "G_NONE\t\t\t",
]) + "\n"


class ReverseLiteratureIdentifierTests(unittest.TestCase):
    def _run(self, root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        lookup = root / "interval_reverse_lookup.tsv"
        lookup.write_text(LOOKUP, encoding="utf-8")

        # Cache keyed on the query the fixed code should build: the symbol, not the defline.
        terms = _phenotype_terms(PHENOTYPE)
        cache_dir = root / "cache"
        cache_dir.mkdir()
        (cache_dir / "tfl1.json").write_text(json.dumps({
            "query": _query_for("TFL1", terms),
            "results": [{"pmid": "41635692", "title": "Meristem fate: to terminate, or not?",
                         "abstractText": "determinate inflorescence architecture",
                         "journalTitle": "J Test", "pubYear": "2026", "authorString": "Tester A"}],
        }), encoding="utf-8")
        (cache_dir / "trembl.json").write_text(json.dumps({
            "query": _query_for("A0A1I9LLD6", terms),
            "results": [{"pmid": "1", "title": "Accession fallback", "abstractText": "",
                         "journalTitle": "J Test", "pubYear": "2026", "authorString": "Tester B"}],
        }), encoding="utf-8")

        out = root / "revlit"
        # No --allow-external-api: everything below is served from the cache.
        self.assertEqual(0, main(["reverse-literature", "--lookup", str(lookup), "--phenotype", PHENOTYPE,
                                  "--cache", str(cache_dir), "--out-dir", str(out), "--case-id", "gate"]))
        status = pd.read_csv(out / "reverse_literature_status.tsv", sep="\t", dtype=str).fillna("")
        literature = pd.read_csv(out / "reverse_literature.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "reverse_literature_manifest.json").read_text(encoding="utf-8"))
        return status, literature, manifest

    def test_bare_accession_uses_its_separately_supplied_symbol(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            status, literature, _ = self._run(root)
            row = status[status["reference_gene"].eq("P93003")].iloc[0]
            self.assertEqual("TFL1", row["query_term"])
            self.assertEqual("TFL1", row["reference_gene_symbol"])
            self.assertEqual("reference_gene_symbol", row["identifier_source"])
            self.assertEqual("retrieved", row["query_status"])
            self.assertEqual("1", row["article_count"])
            # The provenance keeps the original accession while the query uses
            # the separately supplied symbol.
            self.assertEqual("P93003", literature.iloc[0]["reference_gene"])
            self.assertIn('"TFL1"', literature.iloc[0]["query"])
            self.assertIn("PMID:41635692", set(literature["source_id"]))

            # The real B semantic-bundle consumer accepts the CLI's link rows
            # and attaches them to the interval gene without upgrading support.
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "G_TFL", "family_id": "F_TFL"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            bundle_result = build_route_b_semantic_bundle(
                str(atlas), str(root / "bundle"), phenotype=PHENOTYPE,
                literature_arg=str(root / "revlit" / "reverse_literature.tsv"),
                literature_status_arg=str(root / "revlit" / "reverse_literature_status.tsv"),
            )
            units = pd.read_csv(bundle_result["units"], sep="\t", dtype=str).fillna("")
            self.assertEqual("1", units.iloc[0]["literature_link_count"])
            bundle = json.loads(Path(bundle_result["bundle"]).read_text(encoding="utf-8"))
            self.assertIn("PMID:41635692", bundle["units"][0]["stable_source_ids"])

    def test_trembl_defline_falls_back_to_the_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, _, _ = self._run(Path(tmp))
            row = status[status["reference_gene"].eq("tr|A0A1I9LLD6|A0A1I9LLD6_ARATH")].iloc[0]
            self.assertEqual("A0A1I9LLD6", row["query_term"])
            self.assertEqual("uniprot_defline_accession", row["identifier_source"])
            self.assertEqual("retrieved", row["query_status"])

    def test_unusable_identifier_is_not_reported_as_absent_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, _, manifest = self._run(Path(tmp))
            row = status[status["reference_gene"].eq("weird|x")].iloc[0]
            self.assertEqual("skipped_unusable_reference_identifier", row["query_status"])
            self.assertEqual("", row["query_term"])
            # Not searched and no hits must not collapse into one number.
            self.assertEqual(1, manifest["skipped_unusable_reference_identifier"])
            self.assertEqual(0, manifest["no_hit"])
            self.assertEqual(2, manifest["query_attempted"])

    def test_missing_reference_keeps_its_own_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, _, _ = self._run(Path(tmp))
            row = status[status["reference_gene"].eq("")].iloc[0]
            self.assertEqual("skipped_no_reference_gene", row["query_status"])
            self.assertEqual("none", row["identifier_source"])

    def test_conflicting_symbols_are_not_arbitrarily_searched(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lookup = root / "conflicting.tsv"
            lookup.write_text(
                "gene_id\treference_gene\treference_gene_symbol\n"
                "G_A\tsp|P93003|TFL1_ARATH\tTFL1\n"
                "G_B\tsp|P93003|TFL1_ARATH\tFT\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["reverse-literature", "--lookup", str(lookup),
                                      "--phenotype", PHENOTYPE, "--out-dir", str(out)]))
            status = pd.read_csv(out / "reverse_literature_status.tsv", sep="\t", dtype=str).fillna("")
            row = status.iloc[0]
            self.assertEqual("sp|P93003|TFL1_ARATH", row["reference_gene"])
            self.assertEqual("FT;TFL1", row["reference_gene_symbol"])
            self.assertEqual("", row["query_term"])
            self.assertEqual("conflicting_reference_gene_symbols", row["identifier_source"])
            self.assertEqual("not_evaluated_conflicting_reference_gene_symbols", row["query_status"])

    def test_searchable_reference_derivation(self) -> None:
        self.assertEqual(("TFL1", "uniprot_defline_symbol"), _searchable_reference("sp|P93003|TFL1_ARATH"))
        self.assertEqual(("A0A1I9LLD6", "uniprot_defline_accession"),
                         _searchable_reference("tr|A0A1I9LLD6|A0A1I9LLD6_ARATH"))
        self.assertEqual(("TFL1", "as_supplied"), _searchable_reference("TFL1"))
        self.assertEqual(("", "none"), _searchable_reference(""))
        self.assertEqual(("", "unusable_identifier"), _searchable_reference("weird|x"))

    def test_cache_key_uses_the_derived_query(self) -> None:
        # Guards the cache contract: a cached defline query must not be what a rerun looks up.
        terms = _phenotype_terms(PHENOTYPE)
        self.assertNotEqual(safe_id(_query_for("sp|P93003|TFL1_ARATH", terms)),
                            safe_id(_query_for("TFL1", terms)))


if __name__ == "__main__":
    unittest.main()
