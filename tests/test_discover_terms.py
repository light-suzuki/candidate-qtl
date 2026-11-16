from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


GENOME = ">chr1\n" + "ATG" * 300 + "\n"


def _gff(*genes: tuple[str, int, int, str]) -> str:
    lines = ["##gff-version 3"]
    for gid, start, end, note in genes:
        lines.append(f"chr1\tfixture\tgene\t{start}\t{end}\t.\t+\t.\tID={gid};Note={note}")
        lines.append(f"chr1\tfixture\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.1;Parent={gid}")
        lines.append(f"chr1\tfixture\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gid}.1")
    return "\n".join(lines) + "\n"


class DiscoverTermTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fasta = self.root / "genome.fa"; self.fasta.write_text(GENOME, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, gff_text: str, phenotype: str, extra: list[str] | None = None) -> pd.DataFrame:
        gff = self.root / "genes.gff3"; gff.write_text(gff_text, encoding="utf-8")
        out = self.root / f"out_{abs(hash(gff_text)) % 10000}"
        args = ["discover", "--fasta", str(self.fasta), "--gff3", str(gff), "--interval", "chr1:1-900",
                "--phenotype", phenotype, "--out-dir", str(out), "--offline", "--no-auto-local-db"]
        for term in extra or []:
            args += ["--approved-terms", term] if term.endswith(".tsv") else ["--phenotype-term", term]
        self.assertEqual(0, main(args))
        return pd.read_csv(out / "04_evidence.tsv", sep="\t").fillna("")

    def test_generic_anatomy_alone_is_context_not_scored(self) -> None:
        evidence = self._run(
            _gff(("StemGene", 100, 200, "stem-associated protein"), ("FascGene", 300, 400, "fasciation regulator")),
            "stem fasciation",
        )
        stem = dict(zip(evidence[evidence.gene_id == "StemGene"].evidence_code, evidence[evidence.gene_id == "StemGene"].direction))
        fasc = set(evidence[evidence.gene_id == "FascGene"].evidence_code)
        # A generic anatomy word ("stem") only produces context, never scored support.
        self.assertEqual(stem.get("trait_context"), "context")
        self.assertNotIn("trait_function_match", stem)
        # A specific trait word ("fasciation") is scored support.
        self.assertIn("trait_function_match", fasc)

    def test_exact_phrase_match_is_high_confidence(self) -> None:
        evidence = self._run(_gff(("PhraseGene", 100, 300, "meristem enlargement factor")), "meristem enlargement")
        row = evidence[(evidence.gene_id == "PhraseGene") & (evidence.evidence_code == "trait_function_match")].iloc[0]
        self.assertEqual(row["confidence"], "high")

    def test_approved_terms_file_controls_vocabulary(self) -> None:
        approved = self.root / "approved.tsv"
        approved.write_text("term\tweight_class\nfasciation\tdirect\n", encoding="utf-8")
        evidence = self._run(
            _gff(("StemGene", 100, 200, "stem-associated protein"), ("FascGene", 300, 400, "fasciation regulator")),
            "stem fasciation",
            extra=[str(approved)],
        )
        # Only the approved term drives a match; the generic "stem" is not in the vocabulary.
        self.assertTrue(evidence[evidence.gene_id == "StemGene"]["evidence_code"].isin({"trait_function_match", "trait_context"}).sum() == 0)
        self.assertIn("trait_function_match", set(evidence[evidence.gene_id == "FascGene"].evidence_code))

    def test_manifest_records_term_classes_and_source(self) -> None:
        gff = self.root / "g.gff3"; gff.write_text(_gff(("FascGene", 300, 400, "fasciation regulator")), encoding="utf-8")
        out = self.root / "m"
        self.assertEqual(0, main(["discover", "--fasta", str(self.fasta), "--gff3", str(gff), "--interval", "chr1:1-900",
                                  "--phenotype", "stem fasciation", "--out-dir", str(out), "--offline", "--no-auto-local-db"]))
        manifest = json.loads((out / "discover_manifest.json").read_text())
        terms = manifest["phenotype_terms"]
        self.assertEqual(terms["source"], "auto_split")
        self.assertIn("stem", terms["broad"])          # generic anatomy set aside
        self.assertIn("fasciation", terms["scored"])   # specific trait token


if __name__ == "__main__":
    unittest.main()
