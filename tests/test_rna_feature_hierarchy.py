"""Gene/transcript/product identity regressions; no phenotype or benchmark truth."""
import copy
import csv
import json
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from candidate_qtl.reference_extract import _normalize_features


def feature(kind, identifier, parent="", derives=""):
    return {"seqid": "chr1", "type": kind, "start": 1, "end": 100,
            "strand": "+", "phase": ".", "attrs": {"ID": identifier, "Parent": parent, "Derives_from": derives}}


class RnaFeatureHierarchyTests(unittest.TestCase):
    def test_primary_and_products_keep_original_identity_in_both_orders(self):
        for primary in ("miRNA_primary_transcript", "primary_transcript", "ncRNA"):
            for relation in ("Parent", "Derives_from"):
                for reverse in (False, True):
                    with self.subTest(primary=primary, relation=relation, reverse=reverse):
                        rows = [feature("gene", "G1"), feature(primary, "G1.1", "G1"),
                                feature("miRNA", "mature5"), feature("miRNA", "mature3")]
                        for row in rows[2:]:
                            row["attrs"][relation] = "G1.1"
                        if reverse:
                            rows.reverse()
                        attributes = {row["attrs"]["ID"]: copy.deepcopy(row["attrs"]) for row in rows}
                        normalized, qc = _normalize_features(rows)
                        self.assertEqual(["G1"], [r["attrs"]["ID"] for r in normalized if r["type"] == "gene"])
                        self.assertEqual(attributes, {r["attrs"]["ID"]: r["attrs"] for r in normalized})
                        self.assertEqual(["mature3", "mature5"], qc["rna_product_ids"])
                        self.assertEqual([], qc["unresolved_rna_product_ids"])
                        self.assertEqual(0, qc["placeholder_gene_count"])

    def test_gene_less_primary_keeps_one_placeholder_without_promoting_products(self):
        rows = [feature("miRNA", "mature5", "T"), feature("miRNA_primary_transcript", "T")]
        normalized, qc = _normalize_features(rows)
        self.assertEqual(["gene:T"], [r["attrs"]["ID"] for r in normalized if r["type"] == "gene"])
        self.assertEqual(1, qc["transcripts_without_gene_parent"])
        self.assertEqual(["gene:T"], rows[1]["_resolved_gene_parents"])
        self.assertEqual([], qc["unresolved_rna_product_ids"])

    def test_nested_transcripts_with_missing_gene_are_order_invariant(self):
        rows = [feature("primary_transcript", "outer", "missing_gene"), feature("transcript", "inner", "outer")]
        rows[0]["end"] = 25
        forward, forward_qc = _normalize_features(copy.deepcopy(rows))
        backward, backward_qc = _normalize_features(copy.deepcopy(rows[::-1]))
        self.assertEqual(forward_qc, backward_qc)
        for normalized in (forward, backward):
            self.assertEqual(["missing_gene"], [r["attrs"]["ID"] for r in normalized if r["type"] == "gene"])
            self.assertEqual(100, next(r for r in normalized if r["type"] == "gene")["end"])
            for row in normalized[:2]:
                self.assertEqual(["missing_gene"], row["_resolved_gene_parents"])

    def test_unresolved_product_never_creates_a_gene(self):
        normalized, qc = _normalize_features([feature("miRNA", "product", derives="missing_primary")])
        self.assertFalse(any(r["type"] == "gene" for r in normalized))
        self.assertEqual(["product"], qc["unresolved_rna_product_ids"])
        self.assertIn("missing_primary", qc["missing_parent_ids"])

    def test_cross_type_collision_and_cycles_remain_explicit_without_placeholder(self):
        examples = (
            ([feature("gene", "same"), feature("miRNA_primary_transcript", "same", "same")], "cross_type_id_collisions", "same"),
            ([feature("primary_transcript", "T1", "T2"), feature("primary_transcript", "T2", "T1")], "cyclic_parent_ids", "T1"),
            ([feature("gene", "G1"), feature("gene", "G2"), feature("miRNA_primary_transcript", "T", "G1,G2")], "ambiguous_hierarchy_ids", "T"),
        )
        for rows, key, identifier in examples:
            with self.subTest(key=key):
                _, qc = _normalize_features(rows)
                self.assertEqual(0, qc["placeholder_gene_count"])
                self.assertIn(identifier, qc[key])
                self.assertTrue(all(not r.get("_resolved_gene_parents") for r in rows if "transcript" in r["type"]))

    def test_cds_parent_product_is_not_recycled_as_a_transcript(self):
        rows = [feature("gene", "G1"), feature("primary_transcript", "T", "G1"),
                feature("miRNA", "product", "T"), feature("CDS", "cds", "product")]
        normalized, qc = _normalize_features(rows)
        self.assertEqual(0, qc["placeholder_gene_count"])
        self.assertEqual(["product"], qc["unresolved_cds_parent_ids"])
        self.assertEqual(len(rows), len(normalized))

    def test_real_extract_cli_keeps_primary_and_products_out_of_gene_roster(self):
        for primary, relation in (("miRNA_primary_transcript", "Parent"),
                                  ("miRNA_primary_transcript", "Derives_from"),
                                  ("primary_transcript", "Parent"), ("ncRNA", "Parent")):
            with self.subTest(primary=primary, relation=relation), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fasta = root / "genome.fa"
                fasta.write_text(">chr1\n" + "A" * 100 + "\n")
                gff = root / "genes.gff3"
                gff.write_text("chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=G1\n"
                               f"chr1\ttest\t{primary}\t1\t100\t.\t+\t.\tID=T;Parent=G1\n"
                               f"chr1\ttest\tmiRNA\t11\t31\t.\t+\t.\tID=mature5;{relation}=T\n"
                               f"chr1\ttest\tmiRNA\t71\t91\t.\t+\t.\tID=mature3;{relation}=T\n")
                before = gff.read_bytes()
                out = root / "out"
                executed = subprocess.run([sys.executable, "-m", "candidate_qtl", "extract", "--fasta", str(fasta),
                                           "--gff3", str(gff), "--interval", "chr1:1-100", "--out-dir", str(out)],
                                          capture_output=True, text=True)
                self.assertEqual(0, executed.returncode, executed.stderr)
                with (out / "candidates.tsv").open() as stream:
                    self.assertEqual(["G1"], [r["gene_id"] for r in csv.DictReader(stream, delimiter="\t")])
                with (out / "transcripts.tsv").open() as stream:
                    self.assertEqual(["T"], [r["transcript_id"] for r in csv.DictReader(stream, delimiter="\t")])
                manifest = json.loads((out / "extract_manifest.json").read_text())
                self.assertEqual(0, manifest["placeholder_gene_count"])
                self.assertEqual(before, gff.read_bytes())

    def test_public_araport_raw_four_lines_via_real_extract_cli(self):
        fixture = Path(__file__).parent / "fixtures/gff3/araport11_mirna_raw.gff3"
        raw = fixture.read_bytes()
        self.assertEqual("75ceba97e881c1cef3baa59e1cd09ee0a8de974d9416fb0ead9157b2eed2b881",
                         hashlib.sha256(raw).hexdigest())
        features = [line.split("\t") for line in raw.decode("utf-8").splitlines() if not line.startswith("#")]
        self.assertEqual(4, len(features))
        chrom = features[0][0]
        end = max(int(row[4]) for row in features)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "synthetic_coordinate_reference.fa"
            # The annotation is real/raw; this N-only sequence tests identity,
            # not Arabidopsis sequence biology or scientific performance.
            fasta.write_text(f">{chrom}\n" + "N" * end + "\n")
            out = root / "out"
            executed = subprocess.run([sys.executable, "-m", "candidate_qtl", "extract", "--fasta", str(fasta),
                                       "--gff3", str(fixture), "--interval", f"{chrom}:1-{end}", "--out-dir", str(out)],
                                      capture_output=True, text=True)
            self.assertEqual(0, executed.returncode, executed.stderr)
            with (out / "candidates.tsv").open() as stream:
                self.assertEqual(["AT3G09285"], [r["gene_id"] for r in csv.DictReader(stream, delimiter="\t")])
            with (out / "transcripts.tsv").open() as stream:
                self.assertEqual(["AT3G09285.1"], [r["transcript_id"] for r in csv.DictReader(stream, delimiter="\t")])
            manifest = json.loads((out / "extract_manifest.json").read_text())
            self.assertEqual(0, manifest["placeholder_gene_count"])
            self.assertEqual([], manifest["gff_dialect_qc"]["unresolved_rna_product_ids"])
            self.assertEqual(raw, fixture.read_bytes())


if __name__ == "__main__":
    unittest.main()
