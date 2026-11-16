from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class VariantEvidenceTests(unittest.TestCase):
    def test_snpeff_vep_consequences_are_compared_as_unvalidated_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "ref.fa"; raw = root / "parents.vcf"; snpeff = root / "snpeff.vcf"; vep = root / "vep.vcf"
            fasta.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")
            raw.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tT,G\t.\tPASS\t.\n", encoding="utf-8")
            ann_t = "T|stop_gained&splice_region_variant|HIGH|GENE1|GENE1|transcript|TX1|protein_coding|1/2|c.10A>T|p.K4*|||||"
            ann_g = "G|missense_variant|MODERATE|GENE1|GENE1|transcript|TX2|protein_coding|1/2|c.10A>G|p.K4R|||||"
            snpeff.write_text("##fileformat=VCFv4.2\n##SnpEffVersion=\"5.2\"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + f"chr1\t10\t.\tA\tT,G\t.\tPASS\tANN={ann_t},{ann_g}\n", encoding="utf-8")
            vep.write_text(
                "##fileformat=VCFv4.2\n##VEP=112\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|HGVSc|HGVSp\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t10\t.\tA\tT,G\t.\tPASS\tCSQ=T|stop_gained|HIGH|GENE1|GENE1|Transcript|TX1|c.10A>T|p.K4*,G|5_prime_UTR_variant&upstream_gene_variant|MODIFIER|GENE1|GENE1|Transcript|TX2||\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["variant-evidence", "--vcf", str(raw), "--fasta", str(fasta), "--snpeff-vcf", str(snpeff), "--vep-vcf", str(vep), "--out-dir", str(out)]))
            concordance = pd.read_csv(out / "variant_annotation_concordance.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("partial", set(concordance["concordance"]))
            self.assertIn("discordant", set(concordance["concordance"]))
            evidence = pd.read_csv(out / "variant_evidence.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"context"}, set(evidence["direction"]))
            self.assertEqual({"low"}, set(evidence["confidence"]))
            manifest = json.loads((out / "variant_evidence_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["normalization_executed"])
            self.assertEqual("5.2", manifest["tool_versions"]["snpeff"]); self.assertEqual("112", manifest["tool_versions"]["vep"])
            self.assertIn("not independent", manifest["independence_rule"])
            self.assertIn("predicted_unvalidated", manifest["wet_gate"])
            script = (out / "run_variant_annotation.sh").read_text(encoding="utf-8")
            self.assertIn("bcftools norm", script); self.assertIn("--check-ref w", script); self.assertIn("-m -any", script)


if __name__ == "__main__": unittest.main()
