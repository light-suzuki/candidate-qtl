from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class StructuralVariantTests(unittest.TestCase):
    def test_syri_and_read_calls_are_reconciled_without_scoring_qtl_overlap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); reference=root/"ref.fa"; query=root/"query.fa"; candidates=root/"candidates.tsv"; syri=root/"syri.out"; reads=root/"sniffles.tsv"; out=root/"out"; decision=root/"decision.json"
            reference.write_text(">chr1\n"+"A"*1000+"\n",encoding="utf-8"); query.write_text(">chr1\n"+"A"*1000+"\n",encoding="utf-8")
            candidates.write_text("gene_id\tchrom\tstart\tend\tstrand\nGene1\tchr1\t100\t200\t+\n",encoding="utf-8")
            syri.write_text(
                "chr1\t120\t180\t-\t-\tchr1\t120\t180\tINV1\t-\tINV\t-\n"
                "chr1\t300\t360\t-\t-\tchr1\t300\t360\tDEL1\tSYN1\tDEL\t-\n",
                encoding="utf-8",
            )
            reads.write_text(
                "chrom\tstart\tend\tsvtype\tcall_id\tsupport\tmapq\tsample\tanalysis_chain\tflags\n"
                "chr1\t125\t185\tINV\tR1\t12\t50\tparentA\treads:ont:sniffles\t\n"
                "chr1\t305\t365\tDEL\tR2\t8\t10\tparentA\treads:ont:sniffles\trepeat_rich\n",
                encoding="utf-8",
            )
            decision.write_text('{"mode":"broaden_to_regulatory_sv_te_cnv"}',encoding="utf-8")
            self.assertEqual(0,main(["sv-plan","--reference",str(reference),"--query",str(query),"--qtl-interval","chr1:1-500","--candidates",str(candidates),"--syri-out",str(syri),"--call-table",f"sniffles={reads}","--decision",str(decision),"--out-dir",str(out)]))
            calls=pd.read_csv(out/"canonical_sv.tsv",sep="\t",dtype=str).fillna("")
            self.assertIn("gene_body",set(calls["overlap_class"])); self.assertIn("intergenic_qtl",set(calls["overlap_class"]))
            concordance=pd.read_csv(out/"sv_concordance.tsv",sep="\t",dtype=str).fillna("")
            self.assertEqual({"assembly_read_concordant"},set(concordance["status"]))
            self.assertEqual({"False"},set(concordance["independent_biological_evidence"]))
            evidence=pd.read_csv(out/"sv_evidence.tsv",sep="\t",dtype=str).fillna("")
            self.assertEqual({"context"},set(evidence["direction"])); self.assertEqual({"structural_variant_context"},set(evidence["evidence_code"]))
            manifest=json.loads((out/"sv_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("broaden_to_regulatory_sv_te_cnv",manifest["trigger"])
            self.assertIn("cannot promote",manifest["ranking_rule"]); self.assertIn("same sample",manifest["independence_rule"])
            script=(out/"run_assembly_sv.sh").read_text(encoding="utf-8")
            self.assertIn("minimap2 -ax asm5",script); self.assertIn("syri -c",script)
            self.assertFalse(pd.read_csv(out/"review_required.tsv",sep="\t").empty)


if __name__=="__main__": unittest.main()
