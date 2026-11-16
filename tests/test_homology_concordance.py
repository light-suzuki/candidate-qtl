from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


def hit(q,s,bits,q_len=100,s_len=100,aln=90,identity=70,desc="protein"):
    return f"{q}\t{s}\t{identity}\t{aln}\t{q_len}\t{s_len}\t1\t{aln}\t1\t{aln}\t1e-30\t{bits}\t{desc}\n"


class HomologyConcordanceTests(unittest.TestCase):
    def test_similarity_domain_and_orthology_paths_are_compared_without_vote_inflation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); blast=root/"blast.tsv"; diamond=root/"diamond.tsv"; mmseqs=root/"mmseqs.tsv"
            blast.write_text(hit("Clear","FAM1|A",200)+hit("Paralog","FAM2|A",200)+hit("Paralog","FAM2|B",190)+hit("Truncated","FAM3|A",100,q_len=100,s_len=300,aln=60),encoding="utf-8")
            diamond.write_text(hit("Clear","FAM1|A",205)+hit("Paralog","FAM2|B",198)+hit("Truncated","FAM4|X",95,q_len=100,s_len=300,aln=60),encoding="utf-8")
            mmseqs.write_text(hit("Clear","FAM1|A",210)+hit("Paralog","FAM2|A",202)+hit("Truncated","FAM5|Y",90,q_len=100,s_len=300,aln=60),encoding="utf-8")
            pfam=root/"pfam.tsv"; interpro=root/"interpro.tsv"
            pfam.write_text("gene_id\tdomain_id\tstart\tend\nClear\tPF1\t1\t50\nParalog\tPF2\t5\t70\n",encoding="utf-8")
            interpro.write_text("gene_id\tdomain_id\tstart\tend\nClear\tPF1\t1\t50\nParalog\tPF3\t5\t70\n",encoding="utf-8")
            oma=root/"oma.tsv"; orthofinder=root/"orthofinder.tsv"
            oma.write_text("gene_id\torthogroup\trelationship\nClear\tOG1\tone_to_one\nParalog\tOG2\tone_to_many\n",encoding="utf-8")
            orthofinder.write_text("gene_id\torthogroup\trelationship\nClear\tOG1\tone_to_one\nParalog\tOG3\tmany_to_many\n",encoding="utf-8")
            candidates=root/"candidates.tsv"; candidates.write_text("gene_id\tdry_group\nClear\ttop\nParalog\trescue\nTruncated\tevidence_insufficient\n",encoding="utf-8")
            meta=root/"blast.json"; meta.write_text(json.dumps({"tool_version":"2.17","database_version":"fixture1","taxon_scope":"plants","command":"blastp ..."}),encoding="utf-8")
            out=root/"out"
            self.assertEqual(0,main(["homology-concordance","--hit-table",f"blast={blast}","--hit-table",f"diamond={diamond}","--hit-table",f"mmseqs2={mmseqs}",
                                     "--domain-table",f"pfam={pfam}","--domain-table",f"interpro={interpro}","--orthology-table",f"omamer={oma}","--orthology-table",f"orthofinder={orthofinder}",
                                     "--provider-metadata",f"blast={meta}","--candidates",str(candidates),"--out-dir",str(out),"--top-n","2"]))
            hom=pd.read_csv(out/"homology_concordance.tsv",sep="\t",dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("exact_hit",hom.loc["Clear","concordance"]); self.assertIn("blast",hom.loc["Paralog","paralog_ambiguous_providers"])
            self.assertEqual("domain_only",hom.loc["Truncated","concordance"])
            domains=pd.read_csv(out/"domain_concordance.tsv",sep="\t",dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("exact_architecture",domains.loc["Clear","concordance"]); self.assertEqual("False",domains.loc["Clear","trait_linked"])
            orth=pd.read_csv(out/"orthology_concordance.tsv",sep="\t",dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("exact_orthogroup",orth.loc["Clear","concordance"]); self.assertEqual("True",orth.loc["Paralog","paralog_ambiguous"])
            escalation=pd.read_csv(out/"homology_escalation_plan.tsv",sep="\t",dtype=str).fillna("")
            self.assertEqual({"Clear","Paralog"},set(escalation["gene_id"])); self.assertNotIn("Truncated",set(escalation["gene_id"]))
            evidence=pd.read_csv(out/"homology_evidence.tsv",sep="\t",dtype=str).fillna("")
            self.assertEqual({"context"},set(evidence["direction"])); self.assertEqual({"comparative_annotation_chain"},set(evidence["independence_group"]))
            manifest=json.loads((out/"homology_concordance_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("2.17",manifest["provider_metadata"]["blast"]["tool_version"]); self.assertIn("not independent",manifest["independence_rule"])


if __name__=="__main__": unittest.main()
