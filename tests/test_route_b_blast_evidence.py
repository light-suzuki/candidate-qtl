"""Synthetic residue/gate regressions; no model or biological accuracy claims."""
from __future__ import annotations

import ast
import csv
import itertools
import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.route_b_blast_evidence import IDENTITY_POLICY, aggregate_hsps, run_local_blastp


def hsp(qseq, sseq, *, gene="g", ref="r", qs=1, ss=1, qlen=None, slen=None, score=100):
    nq, ns = len(qseq.replace("-", "")), len(sseq.replace("-", ""))
    matches = sum(q == s and q != "-" for q, s in zip(qseq, sseq))
    return "\t".join(map(str, [
        gene, ref, 100 * matches / len(qseq), len(qseq), qlen or qs+nq-1,
        slen or ss+ns-1, qs, qs+nq-1, ss, ss+ns-1, 1e-20, score, qseq, sseq, "synthetic",
    ]))


class ResidueEvidenceTests(unittest.TestCase):
    def aggregate(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"hsp.tsv"
            p.write_text("\n".join(rows)+"\n", encoding="utf-8")
            return aggregate_hsps(p, 1e-5)[0]

    def run_local(self, rows, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def run(command, **_):
                self.assertIn("qseq sseq", command[command.index("-outfmt")+1])
                Path(command[command.index("-out")+1]).write_text("\n".join(rows)+"\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")
            result, meta = run_local_blastp(
                root/"q.fa", root/"db", root, 1, kwargs.get("min_identity",30),
                kwargs.get("min_coverage",.5), 1e-5, run_tool=run,
                find_binary=lambda _: "blastp", normalise_subject=lambda s: s,
            )
            with Path(meta["pair_evidence_table"]).open() as handle:
                evidence = list(csv.DictReader(handle, delimiter="\t"))
            legacy = Path(meta["hits_table"]).read_text()
            return result, meta, evidence, legacy

    def test_weighted_28_not_max_40(self):
        rows = [hsp("A"*200,"A"*80+"C"*120,qlen=1000,slen=1000),
                hsp("A"*800,"A"*200+"C"*600,qs=201,ss=201,qlen=1000,slen=1000)]
        result, meta, evidence, _ = self.run_local(rows)
        self.assertEqual(result, [])
        self.assertEqual(float(evidence[0]["percent_identity"]),28)
        self.assertEqual(float(evidence[0]["max_hsp_identity"]),40)
        self.assertEqual(float(evidence[0]["query_coverage"]),1)
        self.assertEqual(meta["identity_policy"],IDENTITY_POLICY)
        self.assertFalse(meta["max_hsp_identity_is_gate_input"])

    def test_duplicate_hsps_count_once(self):
        row=hsp("AAAA","AACA")
        self.assertEqual(self.aggregate([row]),self.aggregate([row,row,row]))

    def test_split_alignment_preserves_gate_measurements(self):
        whole=self.aggregate([hsp("AAAAAA","AACCAA")])[0]
        split=self.aggregate([hsp("AAA","AAC",qlen=6,slen=6),
                              hsp("AAA","CAA",qs=4,ss=4,qlen=6,slen=6)])[0]
        for key in ("percent_identity","query_coverage","covered_query_residues","identical_query_residues"):
            self.assertEqual(whole[key],split[key])

    def test_overlap_matching_positions_count_once(self):
        rows=[hsp("AAAA","AACA",qlen=6,slen=6),
              hsp("AAAA","CAAA",qs=3,ss=3,qlen=6,slen=6)]
        p=self.aggregate(rows)[0]
        self.assertEqual(p["covered_query_residues"],6)
        self.assertEqual(p["identical_query_residues"],5)
        self.assertEqual(p["conflicting_query_residues"],0)

    def test_overlap_disagreement_is_conservative(self):
        p=self.aggregate([hsp("AAAA","AAAA"),hsp("AAAA","ACAA")])[0]
        self.assertEqual(p["percent_identity"],75)
        self.assertEqual(p["conflicting_query_residues"],1)

    def test_alternative_subject_positions_remain_conflicted(self):
        p=self.aggregate([hsp("AAAA","AAAA",slen=8),hsp("AAAA","AAAA",ss=5,slen=8)])[0]
        self.assertEqual(p["percent_identity"],0)
        self.assertEqual(p["query_coverage"],1)

    def test_gap_policy_is_explicit(self):
        p=self.aggregate([hsp("AA-AA","AAC-A")])[0]
        self.assertEqual(p["covered_query_residues"],4)
        self.assertEqual(p["percent_identity"],75)

    def test_row_permutations_are_identical(self):
        rows=[hsp("AAAA","AACA",qlen=6,slen=6),
              hsp("AAAA","CAAA",qs=3,ss=3,qlen=6,slen=6),
              hsp("AAA","AAA",gene="b")]
        baseline=self.aggregate(rows)
        for permutation in itertools.permutations(rows):
            self.assertEqual(baseline,self.aggregate(permutation))

    def test_different_subjects_not_pooled(self):
        rows=[hsp("AAAA","AAAA",qlen=10,slen=10),
              hsp("AAAA","AAAA",ref="other",qs=7,ss=7,qlen=10,slen=10)]
        self.assertEqual(self.run_local(rows)[0],[])

    def test_legacy_family_index_columns_preserved(self):
        _,meta,_,legacy=self.run_local([hsp("AAAA","AAAA")])
        self.assertEqual(len(legacy.rstrip().split("\t")),13)
        self.assertEqual(len(meta["raw_alignment_sha256"]),64)
        self.assertEqual(len(meta["pair_evidence_sha256"]),64)

    def test_legacy_rows_fail_instead_of_guessing_identity(self):
        with self.assertRaisesRegex(ValueError,"expected 15 fields"):
            self.aggregate(["g\tr\t40\t1000\t1000\t1000\t1\t1000\t1\t1000\t0\t100\tx"])

    def test_inconsistent_sequence_length_rejected(self):
        with self.assertRaisesRegex(ValueError,"inconsistent sequence length"):
            self.aggregate([hsp("AAAA","AAAA"),hsp("AAAA","AAAA",qlen=5)])

    def test_nonfinite_rejected(self):
        f=hsp("AAAA","AAAA").split("\t"); f[2]="nan"
        with self.assertRaisesRegex(ValueError,"nonfinite"):
            self.aggregate(["\t".join(f)])

    def test_missing_provider_is_not_negative(self):
        rows,meta=run_local_blastp(Path("q"),Path("db"),Path("out"),1,30,.5,1e-5,
            run_tool=lambda *_: self.fail("must not execute"),find_binary=lambda _:None,
            normalise_subject=str)
        self.assertEqual(rows,[])
        self.assertTrue(meta["missing_provider_is_not_negative_evidence"])

    def test_real_wrapper_delegates_with_patchable_runtime(self):
        # Execute the actual wrapper AST without importing unrelated optional providers.
        source=Path(__file__).resolve().parents[1]/"src/candidate_qtl/route_b_prepare.py"
        tree=ast.parse(source.read_text(encoding="utf-8"))
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="_run_blastp")
        marker=object()
        namespace={"Path":Path,"run_local_blastp":lambda *a,**k:marker,
                   "run_tool":object(),"shutil":__import__("shutil"),"_normalise_subject":str}
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),"exec"),namespace)
        self.assertIs(namespace["_run_blastp"](Path("q"),Path("db"),Path("o"),1,30,.5,1e-5),marker)

if __name__ == "__main__":
    unittest.main()
