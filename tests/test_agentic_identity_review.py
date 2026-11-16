"""The review producer must reach the canonical consumer without inventing answers."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO

from candidate_qtl.core import sha256
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.post_freeze_identity import load_identity_assessment
from candidate_qtl.relevance_tiers import annotate_relevance, NEW_TIERS
from scripts.agentic_identity_review import prepare, compile_review, _claims
from tests import test_post_freeze_identity as identity_fixtures
from tests import test_route_a_review_relevance as route_a


class IdentityReviewTests(unittest.TestCase):
    def fixture(self, root):
        (root / "identity").mkdir()
        (root / "route").mkdir()
        actual = {r.id: str(r.seq) for r in SeqIO.parse(str(route_a.FIXTURE / "interval_proteins.fa"), "fasta")}
        fixture = identity_fixtures.PostFreezeIdentityTests().fixture(root / "identity", sequence=route_a.S,
                                                  targets={"p1": actual["g1"], "p2": actual["g2"]})
        packet = root / "packet"
        packet.mkdir()
        proteins = packet / "candidate_proteins.fa"
        proteins.write_text(f">g1\n{actual['g1']}\n>g2\n{actual['g2']}\n", encoding="utf-8")
        _, members, hypotheses, _ = route_a.RouteAReviewRelevanceTests()._run(
            root / "route", {"source_ids": "PMID:1234567", PAPER_FIELD: normalize_paper_reports([fixture["paper"]])})
        atlas = root / "route" / "genes.tsv"
        stage = {"case_id": "route-a-review", "inputs": {"packet_b": str(packet)},
                 "b0": {"interval_gene_atlas": str(atlas)},
                 "a_sources": {"mechanisms_with_verified_sources": str(hypotheses)},
                 "files": {str(p.relative_to(root)): sha256(p) for p in root.rglob("*") if p.is_file()}}
        anchor = root / "stage.json"
        anchor.write_text(json.dumps(stage), encoding="utf-8")
        args = [str(anchor), sha256(anchor), fixture["bridge"]["manifest"],
                members["route_a_review_candidates"], str(Path(members["route_a_review_candidates"]).with_name("route_a_review_candidates_manifest.json")),
                str(root / "identity" / "paper.json"), str(root / "prepared")]
        result = prepare(*args, provider="fixture", model="explicit-fixture-review", run_id="run-1")
        return args, result, atlas, proteins

    def test_cli_compile_preserves_original_answers_and_derives_bound_t1(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args, result, atlas, proteins = self.fixture(root)
            self.assertEqual(1, result["review_count"])
            template = json.loads((root / "prepared" / "review_template.json").read_text())
            self.assertEqual("not_evaluated", template[0]["decision"])
            self.assertEqual({"unknown"}, set(template[0]["review_axes"].values()))
            template[0].update(decision="same_gene", reason="Synthetic source, official CDS, experimental taxon and unique complete model reviewed.")
            template[0]["review_axes"] = dict.fromkeys(template[0]["review_axes"], "yes")
            review = root / "actual-review.json"
            review.write_text(json.dumps(template, indent=2), encoding="utf-8")
            manifest = root / "prepared" / "prepare_manifest.json"
            completed = subprocess.run([sys.executable, "scripts/agentic_identity_review.py", "compile",
                                        "--manifest", str(manifest), "--sha256", sha256(manifest),
                                        "--review", str(review), "--out", str(root / "compiled")], capture_output=True, text=True)
            self.assertEqual(0, completed.returncode, completed.stderr)
            response = json.loads(completed.stdout)
            self.assertEqual(1, response["resolved_evidence_count"])
            claims = _claims(args[3], args[4], "route-a-review", atlas)
            evidence, _ = load_identity_assessment(response["assessment"], case_id="route-a-review",
                                                   atlas_arg=str(atlas), canonical_proteins_arg=str(proteins), claims=claims)
            self.assertEqual(NEW_TIERS[0], annotate_relevance([{"candidate_gene_id": "g1"}], claims + evidence)[0]["relevance_tier"])
            self.assertEqual("unknown", json.loads(evidence[0]["original_paper_report"]["paper_checklist"])["same_gene_as_focal"])

    def test_template_is_not_a_completed_review_and_endpoint_changes_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            manifest = root / "prepared" / "prepare_manifest.json"
            review = root / "prepared" / "review_template.json"
            with self.assertRaisesRegex(ValueError, "decision/reason"):
                compile_review(str(manifest), sha256(manifest), str(review), str(root / "blank"))
            data = json.loads(review.read_text())
            data[0].update(reason="unknown", candidate_gene_id="g2")
            changed = root / "changed.json"
            changed.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "changed frozen"):
                compile_review(str(manifest), sha256(manifest), str(changed), str(root / "bad"))

    def test_unknown_review_keeps_zero_resolved_and_material_tampering_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            manifest = root / "prepared" / "prepare_manifest.json"
            data = json.loads((root / "prepared" / "review_template.json").read_text())
            data[0].update(decision="unknown", reason="The reviewer could not establish paper-to-CDS identity.")
            review = root / "unknown.json"
            review.write_text(json.dumps(data))
            result = compile_review(str(manifest), sha256(manifest), str(review), str(root / "unknown"))
            self.assertEqual(0, result["resolved_evidence_count"])
            source = next((root / "prepared" / "papers").glob("*.json"))
            source.write_text(source.read_text() + " ")
            with self.assertRaisesRegex(ValueError, "artifact mismatch"):
                compile_review(str(manifest), sha256(manifest), str(review), str(root / "tampered"))


if __name__ == "__main__":
    unittest.main()
