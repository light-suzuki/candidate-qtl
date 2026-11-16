"""A family index that exists but does not carry the family must not swallow the bait hits.

Issue #568 let a verified Route A bait reach its interval members by sequence
when no reference-genome family index existed at all. The far more common state
is an index that exists and simply has nothing for this family, reported as
``not_evaluated_family_unresolved_in_index``. Both mean the family question went
unanswered, so both must fall back to the bait's own blast evidence; admitting
only the first discarded every bait hit in ordinary runs, where an index is
present but incomplete.

The projection stays positive-only: no hit is not absence, the genome-wide
paralog guard stays unevaluated, and every threshold-passing copy is kept.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.family_candidate_bridge import family_candidate_bridge

HYPOTHESES = (
    "hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes\n"
    "H_LEG\t11S globulin seed storage proteins\tLEGB\tP05190\t\n"
)
# The index answered for other families and had nothing for this one.
STATUS = (
    "hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\n"
    "H_LEG\tnot_evaluated_family_unresolved_in_index\t\tPMID:123\n"
)
BAIT = (
    "hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id\n"
    "H_LEG\tPMID:123\tsource_verified\tP05190\n"
)
GENES = "gene_id\nPs1\nPs2\nPs3\nPs4\n"
BLAST = (
    "panel_id\tprovider_status\treview_status\tambiguity_status\ttop_hit_target_id\t"
    "top_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\t"
    "paralog_count\tdb_fingerprint\tinterval_proteins_sha256\tblastp_version\n"
    "P05190\tevaluable\tverified\tambiguous_paralog\tPs2\t61.0\t0.90\t0.9\t1e-95\t1\t"
    "db1\tsha1\t2.12.0+\n"
)


def _long_form(rows: list[tuple[str, str, str]]) -> str:
    header = ("panel_id\ttarget_id\tidentity\tquery_coverage\ttarget_coverage\tevalue\tbitscore\t"
              "review_status\tprovider_status\trank_contribution\teligibility\n")
    body = "".join(
        f"P05190\t{target}\t{identity}\t{coverage}\t0.9\t1e-95\t420\tverified\tevaluable\tnone\teligible\n"
        for target, identity, coverage in rows)
    return header + body


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class FamilyUnresolvedInIndexProjectionTests(unittest.TestCase):
    def _run(self, root: Path, *, long_form_rows):
        files = {
            "hypotheses": _write(root, "hypotheses.tsv", HYPOTHESES),
            "status": _write(root, "status.tsv", STATUS),
            "bait": _write(root, "bait.tsv", BAIT),
            "genes": _write(root, "genes.tsv", GENES),
            "blast": _write(root, "blast_evidence.tsv", BLAST),
            "eligible": _write(root, "eligible_interval_hits.tsv", _long_form(long_form_rows)),
        }
        out = root / "out"
        self.assertEqual(0, main([
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case_unresolved",
            "--route-a-hypotheses", str(files["hypotheses"]),
            "--route-a-family-status", str(files["status"]),
            "--route-a-bait-panel", str(files["bait"]),
            "--genes", str(files["genes"]),
            "--route-a-blast", str(files["blast"]),
            "--route-a-blast-eligible", str(files["eligible"]),
        ]))
        bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
        return bridge, manifest, out

    def test_bait_hits_survive_an_index_that_lacks_the_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge, manifest, out = self._run(
                root, long_form_rows=[("Ps2", "61.0", "0.90"), ("Ps3", "52.0", "0.74")])

            route_a = bridge[bridge["route"] == "a"]
            self.assertEqual(["Ps2", "Ps3"], sorted(route_a["gene_id"]))
            self.assertEqual({"blastp"}, set(route_a["match_type"]))
            # Sequence order is not copy evidence, so every eligible target is
            # kept for downstream copy review rather than ranked here.
            self.assertEqual({"none"}, set(route_a["rank_contribution"]))

            projection = manifest["route_a_sequence_family_projection"]
            self.assertEqual(2, projection["projection_count"])
            self.assertFalse(projection["paralog_guard_evaluated"])
            self.assertEqual({"not_evaluated_family_unresolved_in_index"},
                             {row["family_status"] for row in projection["projections"]})
            # An unanswered family question is not a falsified family.
            self.assertEqual(0, manifest["route_a_blast"]["falsified_family_projection_count"])

            union = root / "union"
            self.assertEqual(0, main([
                "supported-member-union", "--out-dir", str(union), "--case-id", "case_unresolved",
                "--route-a-bridge", str(out / "family_candidate_bridge.tsv"),
            ]))
            supported = pd.read_csv(union / "supported_member_union.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"Ps2", "Ps3"}, set(supported["candidate_gene_id"]))

    def test_a_hit_outside_the_interval_projects_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # The only threshold-passing target lies outside the interval
            # universe, so there is nothing to project and nothing may be
            # invented in its place. This matches, byte for byte, what the
            # already-admitted "no index existed" status does in the same
            # situation -- the fix extends that path, it does not soften it.
            bridge, manifest, _ = self._run(root, long_form_rows=[("OutsidePs", "61.0", "0.90")])

            self.assertTrue(bridge[bridge["route"] == "a"].empty)
            self.assertEqual(0, manifest["route_a_sequence_family_projection"]["projection_count"])
            self.assertEqual(0, manifest["route_a_blast"]["falsified_family_projection_count"])

    def test_identity_verified_unbound_blast_hit_is_sealed_review_only(self):
        """An incomplete source claim may reach Wet review, never support."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = _write(root, "hypotheses.tsv", HYPOTHESES)
            status = _write(root, "status.tsv", "hypothesis_id\n")
            canonical = _write(root, "canonical.tsv", "hypothesis_id\n")
            projection = _write(root, "projection.tsv", (
                "hypothesis_id\tbait_id\tprojection_eligibility\treview_status\tfunction_source_ids\n"
                "H_LEG\tP05190\treview_only\treview_required\tPMID:123\n"
            ))
            genes = _write(root, "genes.tsv", GENES)
            review_blast = _write(root, "review_blast.tsv", BLAST.replace("\tverified\t", "\treview_required\t"))
            result = family_candidate_bridge(
                str(root / "out"), route_a_hypotheses_arg=str(hypotheses),
                route_a_family_status_arg=str(status), route_a_bait_panel_arg=str(canonical),
                route_a_projection_panel_arg=str(projection), route_a_blast_arg=str(review_blast),
                genes_arg=str(genes), case_id="case_unresolved",
            )
            review = pd.read_csv(result["review_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["Ps2"], review["candidate_gene_id"].tolist())
            self.assertEqual({"review_required"}, set(review["support_status"]))
            self.assertEqual({"none"}, set(review["rank_contribution"]))
            # The canonical bridge remains empty; review-only material never
            # enters the artifact consumed by supported_member_union.
            bridge = pd.read_csv(result["bridge_tsv"], sep="\t", dtype=str).fillna("")
            self.assertTrue(bridge.empty)
            sealed = json.loads(Path(result["review_manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("route-review-candidates-1.0.0", sealed["contract"])
            self.assertEqual("route_a_review_candidate_freeze", sealed["stage"])


if __name__ == "__main__":
    unittest.main()
