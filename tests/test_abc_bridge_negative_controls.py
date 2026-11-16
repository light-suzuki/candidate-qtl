"""ABC/bridge integrated negative controls (issue #374 Criterion 11).

Two negative controls are pinned here, both inside the ABC -> bridge ->
family-to-copy -> integrated-narrowing environment:

1. **False convergence.** When Route A / B / C all point at the same generic
   family (a broad naming bucket / generic orthogroup), the convergence is
   never independent biological credit or causal proof. The bridge keeps the
   agreeing routes in one dependency cluster, every row carries
   ``rank_contribution=none`` and there is no score column; a row's mapping
   tier never rises because another route agrees. family-to-copy records the
   convergence as ``route_rescue_state`` but never promotes a copy on it:
   with no discriminating axis the copies stay
   ``family_supported_copy_unresolved`` and the shortlist stays empty.

2. **Generic-family BLASTP.** A strong BLASTP hit against a *generic* family is
   similarity evidence only. Route A ``surviving_family_present_in_interval``
   and Route B ``best_hit_homology`` are mapped to ``match_type=blastp`` and
   capped at the sequence-match tiers (tier3/tier4), never
   ``tier1_source_backed_unique`` / ``tier2_source_backed_ambiguous``; the
   identity / coverage / E-value are transcribed but never promote. In
   family-to-copy a generic BLAST best-hit (``family_level_only``) is never a
   one-to-one ortholog and never discriminates, and a copy whose own logic
   would make it ``preferred_copy`` but whose bridge tier is ``tier5_name_only``
   is gate-blocked (``bridge_promotion_gate=blocked_tier5``) and never rises to
   single Top (N1) in integrated narrowing.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.family_candidate_bridge import BRIDGE_COLUMNS
from candidate_qtl.schemas import FAMILY_CANDIDATE_BRIDGE_SCHEMA, validate

# ---------------------------------------------------------------------------
# Shared synthetic fixtures (no private sequences / genomes).
# ---------------------------------------------------------------------------

# False convergence: Route A / B / C all point at the same generic family
# "LRR_kinase" / generic orthogroup "OG_GEN", all citing the same root source.
GEN_HYPOTHESES = "\n".join([
    "hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes",
    "H_GEN\tLRR_kinase\tREF_GEN\tACC_GEN\tG1;G2",
]) + "\n"

GEN_BAIT_PANEL = "\n".join([
    "hypothesis_id\tfunction_source_ids\tsource_provider",
    "H_GEN\tROOT_GEN\tuniprot_2026",
]) + "\n"

GEN_FAMILY_STATUS = "\n".join([
    "hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id",
    "H_GEN\tsurviving_family_present_in_interval\tG1;G2\tROOT_GEN\tREF_GEN\tACC_GEN",
]) + "\n"

# A generic-family best hit: high identity / full coverage, yet only a best
# hit, never a reciprocal / table-backed orthology claim.
GEN_ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\trelationship\torthogroup\tpercent_identity\tquery_coverage\ttarget_coverage\tevalue\tprovider\tparalog_count\tsource_id",
    "G1\tREF1\tbest_hit_homology\tOG_GEN\t99.0\t1.0\t1.0\t1e-100\tlocal_blastp\t\tROOT_GEN",
    "G2\tREF2\tbest_hit_homology\tOG_GEN\t98.0\t1.0\t1.0\t1e-90\tlocal_blastp\t\tROOT_GEN",
]) + "\n"

# Route C agrees on the same generic family via a frozen-retrieval orthogroup
# match (capped at sequence-match by design for route == c).
GEN_MATCH = "\n".join([
    "reference_gene_id\trank\tmatched_interval_gene\tmatch_basis\tfamily_id\torthogroup_id",
    "R_GEN\t1\tG1;G2\tsame_orthogroup\tLRR_kinase\tOG_GEN",
]) + "\n"

GEN_MERGED = "\n".join([
    "reference_gene_id\trank\tdatabase_version\tmodel_version",
    "R_GEN\t1\tdb2026\tesm3",
]) + "\n"

GEN_CLUSTERS = "\n".join([
    "dependency_cluster\tdependency_status\troot_source_id\tsource_ids",
    "dep_GEN\troot_source_known\tROOT_GEN\tROOT_GEN",
]) + "\n"

# Contrast fixture: the same strong blast hit against a *curated* family
# (concrete_interval_gene_matched) vs a *generic* family
# (surviving_family_present_in_interval). Only curation reaches tier1.
BLAST_HYPOTHESES = "\n".join([
    "hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes",
    "H_CUR\tCLV\tREF_CUR\tACC_CUR\tG1",
    "H_BL\tLRR_kinase\tREF_BL\tACC_BL\tG2",
]) + "\n"

BLAST_BAIT_PANEL = "\n".join([
    "hypothesis_id\tfunction_source_ids\tsource_provider",
    "H_CUR\tROOT_CUR\tuniprot_2026",
    "H_BL\tROOT_BL\tuniprot_2026",
]) + "\n"

BLAST_FAMILY_STATUS = "\n".join([
    "hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id",
    "H_CUR\tconcrete_interval_gene_matched\tG1\tROOT_CUR\tREF_CUR\tACC_CUR",
    "H_BL\tsurviving_family_present_in_interval\tG2\tROOT_BL\tREF_BL\tACC_BL",
]) + "\n"

# blast_evidence.tsv: identical 99% / full-coverage unique hits for both baits.
BLAST_EVIDENCE = "\n".join([
    "panel_id\tprovider_status\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\tparalog_count\tambiguity_status\treview_status\trank_contribution",
    "ACC_CUR\tevaluable\tG1\t99.0\t1.0\t1.0\t0.0\t0\tunique\tverified\tnone",
    "ACC_BL\tevaluable\tG2\t99.0\t1.0\t1.0\t0.0\t0\tunique\tverified\tnone",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _bridge_row(
    gene_id: str, family_id: str, route: str, match_type: str,
    mapping_tier: str, ambiguity: str, dependency: str,
) -> str:
    """One canonical bridge row (positions follow BRIDGE_COLUMNS)."""
    values = ["demo", gene_id, family_id, route, "", "", match_type, "", "", "", "", "1",
              ambiguity, "complete", "", "ff", "ff", "", dependency, mapping_tier, "none", ""]
    return "\t".join(values)


def _write_bridge(root: Path, rows: list[str]) -> Path:
    path = root / "bridge.tsv"
    path.write_text("\t".join(BRIDGE_COLUMNS) + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


class FalseConvergenceBridgeTests(unittest.TestCase):
    """Route A / B / C converging on one generic family is clustered and never
    scored, never fused, and never raises any row's mapping tier."""

    def _fixture(self, root: Path) -> dict[str, Path]:
        return {
            "hypotheses": _write(root, "hypotheses.tsv", GEN_HYPOTHESES),
            "bait_panel": _write(root, "bait_panel.tsv", GEN_BAIT_PANEL),
            "family_status": _write(root, "family_status.tsv", GEN_FAMILY_STATUS),
            "orthology": _write(root, "orthology.tsv", GEN_ORTHOLOGY),
            "match": _write(root, "match.tsv", GEN_MATCH),
            "merged": _write(root, "merged.tsv", GEN_MERGED),
            "clusters": _write(root, "clusters.tsv", GEN_CLUSTERS),
        }

    def _run_bridge(self, root: Path) -> pd.DataFrame:
        fixture = self._fixture(root)
        out = root / "bridge_out"
        self.assertEqual(0, main([
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "demo",
            "--route-a-hypotheses", str(fixture["hypotheses"]),
            "--route-a-family-status", str(fixture["family_status"]),
            "--route-a-bait-panel", str(fixture["bait_panel"]),
            "--route-b-orthology", str(fixture["orthology"]),
            "--route-c-retrieval", str(fixture["merged"]),
            "--route-c-match", str(fixture["match"]),
            "--dependency-clusters", str(fixture["clusters"]),
        ]))
        return _read_tsv(out / "family_candidate_bridge.tsv")

    def test_abc_convergence_on_a_generic_family_is_never_a_score_or_fused_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_bridge(Path(tmp))
            # Reporting-only: no score column, rank_contribution none everywhere.
            self.assertEqual(list(BRIDGE_COLUMNS), list(frame.columns))
            self.assertFalse({"dry_score", "score", "dry_group", "rank"} & set(frame.columns))
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            # Agreement produces one parallel row per route for G1, never a sum.
            g1 = frame[frame["gene_id"].eq("G1")]
            self.assertEqual({"a", "b", "c"}, set(g1["route"]))
            self.assertEqual(3, len(g1))
            # Generic-family convergence is sequence-match / naming only: no
            # route claims curated / rbh orthology for the generic family.
            self.assertTrue(g1["match_type"].isin({"blastp", "orthogroup"}).all())

    def test_generic_family_convergence_never_raises_any_row_above_sequence_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_bridge(Path(tmp))
            rows = frame[frame["gene_id"].isin({"G1", "G2"})]
            self.assertFalse(rows.empty)
            # No row reaches the source-backed tiers, however strongly the
            # routes agree: the generic family never earns tier1/tier2.
            self.assertTrue(rows["mapping_tier"].isin({
                "tier3_sequence_match_unique", "tier4_sequence_match_ambiguous",
            }).all())
            # Route A surviving family -> blastp + 2 members -> ambiguous tier4.
            route_a = rows[(rows["route"].eq("a"))]
            self.assertTrue(route_a["match_type"].eq("blastp").all())
            self.assertTrue(route_a["mapping_tier"].eq("tier4_sequence_match_ambiguous").all())
            # Route B best_hit_homology -> blastp (not rbh) -> tier3, even at
            # 99% identity / full coverage.
            route_b = rows[(rows["route"].eq("b"))]
            self.assertTrue(route_b["match_type"].eq("blastp").all())
            self.assertTrue(route_b["mapping_tier"].eq("tier3_sequence_match_unique").all())
            self.assertEqual("99.0", route_b[route_b["gene_id"].eq("G1")].iloc[0]["identity"])
            # Route C orthogroup on the generic family -> capped at tier4.
            route_c = rows[(rows["route"].eq("c"))]
            self.assertTrue(route_c["match_type"].eq("orthogroup").all())
            self.assertTrue(route_c["mapping_tier"].eq("tier4_sequence_match_ambiguous").all())

    def test_converged_routes_share_one_dependency_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_bridge(Path(tmp))
            g1 = frame[frame["gene_id"].eq("G1")]
            # A and B cite the same root source: they share one cluster, so the
            # agreement is never N independent sources.
            a_row = g1[g1["route"].eq("a")].iloc[0]
            b_row = g1[g1["route"].eq("b")].iloc[0]
            self.assertEqual("dep_GEN", a_row["dependency_cluster_id"])
            self.assertEqual("dep_GEN", b_row["dependency_cluster_id"])
            # Route C roots on its own declared database_version.
            c_row = g1[g1["route"].eq("c")].iloc[0]
            self.assertEqual("db2026", c_row["source_id"])
            self.assertNotEqual("dep_GEN", c_row["dependency_cluster_id"])
            manifest = json.loads((Path(tmp) / "bridge_out" / "family_candidate_bridge_manifest.json")
                                  .read_text(encoding="utf-8"))
            entries = {entry["dependency_cluster_id"]: entry for entry in manifest["dependency_clusters"]}
            self.assertIn("dep_GEN", entries)
            self.assertEqual("ROOT_GEN", entries["dep_GEN"]["root_source_id"])
            self.assertEqual("a;b", entries["dep_GEN"]["routes"])
            # Schema stays valid end-to-end.
            validate(frame, FAMILY_CANDIDATE_BRIDGE_SCHEMA)


class GenericFamilyBlastpBridgeTests(unittest.TestCase):
    """A strong generic-family BLASTP hit is similarity evidence only: it is
    transcribed but never promotes a row to a source-backed tier."""

    def _run_route_a(self, root: Path) -> pd.DataFrame:
        out = root / "bridge_out"
        self.assertEqual(0, main([
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "demo",
            "--route-a-hypotheses", str(_write(root, "hypotheses.tsv", BLAST_HYPOTHESES)),
            "--route-a-family-status", str(_write(root, "family_status.tsv", BLAST_FAMILY_STATUS)),
            "--route-a-bait-panel", str(_write(root, "bait_panel.tsv", BLAST_BAIT_PANEL)),
            "--route-a-blast", str(_write(root, "blast_evidence.tsv", BLAST_EVIDENCE)),
        ]))
        return _read_tsv(out / "family_candidate_bridge.tsv")

    def test_route_a_surviving_generic_family_blastp_is_sequence_match_not_source_backed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_route_a(Path(tmp))
            # Same blast evidence, different family status: curation reaches
            # tier1, the generic (surviving, non-curated) family stays at the
            # sequence-match tier3. The blast hit itself never promotes.
            curated = frame[(frame["gene_id"].eq("G1")) & frame["route"].eq("a")].iloc[0]
            self.assertEqual("curated", curated["match_type"])
            self.assertEqual("tier1_source_backed_unique", curated["mapping_tier"])
            generic = frame[(frame["gene_id"].eq("G2")) & frame["route"].eq("a")].iloc[0]
            self.assertEqual("blastp", generic["match_type"])
            self.assertEqual("tier3_sequence_match_unique", generic["mapping_tier"])
            # The identity / coverage are recorded (evidence kept) yet capped.
            self.assertEqual("99.0", generic["identity"])
            self.assertEqual("1.0", generic["query_coverage"])
            self.assertEqual("resolved_single_copy", generic["ambiguity_status"])
            self.assertTrue(frame["rank_contribution"].eq("none").all())

    def test_route_b_best_hit_homology_is_never_source_backed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "bridge_out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "demo",
                "--route-b-orthology", str(_write(root, "orthology.tsv", GEN_ORTHOLOGY)),
            ]))
            frame = _read_tsv(out / "family_candidate_bridge.tsv")
            row = frame[frame["gene_id"].eq("G1")].iloc[0]
            # A generic best hit is not reciprocal orthology: blastp, tier3,
            # never rbh / orthogroup and never tier1/tier2.
            self.assertEqual("blastp", row["match_type"])
            self.assertEqual("tier3_sequence_match_unique", row["mapping_tier"])
            self.assertEqual("99.0", row["identity"])
            self.assertNotEqual("rbh", row["match_type"])
            self.assertNotEqual("tier1_source_backed_unique", row["mapping_tier"])


class AbcBridgeCopyNegativeControlTests(unittest.TestCase):
    """At the copy / narrowing level: convergence and generic-family BLASTP
    never promote a copy, and a copy the bridge caps at name-only is
    gate-blocked and never rises to single Top."""

    def _copy_fixture(self, root: Path, members: str) -> dict[str, Path]:
        _write(root, "families.tsv", "\n".join([
            "hypothesis_id\tgene_family\tmember_interval_genes",
            f"hyp_GEN\tLRR_kinase\t{members}",
        ]) + "\n")
        atlas_rows = ["gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene"]
        atlas_rows += [f"{gene}\tchr1\t{1000 + i * 200000}\t{3000 + i * 200000}\tLRR\tLRR-Pkinase\tREF_{gene}"
                       for i, gene in enumerate(members.split(";"))]
        _write(root, "atlas.tsv", "\n".join(atlas_rows) + "\n")
        orthology_rows = ["gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup"]
        orthology_rows += [f"{gene}\tfamily_level_only\torthofinder\tOF_chain\tREF_{gene}\tOG_GEN"
                           for gene in members.split(";")]
        _write(root, "orthology.tsv", "\n".join(orthology_rows) + "\n")
        return {
            "families": root / "families.tsv", "atlas": root / "atlas.tsv",
            "orthology": root / "orthology.tsv",
        }

    def _run_family_to_copy(self, root: Path, bridge: Path) -> Path:
        fixture = self._copy_fixture(root, "G1;G2")
        out = root / "copy"
        self.assertEqual(0, main(["family-to-copy", "--families", str(fixture["families"]),
                                  "--out-dir", str(out), "--case-id", "demo",
                                  "--atlas", str(fixture["atlas"]),
                                  "--orthology", str(fixture["orthology"]),
                                  "--bridge", str(bridge)]))
        return out

    def test_abc_false_convergence_does_not_promote_a_copy(self) -> None:
        """The false-convergence bridge (route a/b/c all agreeing) is chained
        into family-to-copy: the convergence is recorded as multi-route support
        yet no copy is promoted -- the generic BLAST best-hit members have no
        discriminating axis and stay unresolved."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = {
                "hypotheses": _write(root, "hypotheses.tsv", GEN_HYPOTHESES),
                "bait_panel": _write(root, "bait_panel.tsv", GEN_BAIT_PANEL),
                "family_status": _write(root, "family_status.tsv", GEN_FAMILY_STATUS),
                "orthology": _write(root, "orthology.tsv", GEN_ORTHOLOGY),
                "match": _write(root, "match.tsv", GEN_MATCH),
                "merged": _write(root, "merged.tsv", GEN_MERGED),
                "clusters": _write(root, "clusters.tsv", GEN_CLUSTERS),
            }
            bridge_out = root / "bridge_out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(bridge_out), "--case-id", "demo",
                "--route-a-hypotheses", str(fixture["hypotheses"]),
                "--route-a-family-status", str(fixture["family_status"]),
                "--route-a-bait-panel", str(fixture["bait_panel"]),
                "--route-b-orthology", str(fixture["orthology"]),
                "--route-c-retrieval", str(fixture["merged"]),
                "--route-c-match", str(fixture["match"]),
                "--dependency-clusters", str(fixture["clusters"]),
            ]))
            out = self._run_family_to_copy(root, bridge_out / "family_candidate_bridge.tsv")
            discrimination = _read_tsv(out / "family_copy_discrimination.tsv").set_index("candidate_gene_id")
            # Multi-route convergence is recorded as support, not credit.
            self.assertEqual({"multi_route_supported"}, set(discrimination["route_rescue_state"]))
            # ... but it never promotes a copy: no discriminating axis exists.
            self.assertTrue(discrimination["disambiguation_status"].eq("family_supported_copy_unresolved").all())
            self.assertTrue(discrimination["rank_contribution"].eq("none").all())
            shortlist = _read_tsv(out / "family_copy_shortlist.tsv")
            self.assertEqual(0, len(shortlist))
            # Not a copy-selection loss either: nothing was ever narrowed.
            self.assertTrue(discrimination["copy_selection_loss"].eq("false").all())

    def test_generic_bucket_copy_is_gate_blocked_and_never_rises_to_single_top(self) -> None:
        """A sole copy whose own homology makes it internally preferred_copy but
        whose bridge tier is tier5_name_only (a generic naming bucket) is
        gate-blocked from the shortlist (blocked_tier5) and never reaches N1 in
        integrated narrowing -- its disambiguation status is preserved in the
        artifact but it is never promoted."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._copy_fixture(root, "G_S")
            bridge = _write_bridge(root, [_bridge_row(
                "G_S", "LRR_kinase", "b", "name_only",
                "tier5_name_only", "paralog_ambiguous", "dep_GEN",
            )])
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(fixture["families"]),
                                      "--out-dir", str(out), "--case-id", "demo",
                                      "--atlas", str(fixture["atlas"]),
                                      "--orthology", str(fixture["orthology"]),
                                      "--bridge", str(bridge)]))
            discrimination = _read_tsv(out / "family_copy_discrimination.tsv").set_index("candidate_gene_id")
            # The internal disambiguation status is preserved, untouched...
            self.assertEqual("preferred_copy", discrimination.loc["G_S", "disambiguation_status"])
            # ... but the generic bucket tier blocks promotion.
            self.assertEqual("tier5_name_only", discrimination.loc["G_S", "bridge_mapping_tier"])
            self.assertEqual("blocked_tier5", discrimination.loc["G_S", "bridge_promotion_gate"])
            self.assertEqual("covered", discrimination.loc["G_S", "bridge_missingness"])
            shortlist = _read_tsv(out / "family_copy_shortlist.tsv")
            self.assertEqual(0, len(shortlist))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["bridge"]["blocked_count"])

            # Integrated narrowing over those artifacts: the gate-blocked
            # preferred copy never claims the single Top N1 tier.
            _write(root, "ranked.tsv", "\n".join([
                "gene_id\tchrom\tstart\tend\tdry_score\tdry_rank\tdry_group\ttop_eligibility\tscreening_rank",
                "G_S\tchr1\t1000\t3000\t5.0\t1\ttop\tstrong_genetic_evidence\t1",
            ]) + "\n")
            _write(root, "recovery.tsv", "\n".join([
                "hypothesis_id\tgene_family\tmember_interval_genes",
                "hyp_GEN\tLRR_kinase\tG_S",
            ]) + "\n")
            narrow_out = root / "narrow"
            self.assertEqual(0, main(["integrated-narrowing", "--ranked", str(root / "ranked.tsv"),
                                      "--out-dir", str(narrow_out), "--case-id", "demo",
                                      "--family-recovery", str(root / "recovery.tsv"),
                                      "--family-copy", str(out / "family_copy_discrimination.tsv"),
                                      "--family-copy-shortlist", str(out / "family_copy_shortlist.tsv")]))
            narrowed = _read_tsv(narrow_out / "integrated_narrowing.tsv").set_index("gene_id")
            self.assertNotEqual("N1", narrowed.loc["G_S", "narrowing_tier"])
            self.assertEqual("N3", narrowed.loc["G_S", "narrowing_tier"])
            self.assertEqual("base_top_reserve_rescue", narrowed.loc["G_S", "narrowing_reason"])
            self.assertEqual("preferred_copy", narrowed.loc["G_S", "copy_state"])
            self.assertEqual("", narrowed.loc["G_S", "first_wet_priority"])
            self.assertEqual("none", narrowed.loc["G_S", "rank_contribution"])


if __name__ == "__main__":
    unittest.main()
