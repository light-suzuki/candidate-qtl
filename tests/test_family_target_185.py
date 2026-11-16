"""Issue #185: target-relative copy discrimination.

The old heuristic treated the family-internal modal subfamily / modal domain
architecture as the implicit answer, which biases plant copy discrimination
toward the majority copy even when the true functional ortholog is a minority
copy (duplication / subfunctionalization / lineage-specific expansion). These
tests verify the target-relative semantics: with a verified target
specification the copies are scored by how closely they match the target; the
modal comparisons are descriptive statistics and outlier QC only; without a
target no copy is promoted on them.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main

TARGET_HEADER = ("target_id\thypothesis_id\treference_gene_id\treference_species\texpected_family_id"
                 "\texpected_orthogroup\texpected_subfamily\texpected_domain_architecture"
                 "\texpected_synteny_anchor\texpected_expression_context\tsource_ids\tverification_status")


def _write(root: Path, name: str, header: str, records: str) -> str:
    (root / name).write_text(header + "\n" + records + "\n", encoding="utf-8")
    return name


def _families(root: Path, records: str) -> str:
    return _write(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes", records)


def _atlas(root: Path, records: str) -> str:
    return _write(root, "atlas.tsv",
                  "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene", records)


def _orthology(root: Path, records: str) -> str:
    return _write(root, "orthology.tsv",
                  "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup", records)


def _subfamily(root: Path, records: str) -> str:
    return _write(root, "subfamily.tsv", "gene_id\tsubfamily", records)


def _position(root: Path, records: str) -> str:
    return _write(root, "position.tsv", "gene_id\tposition_class", records)


def _target(root: Path, row: str) -> str:
    return _write(root, "target.tsv", TARGET_HEADER, row)


def _run(root: Path, *, orthology: str | None = None, subfamily_map: str | None = None,
         target: str | None = None, position: str | None = None, out_name: str = "copy") -> Path:
    out = root / out_name
    argv = ["family-to-copy", "--families", str(root / "families.tsv"),
            "--atlas", str(root / "atlas.tsv"), "--out-dir", str(out)]
    for flag, value in (("--orthology", orthology), ("--subfamily-map", subfamily_map),
                        ("--family-target", target), ("--position", position)):
        if value:
            argv += [flag, str(root / value)]
    assert main(argv) == 0, argv
    return out


class TargetRelativeMinorityOrthologTests(unittest.TestCase):
    """Issue test A: the minority copy that matches the target is the preferred one."""

    def test_minority_true_ortholog_wins_over_the_family_modal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Family modal is Y (copies B and C), but the verified target is X:
            # copy A matches the target on both axes and must be preferred.
            _families(root, "hyp_min\tFAMX\tA;B;C\n")
            _atlas(root, "A\tchr1\t1000\t3000\tFAMX\tX\tREF_A\n"
                         "B\tchr1\t5000\t7000\tFAMX\tY\tREF_B\n"
                         "C\tchr1\t9000\t11000\tFAMX\tY\tREF_C\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_A\tOG-X\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_B\tOG-X\n"
                             "C\tone_to_many\torthofinder\tOF\tREF_C\tOG-X\n")
            _subfamily(root, "A\tX\nB\tY\nC\tY\n")
            _target(root, "tgt_min\thyp_min\tREF_T\tPisum sativum\tFAMX\tOG-X\tX\tX\t\t\tpanel\tverified")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", target="target.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("target_match", frame.loc["A", "subfamily_vs_target"])
            self.assertEqual("target_match", frame.loc["A", "domain_vs_target"])
            self.assertIn("subfamily_match", frame.loc["A", "leading_axes"])
            self.assertIn("domain_vs_target", frame.loc["A", "leading_axes"])
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            # The majority of B/C does not weaken the target-matched minority.
            self.assertEqual("target_mismatch", frame.loc["B", "subfamily_vs_target"])
            self.assertEqual("target_mismatch", frame.loc["C", "domain_vs_target"])
            self.assertEqual("copy_weakened", frame.loc["B", "disambiguation_status"])
            self.assertEqual("copy_weakened", frame.loc["C", "disambiguation_status"])
            # The modal stays descriptive: the family modal is Y even though the
            # preferred copy is the X minority.
            self.assertEqual("Y", frame.loc["A", "subfamily_family_modal"])
            self.assertEqual("architecture_mismatch", frame.loc["A", "domain_vs_family_modal"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("target_relative", manifest["selection_basis"]["FAMX"])
            self.assertEqual(3, manifest["selection_basis_counts"]["target_relative"])


class ModalArchitectureAsQcOnlyTests(unittest.TestCase):
    """Issue test B: the modal comparison is outlier QC, not selection."""

    def test_truncated_copy_is_flagged_as_outlier_and_as_target_mismatch_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Two intact copies and one truncated copy (LRR only). The target
            # architecture is the intact one; the modal comparison and the
            # target comparison must be recorded as separate columns and the
            # modal mismatch must not promote or weaken any copy by itself.
            _families(root, "hyp_qc\tQCFAM\tA;B;C\n")
            _atlas(root, "A\tchr1\t1000\t3000\tQCFAM\tLRR-Pkinase\tREF_A\n"
                         "B\tchr1\t5000\t7000\tQCFAM\tLRR-Pkinase\tREF_B\n"
                         "C\tchr1\t9000\t11000\tQCFAM\tLRR\tREF_C\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_A\tOG-QC\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_B\tOG-QC\n"
                             "C\tone_to_many\torthofinder\tOF\tREF_C\tOG-QC\n")
            _subfamily(root, "A\tSX\nB\tSX\nC\tSX\n")
            _target(root, "tgt_qc\thyp_qc\tREF_T\tPisum sativum\tQCFAM\tOG-QC\tSX\tLRR-Pkinase\t\t\tpanel\tverified")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", target="target.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("target_match", frame.loc["A", "domain_vs_target"])
            self.assertEqual("architecture_match", frame.loc["A", "domain_vs_family_modal"])
            self.assertEqual("architecture_match", frame.loc["A", "domain_architecture_match"])
            self.assertEqual("target_mismatch", frame.loc["C", "domain_vs_target"])
            self.assertEqual("architecture_mismatch", frame.loc["C", "domain_vs_family_modal"])
            self.assertEqual("architecture_mismatch", frame.loc["C", "domain_architecture_match"])
            # The modal axis is outlier QC and never a selection axis.
            axes = set(frame.loc["A", "discriminating_axes"].split(";"))
            self.assertIn("domain_vs_target", axes)
            self.assertNotIn("domain_architecture_match", axes)
            # Both intact copies tie on the target; the outlier is weakened.
            self.assertEqual("co_preferred_tie", frame.loc["A", "disambiguation_status"])
            self.assertEqual("co_preferred_tie", frame.loc["B", "disambiguation_status"])
            self.assertEqual("copy_weakened", frame.loc["C", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A", "B"}, set(shortlist["candidate_gene_id"]))


class NoTargetTests(unittest.TestCase):
    """Issue test C: without a target the modal alone never picks a copy."""

    def test_no_target_means_no_preferred_copy_from_modal_alone(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _families(root, "hyp_not\tNOTFAM\tA;B\n")
            _atlas(root, "A\tchr1\t1000\t3000\tNOTFAM\tX\tREF_A\n"
                         "B\tchr1\t900000\t902000\tNOTFAM\tY\tREF_B\n")
            _orthology(root, "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG-N\n"
                             "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG-N\n")
            _subfamily(root, "A\tX\nB\tY\n")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertTrue(frame["subfamily_vs_target"].eq("subfamily_unresolved").all())
            self.assertTrue(frame["domain_vs_target"].eq("target_not_evaluated").all())
            self.assertTrue(frame["disambiguation_status"].eq("family_supported_copy_unresolved").all())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("copy_target_missing", reasons)
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("unresolved_no_target", manifest["selection_basis"]["NOTFAM"])


class DuplicatedTargetTests(unittest.TestCase):
    """Issue test D: a target duplicated in the focal genome keeps its tie."""

    def test_two_copies_matching_the_target_never_produce_a_single_winner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _families(root, "hyp_dup\tDUPFAM\tA;B\n")
            _atlas(root, "A\tchr1\t1000\t3000\tDUPFAM\tDA-arch\tREF_T\n"
                         "B\tchr1\t900000\t902000\tDUPFAM\tDA-arch\tREF_T\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_T\tOG-D\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_T\tOG-D\n")
            _subfamily(root, "A\tSX\nB\tSX\n")
            _position(root, "A\tbroad_interval_only\nB\tbroad_interval_only\n")
            _target(root, "tgt_dup\thyp_dup\tREF_T\tPisum sativum\tDUPFAM\tOG-D\tSX\tDA-arch\t\t\tpanel\tverified")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv",
                       target="target.tsv", position="position.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertTrue(frame["subfamily_vs_target"].eq("target_match").all())
            self.assertTrue(frame["domain_vs_target"].eq("target_match").all())
            # Both copies match the target but nothing else separates them: the
            # one-to-many duplication is preserved, no single winner is made.
            self.assertTrue(frame["disambiguation_status"].eq("paralog_ambiguous").all())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("copy_paralog_ambiguous", reasons)

    def test_duplicated_target_blocks_a_sole_target_match_winner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # B shares the target subfamily with A, so the reference target is
            # duplicated in the focal genome: A's unique domain match is not
            # enough to make a sole preferred copy without locus evidence.
            _families(root, "hyp_dg\tDGFAM\tA;B;C\n")
            _atlas(root, "A\tchr1\t1000\t3000\tDGFAM\tDX\tREF_A\n"
                         "B\tchr1\t5000\t7000\tDGFAM\tDY\tREF_B\n"
                         "C\tchr1\t9000\t11000\tDGFAM\tDZ\tREF_C\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_A\tOG-DG\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_B\tOG-DG\n"
                             "C\tone_to_many\torthofinder\tOF\tREF_C\tOG-DG\n")
            _subfamily(root, "A\tDX\nB\tDX\nC\tDZ\n")
            _target(root, "tgt_dg\thyp_dg\tREF_T\tPisum sativum\tDGFAM\tOG-DG\tDX\tDX\t\t\tpanel\tverified")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", target="target.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("target_match", frame.loc["A", "subfamily_vs_target"])
            self.assertEqual("target_match", frame.loc["A", "domain_vs_target"])
            self.assertEqual("paralog_ambiguous", frame.loc["A", "disambiguation_status"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("copy_target_duplicated", reasons)


class PositionOverridesMajorityTests(unittest.TestCase):
    """Issue test E: locus-specific position evidence discriminates without a target."""

    def test_recombinant_core_minority_copy_is_preferred_without_a_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _families(root, "hyp_pos\tPOSFAM\tA;B;C\n")
            _atlas(root, "A\tchr1\t1000\t3000\tPOSFAM\tP1\tREF_A\n"
                         "B\tchr1\t5000\t7000\tPOSFAM\tP2\tREF_B\n"
                         "C\tchr1\t9000\t11000\tPOSFAM\tP3\tREF_C\n")
            _orthology(root, "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG-P\n"
                             "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG-P\n"
                             "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG-P\n")
            _subfamily(root, "A\tSA\nB\tSB\nC\tSB\n")
            _position(root, "A\trecombinant_core\nB\tsupport_interval\nC\tsupport_interval\n")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", position="position.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertIn("position_discrimination", frame.loc["A", "leading_axes"])
            self.assertEqual("copy_weakened", frame.loc["B", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("locus_specific_only", manifest["selection_basis"]["POSFAM"])


class PendingReviewTargetTests(unittest.TestCase):
    """An LLM-only / unverified target never becomes the sole basis of a preferred copy."""

    def test_pending_review_target_cannot_promote_a_copy_on_its_own(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _families(root, "hyp_pend\tPENDFAM\tA;B\n")
            _atlas(root, "A\tchr1\t1000\t3000\tPENDFAM\tPX\tREF_A\n"
                         "B\tchr1\t5000\t7000\tPENDFAM\tPY\tREF_B\n")
            _orthology(root, "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG-PD\n"
                             "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG-PD\n")
            _subfamily(root, "A\tPX\nB\tPY\n")
            _target(root, "tgt_pend\thyp_pend\tREF_T\tPisum sativum\tPENDFAM\tOG-PD\tPX\tPX\t\t\tllm_proposal\tpending_review")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", target="target.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("pending_review", frame.loc["A", "target_verification_status"])
            self.assertEqual("target_match", frame.loc["A", "subfamily_vs_target"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["A", "disambiguation_status"])
            self.assertIn("pending_review", frame.loc["A", "review_reason"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("copy_target_pending_review", reasons)
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("PENDFAM", manifest["target_relative"]["pending_review_targets_held_back"])


class TargetMatchingTests(unittest.TestCase):
    """Target matching through hypothesis_id, expected_family_id and expected_orthogroup."""

    def test_target_matched_by_expected_orthogroup_applies_to_the_copy_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # No hypothesis_id / expected_family_id on the target: it links to
            # the copy set through the orthogroup layer.
            _families(root, "hyp_og\tOGFAM\tA;B\n")
            _atlas(root, "A\tchr1\t1000\t3000\tOGFAM\tX\tREF_A\n"
                         "B\tchr1\t5000\t7000\tOGFAM\tY\tREF_B\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_A\tOG-T\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_B\tOG-T\n")
            _subfamily(root, "A\tSX\nB\tSY\n")
            _target(root, "tgt_og\t\tREF_T\tPisum sativum\t\tOG-T\tSX\tX\t\t\tpanel\tverified")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv", target="target.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("tgt_og", frame.loc["A", "target_id"])
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertEqual("target_match", frame.loc["A", "subfamily_vs_target"])
            self.assertEqual("target_mismatch", frame.loc["B", "subfamily_vs_target"])

    def test_backward_compatible_columns_are_present_and_honest_without_a_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _families(root, "hyp_bc\tBCFAM\tA;B\n")
            _atlas(root, "A\tchr1\t1000\t3000\tBCFAM\tX\tREF_A\n"
                         "B\tchr1\t5000\t7000\tBCFAM\tX\tREF_B\n")
            _orthology(root, "A\tone_to_many\torthofinder\tOF\tREF_A\tOG-BC\n"
                             "B\tone_to_many\torthofinder\tOF\tREF_B\tOG-BC\n")
            _subfamily(root, "A\tSX\nB\tSX\n")
            out = _run(root, orthology="orthology.tsv", subfamily_map="subfamily.tsv")
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            for column in ("target_id", "target_verification_status", "target_reference_gene",
                           "target_subfamily", "target_domain_architecture", "subfamily_vs_target",
                           "domain_vs_target", "domain_vs_family_modal", "subfamily_family_modal",
                           "target_specific_axes", "selection_basis"):
                self.assertIn(column, set(frame.columns), column)
            self.assertTrue(frame["subfamily_vs_target"].eq("subfamily_unresolved").all())
            self.assertTrue(frame["domain_vs_target"].eq("target_not_evaluated").all())
            # Legacy columns keep their modal semantics for the family modal.
            self.assertTrue(frame["domain_architecture_match"].eq("architecture_match").all())
            self.assertEqual("SX", frame.iloc[0]["subfamily_family_modal"])


if __name__ == "__main__":
    unittest.main()
