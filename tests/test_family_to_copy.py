"""Family-to-copy disambiguation (issue #170)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import family_to_copy, schemas
from candidate_qtl.cli import main


FAMILIES = "\n".join([
    "hypothesis_id\tgene_family\tmember_interval_genes",
    "rev_fam_CLV\tCLV\tG_A;G_B;G_C",
    "rev_fam_WOX\tWOX\tG_W",
]) + "\n"

ATLAS = "\n".join([
    "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
    "G_A\tchr1\t1000\t3000\tCLV1\tLRR-Pkinase\tCLV1",
    "G_B\tchr1\t5000\t7000\tCLV2\tLRR-Pkinase\tCLV2",
    "G_C\tchr1\t900000\t902000\tCLV3\tLRR\tCLV3",
    "G_W\tchr1\t30000\t32000\tWOX4\tHomeobox\tWOX4",
]) + "\n"


class FamilyToCopyTests(unittest.TestCase):
    def _audit_pair(self, root: Path, relationship: str, *, target: bool = False) -> pd.DataFrame:
        def write(name, rows):
            path = root / name
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            return str(path)

        families = write("families.tsv", [{"hypothesis_id": "h", "family_id": "F", "member_interval_genes": "gA;gB"}])
        atlas = write("atlas.tsv", [{"gene_id": gene, "evolutionary_family_id": "F"} for gene in ("gA", "gB")])
        orthology = write("orthology.tsv", [
            {"gene_id": "gA", "reference_gene_id": "refA", "orthogroup": "OG1",
             "relationship": "one_to_one", "evidence_chain": "provider1"},
            {"gene_id": "gB", "reference_gene_id": "refB" if target else "refA",
             "orthogroup": "OG2" if target else "OG1", "relationship": relationship, "evidence_chain": "provider1"},
        ])
        target_path = write("target.tsv", [{
            "target_id": "t", "hypothesis_id": "h", "reference_gene_id": "refA", "reference_species": "plant",
            "expected_family_id": "F", "expected_orthogroup": "OG1", "expected_subfamily": "OG1",
            "expected_domain_architecture": "", "expected_synteny_anchor": "", "expected_expression_context": "",
            "source_ids": "fixture", "verification_status": "verified",
        }]) if target else None
        family_to_copy.family_to_copy(families, str(root / "copy"), atlas_arg=atlas,
                                      orthology_arg=orthology, target_arg=target_path)
        return pd.read_csv(root / "copy/family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")

    def test_unresolved_orthology_is_not_comparative_counterevidence(self):
        for relationship in ("", "unmapped_ambiguous", "unknown"):
            with self.subTest(relationship=relationship), tempfile.TemporaryDirectory() as tmp:
                frame = self._audit_pair(Path(tmp), relationship)
                self.assertEqual({"gA", "gB"}, set(frame.index))
                self.assertFalse(frame["disambiguation_status"].isin({"preferred_copy", "copy_weakened"}).any())
                self.assertEqual("unmapped_ambiguous", frame.loc["gB", "orthology_class"])
                self.assertEqual("", frame.loc["gB", "copy_specific_counterevidence"])
                self.assertFalse(frame["discriminating_axes"].str.contains("orthology_class").any())

    def test_shared_chain_positive_axis_is_not_counterevidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._audit_pair(Path(tmp), "many_to_many", target=True)
            self.assertEqual("target_match", frame.loc["gA", "subfamily_vs_target"])
            self.assertEqual("orthology_class", frame.loc["gA", "leading_axes"])
            self.assertEqual("", frame.loc["gA", "copy_specific_counterevidence"])
            self.assertIn("orthology_class=many_to_many", frame.loc["gB", "copy_specific_counterevidence"])
            self.assertIn("subfamily_match=", frame.loc["gB", "copy_specific_counterevidence"])
            self.assertEqual({"gA", "gB"}, set(frame.index))

    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "families.tsv", "atlas.tsv", "orthology.tsv", "synteny.tsv", "expression.tsv", "position.tsv", "subfamily.tsv",
            "target.tsv")}
        paths["families.tsv"].write_text(FAMILIES, encoding="utf-8")
        paths["atlas.tsv"].write_text(ATLAS, encoding="utf-8")
        paths["target.tsv"].write_text(
            "target_id\thypothesis_id\treference_gene_id\treference_species\texpected_family_id\texpected_orthogroup\texpected_subfamily\texpected_domain_architecture\texpected_synteny_anchor\texpected_expression_context\tsource_ids\tverification_status\n",
            encoding="utf-8")
        paths["orthology.tsv"].write_text("\n".join([
            "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
            "G_A\tone_to_one\torthofinder\tOF_chain\tCLV1\tOG1",
            "G_A\tone_to_one\tomamer\tOF_chain\tCLV1\tOG1",
            "G_B\tone_to_many\torthofinder\tOF_chain\tCLV2\tOG1",
            "G_C\tmany_to_many\torthofinder\tOF_chain\tCLV3\tOG2",
            "G_W\tone_to_one\torthofinder\tOF_chain\tWOX4\tOG3",
        ]) + "\n", encoding="utf-8")
        paths["synteny.tsv"].write_text("\n".join([
            "gene_id\tsynteny_status\tprovider",
            "G_A\tsynteny_supported\tjcvi",
            "G_B\tsynteny_absent\tjcvi",
            "G_C\tsynteny_absent\tjcvi",
        ]) + "\n", encoding="utf-8")
        paths["expression.tsv"].write_text("\n".join([
            "gene_id\trelevant_expression_status",
            "G_A\texpressed_in_relevant_context",
            "G_B\tunmeasured",
            "G_C\texpressed_elsewhere_only",
        ]) + "\n", encoding="utf-8")
        paths["position.tsv"].write_text("\n".join([
            "gene_id\tposition_class",
            "G_A\tbroad_interval_only",
            "G_B\tbroad_interval_only",
            "G_C\tbroad_interval_only",
        ]) + "\n", encoding="utf-8")
        paths["subfamily.tsv"].write_text("gene_id\tsubfamily\nG_A\tCLV1_BAM_like\nG_B\tCLV1_BAM_like\nG_C\tPXY_TDR_like\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, keys: tuple[str, ...] = ("atlas", "orthology", "synteny", "expression", "position", "subfamily")) -> Path:
        fixture = self._fixture(root)
        flags = {"atlas": "--atlas", "orthology": "--orthology", "synteny": "--synteny",
                 "expression": "--expression", "position": "--position", "subfamily": "--subfamily-map"}
        argv = ["family-to-copy", "--families", str(fixture["families.tsv"]), "--out-dir", str(root / "copy"), "--case-id", "demo"]
        for key in keys:
            argv += [flags[key], str(fixture[f"{key}.tsv"])]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def test_discriminating_axes_pick_a_first_wet_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertIn("orthology_class", frame.loc["G_A", "leading_axes"])
            self.assertIn("synteny_status", frame.loc["G_A", "leading_axes"])
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("G_A", set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["G_A", "first_wet_priority"])

    def test_preferred_copy_is_not_a_causal_claim(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("never a causal claim", manifest["meaning_of_preferred_copy"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(shortlist["rank_contribution"].eq("none").all())
            for path in sorted(out.glob("*.tsv")):
                columns = set(pd.read_csv(path, sep="\t", dtype=str).columns)
                self.assertFalse({"dry_score", "score", "dry_group"} & columns, path.name)

    def test_family_level_support_is_kept_separate_from_copy_specific_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertIn("member of CLV", frame.loc["G_B", "family_level_support"])
            # G_B loses on orthology, synteny and subfamily-free axes; those are counterevidence.
            self.assertIn("orthology_class=one_to_many", frame.loc["G_B", "copy_specific_counterevidence"])
            self.assertNotIn("member of", frame.loc["G_B", "copy_specific_support"])

    def test_unmeasured_expression_is_not_a_negative(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("unmeasured", frame.loc["G_B", "relevant_expression_status"])
            # G_B is weakened by orthology/synteny, never by the absence of a measurement.
            self.assertNotIn("relevant_expression_status=unmeasured", frame.loc["G_B", "copy_specific_support"])

    def test_broad_interval_membership_does_not_discriminate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            axes = set(frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertNotIn("position_discrimination", axes)

    def test_one_evidence_chain_is_counted_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            # No curated subfamily map, so the subfamily is derived from the orthogroup
            # and therefore rests on the same evidence chain as the orthology class.
            # G_C belongs to a different orthogroup (OG2), so it is outside the
            # validated copy set {G_A, G_B} and is never compared as a paralog.
            out = self._run(Path(tmp), keys=("atlas", "orthology", "synteny", "expression", "position"))
            relationships = pd.read_csv(out / "family_copy_relationships.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # Two providers agree on G_A; that is robustness, and the class is still one_to_one.
            self.assertEqual("omamer;orthofinder", relationships.loc["G_A", "orthology_providers"])
            self.assertEqual("one_to_one", relationships.loc["G_A", "orthology_class"])
            self.assertEqual("orthogroup", relationships.loc["G_A", "subfamily_source"])
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # Inside the validated component both members share OG1, so the
            # orthogroup-derived subfamily is identical and can never add a
            # second, orthology-derived vote.
            self.assertNotIn("subfamily_match", frame.loc["G_A", "discriminating_axes"])
            # G_C (OG2) is not compared as a paralog of the OG1 component.
            self.assertEqual("not_in_validated_copy_set", frame.loc["G_C", "disambiguation_status"])
            self.assertNotIn("G_C", set(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)["candidate_gene_id"]))
            leading = frame.loc["G_A", "leading_axes"].split(";")
            self.assertIn("orthology_class", leading)
            self.assertNotIn("subfamily_match", leading)

    def test_a_curated_subfamily_is_an_independent_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            # G_A and G_B are paralogs of the same orthogroup but belong to
            # different curated subfamilies. Since issue #185 the subfamily axis
            # is target-relative: without a target the curated map alone cannot
            # drive selection, so a verified target supplies the expected
            # subfamily and the curated map discriminates within the validated
            # copy set on its own evidence chain.
            fixture["subfamily.tsv"].write_text(
                "gene_id\tsubfamily\nG_A\tCLV1_BAM_like\nG_B\tPXY_TDR_like\nG_C\tCLV3_peptide\n", encoding="utf-8")
            fixture["target.tsv"].write_text("\n".join([
                "target_id\thypothesis_id\treference_gene_id\treference_species\texpected_family_id\texpected_orthogroup\texpected_subfamily\texpected_domain_architecture\texpected_synteny_anchor\texpected_expression_context\tsource_ids\tverification_status",
                "tgt_clv\trev_fam_CLV\tCLV1\tArabidopsis\tCLV\tOG1\tCLV1_BAM_like\tLRR-Pkinase\t\t\tpaper:X\tverified",
            ]) + "\n", encoding="utf-8")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(fixture["families.tsv"]),
                                      "--atlas", str(fixture["atlas.tsv"]), "--orthology", str(fixture["orthology.tsv"]),
                                      "--subfamily-map", str(fixture["subfamily.tsv"]),
                                      "--family-target", str(fixture["target.tsv"]),
                                      "--out-dir", str(out), "--case-id", "demo"]))
            relationships = pd.read_csv(out / "family_copy_relationships.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("curated_map", relationships.loc["G_A", "subfamily_source"])
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            leading = frame.loc["G_A", "leading_axes"].split(";")
            # A curated map is not derived from the orthology chain, so it counts on its own.
            self.assertIn("subfamily_match", leading)
            self.assertEqual(len(leading), len(set(leading)))
            self.assertEqual("target_match", frame.loc["G_A", "subfamily_vs_target"])
            self.assertEqual("target_mismatch", frame.loc["G_B", "subfamily_vs_target"])

    def test_family_level_only_without_axes_never_promotes_a_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, keys=("atlas",))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # The atlas alone carries only a symbol-stem naming bucket (no
            # orthogroup, no curated family, no domain provenance): issue #184
            # gates the family as provisional, so no member is compared as a
            # copy and none is promoted.
            self.assertTrue(frame["disambiguation_status"].eq("family_resolution_required").all())
            self.assertTrue(frame["copy_set_status"].eq("provisional").all())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("family_resolution_required", reasons)

    def test_no_discriminating_input_leaves_everything_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, keys=())
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frame["disambiguation_status"].eq("insufficient_discriminating_evidence").all())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("no_discriminating_input_supplied", reasons)

    def test_tandem_ambiguity_is_preserved_instead_of_forcing_a_winner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            # Two tandem copies, identical on every supplied axis.
            fixture["families.tsv"].write_text("hypothesis_id\tgene_family\tmember_interval_genes\nfam\tCLV\tG_A;G_B\n", encoding="utf-8")
            fixture["atlas.tsv"].write_text("\n".join([
                "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                "G_A\tchr1\t1000\t3000\tCLV1\tLRR-Pkinase\tCLV1",
                "G_B\tchr1\t4000\t6000\tCLV1\tLRR-Pkinase\tCLV1",
            ]) + "\n", encoding="utf-8")
            fixture["orthology.tsv"].write_text(
                "gene_id\trelationship\tprovider\treference_gene_id\torthogroup\n"
                "G_A\tclose_paralog\torthofinder\tCLV1\tOG1\n"
                "G_B\tclose_paralog\torthofinder\tCLV1\tOG1\n", encoding="utf-8")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(fixture["families.tsv"]),
                                      "--atlas", str(fixture["atlas.tsv"]), "--orthology", str(fixture["orthology.tsv"]),
                                      "--out-dir", str(out)]))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frame["disambiguation_status"].eq("tandem_ambiguous").all())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("copy_tandem_ambiguous", reasons)

    def test_subfamily_is_derived_or_honestly_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, keys=("atlas", "subfamily"))
            relationships = pd.read_csv(out / "family_copy_relationships.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("CLV1_BAM_like", relationships.loc["G_A", "subfamily"])
            self.assertEqual("curated_map", relationships.loc["G_A", "subfamily_source"])
            # No curated entry and no orthogroup for G_W: unresolved, never invented.
            self.assertEqual("", relationships.loc["G_W", "subfamily"])
            self.assertEqual("unresolved", relationships.loc["G_W", "subfamily_source"])

    def test_single_member_family_and_reduction_ratio(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", frame.loc["G_W", "disambiguation_status"])
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(4, manifest["family_member_count"])
            self.assertEqual(2, manifest["first_wet_candidate_count"])
            self.assertEqual(2.0, manifest["family_to_first_wet_reduction"])


def _write_table(root: Path, name: str, header: str, records: str) -> Path:
    path = root / name
    path.write_text(header + "\n" + records + "\n", encoding="utf-8")
    return path


class MissingnessIsNotComparativeTests(unittest.TestCase):
    """Issue #281: unmeasured / not-evaluated values are never comparative loss votes."""

    def _run(self, root: Path, *, synteny: str | None = None, expression: str | None = None,
             position: str | None = None, subfamily: str | None = None, target: str | None = None) -> Path:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     "hyp281\tFAM281\tA;B")
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                     "A\tchr1\t1000\t3000\tFAM281\tDA\tREF_A\nB\tchr1\t200000\t202000\tFAM281\t\tREF_B")
        _write_table(root, "orthology.tsv",
                     "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                     "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG281\n"
                     "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG281")
        argv = ["family-to-copy", "--families", str(root / "families.tsv"), "--atlas", str(root / "atlas.tsv"),
                "--orthology", str(root / "orthology.tsv"), "--out-dir", str(root / "copy"), "--case-id", "demo"]
        for flag, content in (("--synteny", synteny), ("--expression", expression), ("--position", position),
                              ("--subfamily-map", subfamily), ("--family-target", target)):
            if content is not None:
                argv += [flag, str(_write_table(root, flag[2:] + ".tsv", content[0], content[1]))]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def _frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")

    def test_synteny_missingness_does_not_discriminate_or_weaken(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A is measured adverse (synteny_absent), B is simply not evaluated.
            out = self._run(root, synteny=("gene_id\tsynteny_status\tprovider", "A\tsynteny_absent\tjcvi"))
            frame = self._frame(out)
            self.assertNotIn("synteny_status", frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertTrue(frame["disambiguation_status"].eq("family_supported_copy_unresolved").all())
            self.assertFalse(frame["disambiguation_status"].eq("preferred_copy").any())
            self.assertFalse(frame["disambiguation_status"].eq("copy_weakened").any())
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_expression_missingness_does_not_discriminate_or_weaken(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A is measured negative, B is unmeasured.
            out = self._run(root, expression=("gene_id\trelevant_expression_status",
                                              "A\tnot_detected_in_measured_context"))
            frame = self._frame(out)
            self.assertNotIn("relevant_expression_status", frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertTrue(frame["disambiguation_status"].eq("family_supported_copy_unresolved").all())
            self.assertFalse(frame["disambiguation_status"].eq("preferred_copy").any())
            self.assertFalse(frame["disambiguation_status"].eq("copy_weakened").any())
            self.assertEqual("", frame.iloc[0]["copy_specific_counterevidence"])
            self.assertEqual("", frame.iloc[1]["copy_specific_counterevidence"])

    def test_target_missingness_does_not_discriminate_or_weaken(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A mismatches the verified target on the domain axis; B is not
            # evaluated (no architecture), so the target axes cannot compare.
            out = self._run(root,
                            target=("target_id\thypothesis_id\treference_gene_id\treference_species"
                                    "\texpected_family_id\texpected_orthogroup\texpected_subfamily"
                                    "\texpected_domain_architecture\texpected_synteny_anchor"
                                    "\texpected_expression_context\tsource_ids\tverification_status",
                                    "tgt281\thyp281\tREF_T\tPisum sativum\tFAM281\tOG281\tSX\tDX\t\t\tpanel\tverified"))
            frame = self._frame(out)
            self.assertEqual("target_mismatch", frame.loc["A", "domain_vs_target"])
            self.assertEqual("target_not_evaluated", frame.loc["B", "domain_vs_target"])
            axes = set(frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertNotIn("subfamily_match", axes)
            self.assertNotIn("domain_vs_target", axes)
            self.assertFalse(frame["disambiguation_status"].eq("preferred_copy").any())
            self.assertFalse(frame["disambiguation_status"].eq("copy_weakened").any())
            self.assertEqual("", frame.loc["B", "copy_specific_counterevidence"])

    def test_two_evaluable_members_keep_the_existing_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Both members evaluable and adverse is still a real comparison.
            out = self._run(root, synteny=("gene_id\tsynteny_status\tprovider",
                                           "A\tsynteny_supported\tjcvi\nB\tsynteny_absent\tjcvi"))
            frame = self._frame(out)
            self.assertIn("synteny_status", frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertEqual("copy_weakened", frame.loc["B", "disambiguation_status"])

    def test_uninformative_bystander_is_unresolved_not_weakened(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A and B discriminate on synteny; C is a not-evaluated bystander.
            # C must be unresolved (coverage gap), never weakened by its own
            # missingness, and its missing value must not appear as counterevidence.
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp281b\tFAM281B\tA;B;C")
            _write_table(root, "atlas.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                         "A\tchr1\t1000\t3000\tFAM281B\tDA\tREF_A\nB\tchr1\t5000\t7000\tFAM281B\tDA\tREF_B\n"
                         "C\tchr1\t9000\t11000\tFAM281B\t\tREF_C")
            _write_table(root, "orthology.tsv",
                         "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                         "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG281\n"
                         "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG281\n"
                         "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG281")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                      "--atlas", str(root / "atlas.tsv"), "--orthology", str(root / "orthology.tsv"),
                                      "--synteny", str(_write_table(root, "synteny.tsv",
                                                                    "gene_id\tsynteny_status\tprovider",
                                                                    "A\tsynteny_supported\tjcvi\n"
                                                                    "B\tsynteny_absent\tjcvi\n"
                                                                    "C\tsynteny_not_evaluated\tjcvi")),
                                      "--out-dir", str(out)]))
            frame = self._frame(out)
            self.assertIn("synteny_status", frame.iloc[0]["discriminating_axes"].split(";"))
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertEqual("copy_weakened", frame.loc["B", "disambiguation_status"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["C", "disambiguation_status"])
            self.assertEqual("", frame.loc["C", "copy_specific_counterevidence"])


class RecombinantCoreGateTests(unittest.TestCase):
    """Issue #282: recombinant_core is a hard positional gate over preferred copies."""

    def _fixture(self, root: Path, members: str, *, orthology: str, position: str,
                 synteny: str | None = None, expression: str | None = None) -> Path:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     f"hyp282\tFAM282\t{members}")
        atlas_rows = "\n".join(
            f"{gene}\tchr1\t{1000 + i * 4000}\t{3000 + i * 4000}\tFAM282\tDA\tREF_{gene}"
            for i, gene in enumerate(members.split(";")))
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                     atlas_rows)
        _write_table(root, "orthology.tsv",
                     "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                     orthology)
        _write_table(root, "position.tsv", "gene_id\tposition_class", position)
        argv = ["family-to-copy", "--families", str(root / "families.tsv"), "--atlas", str(root / "atlas.tsv"),
                "--orthology", str(root / "orthology.tsv"), "--position", str(root / "position.tsv"),
                "--out-dir", str(root / "copy"), "--case-id", "demo"]
        for flag, content in (("--synteny", synteny), ("--expression", expression)):
            if content is not None:
                argv += [flag, str(_write_table(root, flag[2:] + ".tsv", content[0], content[1]))]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def _frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")

    def test_recombinant_core_copy_is_not_overtaken_by_comparative_axes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A lies in the recombinant-defined core (position only); B leads
            # orthology, synteny and expression. B must not become preferred.
            out = self._fixture(
                root, "A;B",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_one\torthofinder\tOF\tREF_B\tOG282",
                position="A\trecombinant_core\nB\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_absent\tjcvi\nB\tsynteny_supported\tjcvi"),
                expression=("gene_id\trelevant_expression_status",
                            "A\texpressed_elsewhere_only\nB\texpressed_in_relevant_context"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertIn("position_discrimination", frame.loc["A", "leading_axes"])
            self.assertNotIn("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["B", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])

    def test_multiple_cores_tie_break_among_core_members_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A and B are both core; C is a support_interval copy leading
            # orthology. The comparison must happen inside the core (A beats B on
            # synteny); C can never win.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_many\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tone_to_one\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core\nC\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_supported\tjcvi\nB\tsynteny_absent\tjcvi\nC\tsynteny_not_evaluated\tjcvi"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["C", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))

    def test_two_indistinguishable_cores_keep_a_tie(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Both cores are identical on every supplied axis: the tie must be
            # preserved inside the core rather than handed to the non-core copy.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core\nC\tsupport_interval")
            frame = self._frame(out)
            self.assertEqual({"co_preferred_tie"}, set(frame.loc[["A", "B"], "disambiguation_status"]))
            self.assertNotIn("preferred_copy", set(frame["disambiguation_status"]))
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A", "B"}, set(shortlist["candidate_gene_id"]))
            self.assertTrue(shortlist.set_index("candidate_gene_id").loc[["A", "B"], "first_wet_priority"].eq("co_first").all())

    def test_indistinguishable_cores_keep_the_tie_even_when_a_non_core_leads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A and B are both core and identical (synteny_absent); C is not
            # position-evaluated and carries the only synteny value of the
            # family. Since issue #282 the comparison runs inside the core
            # only, where nothing discriminates: the cores keep their
            # co-preferred tie, C can never win, and C's row stays in the
            # artifact with the core gate recorded in the review_reason.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_absent\tjcvi\nB\tsynteny_absent\tjcvi\nC\tsynteny_supported\tjcvi"))
            frame = self._frame(out)
            self.assertEqual({"co_preferred_tie"}, set(frame.loc[["A", "B"], "disambiguation_status"]))
            self.assertNotIn("preferred_copy", set(frame["disambiguation_status"]))
            self.assertEqual("tandem_ambiguous", frame.loc["C", "disambiguation_status"])
            self.assertIn("recombinant_core gate", frame.loc["C", "review_reason"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A", "B"}, set(shortlist["candidate_gene_id"]))

    def test_without_a_core_the_comparative_axes_keep_today_s_behavior(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # No recombinant_core anywhere: the copy leading three axes wins as before.
            out = self._fixture(
                root, "A;B",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_one\torthofinder\tOF\tREF_B\tOG282",
                position="A\tpeak_overlap\nB\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_absent\tjcvi\nB\tsynteny_supported\tjcvi"),
                expression=("gene_id\trelevant_expression_status",
                            "A\texpressed_elsewhere_only\nB\texpressed_in_relevant_context"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["A", "disambiguation_status"])

    def test_sole_core_member_is_kept_as_first_wet_priority(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # The core contains exactly one copy (A) while the non-core copy (B)
            # is not position-evaluated and leads every comparative axis. A must
            # stay the first-Wet priority; B must not become preferred.
            out = self._fixture(
                root, "A;B",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_one\torthofinder\tOF\tREF_B\tOG282",
                position="A\trecombinant_core",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_absent\tjcvi\nB\tsynteny_supported\tjcvi"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["B", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])

    def test_core_copy_with_partial_synteny_is_not_overtaken_by_a_supported_non_core(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Requirement A: the core copy (A, synteny_partial) must stay
            # first-Wet even though the non-core copy (B, support_interval)
            # carries the best synteny tier of the family. The comparative axes
            # rank inside the core; B can never become preferred.
            out = self._fixture(
                root, "A;B",
                orthology="A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG282",
                position="A\trecombinant_core\nB\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_partial\tjcvi\nB\tsynteny_supported\tjcvi"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["B", "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])

    def test_two_cores_compare_inside_the_pool_only_and_a_leading_non_core_loses(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Requirement B: A and B are both core; C (support_interval,
            # one_to_one, expression_supported) leads every family-wide axis.
            # The comparison must run inside the core only, so A
            # (synteny_supported) beats B (synteny_absent); C can never become
            # preferred, and C's gate is auditable in its review_reason.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_many\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tone_to_one\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core\nC\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_supported\tjcvi\nB\tsynteny_absent\tjcvi"),
                expression=("gene_id\trelevant_expression_status",
                            "C\texpressed_in_relevant_context"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertNotIn("preferred_copy", frame.loc["C", "disambiguation_status"])
            self.assertEqual("family_supported_copy_unresolved", frame.loc["C", "disambiguation_status"])
            self.assertIn("recombinant_core gate", frame.loc["C", "review_reason"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])

    def test_sole_core_with_no_comparative_evidence_stays_first_wet(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Requirement C: the sole recombinant_core member is held as the
            # first-Wet candidate even when no axis discriminates at all (the
            # non-core copies are not even position-evaluated). The review_reason
            # records the core gate; this is a positional hold, not fabricated
            # comparative evidence.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core")
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["A", "disambiguation_status"])
            self.assertIn("recombinant_core gate", frame.loc["A", "review_reason"])
            self.assertNotIn("preferred_copy", frame.loc[["B", "C"], "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])

    def test_non_core_members_stay_in_the_artifact_with_an_auditable_gate_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Requirement E: non-core copies are never removed from the
            # artifact; their discrimination rows record the recombinant_core
            # gate so the hold is auditable without a new column.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_many\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tone_to_one\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core\nC\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_supported\tjcvi\nB\tsynteny_absent\tjcvi"))
            frame = self._frame(out)
            for gene_id in ("A", "B", "C"):
                self.assertIn(gene_id, set(frame.index))
            self.assertIn("recombinant_core gate", frame.loc["C", "review_reason"])
            self.assertEqual("", frame.loc["A", "review_reason"])
            # The schema is unchanged by the gate: the same column set as a
            # run with no core member anywhere.
            no_core = self._fixture(
                root, "A;B;C",
                orthology="A\tone_to_many\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tone_to_many\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tone_to_one\torthofinder\tOF\tREF_C\tOG282",
                position="A\tpeak_overlap\nB\tpeak_overlap\nC\tpeak_overlap",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_supported\tjcvi\nB\tsynteny_absent\tjcvi"))
            columns_core = set(pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).columns)
            columns_no_core = set(pd.read_csv(no_core / "family_copy_discrimination.tsv", sep="\t", dtype=str).columns)
            self.assertEqual(columns_no_core, columns_core)

    def test_pool_comparison_picks_the_partial_core_over_a_supported_non_core(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Requirement F (real-data probe): two cores (A=synteny_absent,
            # B=synteny_partial) plus a non-core (C=support_interval,
            # synteny_supported). Inside the pool synteny discriminates with
            # best=partial, so B must become preferred_copy -- C's supported
            # tier neither steals the axis nor turns position discriminating
            # against the cores.
            out = self._fixture(
                root, "A;B;C",
                orthology="A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG282\n"
                          "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG282\n"
                          "C\tfamily_level_only\torthofinder\tOF\tREF_C\tOG282",
                position="A\trecombinant_core\nB\trecombinant_core\nC\tsupport_interval",
                synteny=("gene_id\tsynteny_status\tprovider",
                         "A\tsynteny_absent\tjcvi\nB\tsynteny_partial\tjcvi\nC\tsynteny_supported\tjcvi"))
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["B", "disambiguation_status"])
            self.assertIn("synteny_status", frame.loc["B", "leading_axes"])
            self.assertNotIn("preferred_copy", frame.loc[["A", "C"], "disambiguation_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"B"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.set_index("candidate_gene_id").loc["B", "first_wet_priority"])


BRIDGE_HEADER = ("case_id\tgene_id\tfamily_id\troute\treference_accession\tsource_id\tmatch_type\tidentity"
                 "\tquery_coverage\ttarget_coverage\tevalue\tparalog_count\tambiguity_status\tprovider_status"
                 "\tprovider_version\tinput_fingerprint\toutput_fingerprint\treview_status"
                 "\tdependency_cluster_id\tmapping_tier\trank_contribution\treason")


def _bridge_rows(records: list[tuple[str, str, str, str, str, str]]) -> str:
    """Compact bridge record writer: (gene_id, family_id, route, mapping_tier, ambiguity_status, dependency_cluster_id)."""
    lines = [BRIDGE_HEADER]
    for gene_id, family_id, route, mapping_tier, ambiguity, dependency in records:
        lines.append("\t".join(["demo", gene_id, family_id, route, "", "", "", "", "", "", "", "", ambiguity,
                                "complete", "", "ff", "ff", "", dependency, mapping_tier, "none", ""]))
    return "\n".join(lines) + "\n"


class FamilyToCopyBridgeTests(unittest.TestCase):
    """Issue #374: family_candidate_bridge wiring into family-to-copy.

    The bridge is a parallel, record-only layer: mapping tier / dependency /
    ambiguity / route support are recorded in new columns, tier5/tier6 members
    are blocked from the shortlist without touching their disambiguation
    status, and without a bridge the legacy schema stays byte-identical.
    """

    def _fixture(self, root: Path, *, bridge_records: list[tuple[str, str, str, str, str, str]]) -> Path:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     "hyp374\tCLV\tG_A;G_B;G_C\nhyp374b\tWOX\tG_W\nhyp374c\tRES\tG_R1;G_R2")
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                     "G_A\tchr1\t1000\t3000\tCLV\tLRR-Pkinase\tCLV1\n"
                     "G_B\tchr1\t5000\t7000\tCLV\tLRR-Pkinase\tCLV2\n"
                     "G_C\tchr1\t900000\t902000\tCLV\tLRR\tCLV3\n"
                     "G_W\tchr1\t30000\t32000\tWOX\tHomeobox\tWOX4\n"
                     "G_R1\tchr2\t1000\t3000\tRES\t\t\n"
                     "G_R2\tchr2\t5000\t7000\tRES\t\t")
        _write_table(root, "orthology.tsv",
                     "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                     "G_A\tone_to_one\torthofinder\tOF_chain\tCLV1\tOG1\n"
                     "G_B\tone_to_many\torthofinder\tOF_chain\tCLV2\tOG1\n"
                     "G_C\tmany_to_many\torthofinder\tOF_chain\tCLV3\tOG2\n"
                     "G_W\tone_to_one\torthofinder\tOF_chain\tWOX4\tOG3")
        _write_table(root, "synteny.tsv", "gene_id\tsynteny_status\tprovider",
                     "G_A\tsynteny_supported\tjcvi\nG_B\tsynteny_absent\tjcvi\nG_C\tsynteny_absent\tjcvi")
        _write_table(root, "expression.tsv", "gene_id\trelevant_expression_status",
                     "G_A\texpressed_in_relevant_context\nG_B\tunmeasured\nG_C\texpressed_elsewhere_only")
        _write_table(root, "position.tsv", "gene_id\tposition_class",
                     "G_A\tbroad_interval_only\nG_B\tbroad_interval_only\nG_C\tbroad_interval_only")
        _write_table(root, "subfamily.tsv", "gene_id\tsubfamily",
                     "G_A\tCLV1_BAM_like\nG_B\tCLV1_BAM_like\nG_C\tPXY_TDR_like")
        bridge = root / "bridge.tsv"
        bridge.write_text(_bridge_rows(bridge_records), encoding="utf-8")
        argv = ["family-to-copy", "--families", str(root / "families.tsv"), "--out-dir", str(root / "copy"),
                "--case-id", "demo", "--atlas", str(root / "atlas.tsv"), "--orthology", str(root / "orthology.tsv"),
                "--synteny", str(root / "synteny.tsv"), "--expression", str(root / "expression.tsv"),
                "--position", str(root / "position.tsv"), "--subfamily-map", str(root / "subfamily.tsv"),
                "--bridge", str(bridge)]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def _frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")

    @staticmethod
    def _main_bridge() -> list[tuple[str, str, str, str, str, str]]:
        return [
            ("G_A", "CLV", "a", "tier1_source_backed_unique", "resolved_single_copy", "dep_CLV"),
            ("G_B", "CLV", "a", "tier2_source_backed_ambiguous", "resolved_single_copy", "dep_CLV"),
            ("G_C", "CLV", "b", "tier5_name_only", "paralog_ambiguous", "dep_CLV"),
            ("G_W", "WOX", "b", "tier5_name_only", "paralog_ambiguous", "dep_WOX"),
            ("G_R1", "RES", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep_RES"),
            ("G_R2", "RES", "a", "tier4_sequence_match_ambiguous", "coverage_insufficient", "dep_RES"),
        ]

    def test_without_a_bridge_the_legacy_schema_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run_without_bridge(Path(tmp))
            columns = set(pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).columns)
            self.assertFalse({"bridge_mapping_tier", "bridge_promotion_gate", "bridge_missingness",
                              "route_rescue_state", "copy_selection_loss"} & columns)
            shortlist_columns = set(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).columns)
            self.assertNotIn("bridge_mapping_tier", shortlist_columns)
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["bridge"]["supplied"])
            self.assertIsNone(manifest["bridge"]["bridge_tsv_sha256"])
            self.assertEqual(0, manifest["bridge"]["blocked_count"])
            self.assertIsNone(manifest["inputs"]["bridge"])

    def _run_without_bridge(self, root: Path) -> Path:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     "hyp374\tCLV\tG_A;G_B;G_C\nhyp374b\tWOX\tG_W\nhyp374c\tRES\tG_R1;G_R2")
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                     "G_A\tchr1\t1000\t3000\tCLV\tLRR-Pkinase\tCLV1\n"
                     "G_B\tchr1\t5000\t7000\tCLV\tLRR-Pkinase\tCLV2\n"
                     "G_C\tchr1\t900000\t902000\tCLV\tLRR\tCLV3\n"
                     "G_W\tchr1\t30000\t32000\tWOX\tHomeobox\tWOX4\n"
                     "G_R1\tchr2\t1000\t3000\tRES\t\t\n"
                     "G_R2\tchr2\t5000\t7000\tRES\t\t")
        _write_table(root, "orthology.tsv",
                     "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                     "G_A\tone_to_one\torthofinder\tOF_chain\tCLV1\tOG1\n"
                     "G_B\tone_to_many\torthofinder\tOF_chain\tCLV2\tOG1\n"
                     "G_C\tmany_to_many\torthofinder\tOF_chain\tCLV3\tOG2\n"
                     "G_W\tone_to_one\torthofinder\tOF_chain\tWOX4\tOG3")
        _write_table(root, "synteny.tsv", "gene_id\tsynteny_status\tprovider",
                     "G_A\tsynteny_supported\tjcvi\nG_B\tsynteny_absent\tjcvi\nG_C\tsynteny_absent\tjcvi")
        _write_table(root, "expression.tsv", "gene_id\trelevant_expression_status",
                     "G_A\texpressed_in_relevant_context\nG_B\tunmeasured\nG_C\texpressed_elsewhere_only")
        _write_table(root, "position.tsv", "gene_id\tposition_class",
                     "G_A\tbroad_interval_only\nG_B\tbroad_interval_only\nG_C\tbroad_interval_only")
        _write_table(root, "subfamily.tsv", "gene_id\tsubfamily",
                     "G_A\tCLV1_BAM_like\nG_B\tCLV1_BAM_like\nG_C\tPXY_TDR_like")
        argv = ["family-to-copy", "--families", str(root / "families.tsv"), "--out-dir", str(root / "copy"),
                "--case-id", "demo", "--atlas", str(root / "atlas.tsv"), "--orthology", str(root / "orthology.tsv"),
                "--synteny", str(root / "synteny.tsv"), "--expression", str(root / "expression.tsv"),
                "--position", str(root / "position.tsv"), "--subfamily-map", str(root / "subfamily.tsv")]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def test_tier1_member_is_promoted_and_the_gate_records_passed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), bridge_records=self._main_bridge())
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("tier1_source_backed_unique", frame.loc["G_A", "bridge_mapping_tier"])
            self.assertEqual("passed", frame.loc["G_A", "bridge_promotion_gate"])
            self.assertEqual("covered", frame.loc["G_A", "bridge_missingness"])
            self.assertEqual("dep_CLV", frame.loc["G_A", "bridge_dependency_cluster_id"])
            self.assertEqual("resolved_single_copy", frame.loc["G_A", "bridge_ambiguity_status"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"G_A"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual("tier1_source_backed_unique",
                             shortlist.set_index("candidate_gene_id").loc["G_A", "bridge_mapping_tier"])
            self.assertTrue(shortlist["rank_contribution"].eq("none").all())
            self.assertTrue(frame["rank_contribution"].eq("none").all())

    def test_tier5_and_tier6_members_never_enter_the_shortlist_but_keep_their_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # G_W (WOX) is a sole-copy family whose internal logic makes it
            # preferred_copy; the bridge tier5 gate must keep it out of the
            # shortlist while the disambiguation status itself stays untouched.
            out = self._fixture(root, bridge_records=self._main_bridge())
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_W", "disambiguation_status"])
            self.assertEqual("blocked_tier5", frame.loc["G_W", "bridge_promotion_gate"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("G_W", set(shortlist["candidate_gene_id"]))
            # A co-preferred tie inside a component is gated identically.
            root2 = Path(tmp) / "tier6"
            root2.mkdir()
            _write_table(root2, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes", "hyp374d\tTIE\tA;B")
            _write_table(root2, "atlas.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                         "A\tchr1\t1000\t3000\tTIE\tDA\tREF_A\nB\tchr1\t5000\t7000\tTIE\tDA\tREF_B")
            _write_table(root2, "orthology.tsv",
                         "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                         "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG_TIE\nB\tfamily_level_only\torthofinder\tOF\tREF_B\tOG_TIE")
            _write_table(root2, "position.tsv", "gene_id\tposition_class",
                         "A\trecombinant_core\nB\trecombinant_core")
            bridge2 = root2 / "bridge.tsv"
            bridge2.write_text(_bridge_rows([
                ("A", "TIE", "a", "tier6_not_evaluated", "not_evaluated", "dep_TIE"),
                ("B", "TIE", "a", "tier6_not_evaluated", "not_evaluated", "dep_TIE"),
            ]), encoding="utf-8")
            argv = ["family-to-copy", "--families", str(root2 / "families.tsv"), "--out-dir", str(root2 / "copy"),
                    "--case-id", "demo", "--atlas", str(root2 / "atlas.tsv"), "--orthology", str(root2 / "orthology.tsv"),
                    "--position", str(root2 / "position.tsv"), "--bridge", str(bridge2)]
            self.assertEqual(0, main(argv))
            tie_frame = pd.read_csv(root2 / "copy" / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual({"co_preferred_tie"}, set(tie_frame["disambiguation_status"]))
            self.assertEqual({"blocked_tier6"}, set(tie_frame["bridge_promotion_gate"]))
            tie_shortlist = pd.read_csv(root2 / "copy" / "family_copy_shortlist.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(tie_shortlist))
            self.assertEqual(2, json.loads((root2 / "copy" / "family_copy_manifest.json").read_text(encoding="utf-8"))["bridge"]["blocked_count"])

    def test_a_member_without_a_bridge_row_is_reported_not_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            # The bridge covers the CLV family only: G_W has no bridge row, so
            # its promotion is not blocked (absence is reported, never used as
            # a block) and its missingness is explicit.
            out = self._fixture(Path(tmp), bridge_records=[
                ("G_A", "CLV", "a", "tier1_source_backed_unique", "resolved_single_copy", "dep_CLV"),
                ("G_B", "CLV", "a", "tier2_source_backed_ambiguous", "resolved_single_copy", "dep_CLV"),
            ])
            frame = self._frame(out)
            self.assertEqual("no_bridge_row", frame.loc["G_W", "bridge_missingness"])
            self.assertEqual("", frame.loc["G_W", "bridge_mapping_tier"])
            self.assertEqual("passed", frame.loc["G_W", "bridge_promotion_gate"])
            self.assertEqual("covered", frame.loc["G_A", "bridge_missingness"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("G_W", set(shortlist["candidate_gene_id"]))

    def test_bridge_with_missing_required_columns_fails_closed(self) -> None:
        """Issue #374 review: a bridge without the required gene_id / family_id
        / route / mapping_tier columns is rejected before any row is read."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp_vocab\tCLV\tG_A")
            bridge = root / "bridge.tsv"
            bridge.write_text("case_id\tgene_id\tfamily_id\nx\tG_A\tCLV\n", encoding="utf-8")
            with self.assertRaises(ValueError) as raised:
                main(["family-to-copy", "--families", str(root / "families.tsv"),
                      "--out-dir", str(root / "copy"), "--case-id", "demo", "--bridge", str(bridge)])
            self.assertIn("missing columns", str(raised.exception))
            self.assertIn("mapping_tier", str(raised.exception))
            self.assertIn("route", str(raised.exception))

    def test_bridge_outside_the_canonical_vocabulary_fails_closed(self) -> None:
        """Issue #374 review: a mapping_tier typo / legacy value or a route
        outside a/b/c is rejected (ValueError) so an unknown value can never
        read as ``passed`` and get promoted."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp_vocab\tCLV\tG_A")
            cases = [
                (("G_A", "CLV", "a", "tier5_name_only_typo", "resolved_single_copy", "dep"), "canonical vocabulary"),
                (("G_A", "CLV", "d", "tier1_source_backed_unique", "resolved_single_copy", "dep"), "a/b/c"),
                (("G_A", "CLV", "a", "tier0_legacy", "resolved_single_copy", "dep"), "canonical vocabulary"),
            ]
            for record, expected in cases:
                with self.subTest(record=record):
                    bridge = root / "bridge.tsv"
                    bridge.write_text(_bridge_rows([record]), encoding="utf-8")
                    with self.assertRaises(ValueError) as raised:
                        main(["family-to-copy", "--families", str(root / "families.tsv"),
                              "--out-dir", str(root / "copy"), "--case-id", "demo",
                              "--bridge", str(bridge)])
                    self.assertIn(expected, str(raised.exception))

    def test_co_preferred_tie_with_one_tier5_side_promotes_only_the_source_backed_side(self) -> None:
        """Issue #374 review: inside one component the promotion gate applies
        per member. A co-preferred tie of a tier1 member and a tier5 member
        puts only the tier1 side into the shortlist; the tier5 side keeps its
        tie status but is blocked."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_tie_fixture(root, [
                ("A", "TIE", "a", "tier1_source_backed_unique", "resolved_single_copy", "dep_TIE"),
                ("B", "TIE", "a", "tier5_name_only", "paralog_ambiguous", "dep_TIE"),
            ])
            out = root / "copy"
            tie_frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t",
                                    dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual({"co_preferred_tie"}, set(tie_frame["disambiguation_status"]))
            self.assertEqual("passed", tie_frame.loc["A", "bridge_promotion_gate"])
            self.assertEqual("blocked_tier5", tie_frame.loc["B", "bridge_promotion_gate"])
            tie_shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"A"}, set(tie_shortlist["candidate_gene_id"]))
            self.assertEqual("tier1_source_backed_unique",
                             tie_shortlist.set_index("candidate_gene_id").loc["A", "bridge_mapping_tier"])
            self.assertEqual("co_first", tie_shortlist.set_index("candidate_gene_id").loc["A", "first_wet_priority"])
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["bridge"]["blocked_count"])
            # A promoted side exists, so the family-level basis is not lost.
            self.assertEqual({"TIE": "locus_specific_only"}, manifest["selection_basis"])

    def test_fully_blocked_tie_never_counts_as_a_selection(self) -> None:
        """Issue #374 review: when every promoted member of a family is gate
        blocked, made_preferred must stay False and the family's selection
        basis is unresolved_no_target -- a blocked member is never "selected"."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_tie_fixture(root, [
                ("A", "TIE", "a", "tier5_name_only", "paralog_ambiguous", "dep_TIE"),
                ("B", "TIE", "a", "tier5_name_only", "paralog_ambiguous", "dep_TIE"),
            ])
            out = root / "copy"
            tie_frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t",
                                    dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual({"co_preferred_tie"}, set(tie_frame["disambiguation_status"]))
            self.assertEqual({"blocked_tier5"}, set(tie_frame["bridge_promotion_gate"]))
            tie_shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(tie_shortlist))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["bridge"]["blocked_count"])
            self.assertEqual({"TIE": "unresolved_no_target"}, manifest["selection_basis"])

    def _write_tie_fixture(self, root: Path, bridge_records: list[tuple[str, str, str, str, str, str]]) -> None:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     "hyp_tie\tTIE\tA;B")
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                     "A\tchr1\t1000\t3000\tTIE\tDA\tREF_A\nB\tchr1\t5000\t7000\tTIE\tDA\tREF_B")
        _write_table(root, "orthology.tsv",
                     "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                     "A\tfamily_level_only\torthofinder\tOF\tREF_A\tOG_TIE\n"
                     "B\tfamily_level_only\torthofinder\tOF\tREF_B\tOG_TIE")
        _write_table(root, "position.tsv", "gene_id\tposition_class",
                     "A\trecombinant_core\nB\trecombinant_core")
        bridge = root / "bridge.tsv"
        bridge.write_text(_bridge_rows(bridge_records), encoding="utf-8")
        out = root / "copy"
        self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                  "--out-dir", str(out), "--case-id", "demo",
                                  "--atlas", str(root / "atlas.tsv"),
                                  "--orthology", str(root / "orthology.tsv"),
                                  "--position", str(root / "position.tsv"),
                                  "--bridge", str(bridge)]))

    def test_rescue_state_and_copy_selection_loss_are_separate_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), bridge_records=self._main_bridge())
            frame = self._frame(out)
            # RES: both members ride on route a only -> rescue_a_only; the bridge
            # has tier1-4 hits yet family_to_copy could not resolve the family
            # (family_resolution_required), so copy_selection_loss=true -- a
            # different column and vocabulary from route rescue.
            self.assertEqual("rescue_a_only", frame.loc["G_R1", "route_rescue_state"])
            self.assertEqual("rescue_a_only", frame.loc["G_R2", "route_rescue_state"])
            self.assertEqual("true", frame.loc["G_R1", "copy_selection_loss"])
            self.assertEqual("true", frame.loc["G_R2", "copy_selection_loss"])
            self.assertEqual("family_resolution_required", frame.loc["G_R1", "disambiguation_status"])
            # CLV: routes a+b -> multi_route_supported, and a copy WAS selected, so no loss.
            self.assertEqual("multi_route_supported", frame.loc["G_A", "route_rescue_state"])
            self.assertEqual("false", frame.loc["G_A", "copy_selection_loss"])
            # WOX: route b only -> rescue_b_only; its sole member is tier5, so the
            # family has no tier1-4 bridge hit and no loss is reported.
            self.assertEqual("rescue_b_only", frame.loc["G_W", "route_rescue_state"])
            self.assertEqual("false", frame.loc["G_W", "copy_selection_loss"])
            self.assertNotEqual(frame.loc["G_W", "route_rescue_state"], frame.loc["G_W", "copy_selection_loss"])

    def test_manifest_reports_bridge_totals(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), bridge_records=self._main_bridge())
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["bridge"]["supplied"])
            self.assertIsNotNone(manifest["bridge"]["bridge_tsv_sha256"])
            # Only G_W (sole preferred copy on tier5) is blocked.
            self.assertEqual(1, manifest["bridge"]["blocked_count"])
            self.assertEqual(1, manifest["bridge"]["copy_selection_loss_count"])
            self.assertEqual({
                "rescue_a_only": 1, "rescue_b_only": 1, "rescue_c_only": 0, "multi_route_supported": 1,
                "not_evaluated": 0,
            }, manifest["bridge"]["bridge_rescue_counts"])
            self.assertEqual(1, manifest["first_wet_candidate_count"])
            self.assertIsNotNone(manifest["inputs"]["bridge"]["sha256"])
            self.assertIn("blocked_tier5", manifest["bridge"]["basis"])

    def test_cross_family_bridge_rows_never_leak_between_families(self) -> None:
        """Issue #374 T3 F1: the same gene carries a strong (tier1) bridge row
        under one family and a weak (tier5 / tier6) row under another. The
        strong family's tier must never promote the gene under the weak
        family -- the bridge is looked up by (gene_id, family_id), not by
        gene_id alone."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp_ax\tFAMILY_AX\tG_X\nhyp_bx\tFAMILY_BX\tG_X\n"
                         "hyp_ay\tFAMILY_AY\tG_Y\nhyp_by\tFAMILY_BY\tG_Y")
            _write_table(root, "atlas.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                         "G_X\tchr1\t1000\t3000\tX\tLRR-Pkinase\tREF_X\n"
                         "G_Y\tchr1\t5000\t7000\tY\tLRR-Pkinase\tREF_Y")
            _write_table(root, "orthology.tsv",
                         "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                         "G_X\tone_to_one\torthofinder\tOF_chain\tREF_X\tOG_X\n"
                         "G_Y\tone_to_one\torthofinder\tOF_chain\tREF_Y\tOG_Y")
            bridge = root / "bridge.tsv"
            bridge.write_text(_bridge_rows([
                ("G_X", "FAMILY_AX", "a", "tier1_source_backed_unique", "resolved_single_copy", "dep_X_strong"),
                ("G_X", "FAMILY_BX", "a", "tier5_name_only", "paralog_ambiguous", "dep_X_weak"),
                ("G_Y", "FAMILY_AY", "a", "tier1_source_backed_unique", "resolved_single_copy", "dep_Y_strong"),
                ("G_Y", "FAMILY_BY", "a", "tier6_not_evaluated", "not_evaluated", "dep_Y_weak"),
            ]), encoding="utf-8")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                      "--out-dir", str(out), "--case-id", "demo",
                                      "--atlas", str(root / "atlas.tsv"),
                                      "--orthology", str(root / "orthology.tsv"),
                                      "--bridge", str(bridge)]))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("") \
                .set_index(["family_id", "candidate_gene_id"])

            # Both genes are their family's sole member with their own homology
            # layer, so each family's own logic makes them preferred_copy --
            # the bridge gate is the only thing that can still block them.
            self.assertEqual("preferred_copy", frame.loc[("FAMILY_AX", "G_X"), "disambiguation_status"])
            self.assertEqual("preferred_copy", frame.loc[("FAMILY_BX", "G_X"), "disambiguation_status"])
            self.assertEqual("preferred_copy", frame.loc[("FAMILY_AY", "G_Y"), "disambiguation_status"])
            self.assertEqual("preferred_copy", frame.loc[("FAMILY_BY", "G_Y"), "disambiguation_status"])

            # FAMILY_AX / FAMILY_AY: the gene's own (strong) bridge row -> passed.
            self.assertEqual("tier1_source_backed_unique", frame.loc[("FAMILY_AX", "G_X"), "bridge_mapping_tier"])
            self.assertEqual("passed", frame.loc[("FAMILY_AX", "G_X"), "bridge_promotion_gate"])
            self.assertEqual("dep_X_strong", frame.loc[("FAMILY_AX", "G_X"), "bridge_dependency_cluster_id"])
            self.assertEqual("tier1_source_backed_unique", frame.loc[("FAMILY_AY", "G_Y"), "bridge_mapping_tier"])
            self.assertEqual("passed", frame.loc[("FAMILY_AY", "G_Y"), "bridge_promotion_gate"])

            # FAMILY_BX / FAMILY_BY: must read this family's OWN (weak) row --
            # never FAMILY_AX/FAMILY_AY's tier1 leaking across. Before the F1
            # fix these rows read bridge_mapping_tier=tier1_source_backed_unique
            # and bridge_promotion_gate=passed, wrongly promoting them.
            self.assertEqual("tier5_name_only", frame.loc[("FAMILY_BX", "G_X"), "bridge_mapping_tier"])
            self.assertEqual("blocked_tier5", frame.loc[("FAMILY_BX", "G_X"), "bridge_promotion_gate"])
            self.assertEqual("dep_X_weak", frame.loc[("FAMILY_BX", "G_X"), "bridge_dependency_cluster_id"])
            self.assertEqual("tier6_not_evaluated", frame.loc[("FAMILY_BY", "G_Y"), "bridge_mapping_tier"])
            self.assertEqual("blocked_tier6", frame.loc[("FAMILY_BY", "G_Y"), "bridge_promotion_gate"])

            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            promoted = set(zip(shortlist["family_id"], shortlist["candidate_gene_id"]))
            self.assertEqual({("FAMILY_AX", "G_X"), ("FAMILY_AY", "G_Y")}, promoted)

    def test_same_family_multiple_bridge_rows_still_resolve_to_the_best_tier(self) -> None:
        """Regression guard for the F1 fix: a gene with more than one bridge
        row under its OWN family (different routes) still resolves to the
        strongest (lowest) tier claimed within that family, exactly as
        before -- the keyed lookup narrows the scope to (gene_id, family_id),
        it does not stop the existing best-tier-wins tie-break inside that
        scope."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp_z\tFAMILY_Z\tG_Z")
            _write_table(root, "atlas.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                         "G_Z\tchr1\t1000\t3000\tZ\tLRR-Pkinase\tREF_Z")
            _write_table(root, "orthology.tsv",
                         "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                         "G_Z\tone_to_one\torthofinder\tOF_chain\tREF_Z\tOG_Z")
            bridge = root / "bridge.tsv"
            bridge.write_text(_bridge_rows([
                ("G_Z", "FAMILY_Z", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep_Z"),
                ("G_Z", "FAMILY_Z", "b", "tier1_source_backed_unique", "resolved_single_copy", "dep_Z"),
            ]), encoding="utf-8")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                      "--out-dir", str(out), "--case-id", "demo",
                                      "--atlas", str(root / "atlas.tsv"),
                                      "--orthology", str(root / "orthology.tsv"),
                                      "--bridge", str(bridge)]))
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t",
                                dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("tier1_source_backed_unique", frame.loc["G_Z", "bridge_mapping_tier"])
            self.assertEqual("passed", frame.loc["G_Z", "bridge_promotion_gate"])
            self.assertEqual("multi_route_supported", frame.loc["G_Z", "route_rescue_state"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"G_Z"}, set(shortlist["candidate_gene_id"]))

    def test_route_rescue_state_is_not_evaluated_when_every_bridge_row_is_tier6(self) -> None:
        """Issue #374 T3 F2: a family whose only bridge coverage is
        tier6_not_evaluated must not read as rescue_<route>_only /
        multi_route_supported -- that would convert a provider gap into route
        support. It must read as the explicit "not_evaluated" gap state, kept
        apart from true absence (no bridge row at all, which stays "")."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_tie_fixture(root, [
                ("A", "TIE", "a", "tier6_not_evaluated", "not_evaluated", "dep_TIE"),
                ("B", "TIE", "a", "tier6_not_evaluated", "not_evaluated", "dep_TIE"),
            ])
            out = root / "copy"
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t",
                                dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual({"not_evaluated"}, set(frame["route_rescue_state"]))
            self.assertEqual({"co_preferred_tie"}, set(frame["disambiguation_status"]))
            self.assertEqual({"blocked_tier6"}, set(frame["bridge_promotion_gate"]))
            # copy_selection_loss stays false: tier6 is never a source tier, so
            # the gap must not additionally register as a copy-selection loss.
            self.assertEqual({"false"}, set(frame["copy_selection_loss"]))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual({
                "rescue_a_only": 0, "rescue_b_only": 0, "rescue_c_only": 0,
                "multi_route_supported": 0, "not_evaluated": 1,
            }, manifest["bridge"]["bridge_rescue_counts"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(shortlist))


ROUTE_A_BLAST_HEADER = ("panel_id\treview_status\tprovider_status\tquery_sequence_sha256\tquery_length"
                        "\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage"
                        "\ttop_hit_evalue\ttop_hit_bitscore\thit_count\tparalog_count\tambiguity_status"
                        "\tblastp_version\tdb_fingerprint\tinterval_proteins_sha256\trank_contribution\tnote")


def _route_a_blast_rows(records: list[tuple[str, str, str, str, str, str, str, str]]) -> str:
    """Compact blast_evidence writer: (panel_id, review_status, provider_status,
    top_hit_target_id, ambiguity_status, identity, query_coverage, target_coverage)."""
    lines = [ROUTE_A_BLAST_HEADER]
    for panel_id, review_status, provider_status, target_id, ambiguity, identity, query_cov, target_cov in records:
        lines.append("\t".join([panel_id, review_status, provider_status, "sha_query", "300", target_id,
                                identity, query_cov, target_cov, "1e-50", "800", "1", "0", ambiguity,
                                "blastp 2.15.0", "db_fp", "proteins_fp", "none", ""]))
    return "\n".join(lines) + "\n"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FamilyToCopyRouteABlastTests(unittest.TestCase):
    """Issue #374 route A layer: verified unique blast_evidence rows inject a
    route_a_blastp homology layer for the hit member gene; provenance is
    recorded in the copy_set_evidence_source discrimination column (present
    only when the input is supplied)."""

    def _fixture(self, root: Path, *, members: str = "G_A",
                 blast_records: list[tuple[str, str, str, str, str, str, str, str]],
                 bridge_records: list[tuple[str, str, str, str, str, str]] | None = None,
                 atlas_domain: bool = False,
                 atlas_extra_rows: list[str] | None = None) -> Path:
        _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                     f"hyp374ra\tCLV\t{members}")
        atlas_rows = ["G_A\tchr1\t1000\t3000\tCLV" + ("\tDA\tdom1:2;3\tsource_x\t" if atlas_domain else "\t\t\t\t")]
        if "G_B" in members:
            atlas_rows.append("G_B\tchr1\t5000\t7000\tCLV\t\t\t\t")
        if atlas_extra_rows:
            atlas_rows.extend(atlas_extra_rows)
        _write_table(root, "atlas.tsv",
                     "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\tdomain_ids\tdomain_source\torthology_reference_gene",
                     "\n".join(atlas_rows))
        root.joinpath("blast.tsv").write_text(_route_a_blast_rows(blast_records), encoding="utf-8")
        argv = ["family-to-copy", "--families", str(root / "families.tsv"), "--out-dir", str(root / "copy"),
                "--case-id", "demo", "--atlas", str(root / "atlas.tsv"),
                "--route-a-blast", str(root / "blast.tsv")]
        if bridge_records is not None:
            bridge = root / "bridge.tsv"
            bridge.write_text(_bridge_rows(bridge_records), encoding="utf-8")
            argv += ["--bridge", str(bridge)]
        self.assertEqual(0, main(argv))
        return root / "copy"

    def _frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")

    @staticmethod
    def _unique_hit() -> list[tuple[str, str, str, str, str, str, str, str]]:
        return [("B1", "verified", "evaluable", "G_A", "unique", "66.8", "0.9", "0.85")]

    def test_unique_hit_promotes_single_member_family_to_preferred_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("validated", frame.loc["G_A", "copy_set_status"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertIn("route_a_blast:B1", frame.loc["G_A", "copy_set_id"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(shortlist))
            self.assertEqual("G_A", shortlist.iloc[0]["candidate_gene_id"])
            self.assertTrue(shortlist["rank_contribution"].eq("none").all())

    def test_ambiguous_paralog_never_injects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "ambiguous_paralog", "66.8", "0.9", "0.85")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("provisional", frame.loc["G_A", "copy_set_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_no_threshold_hit_never_injects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "no_threshold_hit", "", "", "")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_coverage_below_threshold_never_injects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "unique", "66.8", "0.4", "0.85")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_coverage_fraction_one_point_zero_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "unique", "65.58", "0.966912", "1.0")],
                bridge_records=[("G_A", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertIn("route_a_blast:B1", frame.loc["G_A", "copy_set_id"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(shortlist))

    def test_identity_percent_below_thirty_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "unique", "29.9", "0.966912", "0.85")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_coverage_percent_form_still_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "unique", "66.8", "96.5", "0.85")],
                bridge_records=[("G_A", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertIn("route_a_blast:B1", frame.loc["G_A", "copy_set_id"])

    def test_non_numeric_identity_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "verified", "evaluable", "G_A", "unique", "", "0.966912", "0.85")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_skipped_review_required_row_never_injects(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=[
                ("B1", "skipped_review_required", "evaluable", "G_A", "unique", "66.8", "0.9", "0.85")])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_bridge_tier6_gate_wins_over_route_a_layer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "CLV", "a", "tier6_not_evaluated", "not_evaluated", "dep_CLV")])
            frame = self._frame(out)
            # Issue #398: a tier6-only bridge never counts as guard evaluation,
            # so the homology layer is held back fail-closed before promotion.
            self.assertNotIn(frame.loc["G_A", "disambiguation_status"], {"preferred_copy", "co_preferred_tie"})
            self.assertEqual("paralog_guard_not_evaluated", frame.loc["G_A", "route_a_blast_guard"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_bridge_tier5_gate_wins_over_route_a_layer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "CLV", "a", "tier5_name_only", "paralog_ambiguous", "dep_CLV")])
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("blocked_tier5", frame.loc["G_A", "bridge_promotion_gate"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_in_interval_paralog_holds_unique_hit(self) -> None:
        """Issue #398 acceptance 1: a unique hit whose member has an in-interval
        relative under the same bridge family (route, family_id) is NOT promoted;
        the machine-readable reason is recorded in route_a_blast_guard."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "LYK", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep"),
                ("G_X", "LYK", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep")],
                atlas_extra_rows=["G_X\tchr1\t9000\t11000\tCLV\t\t\t\t"])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("provisional", frame.loc["G_A", "copy_set_status"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual("interval_paralog_present", frame.loc["G_A", "route_a_blast_guard"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["route_a_blast"]["held_interval_paralog_gene_count"])
            self.assertEqual(0, manifest["route_a_blast"]["injected_member_count"])

    def test_tier6_bridge_row_alone_never_counts_as_a_family_assignment(self) -> None:
        """Issue #398 input rule: a tier6_not_evaluated bridge row is a provider
        gap, never a family assignment, so the guard cannot be evaluated and fails closed."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "LYK", "b", "tier6_not_evaluated", "not_evaluated", "dep"),
                ("G_X", "LYK", "b", "tier6_not_evaluated", "not_evaluated", "dep")],
                atlas_extra_rows=["G_X\tchr1\t9000\t11000\tCLV\t\t\t\t"])
            frame = self._frame(out)
            self.assertNotIn(frame.loc["G_A", "disambiguation_status"], {"preferred_copy", "co_preferred_tie"})
            self.assertEqual("paralog_guard_not_evaluated", frame.loc["G_A", "route_a_blast_guard"])
            self.assertEqual("", frame.loc["G_A", "copy_set_evidence_source"])

    def test_unique_hit_without_in_interval_paralog_still_promotes(self) -> None:
        """Issue #398 acceptance 2: a single hit with no in-interval relative
        under the same bridge family keeps promoting (pea_sgr_i/peam1/a-like)."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "SGR", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("validated", frame.loc["G_A", "copy_set_status"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual("", frame.loc["G_A", "route_a_blast_guard"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(shortlist))
            self.assertEqual("G_A", shortlist.iloc[0]["candidate_gene_id"])

    def test_manifest_marks_paralog_guard_evaluated_when_bridge_supplied(self) -> None:
        """Issue #398: with a bridge the guard actually runs, so a reader can
        tell "no in-interval paralog found" apart from "none looked for". The
        no-bridge case is asserted in test_manifest_records_route_a_blast_section."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "CLV", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            section = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))["route_a_blast"]
            self.assertTrue(section["paralog_guard_evaluated"])
            self.assertEqual(1, section["injected_member_count"])

    def test_paralog_guard_sits_below_bridge_gate(self) -> None:
        """Issue #398 acceptance 5: the guard is a layer-internal condition. When
        it fires, a passing bridge gate cannot rescue the row into promotion; the
        gate only ever blocks a promotion the layer actually reached."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "LYK", "b", "tier1_source_backed_unique", "resolved_single_copy", "dep"),
                ("G_X", "LYK", "b", "tier1_source_backed_unique", "resolved_single_copy", "dep")],
                atlas_extra_rows=["G_X\tchr1\t9000\t11000\tCLV\t\t\t\t"])
            frame = self._frame(out)
            self.assertEqual("family_resolution_required", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("interval_paralog_present", frame.loc["G_A", "route_a_blast_guard"])
            # No promotion point was reached, so there is nothing for the gate to gate.
            self.assertEqual("", frame.loc["G_A", "bridge_promotion_gate"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_paralog_guard_ignores_relatives_inside_the_family(self) -> None:
        """Issue #398: a relative that IS a family member is already covered by
        the component logic (a singleton can never promote a multi-member
        family), so the guard leaves the layer and its provenance intact there."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), members="G_A;G_X", blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "LYK", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep"),
                ("G_X", "LYK", "b", "tier3_sequence_match_unique", "resolved_single_copy", "dep")],
                atlas_extra_rows=["G_X\tchr1\t9000\t11000\tCLV\t\t\t\t"])
            frame = self._frame(out)
            self.assertEqual("not_in_validated_copy_set", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual("", frame.loc["G_A", "route_a_blast_guard"])

    def test_multi_member_family_singleton_hit_stays_not_in_validated_copy_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), members="G_A;G_B", blast_records=self._unique_hit(), bridge_records=[
                ("G_A", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep"),
                ("G_B", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep")])
            frame = self._frame(out)
            self.assertEqual("not_in_validated_copy_set", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("not_applicable", frame.loc["G_A", "copy_set_status"])
            self.assertIn("route_a_blast:B1", frame.loc["G_A", "copy_set_id"])
            self.assertEqual("route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])
            self.assertEqual("family_resolution_required", frame.loc["G_B", "disambiguation_status"])
            self.assertEqual(0, len(pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str)))

    def test_advisory_rows_only_for_families_containing_hit_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_table(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp374a\tFAM1\tG_A;G_B\n"
                         "hyp374b\tFAM2\tG_C;G_D")
            _write_table(root, "atlas.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\tdomain_ids\tdomain_source\torthology_reference_gene",
                         "\n".join([
                             "G_A\tchr1\t1000\t3000\tFAM1\t\t\t\t",
                             "G_B\tchr1\t5000\t7000\tFAM1\t\t\t\t",
                             "G_C\tchr2\t1000\t3000\tFAM2\t\t\t\t",
                             "G_D\tchr2\t5000\t7000\tFAM2\t\t\t\t",
                         ]))
            _write_table(root, "bridge.tsv", "gene_id\tfamily_id\troute\tmapping_tier",
                         "G_A\tFAM1\ta\ttier3_sequence_match_unique\n"
                         "G_B\tFAM1\ta\ttier3_sequence_match_unique\n"
                         "G_C\tFAM2\ta\ttier3_sequence_match_unique\n"
                         "G_D\tFAM2\ta\ttier3_sequence_match_unique")
            root.joinpath("blast.tsv").write_text(_route_a_blast_rows(self._unique_hit()), encoding="utf-8")
            out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                      "--out-dir", str(out), "--case-id", "demo",
                                      "--atlas", str(root / "atlas.tsv"),
                                      "--bridge", str(root / "bridge.tsv"),
                                      "--route-a-blast", str(root / "blast.tsv")]))
            reviews = pd.read_csv(out / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")
            mixed = reviews[reviews["reason"].eq("family_bucket_mixes_homologous_families")]
            self.assertEqual(["FAM1"], sorted(set(mixed["target"])))
            self.assertNotIn("FAM2", set(mixed["target"]))
            self.assertEqual("FAM1", mixed.iloc[0]["target"])

    def test_atlas_and_route_a_sources_are_joined_and_sorted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), blast_records=self._unique_hit(),
                                bridge_records=[("G_A", "CLV", "a", "tier3_sequence_match_unique", "resolved_single_copy", "dep")],
                                atlas_domain=True)
            frame = self._frame(out)
            self.assertEqual("preferred_copy", frame.loc["G_A", "disambiguation_status"])
            self.assertEqual("atlas_domain;route_a_blastp", frame.loc["G_A", "copy_set_evidence_source"])

    def test_without_route_a_blast_legacy_schema_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = self._fixture(root, blast_records=[("B1", "verified", "evaluable", "G_A", "no_threshold_hit", "", "", "")])
            _write_table(root, "families2.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp374ra\tCLV\tG_A")
            _write_table(root, "atlas2.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\tdomain_ids\tdomain_source\torthology_reference_gene",
                         "G_A\tchr1\t1000\t3000\tCLV\t\t\t\t")
            legacy = root / "legacy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families2.tsv"),
                                      "--out-dir", str(legacy), "--case-id", "demo",
                                      "--atlas", str(root / "atlas2.tsv")]))
            for name in ("family_copy_members.tsv", "family_copy_relationships.tsv",
                         "family_copy_shortlist.tsv", "family_copy_review_required.tsv"):
                self.assertEqual(_sha256(legacy / name), _sha256(base / name), name)
            legacy_discrimination = pd.read_csv(legacy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("copy_set_evidence_source", legacy_discrimination.columns)
            self.assertNotIn("route_a_blast_guard", legacy_discrimination.columns)
            supplied_discrimination = pd.read_csv(base / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("copy_set_evidence_source", supplied_discrimination.columns)
            self.assertIn("route_a_blast_guard", supplied_discrimination.columns)
            self.assertTrue(supplied_discrimination["copy_set_evidence_source"].eq("").all())
            self.assertTrue(supplied_discrimination["route_a_blast_guard"].eq("").all())
            pd.testing.assert_frame_equal(
                supplied_discrimination.drop(columns=["copy_set_evidence_source", "route_a_blast_guard"]),
                legacy_discrimination)

    def test_copy_set_evidence_source_values_on_canonical_vocabulary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._fixture(Path(tmp), members="G_A;G_B",
                                blast_records=self._unique_hit(), atlas_domain=True)
            frame = pd.read_csv(out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            for cell in frame["copy_set_evidence_source"]:
                for token in cell.split(";"):
                    self.assertIn(token, schemas.COPY_SET_EVIDENCE_SOURCES | {""})
            # The fail-closed gate itself rejects an off-vocabulary token.
            with self.assertRaises(schemas.SchemaValidationError):
                family_to_copy._validate_copy_set_evidence_sources(
                    [{"copy_set_evidence_source": "route_a_blastp;bogus_source"}])
            # The issue #398 guard column is on the canonical vocabulary too.
            for cell in frame["route_a_blast_guard"]:
                self.assertIn(cell, schemas.ROUTE_A_BLAST_GUARD_REASONS | {""})
            with self.assertRaises(schemas.SchemaValidationError):
                family_to_copy._validate_route_a_blast_guard(
                    [{"route_a_blast_guard": "interval_paralog_present;bogus"}])

    def test_manifest_records_route_a_blast_section(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._fixture(root, blast_records=self._unique_hit())
            manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
            section = manifest["route_a_blast"]
            self.assertTrue(section["supplied"])
            self.assertIsNotNone(section["blast_evidence_tsv_sha256"])
            self.assertEqual(1, section["evidence_rows"])
            self.assertEqual(0, section["eligible_unique_rows"])
            self.assertEqual(0, section["injected_member_count"])
            self.assertEqual(0, section["held_interval_paralog_gene_count"])
            self.assertEqual(1, section["held_guard_not_evaluated_gene_count"])
            # Issue #398: this run has no bridge, so the paralog guard could not
            # run and the promotion was held back fail-closed.
            self.assertFalse(section["paralog_guard_evaluated"])
            self.assertIsNotNone(manifest["inputs"]["route_a_blast"])
            self.assertIn("route_a_blastp", section["basis"])
            legacy = root / "legacy"
            _write_table(root, "families2.tsv", "hypothesis_id\tgene_family\tmember_interval_genes",
                         "hyp374ra\tCLV\tG_A")
            _write_table(root, "atlas2.tsv",
                         "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\tdomain_ids\tdomain_source\torthology_reference_gene",
                         "G_A\tchr1\t1000\t3000\tCLV\t\t\t\t")
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families2.tsv"),
                                      "--out-dir", str(legacy), "--case-id", "demo",
                                      "--atlas", str(root / "atlas2.tsv")]))
            legacy_manifest = json.loads((legacy / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(legacy_manifest["route_a_blast"]["supplied"])
            self.assertEqual(0, legacy_manifest["route_a_blast"]["evidence_rows"])
            self.assertEqual(0, legacy_manifest["route_a_blast"]["eligible_unique_rows"])
            self.assertEqual(0, legacy_manifest["route_a_blast"]["injected_member_count"])
            self.assertEqual(0, legacy_manifest["route_a_blast"]["held_interval_paralog_gene_count"])
            self.assertIsNone(legacy_manifest["route_a_blast"]["blast_evidence_tsv_sha256"])
            self.assertIsNone(legacy_manifest["inputs"]["route_a_blast"])


if __name__ == "__main__":
    unittest.main()
