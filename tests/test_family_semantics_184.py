"""Issue #184: discovery bucket / homologous family / copy-set separation.

The family-index's coarse family value is a *discovery bucket* (naming /
broad functional search). It must never, on its own, define the copy set that
family-to-copy discriminates over: symbol stems and pathway component names do
not guarantee a homologous family. These tests run the chain
``family-index -> interval-atlas -> family-to-copy`` on the issue's fixtures.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


def _hit(query: str, subject: str, identity: float, bitscore: float, *, qlen: int = 300) -> str:
    # outfmt 6: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
    return f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}"


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def _genes(root: Path, records: str) -> Path:
    return _write(root, "genes.tsv", "gene_id\tchrom\tstart\tend\n" + records)


def _blast(root: Path, hits: list[str]) -> Path:
    return _write(root, "hits.tsv", "\n".join(hits) + "\n")


def _candidates(root: Path, records: str) -> Path:
    return _write(root, "candidates.tsv", "gene_id\tchrom\tstart\tend\tprotein_length\n" + records)


def _run_family_index(root: Path, family_map: str | None = None, out_name: str = "index") -> Path:
    out = root / out_name
    argv = ["family-index", "--blast-hits", str(root / "hits.tsv"), "--genes", str(root / "genes.tsv"),
            "--out-dir", str(out), "--reference-label", "TestRef"]
    if family_map:
        argv += ["--family-map", str(root / family_map)]
    assert main(argv) == 0
    return out


def _run_interval_atlas(root: Path, index: Path, out_name: str = "atlas",
                        domains: Path | None = None, orthology: Path | None = None) -> Path:
    out = root / out_name
    argv = ["interval-atlas", "--candidates", str(root / "candidates.tsv"),
            "--family-index", str(index / "family_index.tsv"), "--out-dir", str(out)]
    if domains:
        argv += ["--domains", str(domains)]
    if orthology:
        argv += ["--orthology", str(orthology)]
    assert main(argv) == 0
    return out


def _run_family_to_copy(root: Path, families: str, atlas: Path, out_name: str = "copy",
                        orthology: Path | None = None) -> Path:
    out = root / out_name
    argv = ["family-to-copy", "--families", str(root / families), "--atlas", str(atlas / "interval_gene_atlas.tsv"),
            "--out-dir", str(out)]
    if orthology:
        argv += ["--orthology", str(orthology)]
    assert main(argv) == 0
    return out


def _families(root: Path, records: str) -> Path:
    return _write(root, "families.tsv",
                  "hypothesis_id\tgene_family\tmember_interval_genes\n" + records)


def _orthology(root: Path, records: str) -> Path:
    return _write(root, "orthology.tsv",
                  "gene_id\torthogroup\trelationship\treference_gene_id\tprovider\tevidence_chain\n" + records)


def _domains(root: Path, records: str) -> Path:
    return _write(root, "domains.tsv",
                  "gene_id\tdomain_id\tdomain_name\tstart\tsource\n" + records)


class SameStemNonHomologousTests(unittest.TestCase):
    """Issue fixture 1: same symbol stem, different domains / orthogroups."""

    def test_same_stem_genes_are_split_by_homology_not_compared_as_paralogs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _blast(root, [
                _hit("G1", "sp|P0CLV1|CLV1_ARATH", 62.0, 400),
                _hit("G2", "sp|P0CLV2|CLV2_ARATH", 60.0, 380),
                _hit("G3", "sp|P0CLV3|CLV3_ARATH", 58.0, 360),
            ])
            _genes(root, "G1\tchr1\t1000\t3000\nG2\tchr1\t4000\t6000\nG3\tchr1\t7000\t9000\n")
            _candidates(root, "G1\tchr1\t1000\t3000\t900\nG2\tchr1\t4000\t6000\t850\nG3\tchr1\t7000\t9000\t300\n")
            _domains(root, "G1\tPF00560\tLRR\t50\tPfam\nG1\tPF00069\tPkinase\t600\tPfam\n"
                           "G2\tPF00400\tWD40\t100\tPfam\nG3\tPF00557\tPeptidase\t30\tPfam\n")
            _orthology(root, "G1\tOG10\tone_to_one\tCLV1\torthofinder\tOF_chain\n"
                             "G2\tOG11\tone_to_one\tCLV2\torthofinder\tOF_chain\n"
                             "G3\tOG12\tone_to_one\tCLV3\torthofinder\tOF_chain\n")
            _families(root, "hyp_CLV\tCLV\tG1;G2;G3\n")

            index = _run_family_index(root)
            family = pd.read_csv(index / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # The naming bucket groups them together...
            self.assertEqual({"CLV", "CLV", "CLV"}, set(family["discovery_bucket_id"]))
            self.assertEqual({"symbol_stem"}, set(family["discovery_bucket_method"]))
            # ...while the evolutionary family layer is unresolved for all of them.
            self.assertTrue((family["evolutionary_family_id"] == "").all())
            self.assertTrue((family["evolutionary_family_status"] == "unresolved").all())

            atlas_dir = _run_interval_atlas(root, index, domains=root / "domains.tsv", orthology=root / "orthology.tsv")
            atlas = pd.read_csv(atlas_dir / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # Same discovery bucket, but each gene has its own orthogroup / domain family.
            self.assertEqual({"CLV"}, set(atlas["discovery_bucket_id"]))
            self.assertEqual({"OG10", "OG11", "OG12"}, set(atlas["orthogroup"]))
            self.assertEqual(3, atlas["copy_set_id"].nunique())

            copy = _run_family_to_copy(root, "families.tsv", atlas_dir, orthology=root / "orthology.tsv")
            frame = pd.read_csv(copy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # No member is inside a validated copy set with the others: they are
            # never compared as paralogs and no preferred copy exists.
            self.assertTrue(frame["disambiguation_status"].eq("not_in_validated_copy_set").all())
            self.assertTrue(frame["copy_set_status"].eq("not_applicable").all())
            self.assertEqual(0, len(pd.read_csv(copy / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(copy / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("family_bucket_mixes_homologous_families", reasons)
            manifest = json.loads((copy / "family_copy_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["CLV"], manifest["homology_gate"]["families_mixing_homologous_families"])


class SameFamilyDifferentSymbolsTests(unittest.TestCase):
    """Issue fixture 2: different symbols, same orthogroup -> one copy set."""

    def test_different_symbol_stems_join_the_same_validated_copy_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _blast(root, [
                _hit("G4", "sp|P0NTH1|NTH1_ARATH", 65.0, 420),
                _hit("G5", "sp|P0ATHB|ATHB15_ARATH", 63.0, 400),
            ])
            _genes(root, "G4\tchr1\t1000\t3000\nG5\tchr1\t4000\t6000\n")
            _candidates(root, "G4\tchr1\t1000\t3000\t400\nG5\tchr1\t4000\t6000\t420\n")
            _domains(root, "G4\tPF00046\tHomeobox\t80\tPfam\nG5\tPF00046\tHomeobox\t80\tPfam\n")
            _orthology(root, "G4\tOG777\tone_to_one\tNTH1\torthofinder\tOF_chain\n"
                             "G5\tOG777\tone_to_many\tATHB15\torthofinder\tOF_chain\n")
            _families(root, "hyp_HDZIP\tHD-ZIP\tG4;G5\n")

            index = _run_family_index(root)
            family = pd.read_csv(index / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # Symbol stems differ, so the discovery buckets differ...
            self.assertEqual("NTH", family.loc["G4", "discovery_bucket_id"])
            self.assertEqual("ATHB", family.loc["G5", "discovery_bucket_id"])
            atlas_dir = _run_interval_atlas(root, index, domains=root / "domains.tsv", orthology=root / "orthology.tsv")
            atlas = pd.read_csv(atlas_dir / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # ...but the orthogroup is the same, so the copy set is one.
            self.assertEqual({"OG777"}, set(atlas["orthogroup"]))
            self.assertEqual({"OG777"}, set(atlas["copy_set_id"]))

            copy = _run_family_to_copy(root, "families.tsv", atlas_dir, orthology=root / "orthology.tsv")
            frame = pd.read_csv(copy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # The shared homology layers (orthogroup, shared domain architecture)
            # identify one copy set, regardless of the differing symbol stems.
            self.assertIn("OG777", frame.loc["G4", "copy_set_id"].split(";"))
            self.assertIn("OG777", frame.loc["G5", "copy_set_id"].split(";"))
            self.assertEqual(frame.loc["G4", "copy_set_id"], frame.loc["G5", "copy_set_id"])
            self.assertTrue(frame["copy_set_status"].eq("validated").all())
            # Copy discrimination runs: the one_to_one copy is preferred.
            self.assertEqual("preferred_copy", frame.loc["G4", "disambiguation_status"])
            shortlist = pd.read_csv(copy / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"G4"}, set(shortlist["candidate_gene_id"]))


class FamilyMapAbsentTests(unittest.TestCase):
    """Issue fixture 3: symbol stem only -> provisional, no preferred copy."""

    def test_bucket_only_family_stays_provisional_and_makes_no_preferred_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _blast(root, [_hit("G6", "sp|P0CLV1|CLV1_ARATH", 62.0, 400)])
            _genes(root, "G6\tchr1\t1000\t3000\n")
            _candidates(root, "G6\tchr1\t1000\t3000\t900\n")
            _families(root, "hyp_CLV\tCLV\tG6\n")

            index = _run_family_index(root)
            family = pd.read_csv(index / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("CLV", family.loc["G6", "discovery_bucket_id"])
            self.assertEqual("symbol_stem", family.loc["G6", "discovery_bucket_method"])
            self.assertEqual("unresolved", family.loc["G6", "evolutionary_family_status"])

            atlas_dir = _run_interval_atlas(root, index)
            atlas = pd.read_csv(atlas_dir / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("", atlas.loc["G6", "evolutionary_family_id"])
            self.assertEqual("CLV", atlas.loc["G6", "copy_set_id"])
            self.assertEqual("provisional", atlas.loc["G6", "copy_set_status"])

            copy = _run_family_to_copy(root, "families.tsv", atlas_dir)
            frame = pd.read_csv(copy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # A single-member bucket-only family must not produce a preferred copy.
            self.assertEqual("family_resolution_required", frame.loc["G6", "disambiguation_status"])
            self.assertEqual("provisional", frame.loc["G6", "copy_set_status"])
            self.assertEqual(0, len(pd.read_csv(copy / "family_copy_shortlist.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(copy / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("family_resolution_required", reasons)


class FamilyMapPresentTests(unittest.TestCase):
    """Issue fixture 4: curated family map -> copy discrimination may proceed."""

    def test_curated_family_enables_copy_discrimination(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _blast(root, [
                _hit("G7", "sp|P0CLV1|CLV1_ARATH", 62.0, 400),
                _hit("G8", "sp|P0CLV2|CLV2_ARATH", 60.0, 380),
            ])
            _genes(root, "G7\tchr1\t1000\t3000\nG8\tchr1\t4000\t6000\n")
            _candidates(root, "G7\tchr1\t1000\t3000\t900\nG8\tchr1\t4000\t6000\t850\n")
            _write(root, "family_map.tsv",
                   "subject_id\tfamily_id\tsymbol\n"
                   "sp|P0CLV1|CLV1_ARATH\tLRR-RLK_CLV\tCLV1\n"
                   "sp|P0CLV2|CLV2_ARATH\tLRR-RLK_CLV\tCLV2\n")
            _orthology(root, "G7\tOG1\tone_to_one\tCLV1\torthofinder\tOF_chain\n"
                             "G8\tOG1\tone_to_many\tCLV2\torthofinder\tOF_chain\n")
            _families(root, "hyp_LRR\tLRR-RLK_CLV\tG7;G8\n")

            index = _run_family_index(root, family_map="family_map.tsv")
            family = pd.read_csv(index / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("curated_keyword", family.loc["G7", "discovery_bucket_method"])
            self.assertEqual("LRR-RLK_CLV", family.loc["G7", "evolutionary_family_id"])
            self.assertEqual("curated_family", family.loc["G7", "evolutionary_family_source"])
            self.assertEqual("supported", family.loc["G7", "evolutionary_family_status"])

            atlas_dir = _run_interval_atlas(root, index, orthology=root / "orthology.tsv")
            atlas = pd.read_csv(atlas_dir / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual({"LRR-RLK_CLV"}, set(atlas["evolutionary_family_id"]))
            self.assertEqual({"validated"}, set(atlas["copy_set_status"]))

            copy = _run_family_to_copy(root, "families.tsv", atlas_dir, orthology=root / "orthology.tsv")
            frame = pd.read_csv(copy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # The curated family validates the copy set: discrimination runs and
            # the orthology-leading copy can be preferred.
            self.assertIn("LRR-RLK_CLV", frame.loc["G7", "copy_set_id"].split(";"))
            self.assertEqual(frame.loc["G7", "copy_set_id"], frame.loc["G8", "copy_set_id"])
            self.assertTrue(frame["copy_set_status"].eq("validated").all())
            self.assertEqual("preferred_copy", frame.loc["G7", "disambiguation_status"])


class ClvStyleCollisionTests(unittest.TestCase):
    """Issue fixture 5: same CLV prefix, different protein classes -> not copies."""

    def test_clv3_peptide_is_not_a_paralog_of_clv1_clv2_receptor_kinases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _blast(root, [
                _hit("G9", "sp|P0CLV1|CLV1_ARATH", 62.0, 400),
                _hit("G10", "sp|P0CLV2|CLV2_ARATH", 60.0, 380),
                _hit("G11", "sp|P0CLV3|CLV3_ARATH", 58.0, 360),
            ])
            _genes(root, "G9\tchr1\t1000\t3000\nG10\tchr1\t4000\t6000\nG11\tchr1\t7000\t9000\n")
            _candidates(root, "G9\tchr1\t1000\t3000\t900\nG10\tchr1\t4000\t6000\t850\nG11\tchr1\t7000\t9000\t100\n")
            # G9/G10 are LRR receptor kinases; G11 is a small peptide without a kinase.
            _domains(root, "G9\tPF00560\tLRR\t50\tPfam\nG9\tPF00069\tPkinase\t600\tPfam\n"
                           "G10\tPF00560\tLRR\t50\tPfam\nG10\tPF00069\tPkinase\t600\tPfam\n")
            _orthology(root, "G9\tOG1\tone_to_one\tCLV1\torthofinder\tOF_chain\n"
                             "G10\tOG1\tone_to_many\tCLV2\torthofinder\tOF_chain\n"
                             "G11\tOG2\tmany_to_many\tCLV3\torthofinder\tOF_chain\n")
            _families(root, "hyp_CLV\tCLV\tG9;G10;G11\n")

            index = _run_family_index(root)
            family = pd.read_csv(index / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # The naming collision is real at the bucket level: one shared stem.
            self.assertEqual({"CLV"}, set(family["discovery_bucket_id"]))

            atlas_dir = _run_interval_atlas(root, index, domains=root / "domains.tsv", orthology=root / "orthology.tsv")
            copy = _run_family_to_copy(root, "families.tsv", atlas_dir, orthology=root / "orthology.tsv")
            frame = pd.read_csv(copy / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # CLV1/CLV2 (same orthogroup) are compared as copies; CLV3 is not.
            self.assertEqual("preferred_copy", frame.loc["G9", "disambiguation_status"])
            self.assertIn("OG1", frame.loc["G9", "copy_set_id"].split(";"))
            self.assertNotIn("OG2", frame.loc["G9", "copy_set_id"].split(";"))
            self.assertEqual("not_in_validated_copy_set", frame.loc["G11", "disambiguation_status"])
            self.assertEqual("not_applicable", frame.loc["G11", "copy_set_status"])
            shortlist = pd.read_csv(copy / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("G11", set(shortlist["candidate_gene_id"]))
            reasons = set(pd.read_csv(copy / "family_copy_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("family_bucket_mixes_homologous_families", reasons)


if __name__ == "__main__":
    unittest.main()
