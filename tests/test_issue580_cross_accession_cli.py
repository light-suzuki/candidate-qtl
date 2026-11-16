"""Issue #580: the #524 current-A17 -> R108 bridge must be CLI-reachable."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


A17 = "medtr.A17.gnm5.ann1_6.CUR1"
R108 = "medtr.R108.gnmHiC_1.ann1.R1"


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class CrossAccessionAssemblyBridgeCliTests(unittest.TestCase):
    def test_cross_accession_scope_is_reachable_without_changing_identity_rules(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            functions = _write(
                root,
                "functions.tsv",
                (
                    "native_gene_id\tnative_species\tnative_accession\tnative_assembly\t"
                    "native_annotation\tnative_gene_id_namespace\tphenotype_or_function\t"
                    "source_ids\n"
                    f"{A17}\tMedicago truncatula\tA17\tgnm5\tann1_6\t"
                    "medtr.A17.gnm5.ann1_6\troot function\tPMID:1\n"
                ),
            )
            source_fam3 = _write(
                root,
                "source.tsv",
                f"native_gene_id\tfamily_id\n{A17}\tF1\n",
            )
            target_fam3 = _write(
                root,
                "target.tsv",
                f"native_gene_id\tfamily_id\n{R108}\tF1\n",
            )
            sequence = _write(
                root,
                "sequence.tsv",
                (
                    "source_full_id\ttarget_full_id\treciprocal_sequence_state\t"
                    "sequence_identity\tquery_coverage\tsubject_coverage\n"
                    f"{A17}\t{R108}\trbh\t98.0\t1.0\t1.0\n"
                ),
            )
            out = root / "bridge"

            result = main([
                "assembly-bridge",
                "--functional-evidence", str(functions),
                "--source-fam3", str(source_fam3),
                "--target-fam3", str(target_fam3),
                "--sequence-evidence", str(sequence),
                "--source-namespace", "medtr.A17.gnm5.ann1_6",
                "--target-namespace", "medtr.R108.gnmHiC_1.ann1",
                "--mapping-scope", "cross_accession",
                "--case-id", "issue580",
                "--out-dir", str(out),
            ])
            self.assertEqual(0, result)

            manifest = json.loads((out / "assembly_bridge_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("cross_accession", manifest["mapping_scope"])
            self.assertEqual(1, manifest["counts"]["mapping_verified_one_to_one"])
            self.assertEqual("none", manifest["rank_contribution"])

            mapping = pd.read_csv(out / "assembly_bridge_mapping.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(mapping))
            row = mapping.iloc[0]
            self.assertEqual(A17, row.source_full_id)
            self.assertEqual(R108, row.target_full_id)
            self.assertEqual("cross_accession", row.mapping_scope)
            self.assertEqual("one_to_one", row.relationship)
            self.assertEqual("verified", row.mapping_state)
            self.assertIn("fam3_concordance", row.mapping_basis)
            self.assertIn("pairwise_sequence", row.mapping_basis)


if __name__ == "__main__":
    unittest.main()
