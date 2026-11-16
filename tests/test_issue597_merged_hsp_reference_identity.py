"""Issue #597: the Route B homology link is judged on merged HSPs.

A homolog reported as several HSPs otherwise fails a single-HSP coverage bar,
the interval gene reaches Route B with no reference identity, and the curated
function of its best homolog never becomes semantic review material -- the
semantic stage then sees a generic annotation borrowed from some other member
of the coarse family bucket. This is the Route B twin of #567.
"""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from candidate_qtl.route_b_prepare import _run_blastp

# Extended outfmt includes qseq/sseq before stitle.
# Synthetic alignments below test coverage mechanics, not biological sequences.
LYK3 = "sp|Q6UD73|LYK3_MEDTR LysM domain receptor-like kinase 3 OS=Medicago truncatula GN=LYK3"
OTHER = "tr|G7K7M2|G7K7M2_MEDTR non-specific serine/threonine protein kinase"


def _row(query, subject, identity, qstart, qend, bitscore, title, qlen=1192):
    span = qend - qstart + 1
    matches = round(span * identity / 100)
    qseq = "A" * span
    sseq = "A" * matches + "C" * (span - matches)
    return "\t".join(str(value) for value in (
        query, subject, 100 * matches / span, span, qlen, qlen, qstart, qend, 1, span,
        "0.0", bitscore, qseq, sseq, title))


class MergedHspReferenceIdentityTests(unittest.TestCase):
    def _run(self, rows, *, min_identity=30.0, min_coverage=0.5):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)

            def fake_run(command, **kwargs):
                Path(command[command.index("-out") + 1]).write_text(
                    "\n".join(rows) + "\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")

            with patch("candidate_qtl.route_b_prepare.shutil.which", return_value="/usr/bin/blastp"), \
                    patch("candidate_qtl.route_b_prepare.run_tool", side_effect=fake_run):
                records, meta = _run_blastp(
                    out_dir / "proteins.fa", out_dir / "db", out_dir, 1,
                    min_identity, min_coverage, 1e-5)
            self.all_records = records
            return {record["gene_id"]: record for record in records}, meta

    def test_split_hsps_recover_the_reference_identity(self):
        """The pea K1 witness: 2 HSPs of 0.49 each against LYK3."""
        rows = [
            _row("Psat.cameor.v2.2g244500", "sp|Q6UD73|LYK3_MEDTR", 84.274, 1, 584, 900, LYK3),
            _row("Psat.cameor.v2.2g244500", "sp|Q6UD73|LYK3_MEDTR", 66.277, 585, 1173, 700, LYK3),
        ]
        records, meta = self._run(rows)
        record = records["Psat.cameor.v2.2g244500"]
        self.assertEqual("Q6UD73", record["reference_gene_id"])
        self.assertGreaterEqual(float(record["query_coverage"]), 0.98)
        self.assertEqual("best_hit_homology", record["relationship"])
        self.assertEqual("merged_hsps_per_query_subject_pair", meta["query_coverage_basis"])

    def test_overlapping_hsps_do_not_manufacture_coverage(self):
        rows = [
            _row("gene_overlap", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 500, 900, LYK3),
            _row("gene_overlap", "sp|Q6UD73|LYK3_MEDTR", 80.0, 20, 520, 800, LYK3),
        ]
        records, _ = self._run(rows)
        self.assertEqual({}, records)

    def test_hsps_of_different_subjects_are_not_pooled(self):
        rows = [
            _row("gene_two_subjects", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 500, 900, LYK3),
            _row("gene_two_subjects", "tr|G7K7M2|G7K7M2_MEDTR", 84.0, 600, 1100, 800, OTHER),
        ]
        records, _ = self._run(rows)
        self.assertEqual({}, records)

    def test_eligible_references_are_all_preserved(self):
        """Both references remain homology material, never sole-ortholog proof."""
        rows = [
            _row("gene_two_hits", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 900, 900, LYK3),
            _row("gene_two_hits", "tr|G7K7M2|G7K7M2_MEDTR", 40.0, 1, 900, 300, OTHER),
        ]
        records, _ = self._run(rows)
        self.assertEqual({"Q6UD73", "G7K7M2"},
                         {row["reference_gene_id"] for row in self.all_records})

    def test_single_hsp_behaviour_is_unchanged(self):
        rows = [_row("gene_single", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 700, 900, LYK3)]
        records, _ = self._run(rows)
        self.assertAlmostEqual(float(records["gene_single"]["query_coverage"]), 700 / 1192, places=4)

    def test_identity_threshold_is_unchanged(self):
        """Merging coverage must not relax the identity bar."""
        rows = [
            _row("gene_low_identity", "sp|Q6UD73|LYK3_MEDTR", 20.0, 1, 584, 900, LYK3),
            _row("gene_low_identity", "sp|Q6UD73|LYK3_MEDTR", 19.0, 585, 1173, 700, LYK3),
        ]
        records, _ = self._run(rows)
        self.assertEqual({}, records)


if __name__ == "__main__":
    unittest.main()
