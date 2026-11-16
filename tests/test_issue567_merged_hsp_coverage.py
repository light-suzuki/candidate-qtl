"""Issue #567: family-index coverage is judged on the union of a pair's HSPs.

A homolog reported as several HSPs (duplicated domains, a gap in the gene
model, a segmented alignment) must not fall out of the index because no single
HSP clears the coverage threshold. Merging must not invent coverage either: two
HSPs over the same region stay below the threshold.
"""
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.family_index import family_index

# outfmt 6 + stitle: qseqid sseqid pident length qlen slen qstart qend sstart
#                    send evalue bitscore stitle
HEADER = "{query}\t{subject}\t{identity}\t{length}\t1192\t620\t{qs}\t{qe}\t1\t584\t0.0\t900\t{title}"
LYK3 = "sp|Q6UD73|LYK3_MEDTR LysM domain receptor-like kinase 3 OS=Medicago truncatula GN=LYK3"


def _row(query, subject, identity, qstart, qend, title):
    return HEADER.format(query=query, subject=subject, identity=identity,
                         length=qend - qstart + 1, qs=qstart, qe=qend, title=title)


class MergedHspCoverageTest(unittest.TestCase):
    def _build(self, hit_lines, genes):
        directory = Path(self.tmp.name)
        hits = directory / "hits.tsv"
        hits.write_text("\n".join(hit_lines) + "\n", encoding="utf-8")
        gene_table = directory / "genes.tsv"
        gene_table.write_text(
            "gene_id\tchrom\tstart\tend\n" + "".join(f"{gene}\tchr2\t1\t100\n" for gene in genes),
            encoding="utf-8")
        out = directory / "out"
        result = family_index(str(hits), str(out), genes_arg=str(gene_table),
                              min_identity=30.0, min_query_coverage=0.5, case_id="issue567")
        index = (out / "family_index.tsv").read_text(encoding="utf-8").splitlines()
        manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
        return result, index, manifest

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_split_hsps_are_merged_and_the_gene_is_assigned(self):
        """The pea K1 witness: 2 HSPs of 0.49 each, 84% identity, one family."""
        rows = [
            _row("Psat.cameor.v2.2g244500", "sp|Q6UD73|LYK3_MEDTR", 84.274, 1, 584, LYK3),
            _row("Psat.cameor.v2.2g244500", "sp|Q6UD73|LYK3_MEDTR", 66.277, 585, 1173, LYK3),
        ]
        _, index, manifest = self._build(rows, ["Psat.cameor.v2.2g244500"])
        header, line = index[0].split("\t"), index[1].split("\t")
        record = dict(zip(header, line))
        self.assertEqual(record["gene_id"], "Psat.cameor.v2.2g244500")
        self.assertEqual(record["assignment_status"], "assigned")
        self.assertEqual(record["best_subject"], "sp|Q6UD73|LYK3_MEDTR")
        self.assertGreaterEqual(float(record["best_query_coverage"]), 0.98)
        self.assertEqual(manifest["dark_count"], 0)
        self.assertEqual(manifest["exclusion_breakdown"]["genes_below_min_query_coverage"], 0)
        self.assertEqual(
            manifest["exclusion_breakdown"]["genes_recovered_by_merged_hsp_coverage"],
            ["Psat.cameor.v2.2g244500"])
        self.assertEqual(manifest["quality_thresholds"]["query_coverage_basis"],
                         "merged_hsps_per_query_subject_pair")

    def test_overlapping_hsps_do_not_manufacture_coverage(self):
        """Two HSPs over the same region stay below the threshold."""
        rows = [
            _row("gene_overlap", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 500, LYK3),
            _row("gene_overlap", "sp|Q6UD73|LYK3_MEDTR", 80.0, 20, 520, LYK3),
        ]
        _, index, manifest = self._build(rows, ["gene_overlap"])
        record = dict(zip(index[0].split("\t"), index[1].split("\t")))
        self.assertEqual(record["assignment_status"], "dark_no_qualifying_hit")
        self.assertEqual(manifest["exclusion_breakdown"]["genes_below_min_query_coverage"], 1)
        self.assertEqual(manifest["exclusion_breakdown"]["genes_recovered_by_merged_hsp_coverage"], [])

    def test_hsps_of_different_subjects_are_not_pooled(self):
        """Coverage is per (query, subject); unrelated subjects never combine."""
        other = "sp|P00000|OTHER_MEDTR Unrelated protein OS=Medicago truncatula GN=OTH1"
        rows = [
            _row("gene_two_subjects", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 500, LYK3),
            _row("gene_two_subjects", "sp|P00000|OTHER_MEDTR", 84.0, 600, 1100, other),
        ]
        _, index, manifest = self._build(rows, ["gene_two_subjects"])
        record = dict(zip(index[0].split("\t"), index[1].split("\t")))
        self.assertEqual(record["assignment_status"], "dark_no_qualifying_hit")
        self.assertEqual(manifest["exclusion_breakdown"]["genes_recovered_by_merged_hsp_coverage"], [])

    def test_single_hsp_coverage_is_unchanged(self):
        """A one-HSP hit keeps exactly the coverage it had before the change."""
        rows = [_row("gene_single", "sp|Q6UD73|LYK3_MEDTR", 84.0, 1, 700, LYK3)]
        _, index, manifest = self._build(rows, ["gene_single"])
        record = dict(zip(index[0].split("\t"), index[1].split("\t")))
        self.assertEqual(record["assignment_status"], "assigned")
        self.assertAlmostEqual(float(record["best_query_coverage"]), 700 / 1192, places=4)
        self.assertEqual(manifest["exclusion_breakdown"]["genes_recovered_by_merged_hsp_coverage"], [])


if __name__ == "__main__":
    unittest.main()
