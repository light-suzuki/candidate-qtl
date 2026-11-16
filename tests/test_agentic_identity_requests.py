"""Post-freeze focal identity requests are built from the sealed stage only."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from scripts.agentic_identity_requests import build_requests, protein_map

COLUMNS = ("record_id", "record_kind", "source_id", "reported_gene", "reported_species",
           "reference_gene_id", "focal_species", "gene_correspondence")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(accession: str, species: str, taxid: int, xrefs: list[tuple[str, str]]) -> dict:
    return {"primaryAccession": accession,
            "organism": {"scientificName": species, "taxonId": taxid},
            "uniProtKBCrossReferences": [
                {"database": "EMBL", "id": nucleotide,
                 "properties": [{"key": "ProteinId", "value": protein}]}
                for nucleotide, protein in xrefs]}


class IdentityRequestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        stage_dir = self.root / "prepared"
        records = stage_dir / "inputs" / "record_cache"
        records.mkdir(parents=True)
        (records / "Q11111.json").write_text(json.dumps(
            _record("Q11111", "Synthetic plant", 9001, [("AB000001", "BAA00001.1")])), encoding="utf-8")
        (records / "P22222.json").write_text(json.dumps(
            _record("P22222", "Model plant", 9002, [("AB000002", "BAA00002.1")])), encoding="utf-8")
        table = stage_dir / "a_sources" / "paper_evidence_table.tsv"
        table.parent.mkdir()
        rows = [
            {"record_id": "p1", "record_kind": "paper_report", "source_id": "PMID:1", "reported_gene": "GENA",
             "reported_species": "Synthetic plant", "reference_gene_id": "Q11111", "focal_species": "yes"},
            {"record_id": "p2", "record_kind": "paper_report", "source_id": "PMID:2", "reported_gene": "GENB",
             "reported_species": "Model plant", "reference_gene_id": "P22222", "focal_species": "no"},
            {"record_id": "p3", "record_kind": "paper_report", "source_id": "PMID:3", "reported_gene": "GENC",
             "reported_species": "Synthetic plant", "reference_gene_id": "unknown", "focal_species": "yes"},
        ]
        with table.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        files = {str(p.relative_to(stage_dir)): _sha(p) for p in sorted(stage_dir.rglob("*")) if p.is_file()}
        stage = {"case_id": "c1", "focal_species": "Synthetic plant",
                 "inputs": {"record_cache": str(records)},
                 "a_sources": {"paper_evidence_table": str(table)}, "files": files}
        self.anchor = stage_dir / "stage.json"
        self.anchor.write_text(json.dumps(stage), encoding="utf-8")

    def test_focal_uniprot_reference_becomes_insdc_request_and_others_are_skipped(self):
        result = build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))
        self.assertEqual(1, result["request_count"])
        with Path(result["requests"]).open(encoding="utf-8", newline="") as stream:
            request = next(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(("AB000001", "BAA00001.1", "9001"),
                         (request["reference_accession"], request["cds_protein_id"], request["focal_taxid"]))
        self.assertEqual("Q11111", result["request_provenance"][0]["official_primary_accession"])
        reasons = {row["record_id"]: row["reason"] for row in result["skipped"]}
        self.assertEqual({"p2": "reported_gene_not_focal_species",
                          "p3": "reference_not_uniprot_accession"}, reasons)

    def test_changed_stage_artifact_is_rejected(self):
        table = self.root / "prepared" / "a_sources" / "paper_evidence_table.tsv"
        table.write_text(table.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "stage artifact mismatch"):
            build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))

    def test_known_different_gene_proxy_is_not_labelled_as_reported_gene_sequence(self):
        stage = json.loads(self.anchor.read_text())
        table = Path(stage["a_sources"]["paper_evidence_table"])
        with table.open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        rows[0]["gene_correspondence"] = "different_gene"
        with table.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        stage["files"][str(table.relative_to(self.anchor.parent))] = _sha(table)
        self.anchor.write_text(json.dumps(stage))
        result = build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))
        self.assertEqual(0, result["request_count"])
        self.assertEqual("reference_is_not_the_reported_gene", result["skipped"][0]["reason"])

    def test_official_secondary_accession_is_resolved_but_conflicts_fail(self):
        from scripts.agentic_identity_requests import _records
        path = self.root / "aliases.json"
        record = _record("Q11111", "Synthetic plant", 9001, [("AB000001", "BAA00001.1")])
        record["secondaryAccessions"] = ["Q33333"]
        path.write_text(json.dumps([record]), encoding="utf-8")
        self.assertEqual("Q11111", _records(path)["Q33333"]["primaryAccession"])
        other = _record("Q22222", "Synthetic plant", 9001, [("AB000002", "BAA00002.1")])
        other["secondaryAccessions"] = ["Q33333"]
        path.write_text(json.dumps([record, other]), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "conflicting sealed UniProt"):
            _records(path)

    def test_secondary_accession_reaches_requests_with_official_provenance(self):
        stage = json.loads(self.anchor.read_text())
        cached = Path(stage["inputs"]["record_cache"]) / "Q11111.json"
        record = json.loads(cached.read_text())
        record["primaryAccession"] = "Q33333"
        record["secondaryAccessions"] = ["Q11111"]
        cached.write_text(json.dumps(record))
        stage["files"][str(cached.relative_to(self.anchor.parent))] = _sha(cached)
        self.anchor.write_text(json.dumps(stage))
        result = build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))
        self.assertEqual(1, result["request_count"])
        self.assertEqual(("Q11111", "Q33333", "official_secondary"), tuple(
            result["request_provenance"][0][k] for k in ("reference_gene_id", "official_primary_accession", "accession_binding")))

    def test_unsealed_record_added_after_freeze_is_rejected(self):
        records = self.root / "prepared" / "inputs" / "record_cache"
        (records / "extra.json").write_text(json.dumps(
            _record("Q11111", "Synthetic plant", 9001, [("OTHER", "OTHER.1")])), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not a sealed stage artifact"):
            build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))

    def test_table_reference_cannot_escape_the_sealed_file_inventory(self):
        stage = json.loads(self.anchor.read_text())
        table = Path(stage["a_sources"]["paper_evidence_table"])
        external = self.root / "external.tsv"
        external.write_bytes(table.read_bytes())
        stage["a_sources"]["paper_evidence_table"] = str(external)
        self.anchor.write_text(json.dumps(stage), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not inside the sealed stage"):
            build_requests(str(self.anchor), _sha(self.anchor), str(self.root / "out"))

    def test_protein_map_requires_explicit_gene_tokens(self):
        fasta = self.root / "p.fa"
        fasta.write_text(">P1.1 gene=P1\nMA\n>P2.1 gene=P2\nMB\n", encoding="utf-8")
        result = protein_map(str(fasta), str(self.root / "map.tsv"))
        self.assertEqual(2, result["gene_count"])
        fasta.write_text(">P1.1\nMA\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "gene= token"):
            protein_map(str(fasta), str(self.root / "map2.tsv"))


if __name__ == "__main__":
    unittest.main()
