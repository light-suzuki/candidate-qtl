"""Public sequence / local target correspondence, without causal promotion."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from candidate_qtl.reference_sequence_bridge import build_reference_sequence_bridge, load_reference_sequence_bridge, _run_search
from candidate_qtl.core import sha256

S = "ACDEFGHIKLMNPQRSTVWY" * 5


def _record(path, accession="AA123456.1", genes=None, taxid="3888"):
    genes = genes or [("NODE", "ABC12345.1", S)]
    record = SeqRecord(Seq("ATG" * 1000), id=accession, name=accession.split(".")[0], description="Synthetic official-format record")
    record.annotations = {"molecule_type": "DNA", "organism": "Plant alpha"}
    record.features = [SeqFeature(FeatureLocation(0, len(record)), type="source", qualifiers={"db_xref": [f"taxon:{taxid}"]})]
    cursor = 0
    for gene, protein, sequence in genes:
        record.features.append(SeqFeature(FeatureLocation(cursor, cursor + len(sequence) * 3), type="CDS",
            qualifiers={"gene": [gene], "protein_id": [protein], "translation": [sequence]}))
        cursor += len(sequence) * 3 + 3
    with Path(path).open("w", encoding="utf-8") as handle:
        SeqIO.write(record, handle, "genbank")


def _inputs(root, *, requests=None, targets=None, genes=None):
    cache = root / "records"
    cache.mkdir()
    _record(cache / "AA123456.1.gb", genes=genes)
    requests = requests or [{"request_id": "r1", "source_id": "PMID:1234567", "reported_gene": "NODE",
        "reference_accession": "AA123456.1", "reference_taxid": "3888", "focal_taxid": "3888"}]
    req, proteins, mapping = root / "requests.tsv", root / "proteins.fa", root / "mapping.tsv"
    pd.DataFrame(requests).to_csv(req, sep="\t", index=False)
    targets = targets or {"p1": S + "*"}
    proteins.write_text("".join(f">{protein}\n{sequence}\n" for protein, sequence in targets.items()), encoding="utf-8")
    pd.DataFrame([{"protein_id": protein, "gene_id": f"g{index}"} for index, protein in enumerate(targets, 1)]).to_csv(mapping, sep="\t", index=False)
    return req, proteins, mapping, cache


def _hit(query="q000001", target="p1", *, qlen=100, slen=100, qstart=1, qend=100, sstart=1, send=100, identity=100):
    return {"qseqid": query, "sseqid": target, "pident": identity, "qlen": qlen, "slen": slen,
        "qstart": qstart, "qend": qend, "sstart": sstart, "send": send, "evalue": 1e-40, "bitscore": 200,
        "hsp_count": 1, "query_union_segments": f"{qstart}-{qend}", "target_union_segments": f"{sstart}-{send}",
        "query_coverage": (qend - qstart + 1) / qlen, "target_coverage": (send - sstart + 1) / slen}


class ReferenceSequenceBridgeTests(unittest.TestCase):
    def build(self, root, hits=None, *, requests=None, targets=None, genes=None, state="evaluated"):
        req, proteins, mapping, cache = _inputs(root, requests=requests, targets=targets, genes=genes)
        with patch("candidate_qtl.reference_sequence_bridge._run_search", return_value={"state": state, "hits": hits or []}):
            result = build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache)
        manifest, frame, eligible = load_reference_sequence_bridge(result["manifest"])
        return result, manifest, frame, eligible

    def test_full_correspondence_freezes_sources_and_strips_terminal_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, manifest, frame, hits = self.build(Path(tmp), [_hit()])
            self.assertEqual({"g1"}, set(hits.target_id))
            self.assertEqual("requires_researcher_assessment", frame.iloc[0].locus_identity)
            self.assertEqual("model_coverage_passes_review_threshold", frame.iloc[0].whole_annotation_model_identity)
            self.assertEqual(["ABC12345.1", "NODE"], json.loads(frame.iloc[0].source_gene_names_json))
            self.assertEqual("ABC12345.1", hits.iloc[0].cds_protein_id)
            self.assertEqual("Plant alpha", frame.iloc[0].source_organism)
            self.assertEqual("Plant alpha", hits.iloc[0].source_organism)
            self.assertNotIn("*", (Path(tmp) / "out/target_proteins.fa").read_text())
            self.assertEqual("none", manifest["tier_effect"])
            self.assertEqual("none", manifest["support_effect"])
            self.assertFalse(json.loads((Path(tmp) / "out/source_freeze.json").read_text())["target_accessed"])
            with self.assertRaisesRegex(ValueError, "new output generation"):
                build_reference_sequence_bridge(Path(tmp) / "requests.tsv", Path(tmp) / "proteins.fa", Path(tmp) / "mapping.tsv", Path(tmp) / "out")
            self.assertEqual(1, result["request_count"])

    def test_partial_model_does_not_erase_locus_correspondence(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, frame, hits = self.build(Path(tmp), [_hit(slen=200)], targets={"p1": S * 2})
            self.assertEqual(1, len(hits))
            self.assertEqual("partial_model_correspondence", frame.iloc[0].whole_annotation_model_identity)
            self.assertEqual("1.0", frame.iloc[0].query_coverage)
            self.assertEqual("0.5", frame.iloc[0].target_coverage)

    def test_short_domain_matches_remain_visible_but_are_not_eligible(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, frame, hits = self.build(Path(tmp), [_hit(qend=30, send=30)])
            self.assertTrue(hits.empty)
            self.assertEqual("below_review_threshold", frame.iloc[0].eligibility)
            self.assertEqual("evaluated_with_hit", frame.iloc[0].evaluation_state)

    def test_equal_competitors_and_unplaced_scaffold_are_all_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, frame, hits = self.build(Path(tmp), [_hit(target="chr_protein"), _hit(target="unplaced_scaffold_protein")],
                targets={"chr_protein": S, "unplaced_scaffold_protein": S})
            self.assertEqual({"g1", "g2"}, set(hits.target_id))
            self.assertEqual(2, len(frame))

    def test_two_source_cds_in_disjoint_halves_mark_annotation_model(self):
        request = {"source_id": "PMID:1234567", "reference_accession": "AA123456.1", "reference_taxid": "3888", "focal_taxid": "3888"}
        requests = [{**request, "request_id": "r1", "reported_gene": "NODE_A", "cds_protein_id": "ABC12345.1"},
                    {**request, "request_id": "r2", "reported_gene": "NODE_B", "cds_protein_id": "ABC12346.1"}]
        with tempfile.TemporaryDirectory() as tmp:
            _, _, frame, hits = self.build(Path(tmp), [_hit(slen=200), _hit(query="q000002", slen=200, sstart=101, send=200)],
                requests=requests, targets={"p1": S * 2}, genes=[("NODE_A", "ABC12345.1", S), ("NODE_B", "ABC12346.1", S)])
            self.assertEqual(2, len(hits))
            self.assertEqual({"multiple_reported_genes_in_model"}, set(frame.whole_annotation_model_identity))

    def test_overlapping_query_hsps_do_not_manufacture_full_model_coverage(self):
        request = {"source_id": "PMID:1234567", "reference_accession": "AA123456.1", "reference_taxid": "3888", "focal_taxid": "3888"}
        requests = [{**request, "request_id": "r1", "reported_gene": "NODE_A", "cds_protein_id": "ABC12345.1"},
                    {**request, "request_id": "r2", "reported_gene": "NODE_B", "cds_protein_id": "ABC12346.1"}]
        raw = ("q000001\tp1\t100\t100\t100\t200\t1\t100\t1\t100\t1e-60\t240\n"
               "q000001\tp1\t67\t100\t100\t200\t1\t100\t101\t200\t1e-35\t160\n"
               "q000002\tp1\t99\t100\t100\t200\t1\t100\t101\t200\t1e-60\t239\n"
               "q000002\tp1\t66\t100\t100\t200\t1\t100\t1\t100\t1e-35\t159\n")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            req, proteins, mapping, cache = _inputs(root, requests=requests, targets={"p1": S * 2},
                genes=[("NODE_A", "ABC12345.1", S), ("NODE_B", "ABC12346.1", S)])
            def run(command, **kwargs):
                Path(command[command.index("-out") + 1]).write_text(raw)
                from subprocess import CompletedProcess
                return CompletedProcess(command, 0, "", "")
            with patch("candidate_qtl.reference_sequence_bridge._resolve_runner", return_value=("native", None)), \
                    patch("candidate_qtl.reference_sequence_bridge.build_protein_db", return_value={"available": True}), \
                    patch("candidate_qtl.reference_sequence_bridge._tool_version", return_value="test"), \
                    patch("candidate_qtl.reference_sequence_bridge.run_tool", side_effect=run):
                result = build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache)
            _, frame, hits = load_reference_sequence_bridge(result["manifest"])
            self.assertEqual(4, len(frame))
            self.assertEqual(4, len(hits))
            self.assertEqual({"single_hsp"}, set(frame.alignment_scope))
            self.assertEqual({"1"}, set(frame.hsp_count))
            self.assertEqual({"0.5"}, set(frame.target_coverage))
            self.assertEqual({100, 99, 67, 66}, set(frame.pident.astype(float)))
            self.assertEqual({"multiple_reported_genes_in_model"}, set(frame.whole_annotation_model_identity))
            self.assertEqual(raw, (root / "out/all_hsps.tsv").read_text())

    def test_bac_requires_explicit_cds_even_when_reported_name_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, frame, hits = self.build(Path(tmp), genes=[("NODE", "ABC12345.1", S), ("OTHER", "ABC12346.1", S)])
            self.assertTrue(hits.empty)
            self.assertEqual("not_evaluated_ambiguous_cds", frame.iloc[0].sequence_state)

    def test_reference_taxid_and_cross_species_cannot_authorize_identity(self):
        for taxid, expected in (("1", "not_evaluated_reference_taxid_mismatch"), ("3888", "not_evaluated_cross_species")):
            with self.subTest(taxid=taxid), tempfile.TemporaryDirectory() as tmp:
                request = {"request_id": "r1", "source_id": "PMID:1234567", "reported_gene": "NODE", "reference_accession": "AA123456.1",
                    "reference_taxid": taxid, "focal_taxid": "2"}
                _, _, frame, hits = self.build(Path(tmp), requests=[request])
                self.assertTrue(hits.empty)
                self.assertEqual(expected, frame.iloc[0].sequence_state)

    def test_version_mismatch_is_not_silently_replaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            req, proteins, mapping, cache = _inputs(root)
            _record(cache / "AA123456.1.gb", accession="AA123456.2")
            result = build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache)
            _, frame, hits = load_reference_sequence_bridge(result["manifest"])
            self.assertTrue(hits.empty)
            self.assertEqual("not_evaluated_accession_version_mismatch", frame.iloc[0].sequence_state)

    def test_tool_failure_and_successful_no_hit_stay_distinct(self):
        for state, expected in (("evaluated", "evaluated_no_hit"), ("not_evaluated_blast_unavailable", "not_evaluated_blast_unavailable")):
            with self.subTest(state=state), tempfile.TemporaryDirectory() as tmp:
                _, _, frame, hits = self.build(Path(tmp), state=state)
                self.assertTrue(hits.empty)
                self.assertEqual(expected, frame.iloc[0].evaluation_state)

    def test_record_failure_is_local_and_network_requires_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            req, proteins, mapping, cache = _inputs(root)
            (cache / "AA123456.1.gb").unlink()
            with patch("candidate_qtl.reference_sequence_bridge.urllib.request.urlopen") as network:
                result = build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache)
                network.assert_not_called()
            _, frame, _ = load_reference_sequence_bridge(result["manifest"])
            self.assertEqual("not_evaluated_record_unavailable", frame.iloc[0].sequence_state)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            req, proteins, mapping, cache = _inputs(root)
            (cache / "AA123456.1.gb").unlink()
            with patch("candidate_qtl.reference_sequence_bridge.urllib.request.urlopen", side_effect=OSError("offline")):
                result = build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache, allow_external_api=True)
            _, frame, _ = load_reference_sequence_bridge(result["manifest"])
            self.assertEqual("not_evaluated_provider_failure", frame.iloc[0].sequence_state)

    def test_hash_drift_of_saved_artifact_or_original_input_is_rejected(self):
        for relative in ("proteins.fa", "out/raw/record_000001.txt", "out/eligible_hits.tsv"):
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result, _, _, _ = self.build(root, [_hit()])
                with (root / relative).open("a") as handle:
                    handle.write("changed")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    load_reference_sequence_bridge(result["manifest"])

    def test_protein_gene_mapping_must_be_complete_and_unambiguous(self):
        for kind in ("missing", "duplicate"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                req, proteins, mapping, cache = _inputs(root)
                mapping.write_text("protein_id\tgene_id\nother\tg1\n" if kind == "missing" else "protein_id\tgene_id\np1\tg1\np1\tg2\n")
                with self.assertRaises(ValueError):
                    build_reference_sequence_bridge(req, proteins, mapping, root / "out", records_arg=cache)

    def test_resealed_internally_inconsistent_endpoints_or_hits_are_rejected(self):
        for filename, column in (("reference_sequence_correspondence.tsv", "source_id"), ("eligible_hits.tsv", "target_id")):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result, manifest, _, _ = self.build(root, [_hit()])
                output = root / "out" / filename
                frame = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
                frame.loc[0, column] = "a different endpoint"
                frame.to_csv(output, sep="\t", index=False, lineterminator="\n")
                manifest["artifacts"][filename] = sha256(output)
                Path(result["manifest"]).write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, "endpoint mismatch|contradict correspondence"):
                    load_reference_sequence_bridge(result["manifest"])

    def test_tool_search_keeps_all_targets_without_hsp_or_top_one_cap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def run(command, **kwargs):
                (root / "all_hsps.tsv").write_text("q000001\tp1\t100\t100\t100\t100\t1\t100\t1\t100\t1e-40\t200\n")
                from subprocess import CompletedProcess
                return CompletedProcess(command, 0, "", "")
            with patch("candidate_qtl.reference_sequence_bridge._resolve_runner", return_value=("native", None)), \
                    patch("candidate_qtl.reference_sequence_bridge.build_protein_db", return_value={"available": True}), \
                    patch("candidate_qtl.reference_sequence_bridge._tool_version", return_value="test"), \
                    patch("candidate_qtl.reference_sequence_bridge.run_tool", side_effect=run):
                result = _run_search(root / "queries.fa", root / "proteins.fa", root, 987, 4, "native")
            command = result["command"]
            self.assertEqual("987", command[command.index("-max_target_seqs") + 1])
            self.assertNotIn("-max_hsps", command)
            self.assertEqual(1, len(result["hits"]))
            self.assertEqual("single_hsp", result["hits"][0]["alignment_scope"])


if __name__ == "__main__":
    unittest.main()
