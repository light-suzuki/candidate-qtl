"""Issue #578: chained assembly bridges must retain root native evidence provenance."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.assembly_bridge import build_assembly_bridge
from candidate_qtl.route_b_prepare import _reference_function_table


OLD = "medtr.A17_HM341.gnm4.ann2.ROOT1"
CURRENT = "medtr.A17.gnm5.ann1_6.CUR1"
R108 = "medtr.R108.gnmHiC_1.ann1.R1"


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def _fam3(gene: str) -> str:
    return f"native_gene_id\tfamily_id\n{gene}\tF1\n"


def _sequence(source: str, target: str) -> str:
    return (
        "source_full_id\ttarget_full_id\treciprocal_sequence_state\t"
        "sequence_identity\tquery_coverage\tsubject_coverage\n"
        f"{source}\t{target}\trbh\t98.0\t1.0\t1.0\n"
    )


class MultiHopBridgeProvenanceTests(unittest.TestCase):
    def _first_bridge(self, root: Path) -> tuple[dict[str, object], Path]:
        functions = _write(
            root,
            "functions.tsv",
            (
                "native_gene_id\tnative_species\tnative_accession\tnative_assembly\t"
                "native_annotation\tnative_gene_id_namespace\tsource_record_id\t"
                "phenotype_or_function\tontology_ids\tsource_ids\tsource_confidence\t"
                "gene_symbol\tevidence_type\tmutant_or_perturbation_background\t"
                "provenance_detail\n"
                f"{OLD}\tMedicago truncatula\tA17_HM341\tgnm4\tann2\t"
                "medtr.A17_HM341.gnm4.ann2\tgfr-record-1\tcurated root phenotype\t"
                "TO:0000001;PO:0009010\tPMID:12345\t5\tMtROOT\t"
                "lis_gene_function_registry\tA17\tcurator note: reviewed\n"
            ),
        )
        source_fam3 = _write(root, "old_fam3.tsv", _fam3(OLD))
        target_fam3 = _write(root, "current_fam3.tsv", _fam3(CURRENT))
        sequence = _write(root, "old_current_sequence.tsv", _sequence(OLD, CURRENT))
        out = root / "first"
        result = build_assembly_bridge(
            functions,
            source_fam3,
            target_fam3,
            sequence,
            out,
            source_namespace="medtr.A17_HM341.gnm4.ann2",
            target_namespace="medtr.A17.gnm5.ann1_6",
            mapping_scope="same_accession_annotation_generation",
            case_id="issue578-first",
        )
        return result, out

    def _second_bridge(self, root: Path, first_reference: Path) -> tuple[dict[str, object], Path]:
        source_fam3 = _write(root, "current2_fam3.tsv", _fam3(CURRENT))
        target_fam3 = _write(root, "r108_fam3.tsv", _fam3(R108))
        sequence = _write(root, "current_r108_sequence.tsv", _sequence(CURRENT, R108))
        out = root / "second"
        result = build_assembly_bridge(
            first_reference,
            source_fam3,
            target_fam3,
            sequence,
            out,
            source_namespace="medtr.A17.gnm5.ann1_6",
            target_namespace="medtr.R108.gnmHiC_1.ann1",
            mapping_scope="cross_accession",
            case_id="issue578-second",
        )
        return result, out

    def test_one_hop_retains_root_native_card_and_upstream_detail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, _out = self._first_bridge(root)
            frame = pd.read_csv(result["route_b_reference_function"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(frame))
            row = frame.iloc[0]

            self.assertEqual(OLD, row.reference_gene_id)
            self.assertEqual(CURRENT, row.gene_id)
            self.assertEqual("Medicago truncatula", row.reference_species)
            self.assertEqual("MtROOT", row.gene_symbol)
            self.assertEqual("lis_gene_function_registry", row.evidence_type)

            provenance = json.loads(row.provenance_detail)
            self.assertEqual("curator note: reviewed", provenance["upstream_provenance_detail"])
            bridge = provenance["assembly_bridge"]
            native = bridge["evidence_native"]
            self.assertEqual(OLD, native["gene_id"])
            self.assertEqual("A17_HM341", native["accession"])
            self.assertEqual("gnm4", native["assembly"])
            self.assertEqual("ann2", native["annotation"])
            self.assertEqual("medtr.A17_HM341.gnm4.ann2", native["gene_id_namespace"])
            self.assertEqual("gfr-record-1", native["source_record_id"])
            self.assertEqual("TO:0000001;PO:0009010", native["ontology_ids"])
            self.assertEqual("5", native["source_confidence"])
            self.assertEqual("MtROOT", native["gene_symbol"])
            self.assertEqual("A17", native["mutant_background"])
            self.assertEqual([OLD, CURRENT], bridge["mapping_path"])
            self.assertEqual(1, len(bridge["mapping_hops"]))
            hop = bridge["mapping_hops"][0]
            self.assertEqual(OLD, hop["source_full_id"])
            self.assertEqual(CURRENT, hop["target_full_id"])
            self.assertEqual("same_accession_annotation_generation", hop["mapping_scope"])
            self.assertEqual("verified", hop["mapping_state"])
            self.assertEqual("one_to_one", hop["relationship"])
            self.assertIn("fam3_concordance", hop["mapping_basis"])
            self.assertIn("pairwise_sequence", hop["mapping_basis"])

    def test_two_hops_keep_root_native_evidence_and_ordered_mapping_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first, _first_out = self._first_bridge(root)
            second, _second_out = self._second_bridge(
                root, Path(first["route_b_reference_function"])
            )
            frame = pd.read_csv(second["route_b_reference_function"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(frame))
            row = frame.iloc[0]

            # Existing per-edge semantics remain unchanged: current A17 is the
            # immediate reference for the current-A17 -> R108 edge.
            self.assertEqual(CURRENT, row.reference_gene_id)
            self.assertEqual(R108, row.gene_id)
            self.assertEqual("Medicago truncatula", row.reference_species)
            self.assertEqual("PMID:12345", row.source_id)

            provenance = json.loads(row.provenance_detail)
            self.assertEqual("curator note: reviewed", provenance["upstream_provenance_detail"])
            bridge = provenance["assembly_bridge"]
            native = bridge["evidence_native"]
            self.assertEqual(OLD, native["gene_id"])
            self.assertEqual("gfr-record-1", native["source_record_id"])
            self.assertEqual("TO:0000001;PO:0009010", native["ontology_ids"])
            self.assertEqual("5", native["source_confidence"])
            self.assertEqual([OLD, CURRENT, R108], bridge["mapping_path"])

            hops = bridge["mapping_hops"]
            self.assertEqual(2, len(hops))
            self.assertEqual((OLD, CURRENT), (hops[0]["source_full_id"], hops[0]["target_full_id"]))
            self.assertEqual(
                (CURRENT, R108),
                (hops[1]["source_full_id"], hops[1]["target_full_id"]),
            )
            self.assertEqual("same_accession_annotation_generation", hops[0]["mapping_scope"])
            self.assertEqual("cross_accession", hops[1]["mapping_scope"])
            self.assertTrue(all(hop["mapping_state"] == "verified" for hop in hops))

            # route_b_prepare already has provenance_detail in its normalized
            # schema. Ensure this bridge-owned chain survives that normalization
            # byte-for-byte rather than being discarded or reinterpreted.
            normalized_dir = root / "normalized"
            normalized_dir.mkdir()
            normalized_path, _meta = _reference_function_table(
                [Path(second["route_b_reference_function"])], normalized_dir
            )
            self.assertIsNotNone(normalized_path)
            normalized = pd.read_csv(normalized_path, sep="\t", dtype=str).fillna("")
            self.assertEqual(row.provenance_detail, normalized.iloc[0].provenance_detail)
            self.assertEqual("Medicago truncatula", normalized.iloc[0].reference_species)
            self.assertEqual("MtROOT", normalized.iloc[0].gene_symbol)
            self.assertEqual("lis_gene_function_registry", normalized.iloc[0].evidence_type)

    def test_malformed_bridge_owned_provenance_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            functions = _write(
                root,
                "bad_functions.tsv",
                (
                    "gene_id\tfunction_description\tsource_id\tprovenance_detail\n"
                    f"{CURRENT}\tcurated root phenotype\tPMID:12345\t"
                    "{\"assembly_bridge\":[] }\n"
                ),
            )
            source_fam3 = _write(root, "bad_source_fam3.tsv", _fam3(CURRENT))
            target_fam3 = _write(root, "bad_target_fam3.tsv", _fam3(R108))
            sequence = _write(root, "bad_sequence.tsv", _sequence(CURRENT, R108))
            with self.assertRaisesRegex(ValueError, "assembly_bridge provenance_detail"):
                build_assembly_bridge(
                    functions,
                    source_fam3,
                    target_fam3,
                    sequence,
                    root / "bad_out",
                    mapping_scope="cross_accession",
                )


if __name__ == "__main__":
    unittest.main()
