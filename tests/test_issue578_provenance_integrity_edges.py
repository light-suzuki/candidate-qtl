"""Adversarial provenance checks for Issue #578 assembly-bridge chaining."""

from __future__ import annotations

import json
import unittest

from candidate_qtl.assembly_bridge import _bridge_provenance


OLD = "medtr.A17_HM341.gnm4.ann2.ROOT1"
CURRENT = "medtr.A17.gnm5.ann1_6.CUR1"
R108 = "medtr.R108.gnmHiC_1.ann1.R1"


def _map(source: str, target: str, *, scope: str = "same_accession_annotation_generation") -> dict[str, str]:
    return {
        "source_full_id": source,
        "target_full_id": target,
        "mapping_scope": scope,
        "mapping_basis": "fam3_concordance+pairwise_sequence",
        "relationship": "one_to_one",
        "mapping_state": "verified",
        "mapping_reason": "unique_reciprocal_sequence_and_family_concordance",
        "pangene_dataset": "",
    }


class AssemblyBridgeProvenanceIntegrityEdgeTests(unittest.TestCase):
    def test_parseable_full_id_does_not_fabricate_unsupplied_native_metadata(self) -> None:
        detail = _bridge_provenance(
            {"native_gene_id": OLD, "native_species": "Medicago truncatula"},
            _map(OLD, CURRENT),
            OLD,
            CURRENT,
        )
        native = json.loads(detail)["assembly_bridge"]["evidence_native"]
        self.assertEqual(OLD, native["gene_id"])
        self.assertEqual("Medicago truncatula", native["species"])
        for field in ("species_code", "accession", "assembly", "annotation", "gene_id_namespace"):
            self.assertEqual("", native[field], field)

    def test_existing_root_is_immutable_and_second_hop_appends_only_one_edge(self) -> None:
        first = _bridge_provenance(
            {
                "native_gene_id": OLD,
                "native_accession": "A17_HM341",
                "native_assembly": "gnm4",
                "native_annotation": "ann2",
                "source_record_id": "gfr-1",
            },
            _map(OLD, CURRENT),
            OLD,
            CURRENT,
        )
        second = _bridge_provenance(
            {
                "gene_id": CURRENT,
                "native_accession": "SHOULD_NOT_REWRITE_ROOT",
                "provenance_detail": first,
            },
            _map(CURRENT, R108, scope="cross_accession"),
            CURRENT,
            R108,
        )
        bridge = json.loads(second)["assembly_bridge"]
        self.assertEqual("A17_HM341", bridge["evidence_native"]["accession"])
        self.assertEqual("gfr-1", bridge["evidence_native"]["source_record_id"])
        self.assertEqual([OLD, CURRENT, R108], bridge["mapping_path"])
        self.assertEqual(2, len(bridge["mapping_hops"]))

    def test_malformed_evidence_native_fails_closed(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": [],
                "mapping_path": [OLD, CURRENT],
                "mapping_hops": [{"source_full_id": OLD, "target_full_id": CURRENT}],
            }
        })
        with self.assertRaisesRegex(ValueError, "evidence_native"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)

    def test_empty_root_gene_fails_closed(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": {"gene_id": ""},
                "mapping_path": [CURRENT],
                "mapping_hops": [],
            }
        })
        with self.assertRaisesRegex(ValueError, "non-empty gene_id"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)

    def test_falsy_wrong_type_mapping_hops_fails_closed(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": {"gene_id": OLD},
                "mapping_path": [OLD, CURRENT],
                "mapping_hops": {},
            }
        })
        with self.assertRaisesRegex(ValueError, "mapping_hops must be a list"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)

    def test_non_object_hop_entry_fails_closed(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": {"gene_id": OLD},
                "mapping_path": [OLD, CURRENT],
                "mapping_hops": ["not-an-object"],
            }
        })
        with self.assertRaisesRegex(ValueError, "entries must be objects"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)

    def test_path_root_mismatch_fails_closed(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": {"gene_id": OLD},
                "mapping_path": ["OTHER", CURRENT],
                "mapping_hops": [{"source_full_id": "OTHER", "target_full_id": CURRENT}],
            }
        })
        with self.assertRaisesRegex(ValueError, "start at the root"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)

    def test_hop_ledger_must_exactly_cover_path(self) -> None:
        malformed = json.dumps({
            "assembly_bridge": {
                "evidence_native": {"gene_id": OLD},
                "mapping_path": [OLD, CURRENT],
                "mapping_hops": [],
            }
        })
        with self.assertRaisesRegex(ValueError, "cover every mapping_path edge"):
            _bridge_provenance({"provenance_detail": malformed}, _map(CURRENT, R108), CURRENT, R108)


if __name__ == "__main__":
    unittest.main()
