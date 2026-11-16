"""Projection preflight must fail before the public member builder writes files."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import route_c_direct_members as members


class RouteCMemberPreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "members"
        self.evidence = self.root / "evidence.tsv"
        self.evidence.write_text("record_id\tconflict_state\nr1\tmapping_conflict\n", encoding="utf-8")
        self.manifest = self.root / "projection_manifest.json"
        self.metadata = {
            "projection_evidence_contract": members.PROJECTION_EVIDENCE_CONTRACT,
            "outputs": {"projection_evidence": str(self.evidence)},
            "projection_evidence_sha256": members.sha256(self.evidence),
            "mapping_conflict_record_ids": ["r1"],
            "mapping_conflict_count": 1,
        }
        self.write_manifest()
        # Isolate the writer/preflight boundary, not the projection row schema.
        validator = patch.object(members, "validate_projection_evidence_frame", side_effect=lambda frame: frame)
        validator.start()
        self.addCleanup(validator.stop)

    def write_manifest(self):
        self.manifest.write_text(json.dumps(self.metadata), encoding="utf-8")

    def run_builder(self):
        return members.build_route_c_direct_members(
            "frozen.tsv", "freeze.json", "projection.tsv", str(self.manifest),
            "atlas.tsv", str(self.out), case_id="synthetic",
        )

    def writer(self, *args, **kwargs):
        self.out.mkdir(parents=True, exist_ok=True)
        table = self.out / "route_c_direct_members.tsv"
        manifest = self.out / "route_c_direct_members_manifest.json"
        row = {column: "" for column in members.OUTPUT_COLUMNS}
        row.update(candidate_gene_id="G1", record_id="r1", support_status=members.SUPPORTED_STATUS,
                   membership_effect="supported_union", rank_contribution="none")
        pd.DataFrame([row], columns=members.OUTPUT_COLUMNS).to_csv(table, sep="\t", index=False)
        manifest.write_text(json.dumps({"inputs": {}, "policy": []}), encoding="utf-8")
        return {"members": str(table), "manifest": str(manifest)}

    def test_bad_digest_never_invokes_writer_or_creates_output(self):
        self.metadata["projection_evidence_sha256"] = "0" * 64
        self.write_manifest()
        with patch.object(members, "_legacy_build_route_c_direct_members", side_effect=self.writer) as writer:
            with self.assertRaisesRegex(ValueError, "hash does not match"):
                self.run_builder()
            writer.assert_not_called()
        self.assertFalse(self.out.exists())

    def test_bad_digest_does_not_overwrite_prior_output(self):
        self.writer()
        before = {p.name: p.read_bytes() for p in self.out.iterdir()}
        self.metadata["projection_evidence_sha256"] = "0" * 64
        self.write_manifest()
        with patch.object(members, "_legacy_build_route_c_direct_members", side_effect=self.writer) as writer:
            with self.assertRaises(ValueError):
                self.run_builder()
            writer.assert_not_called()
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.out.iterdir()})

    def test_missing_evidence_never_invokes_writer(self):
        self.evidence.unlink()
        with patch.object(members, "_legacy_build_route_c_direct_members") as writer:
            with self.assertRaisesRegex(ValueError, "artifact or digest is missing"):
                self.run_builder()
            writer.assert_not_called()

    def test_invalid_contract_roster_and_count_fail_before_writer(self):
        changes = [
            {"projection_evidence_contract": "unsupported-contract"},
            {"outputs": []},
            {"mapping_conflict_record_ids": []},
            {"mapping_conflict_record_ids": "r1"},
            {"mapping_conflict_count": 0},
            {"mapping_conflict_count": "not-an-integer"},
        ]
        original = dict(self.metadata)
        for change in changes:
            with self.subTest(change=change):
                self.metadata = {**original, **change}
                self.write_manifest()
                with patch.object(members, "_legacy_build_route_c_direct_members") as writer:
                    with self.assertRaises(ValueError):
                        self.run_builder()
                    writer.assert_not_called()
        self.assertFalse(self.out.exists())

    def test_valid_conflict_still_downgrades_support_and_updates_digest(self):
        with patch.object(members, "_legacy_build_route_c_direct_members", side_effect=self.writer) as writer:
            result = self.run_builder()
            writer.assert_called_once()
        frame = pd.read_csv(result["members"], sep="\t").fillna("")
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual([members.REVIEW_STATUS], frame["support_status"].tolist())
        self.assertEqual(["review_only"], frame["membership_effect"].tolist())
        self.assertEqual(0, result["supported_count"])
        self.assertEqual(1, result["review_count"])
        self.assertEqual(1, result["projection_conflict_support_blocked_count"])
        self.assertEqual(members.sha256(Path(result["members"])), manifest["output_sha256"])

    def test_valid_nonconflict_preserves_support(self):
        self.evidence.write_text("record_id\tconflict_state\nr1\tnone\n", encoding="utf-8")
        self.metadata.update(projection_evidence_sha256=members.sha256(self.evidence),
                             mapping_conflict_record_ids=[], mapping_conflict_count=0)
        self.write_manifest()
        with patch.object(members, "_legacy_build_route_c_direct_members", side_effect=self.writer):
            result = self.run_builder()
        self.assertEqual(1, result["supported_count"])
        self.assertEqual(0, result["projection_conflict_support_blocked_count"])

    def test_legacy_generation_without_ledger_keeps_existing_behavior(self):
        self.metadata = {}
        self.write_manifest()
        expected = {"legacy": True}
        with patch.object(members, "_legacy_build_route_c_direct_members", return_value=expected) as writer:
            self.assertEqual(expected, self.run_builder())
            writer.assert_called_once()


if __name__ == "__main__":
    unittest.main()
