from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.supported_member_union import build_supported_member_union
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _artifact(root: Path, route: str, rows: list[dict], *, case_id: str = "case1") -> tuple[Path, Path]:
    names = {"a": "family_candidate_bridge.tsv", "b": "route_b_verified_candidates.tsv", "c": "route_c_direct_members.tsv"}
    path = root / names[route]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path, manifest_for(path, route=route, case_id=case_id)


class Issue664CanonicalLineageTests(unittest.TestCase):
    def test_bare_unrelated_digest_and_edit_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "b", [{"candidate_gene_id": "G1", "route": "b", "support_status": "verified_semantic_support"}])
            manifest.unlink()
            with self.assertRaisesRegex(ValueError, "producer manifest"):
                build_supported_member_union(str(root / "bare"), route_b_candidates_arg=str(path))
            manifest = manifest_for(path, route="b")
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            payload["unrelated_sha256"] = "0" * 64
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            result = build_supported_member_union(str(root / "ok"), route_b_candidates_arg=str(path), route_b_candidates_manifest_arg=str(manifest))
            self.assertEqual({"G1"}, set(pd.read_csv(result["candidate_set"], sep="\t")["candidate_gene_id"]))
            path.write_text(path.read_text(encoding="utf-8").replace("G1", "G2"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash"):
                build_supported_member_union(str(root / "edited"), route_b_candidates_arg=str(path), route_b_candidates_manifest_arg=str(manifest))

    def test_output_path_and_route_b_claim_ancestry_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "b", [{"candidate_gene_id": "G1", "route": "b", "support_status": "verified_semantic_support"}])
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            payload.pop("output")
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "output path"):
                build_supported_member_union(str(root / "missing-output"), route_b_candidates_arg=str(path), route_b_candidates_manifest_arg=str(manifest))

            manifest = manifest_for(path, route="b")
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            payload.pop("inputs")
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ancestry"):
                build_supported_member_union(str(root / "missing-ancestry"), route_b_candidates_arg=str(path), route_b_candidates_manifest_arg=str(manifest))

            manifest = manifest_for(path, route="b")
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            claims_manifest = Path(payload["inputs"]["claims_manifest"]["path"])
            claims_payload = json.loads(claims_manifest.read_text(encoding="utf-8"))
            claims_payload["artifact_sha256"] = "0" * 64
            claims_payload["outputs"]["claims_sha256"] = "0" * 64
            claims_manifest.write_text(json.dumps(claims_payload), encoding="utf-8")
            payload["inputs"]["claims_manifest"]["sha256"] = hashlib.sha256(claims_manifest.read_bytes()).hexdigest()
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not bound"):
                build_supported_member_union(str(root / "claims-unbound"), route_b_candidates_arg=str(path), route_b_candidates_manifest_arg=str(manifest))

    def test_route_a_manifest_requires_hashed_input_ancestry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "a", [{"gene_id": "G1", "route": "a", "mapping_tier": "tier1"}])
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            payload.pop("inputs")
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "input ancestry"):
                build_supported_member_union(str(root / "missing-a-ancestry"), route_a_bridge_arg=str(path), route_a_manifest_arg=str(manifest))

    def test_exact_route_membership_basis(self) -> None:
        for routes in (("a",), ("b",), ("c",), ("a", "b"), ("a", "c"), ("b", "c"), ("a", "b", "c")):
            with self.subTest(routes=routes), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                args: dict[str, str] = {}
                for route in routes:
                    rows = ([{"gene_id": "G", "route": "a", "mapping_tier": "tier1"}] if route == "a" else
                            [{"candidate_gene_id": "G", "route": "b", "support_status": "verified_semantic_support"}] if route == "b" else
                            [{"candidate_gene_id": "G", "support_status": "verified_direct_support", "projection_state": "inside_current_qtl"}])
                    path, manifest = _artifact(root, route, rows)
                    args[{"a": "route_a_bridge_arg", "b": "route_b_candidates_arg", "c": "route_c_members_arg"}[route]] = str(path)
                    args[{"a": "route_a_manifest_arg", "b": "route_b_candidates_manifest_arg", "c": "route_c_manifest_arg"}[route]] = str(manifest)
                build_supported_member_union(str(root / "union"), case_id="case1", **args)
                wet = wet_candidate_handoff(str(root / "wet"), case_id="case1", wet_budget=20, **args)
                frame = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
                self.assertEqual("verified_route_" + "_".join(routes) + "_support", frame.loc[0, "membership_basis"])
                self.assertEqual({"G"}, set(frame["candidate_gene_id"]))

    def test_route_c_wrong_qtl_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "c", [{"candidate_gene_id": "G", "support_status": "verified_direct_support", "projection_state": "inside_current_qtl"}])
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            payload["qtl_id"] = "wrong"
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "qtl"):
                build_supported_member_union(str(root / "bad"), route_c_members_arg=str(path), route_c_manifest_arg=str(manifest), case_id="case1")

    def test_swapped_case_interval_and_assembly_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a_path, a_manifest = _artifact(root, "a", [{"gene_id": "G1", "route": "a", "mapping_tier": "tier1"}])
            b_path, b_manifest = _artifact(root, "b", [{"candidate_gene_id": "G2", "route": "b", "support_status": "verified_semantic_support"}])
            payload = json.loads(b_manifest.read_text(encoding="utf-8"))
            payload["case_id"] = "other-case"
            b_manifest.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "case_id"):
                build_supported_member_union(str(root / "case-bad"), route_a_bridge_arg=str(a_path), route_a_manifest_arg=str(a_manifest), route_b_candidates_arg=str(b_path), route_b_candidates_manifest_arg=str(b_manifest))

            b_manifest = manifest_for(b_path, route="b", interval_fingerprint="other-interval")
            with self.assertRaisesRegex(ValueError, "interval_fingerprint"):
                build_supported_member_union(str(root / "interval-bad"), route_a_bridge_arg=str(a_path), route_a_manifest_arg=str(a_manifest), route_b_candidates_arg=str(b_path), route_b_candidates_manifest_arg=str(b_manifest))

            b_manifest = manifest_for(b_path, route="b", assembly="other-assembly")
            with self.assertRaisesRegex(ValueError, "assembly"):
                build_supported_member_union(str(root / "assembly-bad"), route_a_bridge_arg=str(a_path), route_a_manifest_arg=str(a_manifest), route_b_candidates_arg=str(b_path), route_b_candidates_manifest_arg=str(b_manifest))

    def test_explicit_missing_optional_wet_input_is_not_silent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "a", [{"gene_id": "G1", "route": "a", "mapping_tier": "tier1"}])
            with self.assertRaisesRegex(ValueError, "explicitly supplied optional input"):
                wet_candidate_handoff(
                    str(root / "wet"),
                    route_a_bridge_arg=str(path),
                    route_a_manifest_arg=str(manifest),
                    case_id="case1",
                    interval_atlas_arg=str(root / "missing-atlas.tsv"),
                )

    def test_explicit_fill_input_from_another_case_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest = _artifact(root, "a", [{"gene_id": "G1", "route": "a", "mapping_tier": "tier1"}])
            atlas = root / "interval_gene_atlas.tsv"
            pd.DataFrame([{"gene_id": "G2", "case_id": "other-case", "family_id": "F2"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(ValueError, "case_id"):
                wet_candidate_handoff(
                    str(root / "wet"), route_a_bridge_arg=str(path),
                    route_a_manifest_arg=str(manifest), case_id="case1",
                    interval_atlas_arg=str(atlas),
                )


if __name__ == "__main__":
    unittest.main()
