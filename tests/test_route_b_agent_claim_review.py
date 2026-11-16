"""Route-B independent-review binding smoke tests (no biological data)."""
from __future__ import annotations

import copy
import json
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

import pandas as pd

from candidate_qtl.agent_claim_review import AgentClaimReview, B_CONTRACT, prepare_review_prompt
from candidate_qtl.core import sha256
from candidate_qtl.cli import main
from candidate_qtl.cli import build_parser
from candidate_qtl.semantic_claims import CLAIM_COLUMNS
from tests.test_route_b_semantic_transfer_verify import _bundle
from tests.test_agent_claim_review import ReviewFixture, spec, write_json


class RouteBAgentClaimReviewTests(unittest.TestCase):
    def test_supported_seed_and_ordinary_claim_complete_canonical_session_without_support_promotion(self):
        from candidate_qtl.claim_review_session import prepare_session
        from candidate_qtl.source_verify import source_verify

        root = Path(self.tmp.name)
        proposals, source, bundle = self._cli_inputs(root)
        ordinary = pd.read_csv(proposals, sep="\t", dtype=str).fillna("").iloc[0].to_dict()
        seed = {**ordinary, "claim_id": "seed-1", "claim_type": "mechanistic_review_seed",
                "mechanism": "hypothesized developmental signaling transfer",
                "phenotype_relation": "could explain focal organ development",
                "uncertainty_note": "reference observation supports one component; whole transfer remains unverified"}
        pd.DataFrame([ordinary, seed], columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False)
        pd.DataFrame([{"claim_id": row["claim_id"], "source_relation_status": "review_required"}
                      for row in (ordinary, seed)]).to_csv(source, sep="\t", index=False)
        captured = source_verify(str(proposals), str(root / "mixed-source-capture"),
                                 cache_arg=str(root / "source-cache.json"), id_column="claim_id", source_column="source_id")
        self.sources = captured["source_review_records"]
        review = json.loads(self.raw.read_text())["reviews"][0]
        reviews = [review, {**copy.deepcopy(review), "claim_id": "seed-1"}]
        write_json(self.raw, {"reviews": reviews})
        originals = {path: path.read_bytes() for path in (proposals, source, bundle, self.raw, Path(self.sources))}
        session = root / "supported-seed-session"
        prepared = prepare_session(str(session), route="b", inputs={
            "hypotheses": str(proposals), "sources": str(self.sources), "proposal": str(self.proposal),
            "records": str(root / "records.json"), "bundle": str(bundle), "source_verification": str(source)},
            proposer_id="fixture-b-proposer", proposer_model="fixture", case_id="case403")
        script = Path(__file__).resolve().parents[1] / "scripts" / "claim_review_session.py"
        command = [sys.executable, str(script), "complete", "--session", str(session),
                   "--request-sha256", prepared["request_sha256"], "--response", str(self.raw),
                   "--reviewer-id", "fixture-b-reviewer", "--reviewer-model", "fixture"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
        completed = json.loads(result.stdout)
        self.assertEqual("completed", completed["state"])
        self.assertEqual({"requested": 2, "returned": 2, "missing": 0}, completed["review_counts"])
        self.assertEqual(1, completed["verification"]["verified_claim_count"])
        final = pd.read_csv(completed["verification"]["claims"], sep="\t", dtype=str).fillna("").set_index("claim_id")
        self.assertEqual("verified", final.loc["h1", "verification_status"])
        self.assertEqual("not_evaluated", final.loc["seed-1", "verification_status"])
        self.assertEqual("mechanistic_review_seed", final.loc["seed-1", "claim_type"])
        self.assertIn("semantic support promotion prohibited", final.loc["seed-1", "verification_reason"])
        candidates = pd.read_csv(completed["candidates"]["candidates"], sep="\t", dtype=str)
        self.assertEqual(["h1"], candidates["claim_ids"].tolist())
        review_candidates = pd.read_csv(completed["candidates"]["review_candidates"], sep="\t", dtype=str)
        self.assertEqual(["seed-1"], review_candidates["claim_id"].tolist())
        manifest = json.loads(Path(completed["verification"]["manifest"]).read_text())
        audit = json.loads(Path(manifest["agent_claim_review"]["audit"]).read_text())
        audited = {row["claim"]["claim_id"]: row for row in audit["claims"]}
        self.assertEqual(reviews[1], audited["seed-1"]["review"])
        self.assertEqual("supported", audited["seed-1"]["review"]["decision"])
        self.assertEqual("verified_complete", audited["seed-1"]["claim"]["binding_state"])
        self.assertEqual(self.raw.read_bytes(), (session / "completed" / "raw_response.json").read_bytes())
        for path, original in originals.items():
            self.assertEqual(original, path.read_bytes(), str(path))
        self.assertEqual(completed, json.loads(subprocess.run(command, capture_output=True, text=True, check=True).stdout))

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        base = ReviewFixture(Path(self.tmp.name))
        frame = pd.read_csv(base.hypotheses, sep="\t", dtype=str).fillna("")
        frame = frame.rename(columns={"hypothesis_id": "claim_id"})
        self.claims = Path(self.tmp.name) / "claims.tsv"
        frame.to_csv(self.claims, sep="\t", index=False)
        self.prompt = Path(self.tmp.name) / "b-prompt.json"
        prepare_review_prompt(str(self.claims), str(base.sources_path), str(self.prompt), route="b")
        review = copy.deepcopy(base.review)
        review["claim_id"] = review.pop("hypothesis_id")
        self.raw = write_json(Path(self.tmp.name) / "b-review.json", {"reviews": [review]})
        self.proposal = write_json(Path(self.tmp.name) / "b-proposal.json", {"fixture": "proposal"})
        packet = {
            "contract": B_CONTRACT, "hypotheses_sha256": sha256(self.claims),
            "source_review_records_sha256": sha256(base.sources_path),
            "proposal_run": {"executor_id": "b-proposer", "model": "fixture", "raw_response": spec(self.proposal)},
            "review_run": {"executor_id": "b-reviewer", "model": "fixture", "raw_response": spec(self.raw), "prompt": spec(self.prompt)},
        }
        self.packet = write_json(Path(self.tmp.name) / "b-packet.json", packet)
        self.record = base.record
        self.sources = base.sources_path

    def test_positive_claim_id_route_b_binding(self):
        checker = AgentClaimReview(str(self.packet), str(self.claims), str(self.sources), route="b")
        row = checker.bind(pd.read_csv(self.claims, sep="\t", dtype=str).iloc[0].to_dict(), "P12345", self.record, ["PMID:12345"])
        self.assertEqual("b", row["route"])
        self.assertEqual("h1", row["claim_id"])
        self.assertEqual("verified_complete", row["binding_state"])

    def test_route_specific_input_id_is_required(self):
        from candidate_qtl.agent_claim_review import CONTRACT

        for route, wrong_id, required_id in (
            ("a", "claim_id", "hypothesis_id"),
            ("b", "hypothesis_id", "claim_id"),
        ):
            with self.subTest(route=route):
                frame = pd.read_csv(self.claims, sep="\t", dtype=str)
                frame = frame.rename(columns={frame.columns[0]: wrong_id})
                wrong_input = Path(self.tmp.name) / f"wrong-{route}.tsv"
                frame.to_csv(wrong_input, sep="\t", index=False)
                prepare_review_prompt(str(wrong_input), str(self.sources), str(self.prompt), route=route)
                packet = json.loads(self.packet.read_text())
                packet["contract"] = B_CONTRACT if route == "b" else CONTRACT
                packet["hypotheses_sha256"] = sha256(wrong_input)
                packet["review_run"]["prompt"] = spec(self.prompt)
                write_json(self.packet, packet)
                with self.assertRaisesRegex(ValueError, f"requires {required_id}"):
                    AgentClaimReview(str(self.packet), str(wrong_input), str(self.sources), route=route)

    def test_wrong_provider_species_stays_review_required(self):
        checker = AgentClaimReview(str(self.packet), str(self.claims), str(self.sources), route="b")
        record = copy.deepcopy(self.record)
        record["organism"]["scientificName"] = "Pisum sativum"
        row = checker.bind(pd.read_csv(self.claims, sep="\t", dtype=str).iloc[0].to_dict(), "P12345", record, ["PMID:12345"])
        self.assertNotEqual("verified_complete", row["binding_state"])

    def _cli_inputs(self, root: Path, citation: str = "PMID:12345"):
        claim = {column: "" for column in CLAIM_COLUMNS}
        claim.update({"claim_id": "h1", "route": "interval_first", "semantic_unit_id": "unit-1",
                      "family_id": "EV1", "member_interval_genes": "G1", "claim_type": "mutant_phenotype",
                      "reference_gene": "P12345", "reference_species": "Arabidopsis thaliana",
                      "evidence_scope": "family", "source_id": citation, "verification_status": "proposed",
                      "verification_reason": "pending", "rank_contribution": "none"})
        proposals = root / "proposals.tsv"
        pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False)
        bundle = root / "bundle.json"
        _bundle(bundle, members=[("G1", "1:1")])
        payload = json.loads(bundle.read_text(encoding="utf-8"))
        payload["units"][0]["literature_link_candidates"][0]["reference_gene"] = "P12345"
        payload["units"][0]["genes"][0]["orthology_reference_gene"] = "P12345"
        bundle.write_text(json.dumps(payload), encoding="utf-8")
        source = root / "source.tsv"
        pd.DataFrame([{"claim_id": "h1", "source_relation_status": "review_required"}]).to_csv(source, sep="\t", index=False)
        source_records = Path(self.sources)
        prepare_review_prompt(str(proposals), str(source_records), str(self.prompt), route="b")
        packet = json.loads(self.packet.read_text(encoding="utf-8"))
        packet["contract"] = B_CONTRACT
        packet["hypotheses_sha256"] = sha256(proposals)
        packet["review_run"]["prompt"] = spec(self.prompt)
        packet["review_run"]["raw_response"] = spec(self.raw)
        write_json(self.packet, packet)
        return proposals, source, bundle

    def test_transfer_cli_independent_positive_and_unrelated_provider_citation_negative(self):
        proposals, source, bundle = self._cli_inputs(Path(self.tmp.name))
        out = Path(self.tmp.name) / "cli-positive"
        args = ["route-b-semantic-transfer-verify", "--proposals", str(proposals), "--source-verification", str(source),
                "--bundle", str(bundle), "--agent-claim-review", str(self.packet), "--source-review-records", str(self.sources),
                "--record-cache", str(self.tmp.name), "--out-dir", str(out), "--case-id", "case403", "--json"]
        result = main(args)
        self.assertEqual(0, result)
        final = pd.read_csv(out / "semantic_claims_verified.tsv", sep="\t", dtype=str).fillna("")
        self.assertEqual("verified", final.iloc[0]["verification_status"])
        record = copy.deepcopy(self.record)
        record["references"] = []
        write_json(Path(self.tmp.name) / "records.json", record)
        out_bad = Path(self.tmp.name) / "cli-negative"
        bad_args = ["route-b-semantic-transfer-verify", "--proposals", str(proposals), "--source-verification", str(source),
                    "--bundle", str(bundle), "--agent-claim-review", str(self.packet), "--source-review-records", str(self.sources),
                    "--record-cache", str(self.tmp.name), "--out-dir", str(out_bad), "--case-id", "case403", "--json"]
        self.assertEqual(0, main(bad_args))
        final_bad = pd.read_csv(out_bad / "semantic_claims_verified.tsv", sep="\t", dtype=str).fillna("")
        self.assertNotEqual("verified", final_bad.iloc[0]["verification_status"])

    def test_guided_cli_keeps_route_b_packet_separate_from_route_a(self):
        parser = build_parser("core")
        action = next(a for a in parser._actions if getattr(a, "dest", "") == "command")
        guided = action.choices["guided-profile"]
        dests = {a.dest for a in guided._actions}
        self.assertIn("agent_claim_review", dests)
        self.assertIn("route_b_agent_claim_review", dests)
        self.assertIn("route_b_source_review_records", dests)
        self.assertNotIn("source_review_records", dests)

    def test_cli_species_resolution_requires_exact_gene_level_taxon(self):
        root = Path(self.tmp.name)
        proposals, source, bundle = self._cli_inputs(root)
        frame = pd.read_csv(proposals, sep="\t", dtype=str).fillna("")
        frame["reference_species"] = ""
        frame.to_csv(proposals, sep="\t", index=False)
        prepare_review_prompt(str(proposals), str(self.sources), str(self.prompt), route="b")
        packet = json.loads(self.packet.read_text())
        packet["hypotheses_sha256"] = sha256(proposals)
        packet["review_run"]["prompt"] = spec(self.prompt)
        write_json(self.packet, packet)
        record = copy.deepcopy(self.record)
        record["organism"]["taxonId"] = 3702
        write_json(root / "records.json", record)
        original_bundle = json.loads(bundle.read_text())
        for name, taxon, endpoint, expected in [
            ("matching", "3702", "P12345", "verified"),
            ("wrong_taxon", "3880", "P12345", "not_evaluated"),
            ("missing_taxon", "", "P12345", "not_evaluated"),
            ("wrong_endpoint", "3702", "Q99999", "not_evaluated"),
        ]:
            with self.subTest(name=name):
                payload = copy.deepcopy(original_bundle)
                gene = payload["units"][0]["genes"][0]
                for key in ("orthology_reference_species", "reference_species", "orthology_species"):
                    gene[key] = ""
                gene["orthology_reference_gene"] = endpoint
                gene["reference_taxon_id"] = taxon
                write_json(bundle, payload)
                out = root / name
                self.assertEqual(0, main([
                    "route-b-semantic-transfer-verify", "--proposals", str(proposals),
                    "--source-verification", str(source), "--bundle", str(bundle),
                    "--agent-claim-review", str(self.packet), "--source-review-records", str(self.sources),
                    "--record-cache", str(root / "records.json"), "--out-dir", str(out), "--json",
                ]))
                final = pd.read_csv(out / "semantic_claims_verified.tsv", sep="\t").iloc[0]
                self.assertEqual(expected, final["verification_status"])
                audit = json.loads((out / "agent_claim_review_audit.json").read_text())
                if expected == "verified":
                    self.assertEqual("frozen_reference_endpoint_and_exact_provider_taxon",
                                     audit["claims"][0]["reference_species_resolution"]["basis"])


if __name__ == "__main__":
    unittest.main()
