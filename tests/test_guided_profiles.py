from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from candidate_qtl import discover as discover_module
from candidate_qtl.cli import main
from candidate_qtl.guided_profiles import resolve_final_manifest


FIXTURES = Path(__file__).parent / "fixtures"

PAYLOAD_SUBSET = "\n".join([
    "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
    "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
    "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
    "meristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
    "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
    "P30924\tArabidopsis thaliana\tCLV1\treviewed receptor kinase with a fasciation phenotype",
    "hormone transport\tPIN\tauxin redistribution changes stem geometry\tPMID:1978859\t"
    "meristem lesions give the same shape\tassay auxin flux\tArabidopsis pin mutant\tfused organs\t"
    "Q9SXK1\tArabidopsis thaliana\tPIN1\tcanonical auxin efflux carrier",
    "cell wall biosynthesis\tCESA\tweakened wall changes stem shape\tPMID:12172025\t"
    "wall phenotypes rarely fasciate\tsection the stem\tArabidopsis irx mutant\tcollapsed xylem\t"
    "Q9SWW6\tArabidopsis thaliana\tCESA9\tsecondary wall cellulose synthase",
]) + "\n"


def _write_screening_inputs(root: Path) -> tuple[Path, Path]:
    fasta, gff = root / "genome.fa", root / "genes.gff3"
    fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\nchr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=Gene1;Note=meristem regulator\n"
        "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=Gene1.1;Parent=Gene1\n"
        "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=Gene1.1\n", encoding="utf-8",
    )
    return fasta, gff


class GuidedProfileTests(unittest.TestCase):
    def test_every_profile_emits_explicit_dry_run_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for profile in ("screening", "full-case", "audit-benchmark"):
                out = root / profile
                self.assertEqual(0, main(["guided-profile", "--profile", profile, "--out-dir", str(out)]))
                manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
                self.assertTrue(manifest["dry_run"])
                self.assertIn("never auto-filled", manifest["review_policy"])
                self.assertIn("no profile-specific", manifest["ranking_policy"])
                for stage in manifest["stages"]:
                    self.assertIn("command", stage)
                    self.assertIn("arguments", stage)
                    self.assertIn("candidate_qtl_version", stage)
                    self.assertIn("external_transmission", stage)
                    self.assertIn("optional_tools", stage)
                    self.assertIn("stop_condition", stage)

    def test_screening_executes_offline_and_promotes_only_to_review(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"
            fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\nchr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=Gene1;Note=meristem regulator\n"
                "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=Gene1.1;Parent=Gene1\n"
                "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=Gene1.1\n", encoding="utf-8",
            )
            screening = root / "screening_profile"
            self.assertEqual(0, main(["guided-profile", "--profile", "screening", "--out-dir", str(screening), "--execute",
                                      "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:90-350", "--phenotype", "meristem enlargement"]))
            discover_manifest = json.loads((screening / "screening" / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(discover_manifest["external_sequence_submission_this_run"])
            promoted = root / "full_case_profile"
            self.assertEqual(0, main(["guided-profile", "--profile", "full-case", "--out-dir", str(promoted), "--execute",
                                      "--screening-dir", str(screening / "screening")]))
            self.assertIn("REVIEW_REQUIRED", (promoted / "promotion" / "case.yaml").read_text(encoding="utf-8"))
            profile_manifest = json.loads((promoted / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("completed", profile_manifest["execution_results"][0]["status"])

    def test_audit_profile_is_answerless_prediction_alias(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tannotation\nG_SBE\tstarch branching enzyme I\n", encoding="utf-8")
            route = root / "route.tsv"
            route.write_text(
                "phenotype_token\tpathway_family\tannotation_keywords\n"
                "starch\tstarch_biosynthesis\tstarch|amylose|branching enzyme\n",
                encoding="utf-8",
            )
            predict = root / "predict.tsv"
            predict.write_text(
                "case_id\tphenotype_text\tinterval\tinterval_genes\n"
                "case_starch\twrinkled seed with reduced starch\tchr1:1-1000\tgenes.tsv\n",
                encoding="utf-8",
            )
            out = root / "audit"
            self.assertEqual(0, main(["guided-profile", "--profile", "audit-benchmark", "--out-dir", str(out), "--execute",
                                      "--predict-manifest", str(predict), "--blind-route", str(route)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            status = {row["command"]: row["status"] for row in manifest["execution_results"]}
            self.assertTrue(manifest["prediction_only"])
            self.assertFalse(manifest["truth_inputs_accepted"])
            self.assertEqual("completed", status["audit-predict"])
            self.assertNotIn("blind-benchmark", status)
            self.assertEqual("not_runnable_missing_input", status["query-eligibility"])
            self.assertTrue((out / "06_prediction_freeze" / "freeze_manifest.json").exists())

    def test_audit_profile_rejects_legacy_answer_carrying_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "audit"
            fixture = FIXTURES / "blind_benchmark"
            self.assertEqual(1, main(["guided-profile", "--profile", "audit-benchmark", "--out-dir", str(out), "--execute",
                                      "--blind-manifest", str(fixture / "manifest.tsv"), "--blind-route", str(fixture / "route.tsv")]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            records = {row["command"]: row for row in manifest["execution_results"]}
            self.assertEqual("failed", records["audit-predict"]["status"])
            self.assertIn("answer_gene_id", records["audit-predict"]["error"])
            self.assertFalse((out / "02_blind_benchmark").exists())


class GuidedProfileProviderParityTests(unittest.TestCase):
    """Guided screening must select the same annotation provider as `discover` (#159)."""

    def _provider(self, manifest_path: Path) -> str:
        return str(json.loads(manifest_path.read_text(encoding="utf-8"))["annotation_provider"])

    def test_guided_screening_uses_the_configured_local_db_like_direct_discover(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _write_screening_inputs(root)
            local_db = root / "local_pep"
            with mock.patch.object(discover_module, "_auto_local_protein_db", return_value=local_db), \
                 mock.patch.object(discover_module, "_best_local_hits", return_value={}):
                direct = root / "direct"
                self.assertEqual(0, main(["discover", "--offline", "--fasta", str(fasta), "--gff3", str(gff),
                                          "--interval", "chr1:90-350", "--phenotype", "meristem enlargement",
                                          "--out-dir", str(direct)]))
                guided = root / "guided"
                self.assertEqual(0, main(["guided-profile", "--profile", "screening", "--out-dir", str(guided), "--execute",
                                          "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:90-350",
                                          "--phenotype", "meristem enlargement"]))
            direct_provider = self._provider(direct / "discover_manifest.json")
            guided_provider = self._provider(guided / "screening" / "discover_manifest.json")
            self.assertEqual("local_blastp", direct_provider)
            self.assertEqual(direct_provider, guided_provider)
            profile_manifest = json.loads((guided / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            record = profile_manifest["execution_results"][0]
            self.assertEqual("local_blastp", record["annotation_provider"])
            self.assertEqual(str(local_db), record["local_protein_db"])
            self.assertIn("auto-discover", profile_manifest["stages"][0]["provider_policy"])

    def test_no_auto_local_db_disables_discovery_explicitly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _write_screening_inputs(root)
            with mock.patch.object(discover_module, "_auto_local_protein_db", return_value=root / "local_pep"), \
                 mock.patch.object(discover_module, "_best_local_hits", return_value={}):
                guided = root / "guided"
                self.assertEqual(0, main(["guided-profile", "--profile", "screening", "--out-dir", str(guided), "--execute",
                                          "--no-auto-local-db", "--fasta", str(fasta), "--gff3", str(gff),
                                          "--interval", "chr1:90-350", "--phenotype", "meristem enlargement"]))
            self.assertEqual("gff3_only", self._provider(guided / "screening" / "discover_manifest.json"))
            profile_manifest = json.loads((guided / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(profile_manifest["stages"][0]["arguments"]["auto_local_db"])
            self.assertIn("disabled", profile_manifest["stages"][0]["provider_policy"])

    def test_explicit_protein_db_is_forwarded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _write_screening_inputs(root)
            explicit = root / "explicit_pep"
            with mock.patch.object(discover_module, "_auto_local_protein_db", return_value=None), \
                 mock.patch.object(discover_module, "_best_local_hits", return_value={}):
                guided = root / "guided"
                self.assertEqual(0, main(["guided-profile", "--profile", "screening", "--out-dir", str(guided), "--execute",
                                          "--protein-db", str(explicit), "--fasta", str(fasta), "--gff3", str(gff),
                                          "--interval", "chr1:90-350", "--phenotype", "meristem enlargement"]))
            manifest = json.loads((guided / "screening" / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("local_blastp", manifest["annotation_provider"])
            self.assertEqual(str(explicit), manifest["local_protein_db"])


class GuidedProfileFailureIsolationTests(unittest.TestCase):
    """A failed stage must not erase the manifest or cancel independent lanes (#160)."""

    def test_failed_audit_lane_does_not_cancel_independent_lanes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            broken = root / "broken_query_manifest.tsv"
            broken.write_text("not_a_valid\tquery_manifest\n", encoding="utf-8")
            fixture = FIXTURES / "blind_benchmark"
            out = root / "audit"
            self.assertEqual(1, main(["guided-profile", "--profile", "audit-benchmark", "--out-dir", str(out), "--execute",
                                      "--query-manifest", str(broken),
                                      "--blind-manifest", str(fixture / "manifest.tsv"), "--blind-route", str(fixture / "route.tsv")]))
            manifest_path = out / "guided_profile_manifest.json"
            self.assertTrue(manifest_path.exists())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            records = {row["command"]: row for row in manifest["execution_results"]}
            self.assertEqual("failed", records["query-eligibility"]["status"])
            self.assertTrue(records["query-eligibility"]["error_type"])
            self.assertTrue(records["query-eligibility"]["started_at"])
            self.assertTrue(records["query-eligibility"]["finished_at"])
            # The legacy answer-carrying fixture is rejected by the independent
            # prediction lane instead of being executed as a benchmark.
            self.assertEqual("failed", records["audit-predict"]["status"])
            self.assertIn("answer_gene_id", records["audit-predict"]["error"])
            self.assertFalse((out / "02_blind_benchmark").exists())
            self.assertEqual("partial", manifest["overall_status"])

    def test_fail_fast_profile_still_writes_a_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _write_screening_inputs(root)
            out = root / "screening_profile"
            with mock.patch.object(discover_module, "_auto_local_protein_db", side_effect=RuntimeError("db probe exploded")):
                self.assertEqual(1, main(["guided-profile", "--profile", "screening", "--out-dir", str(out), "--execute",
                                          "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:90-350",
                                          "--phenotype", "meristem enlargement"]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", manifest["overall_status"])
            record = manifest["execution_results"][0]
            self.assertEqual("failed", record["status"])
            self.assertEqual("RuntimeError", record["error_type"])
            self.assertIn("db probe exploded", record["error"])
            # The planned command and its resolved arguments survive the failure.
            self.assertEqual("discover", manifest["stages"][0]["command"])
            self.assertEqual(str(fasta), manifest["stages"][0]["arguments"]["fasta"])

    def test_dry_run_manifest_reports_overall_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "dry"
            self.assertEqual(0, main(["guided-profile", "--profile", "audit-benchmark", "--out-dir", str(out)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("dry_run", manifest["overall_status"])
            self.assertNotIn("execution_results", manifest)


class GuidedProfileReproducibilityTests(unittest.TestCase):
    """Issue #302: reruns into the same out-dir must leave a machine-determinable
    final manifest instead of silently overwriting the previous record."""

    def _execute_screening(self, root: Path, out: Path, *, fail_probe: bool) -> int:
        fasta, gff = _write_screening_inputs(root)
        args = ["guided-profile", "--profile", "screening", "--out-dir", str(out), "--execute",
                "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:90-350",
                "--phenotype", "meristem enlargement"]
        if fail_probe:
            with mock.patch.object(discover_module, "_auto_local_protein_db",
                                   side_effect=RuntimeError("db probe exploded")):
                return main(args)
        with mock.patch.object(discover_module, "_auto_local_protein_db", return_value=root / "local_pep"), \
             mock.patch.object(discover_module, "_best_local_hits", return_value={}):
            return main(args)

    def test_first_run_manifest_is_final_without_predecessor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(0, self._execute_screening(root, out, fail_probe=False))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["run_id"])
            self.assertIs(True, manifest["final"])
            self.assertNotIn("supersedes", manifest)
            self.assertEqual("completed", manifest["overall_status"])
            self.assertEqual(["guided_profile_manifest.json"],
                             [path.name for path in out.glob("guided_profile_manifest*.json")])
            resolved = resolve_final_manifest(str(out))
            self.assertEqual(str(out / "guided_profile_manifest.json"), resolved["final"])
            self.assertEqual([], resolved["superseded"])
            self.assertEqual([], resolved["unknown"])

    def test_rerun_same_out_dir_marks_previous_manifest_superseded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(0, self._execute_screening(root, out, fail_probe=False))
            first = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(0, self._execute_screening(root, out, fail_probe=False))
            current = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertIs(True, current["final"])
            self.assertNotEqual(first["run_id"], current["run_id"])
            archive = Path(current["supersedes"])
            self.assertTrue(archive.name.startswith("guided_profile_manifest."))
            archived = json.loads(archive.read_text(encoding="utf-8"))
            self.assertEqual(first["run_id"], archived["run_id"])
            self.assertIs(False, archived["final"])
            self.assertEqual(str(out / "guided_profile_manifest.json"), archived["superseded_by"])
            self.assertEqual("completed", archived["overall_status"])
            resolved = resolve_final_manifest(str(out))
            self.assertEqual(str(out / "guided_profile_manifest.json"), resolved["final"])
            self.assertEqual([str(archive)], resolved["superseded"])
            self.assertEqual([], resolved["unknown"])

    def test_failed_attempt_is_not_final_and_rerun_records_supersession(self) -> None:
        """The pea_catalog scenario: a rejected run is replaced by the corrected
        rerun, and the pair is machine-determinable."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(1, self._execute_screening(root, out, fail_probe=True))
            canonical = out / "guided_profile_manifest.json"
            failed = json.loads(canonical.read_text(encoding="utf-8"))
            self.assertIs(False, failed["final"])
            self.assertEqual("failed", failed["overall_status"])
            self.assertEqual(0, self._execute_screening(root, out, fail_probe=False))
            current = json.loads(canonical.read_text(encoding="utf-8"))
            self.assertIs(True, current["final"])
            self.assertNotEqual(failed["run_id"], current["run_id"])
            archive = Path(current["supersedes"])
            self.assertEqual(f"guided_profile_manifest.{failed['run_id']}.json", archive.name)
            archived = json.loads(archive.read_text(encoding="utf-8"))
            self.assertEqual(failed["run_id"], archived["run_id"])
            self.assertIs(False, archived["final"])
            self.assertEqual(str(canonical), archived["superseded_by"])
            resolved = resolve_final_manifest(str(out))
            self.assertEqual(str(canonical), resolved["final"])
            self.assertEqual([str(archive)], resolved["superseded"])
            self.assertEqual([], resolved["unknown"])

    def test_legacy_manifest_without_final_fields_is_unknown_until_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            out.mkdir()
            legacy = {"schema_version": "1.1", "profile": "triple-route", "dry_run": False,
                      "overall_status": "completed", "stages": []}
            legacy_path = out / "guided_profile_manifest.json"
            legacy_path.write_text(json.dumps(legacy), encoding="utf-8")
            legacy_mtime = legacy_path.stat().st_mtime_ns
            resolved = resolve_final_manifest(str(out))
            self.assertIsNone(resolved["final"])
            self.assertEqual([str(legacy_path)], resolved["unknown"])
            self.assertEqual(0, self._execute_screening(root, out, fail_probe=False))
            current = json.loads(legacy_path.read_text(encoding="utf-8"))
            self.assertIs(True, current["final"])
            archive = Path(current["supersedes"])
            self.assertEqual(f"legacy_{legacy_mtime}", archive.stem.split(".", 1)[1])
            archived = json.loads(archive.read_text(encoding="utf-8"))
            self.assertIs(False, archived["final"])
            self.assertEqual(str(legacy_path), archived["superseded_by"])
            resolved = resolve_final_manifest(str(out))
            self.assertEqual(str(legacy_path), resolved["final"])
            self.assertEqual([str(archive)], resolved["superseded"])
            self.assertEqual([], resolved["unknown"])

    def test_dry_run_manifest_is_explicitly_not_final(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "dry"
            self.assertEqual(0, main(["guided-profile", "--profile", "audit-benchmark", "--out-dir", str(out)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertIs(False, manifest["final"])
            self.assertTrue(manifest["run_id"])
            resolved = resolve_final_manifest(str(out))
            self.assertIsNone(resolved["final"])
            self.assertEqual([str(out / "guided_profile_manifest.json")], resolved["not_final"])

    def test_two_final_manifests_resolve_as_ambiguous(self) -> None:
        # Two final=true manifests with no supersession chain are ambiguous: the
        # resolver must not silently pick one.
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "triple"
            out.mkdir()
            first = {"schema_version": "1.2", "profile": "triple-route", "dry_run": False,
                     "run_id": "20260811T000000000000Z", "final": True,
                     "overall_status": "completed", "stages": []}
            second = {"schema_version": "1.2", "profile": "triple-route", "dry_run": False,
                      "run_id": "20260811T000000000001Z", "final": True,
                      "overall_status": "completed", "stages": []}
            (out / "guided_profile_manifest.json").write_text(json.dumps(first), encoding="utf-8")
            (out / "guided_profile_manifest.20260811T000000000001Z.json").write_text(
                json.dumps(second), encoding="utf-8")
            resolved = resolve_final_manifest(str(out))
            self.assertIsNone(resolved["final"])
            self.assertEqual(2, len(resolved["final_candidates"]))


class GuidedProfilePayloadPreflightTests(unittest.TestCase):
    """Issue #380: the mechanism payload is schema-preflighted at plan time and
    a rejection aborts before any stage runs (execute and dry-run alike)."""

    def test_payload_preflight_accepted_recorded_on_dry_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.tsv"
            payload.write_text(PAYLOAD_SUBSET, encoding="utf-8")
            out = root / "triple"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--mechanism-payload", str(payload)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            preflight = manifest["payload_schema_preflight"]
            self.assertTrue(preflight["plan_time_detection"])
            self.assertEqual("accepted", preflight["status"])
            self.assertEqual("canonical_subset", preflight["accepted_schema_version"])
            self.assertEqual(3, preflight["row_count"])
            self.assertEqual(hashlib.sha256(payload.read_bytes()).hexdigest(), preflight["raw_payload_sha256"])
            self.assertEqual([], preflight["errors"])
            self.assertNotIn("freeze_verification", preflight)
            a02 = next(stage for stage in manifest["stages"] if stage["command"] == "mechanism-import")
            self.assertIn("mechanism_freeze_manifest", a02["arguments"])

    def test_payload_rejected_aborts_execute_before_stage_loop(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload_bad.tsv"
            payload.write_text("case_id\tmechanism_id\ttruth_gene\npea_a\tpea_a_M1\tx\n", encoding="utf-8")
            out = root / "triple"
            self.assertEqual(1, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
                                      "--case-id", "pea_a", "--focal-species", "Pisum sativum",
                                      "--phenotype", "anthocyanin pigmentation loss",
                                      "--mechanism-payload", str(payload)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", manifest["overall_status"])
            self.assertNotIn("execution_results", manifest)
            preflight = manifest["payload_schema_preflight"]
            self.assertEqual("rejected", preflight["status"])
            self.assertIn("unknown column(s)", preflight["errors"][0])
            self.assertIn("truth_gene", preflight["errors"][0])
            self.assertFalse((out / "a01_phenotype_investigation").exists())
            self.assertFalse((out / "a02_mechanism_import").exists())
            dry = root / "triple_dry"
            self.assertEqual(1, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(dry),
                                      "--case-id", "pea_a", "--focal-species", "Pisum sativum",
                                      "--phenotype", "anthocyanin pigmentation loss",
                                      "--mechanism-payload", str(payload)]))
            dry_manifest = json.loads((dry / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", dry_manifest["overall_status"])
            self.assertNotIn("execution_results", dry_manifest)
            self.assertEqual("rejected", dry_manifest["payload_schema_preflight"]["status"])

    def test_no_payload_no_preflight_key(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for profile in ("screening", "full-case", "audit-benchmark", "audit-predict", "triple-route"):
                out = root / profile
                self.assertEqual(0, main(["guided-profile", "--profile", profile, "--out-dir", str(out)]))
                manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
                self.assertNotIn("payload_schema_preflight", manifest)
                self.assertTrue(manifest["dry_run"])
                self.assertEqual("dry_run", manifest["overall_status"])


if __name__ == "__main__":
    unittest.main()
