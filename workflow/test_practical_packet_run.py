"""Synthetic staged packet entry-point tests; no model or network calls."""
import csv
import contextlib
import io
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1]), str(HERE.parents[1] / "src")]
import practical_packet_run as run  # noqa: E402
import test_packet_to_canonical as fixtures  # noqa: E402
from tests.test_blast_evidence import MockBlastTools, S  # noqa: E402


class PracticalPacketRunTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        meta = {"case_id": "case900", "focal_species": "Pisum sativum", "phenotype": "fasciation"}
        for lane in ("a", "b"):
            packet = self.root / lane
            packet.mkdir()
            run.write(packet / "phenotype.json", meta)
        b = self.root / "b"
        run.write(b / "case.json", {**meta, "interval": "chr1:1-9999"})
        (b / "candidates.tsv").write_text("gene_id\tchrom\tstart\tend\nG1\tchr1\t1\t400\nG2\tchr1\t500\t900\n")
        (b / "candidate_proteins.fa").write_text(f">G1\n{S}\n>G2\n{S}\n")
        (b / "reference_corpus.tsv").write_text("reference_gene_id\tfunction_description\nREF1\tmeristem receptor\n")
        pd.DataFrame(fixtures.ATLAS_ROWS).to_csv(self.root / "annotations.tsv", sep="\t", index=False)
        pd.DataFrame([{"gene_id": gene, "reference_gene_id": "REF1", "orthogroup": "OG1",
            "relationship": "homolog", "provider": "TEST_ONLY"} for gene in ("G1", "G2")]).to_csv(
                self.root / "orthology.tsv", sep="\t", index=False)
        pd.DataFrame(fixtures.LITERATURE_ROWS).to_csv(self.root / "literature.tsv", sep="\t", index=False)
        run.write(self.root / "payload.json", {"hypotheses": [fixtures._payload_row(proposed_stable_id="Q6UD73")]})
        records = self.root / "records"
        records.mkdir()
        record = fixtures._entry("Q6UD73", "Arabidopsis thaliana")
        record["sequence"]["value"] = S
        run.write(records / "entry.json", record)
        run.write(self.root / "config.json", {
            "candidate_admission_policy": "content_v1",
            "packet_a": "a", "packet_b": "b", "a_payload": "payload.json",
            "record_cache": "records", "annotations": "annotations.tsv", "literature": "literature.tsv",
            "orthology": "orthology.tsv",
            "a_proposer_id": "TEST-a", "a_model": "TEST/model",
            "b_proposer_id": "TEST-b", "b_model": "TEST/model"})

    def prepare(self):
        return run.prepare(self.root / "config.json", self.root / "prepared")

    def test_new_default_initial_candidates_finish_without_independent_review(self):
        config = run.read(self.root / "config.json")
        del config["candidate_admission_policy"]
        config["source_cache"] = "bibliography.json"
        (self.root / "config.json").write_text(json.dumps(config))
        run.write(self.root / "bibliography.json", [{"id": "123", "title": "Synthetic meristem study"},
                                                   {"id": "12345", "title": "Synthetic reference study"}])
        raw = run.read(self.root / "payload.json")
        raw["hypotheses"][0].update(initial_tier="T2", initial_tier_reason="reported similar cellular phenomenon",
            reported_phenotype="fasciated meristem", reported_process="meristem maintenance",
            phenotype_comparison="same_cellular_process")
        (self.root / "payload.json").write_text(json.dumps(raw))
        prepared = self.prepare()
        self.assertEqual("initial_bibliography_v1", prepared["candidate_admission_policy"])
        self.assertIn("initial_tier", Path(prepared["bundle"]["prompts"]).read_text())
        reviewed = self.review(prepared)
        self.assertEqual("ready_for_initial_candidates", reviewed["state"])
        self.assertEqual([], reviewed["independent_review_routes"])
        self.assertEqual({}, reviewed["machine_reviews"])
        before = {p: p.read_bytes() for p in (self.root / "prepared").rglob("*") if p.is_file()}
        stdout = io.StringIO()
        with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"), contextlib.redirect_stdout(stdout):
            self.assertEqual(0, run.main(["finish", "--reviewed", reviewed["manifest"], "--sha256", reviewed["sha256"],
                                          "--out-dir", str(self.root / "initial-finish")]))
        result = json.loads(stdout.getvalue())
        wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t").fillna("")
        self.assertEqual({"G1"}, set(wet.candidate_gene_id))
        self.assertEqual({"G1": "T2_reported_similar_phenomenon"},
                         dict(zip(wet.candidate_gene_id, wet.relevance_tier)))
        wet_manifest = run.read(Path(result["handoff"]["wet"]["wet_candidates_manifest"]))
        trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t").fillna("").set_index("candidate_gene_id")
        self.assertEqual("T4_reported_function_hypothesis", trace.loc["G2", "relevance_tier"])
        self.assertEqual("not_selected_upper_tier_priority", trace.loc["G2", "selection_disposition"])
        self.assertEqual("not_established", trace.loc["G2", "biological_exclusion"])
        scope = run.read(Path(result["handoff"]["wet"]["boundary_scope_manifest"]))
        self.assertEqual({"G1", "G2"}, set(scope["eligible_gene_ids"]))
        self.assertEqual(["G1"], scope["active_gene_ids"])
        with patch("candidate_qtl.boundary_tier_review.DEFAULT_SELECTION_POLICY", "evidence_scoped_v2"), MockBlastTools(
                raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
            legacy = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "initial-v2-finish")
        legacy_wet = pd.read_csv(legacy["handoff"]["wet"]["wet_candidates"], sep="\t")
        self.assertEqual({"G1", "G2"}, set(legacy_wet.candidate_gene_id))
        self.assertEqual("evidence_scoped_v2", legacy["handoff"]["wet"]["first_stage_selection_policy_version"])
        self.assertNotIn("content_check_state", wet.columns)
        for route in ("a", "b"):
            self.assertEqual("not_required_initial_policy", result["reviews"][route]["review_execution_state"])
            self.assertFalse((self.root / "initial-finish" / route / "completed/packet.json").exists())
        self.assertIn("第一段の本候補", Path(result["report"]).read_text())
        self.assertEqual([], run.audit(Path(result["manifest"]), result["sha256"])["independent_review_required_routes"])
        self.assertTrue(all(p.read_bytes() == value for p, value in before.items()))


    def test_prepared_family_index_reaches_interval_match_and_wet(self):
        """Offline packet E2E: a frozen reference-genome index is evaluated, not silently omitted."""
        from candidate_qtl.agent_claim_review import AXES

        source_text = "In Arabidopsis thaliana, loss of Q6UD73 caused fasciation."
        run.write(self.root / "sources.json", {
            "id": "12345", "pmid": "12345", "title": "Q6UD73 meristem",
            "abstractText": source_text,
        })
        config = run.read(self.root / "config.json")
        config["family_index"] = "family_index.tsv"
        config["source_cache"] = "sources.json"
        (self.root / "config.json").write_text(json.dumps(config), encoding="utf-8")
        index = self.root / "family_index.tsv"
        pd.DataFrame([
            {"gene_id": "G1", "assignment_status": "assigned", "best_subject": "Q6UD73",
             "best_subject_symbol": "Q6UD73", "best_subject_stem": "Q6UD73", "family_stems": "LRR-RLK",
             "family_subjects": "Q6UD73", "family_symbols": "LRR-RLK",
             "index_scope": "reference_genome", "chrom": "chr1", "start": "1", "end": "400",
             "best_identity": "100", "top2_bitscore_ratio": "0.2"},
            {"gene_id": "G2", "assignment_status": "assigned", "best_subject": "UNRELATED",
             "best_subject_symbol": "UNRELATED", "best_subject_stem": "UNRELATED", "family_stems": "UNRELATED",
             "family_subjects": "UNRELATED", "family_symbols": "UNRELATED",
             "index_scope": "reference_genome", "chrom": "chr1", "start": "500", "end": "900",
             "best_identity": "90", "top2_bitscore_ratio": "0.2"},
        ]).to_csv(index, sep="\t", index=False)
        manifest = self.root / "family_index_manifest.json"
        run.write(manifest, {"index_scope": "reference_genome",
                             "reference_label": "synthetic_case900",
                             "reference_species": "Pisum sativum"})
        prepared = self.prepare()
        frozen_index = Path(prepared["inputs"]["family_index"])
        frozen_manifest = Path(prepared["inputs"]["family_index_manifest"])
        self.assertEqual(index.read_bytes(), frozen_index.read_bytes())
        self.assertEqual(manifest.read_bytes(), frozen_manifest.read_bytes())
        reviewed = self.review(prepared)
        prompt = run.read(Path(reviewed["sessions"]["a"]["prompt"]))
        source = prompt["source_review_records"]["records"]["PMID:12345"]
        evidence = {
            "hypothesis_id": "H1", "citation": "PMID:12345", "source_record_id": "12345",
            "sections_sha256": source["sections_sha256"], "subject_gene_identity": "Q6UD73",
            "reported_gene_identity": "Q6UD73", "organism": "Arabidopsis thaliana",
            "reported_organism": "Arabidopsis thaliana", "organism_quote": "In Arabidopsis thaliana",
            "section": "abstractText", "quoted_support": source_text,
            "quote_start": 0, "quote_end": len(source_text),
            "experiment_perturbation": "loss of Q6UD73", "reported_phenotype": "fasciation",
            "predicate_relation": "causal_perturbation_supports_target_trait",
            "decision": "supported", "reason": "TEST_ONLY independently reviewed source",
            "axes": dict.fromkeys(AXES, True),
        }
        with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
            finished = self.finish(reviewed, evidence, out_name="family_index_finished")
        handoff = finished["handoff"]
        interval_root = self.root / "family_index_finished" / "handoff" / "a_interval"
        interval_manifest = run.read(interval_root / "route_a_interval_family_match_manifest.json")
        self.assertTrue(interval_manifest["evaluated"])
        self.assertIsNone(interval_manifest["not_evaluated_reason"])
        interval = pd.read_csv(interval_root / "route_a_interval_family_match.tsv",
                               sep="\t", dtype=str).fillna("")
        self.assertEqual("concrete_interval_gene_matched", interval.iloc[0]["family_status"])
        self.assertEqual(["G1"], interval.iloc[0]["interval_gene_ids"].split(";"))
        bridge = pd.read_csv(handoff["a"]["bridge_tsv"], sep="\t", dtype=str).fillna("")
        self.assertIn("G1", set(bridge.loc[bridge["route"].eq("a"), "gene_id"]))
        wet = pd.read_csv(handoff["wet"]["wet_candidates"], sep="\t", dtype=str).fillna("")
        self.assertIn("G1", set(wet["candidate_gene_id"]))
        inputs = run.read(self.root / "family_index_finished" / "handoff" / "inputs.json")
        self.assertEqual(run.digest(frozen_index), inputs["inputs"]["family_index"]["sha256"])
        self.assertEqual(run.digest(frozen_manifest), inputs["inputs"]["family_index_manifest"]["sha256"])
        self.assertEqual("explicit_cli", handoff["execution_context"]["family_index_resolution"]["source"])

    def test_finish_cli_consumes_bound_identity_and_rejects_changed_review(self):
        from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
        from tests.test_post_freeze_identity import PostFreezeIdentityTests, descriptor

        identity_root = self.root / "identity"
        identity_root.mkdir()
        fixture_owner = PostFreezeIdentityTests()
        identity = fixture_owner.fixture(identity_root, sequence=S,
            targets={"p1": S, "p2": "A" * len(S)}, complete_targets=False)
        identity["paper"].update(reference_gene_id="Q6UD73", gene_correspondence="homolog")
        for name in ("b/candidates.tsv", "annotations.tsv", "orthology.tsv", "literature.tsv"):
            path = self.root / name
            path.write_text(path.read_text().replace("G1", "g1").replace("G2", "g2"))
        (self.root / "b/candidate_proteins.fa").write_bytes(identity["canonical"].read_bytes())
        for lane in ("a", "b"):
            path = self.root / lane / "phenotype.json"
            metadata = run.read(path)
            metadata["focal_species"] = "Plant alpha"
            path.write_text(json.dumps(metadata), encoding="utf-8")
        metadata = run.read(self.root / "b/case.json")
        metadata["focal_species"] = "Plant alpha"
        (self.root / "b/case.json").write_text(json.dumps(metadata), encoding="utf-8")
        config = run.read(self.root / "config.json")
        del config["candidate_admission_policy"]
        config["source_cache"] = "bibliography.json"
        (self.root / "config.json").write_text(json.dumps(config), encoding="utf-8")
        run.write(self.root / "bibliography.json", [
            run.read(Path(identity["item"]["paper_source"]["path"])),
            {"id": "123", "title": "Synthetic meristem study"}])
        raw = run.read(self.root / "payload.json")
        raw["hypotheses"][0].update(source_ids="PMID:1234567;UniProtKB:Q6UD73",
            **{PAPER_FIELD: normalize_paper_reports([identity["paper"]])})
        (self.root / "payload.json").write_text(json.dumps(raw), encoding="utf-8")
        prepared = self.prepare()
        reviewed = self.review(prepared)
        identity["item"]["proposal_id"] = "H1"
        identity["data"].update(case_id="case900",
            atlas=descriptor(Path(prepared["b0"]["interval_gene_atlas"])),
            paper_input=descriptor(Path(prepared["a_import"]["mechanisms"]), id_column="hypothesis_id"),
            candidate_proteins=descriptor(Path(prepared["inputs"]["packet_b"]) / "candidate_proteins.fa"))
        fixture_owner.save(identity)
        # The assessment must stay at its original location: its review is relative.
        identity["data"]["review"]["path"] = "review.json"
        identity["path"].write_text(json.dumps(identity["data"]), encoding="utf-8")
        frozen = {p: p.read_bytes() for directory in ("prepared", "reviewed")
                  for p in (self.root / directory).rglob("*") if p.is_file()}

        def finish_cli(name, with_identity):
            argv = ["finish", "--reviewed", reviewed["manifest"], "--sha256", reviewed["sha256"],
                    "--out-dir", str(self.root / name)]
            if with_identity:
                argv.extend(["--identity-assessment", str(identity["path"])])
            stdout = io.StringIO()
            with MockBlastTools(raw_hits="Q6UD73\tg1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"), \
                    contextlib.redirect_stdout(stdout):
                self.assertEqual(0, run.main(argv))
            result = json.loads(stdout.getvalue())
            wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t").fillna("")
            return result, wet.set_index("candidate_gene_id")

        old_result, old = finish_cli("without-identity", False)
        result, wet = finish_cli("with-identity", True)
        # #1002 accepts this same-species report's qualified member projection;
        # the optional strict identity assessment is no longer a T1 prerequisite.
        self.assertEqual("T1_reported_focal_trait", old.loc["g1", "relevance_tier"])
        self.assertEqual("T1_reported_focal_trait", wet.loc["g1", "relevance_tier"])
        self.assertEqual(set(old.index), set(wet.index))
        for column in ("support_status", "route_a_supported", "route_b_supported", "route_c_supported"):
            self.assertEqual(old[column].to_dict(), wet[column].to_dict())
        self.assertEqual(Path(old_result["handoff"]["wet"]["frozen_candidate_set"]).read_bytes(),
                         Path(result["handoff"]["wet"]["frozen_candidate_set"]).read_bytes())
        inputs = run.read(self.root / "with-identity/handoff/inputs.json")
        self.assertEqual(descriptor(identity["path"]), inputs["post_freeze_identity_assessment"])
        self.assertTrue(all(path.read_bytes() == data for path, data in frozen.items()))
        self.assertEqual("unknown", json.loads(json.loads(normalize_paper_reports([identity["paper"]]))[0]
                         ["paper_checklist"])["same_gene_as_focal"])
        identity["review"].write_text("changed after review", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "identity assessment file hash mismatch"):
            finish_cli("tampered-identity", True)
        self.assertFalse((self.root / "tampered-identity/handoff/wet/wet_candidates.tsv").exists())

    def test_prepare_snapshots_and_checks_frozen_record_cache_manifest(self):
        records = self.root / "records"
        (records / "entry.json").rename(records / "Q6UD73.json")
        preimport = run.import_a_payload(self.root / "payload.json", self.root / "preimport",
                                         case_id="case900", model="TEST/model")
        manifest = {
            "contract": "route-a-reference-record-cache-1.0.0",
            "cache_usable_for_prepare": True,
            "records_dir": str(records.resolve()),
            "records_sha256": {"Q6UD73": run.digest(records / "Q6UD73.json")},
            "mechanisms_sha256": run.digest(Path(preimport["mechanisms"])),
        }
        run.write(self.root / "record-cache-manifest.json", manifest)
        config = run.read(self.root / "config.json")
        config["record_cache_manifest"] = "record-cache-manifest.json"
        (self.root / "config.json").unlink()
        run.write(self.root / "config.json", config)
        prepared = self.prepare()
        self.assertIn("record_cache_manifest", prepared["inputs"])
        self.assertTrue(Path(prepared["inputs"]["record_cache_manifest"]).is_file())

    def test_prepare_rejects_incomplete_record_cache_manifest(self):
        records = self.root / "records"
        run.write(self.root / "record-cache-manifest.json", {
            "contract": "route-a-reference-record-cache-1.0.0",
            "cache_usable_for_prepare": False,
            "records_dir": str(records.resolve()),
            "records_sha256": {"entry": run.digest(records / "entry.json")},
        })
        config = run.read(self.root / "config.json")
        config["record_cache_manifest"] = "record-cache-manifest.json"
        (self.root / "config.json").unlink()
        run.write(self.root / "config.json", config)
        with self.assertRaisesRegex(ValueError, "not usable"):
            self.prepare()

    def test_prepare_preserves_healthy_partial_cache_and_failure_receipt(self):
        from candidate_qtl.reference_record_cache import build_reference_record_cache
        from candidate_qtl.reference_gene_resolution import reference_gene_resolve

        payload = {"hypotheses": [fixtures._payload_row(
            proposed_stable_id="Q6UD73", reference_gene_name="REF1"),
            fixtures._payload_row(hypothesis_id="H2", reference_gene_name="REF2")]}
        (self.root / "payload.json").write_text(json.dumps(payload), encoding="utf-8")
        imported = run.import_a_payload(self.root / "payload.json", self.root / "cache-input",
                                       case_id="case900", model="TEST/model")
        record = fixtures._entry("Q6UD73", "Arabidopsis thaliana")
        record["sequence"]["value"] = S
        record["genes"] = [{"geneName": {"value": "REF1"}}]

        def search(name, species, timeout):
            if name == "REF2":
                raise TimeoutError("neutral injected endpoint failure")
            return ["Q6UD73"], "resolved", "fixture query", 1

        with patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates", side_effect=search), \
             patch("candidate_qtl.reference_record_cache._fetch_uniprot", return_value=record):
            cached = build_reference_record_cache(imported["mechanisms"], str(self.root / "cache"),
                                                  allow_external_api=True)
        config = run.read(self.root / "config.json")
        config.update(record_cache=cached["records_dir"], record_cache_manifest=cached["manifest"])
        (self.root / "config.json").write_text(json.dumps(config), encoding="utf-8")
        prepared = self.prepare()
        saved = run.read(Path(prepared["inputs"]["record_cache_manifest"]))
        self.assertFalse(saved["cache_complete"])
        self.assertEqual(saved["failed_query_count"], 1)
        self.assertEqual(Path(saved["query_ledger"]).parent, (self.root / "prepared" / "inputs").resolve())
        self.assertEqual(saved["records_sha256"]["Q6UD73"], run.digest(Path(cached["records_dir"]) / "Q6UD73.json"))
        # The copied cache stays consumable after the original input is gone.
        shutil.rmtree(self.root / "cache")
        shutil.rmtree(self.root / "cache-input")
        resolved = reference_gene_resolve(saved["mechanisms"], str(self.root / "retry-resolution"),
                                         record_cache_arg=saved["records_dir"], allow_external_api=False)
        with Path(resolved["ledger"]).open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        healthy = [row for row in rows if row.get("resolved_stable_id") == "Q6UD73"]
        pending = [row for row in rows if row["hypothesis_id"].startswith("H2")]
        self.assertTrue(healthy)
        self.assertTrue(pending)
        self.assertTrue(all(row["evaluation_status"] == "not_evaluated" for row in pending))

    def test_prepare_rejects_partial_cache_without_query_ledger(self):
        config = run.read(self.root / "config.json")
        manifest = {"contract": "route-a-reference-record-cache-1.0.0", "cache_usable_for_prepare": True,
                    "cache_complete": False, "records_dir": str(self.root / "records"),
                    "records_sha256": {"entry": run.digest(self.root / "records" / "entry.json")}}
        run.write(self.root / "record-cache-manifest.json", manifest)
        config["record_cache_manifest"] = "record-cache-manifest.json"
        (self.root / "config.json").write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "partial record_cache requires"):
            self.prepare()

    def review(self, prepared):
        responses = self.root / "b_response.json"
        fixtures.PacketToCanonicalTests()._responses(prepared["bundle"]["prompts"], [{
            "claim_type": "mutant_phenotype_bridge", "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
            "reference_gene": "REF1", "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup", "source_ids": ["PMID:123"],
            "alternative_hypothesis": "related copy", "uncertainty_note": "transfer unresolved",
        }], responses)
        return run.review(Path(prepared["manifest"]), prepared["sha256"], responses, self.root / "reviewed")

    def test_review_cli_snapshots_post_proposal_materials_without_changing_prepared(self):
        # B's grounded reference is present in B0, but its record and paper
        # were not known when the initial A-only caches were frozen.
        for name in ("annotations.tsv", "orthology.tsv", "literature.tsv", "b/reference_corpus.tsv"):
            path = self.root / name
            path.write_text(path.read_text(encoding="utf-8").replace("REF1", "Q9AAA1"), encoding="utf-8")
        literature = pd.read_csv(self.root / "literature.tsv", sep="\t", dtype=str).fillna("")
        papers = []
        for pmid in ("67890", "99999"):
            paper = literature.iloc[0].to_dict()
            paper["source_id"] = "PMID:" + pmid
            papers.append(paper)
        pd.DataFrame(papers).to_csv(self.root / "literature.tsv", sep="\t", index=False)
        run.write(self.root / "initial-sources.json", {"id": "12345", "pmid": "12345",
                                                    "title": "Initial A citation"})
        config_path = self.root / "config.json"
        config = run.read(config_path)
        config["source_cache"] = "initial-sources.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")
        prepared = self.prepare()
        before = {str(p.relative_to(self.root / "prepared")): p.read_bytes()
                  for p in (self.root / "prepared").rglob("*") if p.is_file()}
        responses = self.root / "fresh-b.json"
        fixtures.PacketToCanonicalTests()._responses(prepared["bundle"]["prompts"], [{
            "claim_type": "mutant_phenotype_bridge", "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
            "reference_gene": "Q9AAA1", "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup", "source_ids": ["PMID:67890", "PMID:99999"],
            "alternative_hypothesis": "related copy", "uncertainty_note": "transfer unresolved",
        }], responses)
        fresh_sources = self.root / "fresh-sources"
        fresh_sources.mkdir()
        a_text = "In Arabidopsis thaliana, loss of Q6UD73 caused fasciation."
        b_text = "In Arabidopsis thaliana, loss of Q9AAA1 caused fasciation."
        for pmid, text in (("12345", a_text), ("67890", b_text)):
            run.write(fresh_sources / f"{pmid}.json", {"id": pmid, "pmid": pmid,
                "title": "Synthetic source fixture", "abstractText": text})
        fresh_records = self.root / "fresh-records"
        shutil.copytree(self.root / "records", fresh_records)
        b_record = fixtures._entry("Q9AAA1", "Arabidopsis thaliana")
        b_record["references"][0]["citation"]["citationCrossReferences"][0]["id"] = "67890"
        run.write(fresh_records / "Q9AAA1.json", b_record)
        result = subprocess.run([
            sys.executable, str(HERE / "practical_packet_run.py"), "review",
            "--prepared", prepared["manifest"], "--sha256", prepared["sha256"],
            "--responses", str(responses), "--source-cache", str(fresh_sources),
            "--record-cache", str(fresh_records), "--out-dir", str(self.root / "supplemented"),
        ], capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
        reviewed = json.loads(result.stdout)
        root = Path(reviewed["manifest"]).parent
        self.assertEqual(["a", "b"], reviewed["independent_review_routes"])
        self.assertEqual({}, reviewed["machine_reviews"])
        self.assertEqual(before, {str(p.relative_to(self.root / "prepared")): p.read_bytes()
                                 for p in (self.root / "prepared").rglob("*") if p.is_file()})
        run.load_stage(Path(prepared["manifest"]), prepared["sha256"], "awaiting_b_response")
        for key, original in (("source_cache", fresh_sources), ("record_cache", fresh_records)):
            material = reviewed["review_materials"][key]
            self.assertEqual(str(original.absolute()), material["source"])
            for name, sha in material["files"].items():
                saved = Path(material["snapshot"]) / name
                self.assertEqual(run.digest(original / name), sha)
                self.assertEqual(sha, run.digest(saved))
                self.assertEqual(sha, reviewed["files"][str(saved.relative_to(root))])
        self.assertEqual(prepared["a_resolved"]["mechanisms_resolved"], reviewed["a_resolved_input"]["path"])
        self.assertEqual(run.digest(Path(reviewed["a_resolved_input"]["path"])),
                         reviewed["a_resolved_input"]["sha256"])
        resolved = pd.read_csv(prepared["a_resolved"]["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
        a_sources = pd.read_csv(reviewed["a_sources"]["mechanisms_with_verified_sources"], sep="\t", dtype=str).fillna("")
        pd.testing.assert_frame_equal(resolved, a_sources[resolved.columns])
        a_prompt = run.read(Path(reviewed["sessions"]["a"]["prompt"]))
        b_prompt = run.read(Path(reviewed["sessions"]["b"]["prompt"]))
        self.assertIn(a_text, str(a_prompt["source_review_records"]["records"]["PMID:12345"]))
        self.assertIn(b_text, str(b_prompt["source_review_records"]["records"]["PMID:67890"]))
        self.assertNotIn("PMID:99999", b_prompt["source_review_records"]["records"])
        verification = pd.read_csv(reviewed["b_sources"]["source_verification"], sep="\t", dtype=str).fillna("")
        self.assertEqual({"cache_hit"}, set(verification.loc[verification.source_id.eq("PMID:67890"), "cache_investigation_state"]))
        self.assertEqual({"not_evaluated"}, set(verification.loc[verification.source_id.eq("PMID:99999"), "source_relation_status"]))
        self.assertEqual({"proposed"}, set(verification["final_semantic_verification_status"]))
        for route in ("a", "b"):
            session_root = root / route
            request = run.read(Path(reviewed["sessions"][route]["request"]))
            saved_record = session_root / request["inputs"]["records"] / "Q9AAA1.json"
            self.assertEqual(b_record, run.read(saved_record))
            self.assertEqual(run.digest(saved_record), request["files"][str(saved_record.relative_to(session_root))])
        a_request = run.read(Path(reviewed["sessions"]["a"]["request"]))
        self.assertEqual(Path(prepared["inputs"]["a_payload"]).read_bytes(),
                         (root / "a" / a_request["inputs"]["proposal"]).read_bytes())
        # The new material is self-contained; original acquisition paths are
        # lineage only and are not needed to finish the reviewed stage.
        shutil.rmtree(fresh_sources)
        shutil.rmtree(fresh_records)
        with MockBlastTools(raw_hits=""):
            finished = self.finish(reviewed, out_name="supplemented-finished")
        self.assertEqual("awaiting_independent_reviews", run.read(Path(reviewed["manifest"]))["state"])
        self.assertTrue(Path(finished["manifest"]).is_file())

    def test_review_source_file_rechecks_exact_resolved_a_without_reimport_or_resolution(self):
        prepared = self.prepare()
        self.review(prepared)
        source = self.root / "later-sources.json"
        run.write(source, {"id": "12345", "pmid": "12345", "title": "Later A source"})
        with patch.object(run, "import_a_payload", side_effect=AssertionError("A proposal must remain frozen")), \
             patch("candidate_qtl.reference_gene_resolution.reference_gene_resolve", side_effect=AssertionError("A IDs must remain frozen")), \
             patch("candidate_qtl.source_verify._fetch_europe_pmc", side_effect=AssertionError("network not allowed")), \
             patch("candidate_qtl.source_verify._fetch_ncbi_eutils", side_effect=AssertionError("network not allowed")):
            reviewed = run.review(Path(prepared["manifest"]), prepared["sha256"], self.root / "b_response.json",
                                  self.root / "source-file-reviewed", source_cache=source)
        self.assertEqual(run.digest(source), reviewed["review_materials"]["source_cache"]["files"]["source_cache.json"])
        request = run.read(Path(reviewed["sessions"]["b"]["request"]))
        self.assertNotIn("Q9AAA1.json", request["files"])

    def test_review_material_guards_and_failure_preserve_prepared(self):
        prepared = self.prepare()
        self.review(prepared)
        anchor = Path(prepared["manifest"])
        invalid = self.root / "record-file.json"
        run.write(invalid, {})
        with self.assertRaisesRegex(ValueError, "regular records directory"):
            run.review(anchor, prepared["sha256"], self.root / "b_response.json", self.root / "invalid-output", record_cache=invalid)
        self.assertFalse((self.root / "invalid-output").exists())
        for cache_key in ("source_cache", "record_cache"):
            with self.subTest(cache=cache_key), self.assertRaisesRegex(ValueError, "separate"):
                run.review(anchor, prepared["sha256"], self.root / "b_response.json", self.root / "records/new-out",
                           **{cache_key: self.root / "records"})
        self.assertFalse((self.root / "records/new-out").exists())
        with self.assertRaises(FileExistsError):
            run.review(anchor, prepared["sha256"], self.root / "b_response.json", self.root / "reviewed", record_cache=self.root / "records")
        bad_source = self.root / "bad-source.json"
        bad_source.write_text("invalid JSON", encoding="utf-8")
        with self.assertRaises(json.JSONDecodeError):
            run.review(anchor, prepared["sha256"], self.root / "b_response.json", self.root / "failed-review", source_cache=bad_source)
        run.load_stage(anchor, prepared["sha256"], "awaiting_b_response")

    def test_review_rejects_symlink_materials(self):
        prepared = self.prepare()
        link = self.root / "records-link"
        try:
            link.symlink_to(self.root / "records", target_is_directory=True)
        except OSError:
            self.skipTest("symlink creation unavailable")
        for cache_key in ("source_cache", "record_cache"):
            with self.subTest(cache=cache_key), self.assertRaisesRegex(ValueError, "regular"):
                run.review(Path(prepared["manifest"]), prepared["sha256"], self.root / "absent-response",
                           self.root / "symlink-output", **{cache_key: link})
        regular = self.root / "nested-symlink"
        regular.mkdir()
        (regular / "entry.json").symlink_to(self.root / "records/entry.json")
        with self.assertRaisesRegex(ValueError, "symlink"):
            run.review(Path(prepared["manifest"]), prepared["sha256"], self.root / "absent-response",
                       self.root / "symlink-output", record_cache=regular)
        self.assertFalse((self.root / "symlink-output").exists())

    def finish(self, reviewed, a_review=None, cli=False, route_c=None, route_c_status=None,
               out_name="finished"):
        config = {}
        for route in ("a", "b"):
            prompt = run.read(Path(reviewed["sessions"][route]["prompt"]))
            key = "hypothesis_id" if route == "a" else "claim_id"
            path = self.root / (route + "_review.json")
            rows = [{key: row[key], "decision": "review_required",
                "reason": "TEST_ONLY unresolved transfer"} for row in prompt["hypotheses"]]
            if path.exists():
                path.unlink()
            run.write(path, {"reviews": [a_review] if route == "a" and a_review else rows})
            config[route] = {"response": path.name, "reviewer_id": "TEST-review-" + route, "reviewer_model": "TEST/model"}
        reviews_path = self.root / "reviews.json"
        if reviews_path.exists():
            reviews_path.unlink()
        run.write(reviews_path, config)
        if cli:
            command = [sys.executable, str(HERE / "practical_packet_run.py"), "finish",
                "--reviewed", reviewed["manifest"], "--sha256", reviewed["sha256"],
                "--reviews", str(self.root / "reviews.json"), "--max-adaptive-rounds", "2",
                "--out-dir", str(self.root / out_name)]
            if route_c:
                command.extend(["--route-c-members", str(route_c[0]), "--route-c-manifest", str(route_c[1])])
            if route_c_status:
                command.extend(["--route-c-status", str(route_c_status)])
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(0, result.returncode, result.stderr)
            return json.loads(result.stdout)
        return run.finish(Path(reviewed["manifest"]), reviewed["sha256"], self.root / "reviews.json", self.root / out_name,
                          *(route_c or (None, None)), route_c_status)

    def agentic_boundary_round(self, request_path, root, *, adopt):
        from scripts.agentic_boundary import merge as merge_boundary
        from scripts.agentic_boundary import prepare as prepare_boundary

        root.mkdir(parents=True)
        response_dir = root / "responses"
        response_dir.mkdir()
        prepared = root / "prepared"
        manifest = prepare_boundary(str(request_path), str(prepared), max_candidates=8)
        remaining_adopt = adopt
        proposal_payloads = []
        for shard in manifest["shards"]:
            records = json.loads((prepared / shard["file"]).read_text(encoding="utf-8"))["candidate_input"]
            decisions = []
            for index, row in enumerate(records):
                decision = "adopt" if remaining_adopt > 0 else "defer"
                if remaining_adopt > 0:
                    remaining_adopt -= 1
                decisions.append({
                    "candidate_gene_id": row["candidate_gene_id"],
                    "decision": decision,
                    "reason": f"ordinary-case fixture comparison {index}",
                })
            proposal_payloads.append({
                "case_id": manifest["case_id"],
                "source_wet_sha256": manifest["source_wet_sha256"],
                "shard_id": shard["shard_id"],
                "status": "completed",
                "executor_id": "TEST-boundary-executor",
                "executor_model": "TEST/model",
                "responses": decisions,
            })
        proposal_paths = []
        for payload in proposal_payloads:
            path = response_dir / f"proposal-{payload['shard_id']}.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            proposal_paths.append(path)
        proposal = root / "proposals.tsv"
        merge_boundary(
            str(prepared), [str(path) for path in proposal_paths], str(proposal),
            role="proposal", executor_id="TEST-boundary-executor",
            executor_model="TEST/model",
        )
        proposal_sha = run.digest(proposal)
        decisions = pd.read_csv(proposal, sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
        review_paths = []
        for payload in proposal_payloads:
            review_payload = {
                "case_id": manifest["case_id"],
                "source_wet_sha256": manifest["source_wet_sha256"],
                "proposal_sha256": proposal_sha,
                "shard_id": payload["shard_id"],
                "status": "completed",
                "reviewer_id": "TEST-boundary-reviewer",
                "reviewer_model": "TEST/model",
                "responses": [
                    {
                        "candidate_gene_id": row["candidate_gene_id"],
                        "decision": decisions.loc[row["candidate_gene_id"], "decision"],
                        "reason": f"independent ordinary-case comparison {row['candidate_gene_id']}",
                    }
                    for row in payload["responses"]
                ],
            }
            path = response_dir / f"review-{payload['shard_id']}.json"
            path.write_text(json.dumps(review_payload), encoding="utf-8")
            review_paths.append(path)
        independent = root / "independent-review.tsv"
        merge_boundary(
            str(prepared), [str(path) for path in review_paths], str(independent),
            role="review", reviewer_id="TEST-boundary-reviewer",
            reviewer_model="TEST/model", proposal_path=str(proposal),
        )
        return proposal, independent

    def test_c_members_flow_through_practical_finish_and_report(self):
        from tests._issue664_helpers import manifest_for
        reviewed = self.review(self.prepare())
        members = self.root / "route_c_members.tsv"
        pd.DataFrame([{"candidate_gene_id": "G1", "support_status": "verified_direct_support",
            "record_id": "RC1", "source_id": "PMID:1", "species_distance": "S1_same_species",
            "evidence_directness": "D1_mapped_gene_with_functional_validation",
            "projection_state": "inside_current_qtl", "projection_basis": "exact_focal_id"}]).to_csv(
                members, sep="\t", index=False)
        manifest = manifest_for(members, route="c", case_id="case900")
        with MockBlastTools(raw_hits=""):
            result = self.finish(reviewed, route_c=(members, manifest))
        wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t")
        self.assertTrue(wet.loc[wet.candidate_gene_id.eq("G1"), "route_c_supported"].iloc[0])
        self.assertIn("route_c_supported=1", Path(result["report"]).read_text(encoding="utf-8"))
        self.assertIn("Direct C execution: `completed_members`", Path(result["report"]).read_text(encoding="utf-8"))
        self.assertEqual(run.digest(members), result["route_c_input"]["members_sha256"])

    def test_c_provider_failure_is_recorded_without_losing_ab_wet(self):
        reviewed = self.review(self.prepare())
        status = self.root / "c_status.json"
        run.write(status, {"state": "completed_members", "reason": "unearned completion"})
        with self.assertRaisesRegex(ValueError, "status does not match"):
            run.finish(Path(reviewed["manifest"]), reviewed["sha256"], self.root / "reviews.json",
                       self.root / "not-created", route_c_status=status)
        status = self.root / "c_status_failure.json"
        run.write(status, {"case_id": "case900", "state": "provider_failure", "reason": "retrieval unavailable"})
        with MockBlastTools(raw_hits=""):
            result = self.finish(reviewed, route_c_status=status)
        self.assertEqual("provider_failure", result["route_c_execution"]["state"])
        self.assertTrue(Path(result["handoff"]["wet"]["wet_candidates"]).is_file())
        report = Path(result["report"]).read_text(encoding="utf-8")
        self.assertIn("Direct C execution: `provider_failure`", report)
        self.assertIn("retrieval unavailable", report)

    def test_b_shard_merge_is_consumed_by_existing_review_stage(self):
        from scripts.agentic_route_b import merge, prepare
        prepared = self.prepare()
        prompts = prepared["bundle"]["prompts"]
        shards = prepare(prompts, str(self.root / "b-shards"), max_prompts=2)
        self.assertEqual(sum(len(shard["prompt_ids"]) for shard in shards["shards"]),
                         shards["prompt_count"])
        baseline = self.root / "baseline-response.json"
        fixtures.PacketToCanonicalTests()._responses(prompts, [], baseline)
        responses = run.read(baseline)["candidates"]
        partial = self.root / "b-partial.json"
        run.write(partial, {"case_id": "case900", "responses": responses})
        merged = self.root / "merged-response.json"
        manifest = merge(prompts, [str(partial)], str(merged), executor="TEST-b", model="TEST/model")
        self.assertEqual("complete", manifest["completion_status"])
        reviewed = run.review(Path(prepared["manifest"]), prepared["sha256"], merged,
                              self.root / "reviewed")
        self.assertEqual("awaiting_independent_reviews", reviewed["state"])

    def test_issue758_agentic_b_execution_manifest_reaches_first_wet(self):
        from scripts.agentic_route_b import merge, prepare

        prepared = self.prepare()
        prompts = prepared["bundle"]["prompts"]
        shards = prepare(prompts, str(self.root / "issue758-b-shards"), max_prompts=1)
        self.assertGreaterEqual(len(shards["shards"]), 1)
        baseline = self.root / "issue758-b-partial.json"
        fixtures.PacketToCanonicalTests()._responses(prompts, [{
            "claim_type": "mutant_phenotype_bridge",
            "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling",
            "phenotype_relation": "could explain fasciation",
            "reference_gene": "REF1",
            "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup",
            "source_ids": ["PMID:123"],
            "alternative_hypothesis": "related receptor copy",
            "uncertainty_note": "transfer unresolved",
        }], baseline)
        partial = self.root / "issue758-b-response.json"
        run.write(partial, {
            "case_id": "case900", "executor": "TEST-b", "model": "TEST/model",
            "responses": run.read(baseline)["candidates"],
        })
        merged = self.root / "issue758-b-merged.json"
        execution = merge(prompts, [str(partial)], str(merged),
                          executor="TEST-b", model="TEST/model")
        execution_manifest = Path(merged).with_name(Path(merged).stem + "_manifest.json")
        self.assertEqual("agentic-route-b-execution-1.1.0", execution["contract"])
        self.assertEqual(execution["completed_count"], execution["response_state_counts"]["completed_claims"])
        forged_counts = self.root / "issue758-b-forged-counts.json"
        forged = run.read(execution_manifest)
        forged["response_state_counts"]["missing_response"] += 1
        run.write(forged_counts, forged)
        with self.assertRaisesRegex(ValueError, "state counts"):
            run.review(Path(prepared["manifest"]), prepared["sha256"], merged,
                       self.root / "issue758-forged-counts-reviewed", b_execution_manifest=forged_counts)
        legacy = self.root / "issue758-b-legacy-manifest.json"
        old_manifest = run.read(execution_manifest)
        old_manifest.pop("response_state_counts")
        run.write(legacy, old_manifest)
        legacy_out = self.root / "issue758-legacy-context"
        legacy_out.mkdir()
        legacy_context = run._route_b_execution_context(legacy, prepared, merged, legacy_out)
        self.assertEqual(execution["response_state_counts"], legacy_context["response_state_counts"])
        self.assertNotIn("response_state_counts", run.read(legacy_out / "b_execution_manifest.json"))
        tampered_manifest = self.root / "issue758-b-tampered-manifest.json"
        tampered = run.read(execution_manifest)
        tampered["response_roster"][0]["raw_text_sha256"] = "0" * 64
        run.write(tampered_manifest, tampered)
        with self.assertRaisesRegex(ValueError, "raw response hash"):
            run.review(
                Path(prepared["manifest"]), prepared["sha256"], merged,
                self.root / "issue758-tampered-reviewed", b_execution_manifest=tampered_manifest,
            )
        reviewed = run.review(
            Path(prepared["manifest"]), prepared["sha256"], merged,
            self.root / "issue758-reviewed", b_execution_manifest=execution_manifest,
        )
        self.assertEqual("complete", reviewed["route_b_execution"]["completion_status"])
        self.assertEqual(execution["response_sha256"], reviewed["route_b_execution"]["response_sha256"])
        with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
            result = self.finish(reviewed, out_name="issue758-finished")
        self.assertEqual("first_stage_frozen", result["state"])
        self.assertEqual("complete", result["route_b_execution"]["completion_status"])
        report = Path(result["report"]).read_text(encoding="utf-8")
        self.assertIn("Route B execution: `complete`", report)
        self.assertIn("Route B response states:", report)
        self.assertIn(f"completed_claims={execution['completed_count']}", report)
        self.assertIn("Direct C execution: `not_evaluated_no_input`", report)
        self.assertTrue(Path(result["handoff"]["wet"]["wet_candidates"]).is_file())

    def test_core_fix_mixed_claims_reach_ordinary_review_and_checkpoint_retry(self):
        """CORE_FIX consumer regression with fixed responses, not a scientific benchmark."""
        from Bio.Data.CodonTable import unambiguous_dna_by_name
        from candidate_qtl.claim_review_session import _validate_route_b_execution
        from scripts.agentic_route_b import merge

        codons = {aa: codon for codon, aa in unambiguous_dna_by_name["Standard"].forward_table.items()}
        cds = "".join(codons[aa] for aa in S) + "TAA"
        stride = len(cds) + 10
        (self.root / "claim-schema.fa").write_text(">chr1\n" + (cds + "N" * 10) * 2 + "\n")
        gff = ["##gff-version 3"]
        for index in (1, 2):
            start, end = (index - 1) * stride + 1, (index - 1) * stride + len(cds)
            gene = f"G{index}"
            gff.extend((f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                        f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
                        f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gene}.t1"))
        (self.root / "claim-schema.gff3").write_text("\n".join(gff) + "\n")
        config = run.read(self.root / "config.json")
        config.pop("packet_a")
        config.pop("packet_b")
        config.update(case_id="case900", focal_species="Pisum sativum", phenotype="fasciation",
                      interval=f"chr1:1-{stride * 2}", reference_fasta="claim-schema.fa",
                      annotation_gff3="claim-schema.gff3", reference_storage="snapshot",
                      reference_corpus="b/reference_corpus.tsv")
        run.write(self.root / "claim-schema-config.json", config)
        prepared = run.prepare(self.root / "claim-schema-config.json", self.root / "claim-schema-prepared")
        self.assertIsNotNone(prepared["ordinary_case"])
        prompts = prepared["bundle"]["prompts"]
        prompt_rows = [json.loads(line) for line in Path(prompts).read_text().splitlines() if line.strip()]
        good = {"claim_type": "mutant_phenotype_bridge", "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
                "reference_gene": "REF1", "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "orthogroup", "source_ids": ["PMID:123"],
                "uncertainty_note": "transfer unresolved"}
        bad = {key: value for key, value in good.items() if key != "claim_type"}
        for name, claims in (("mixed", [good, bad]), ("corrected", [good])):
            run.write(self.root / f"claim-schema-{name}.json", {
                "case_id": "case900", "executor": "TEST-b", "model": "TEST/model",
                "responses": [{"prompt_id": row["prompt_id"], "raw_text": json.dumps({"claims": claims})}
                              for row in prompt_rows],
            })
        checkpoint = self.root / "claim-schema-checkpoint.json"
        initial = merge(prompts, [str(self.root / "claim-schema-mixed.json")], str(checkpoint))
        self.assertEqual(("partial", 0, len(prompt_rows)),
                         (initial["completion_status"], initial["completed_count"], initial["not_evaluated_count"]))
        original_checkpoint = checkpoint.read_bytes()
        resumed = self.root / "claim-schema-resumed.json"
        retried = merge(prompts, [str(checkpoint), str(self.root / "claim-schema-corrected.json")], str(resumed))
        self.assertEqual("complete", retried["completion_status"])
        self.assertEqual(len(prompt_rows), len(retried["retry_history"]))
        self.assertEqual(original_checkpoint, checkpoint.read_bytes())

        for name, response, state, expected_rejects in (("mixed", checkpoint, "partial", len(prompt_rows)),
                                                        ("resumed", resumed, "complete", 0)):
            execution_path = response.with_name(response.stem + "_manifest.json")
            result = subprocess.run([
                sys.executable, str(HERE / "practical_packet_run.py"), "review",
                "--prepared", prepared["manifest"], "--sha256", prepared["sha256"],
                "--responses", str(response), "--b-execution-manifest", str(execution_path),
                "--out-dir", str(self.root / f"claim-schema-{name}-reviewed"),
            ], capture_output=True, text=True)
            self.assertEqual(0, result.returncode, result.stderr)
            reviewed = json.loads(result.stdout)
            validated = _validate_route_b_execution(reviewed["route_b_execution"], "case900",
                                                    bundle_path=prepared["bundle"]["bundle"])
            self.assertEqual(state, validated["state"])
            self.assertEqual(len(prompt_rows), reviewed["b_import"]["claim_row_count"])
            self.assertEqual(expected_rejects, reviewed["b_import"]["rejected_claim_row_count"])
            if name == "mixed":
                with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
                    finished = self.finish(reviewed, out_name="claim-schema-mixed-finished")
                self.assertEqual("partial", finished["route_b_execution"]["completion_status"])
                self.assertIn("Route B execution: `partial`", Path(finished["report"]).read_text())

            # Rehashing self-reported completion flags/counts cannot alter the
            # execution verdict derived from the same raw model claims.
            forged_payload = run.read(response)
            for row in forged_payload["candidates"]:
                row["response_completion_status"] = "complete" if state == "partial" else "partial"
                row["rejected_claims"] = []
            forged_response = self.root / f"claim-schema-{name}-forged.json"
            run.write(forged_response, forged_payload)
            forged = run.read(execution_path)
            forged.update(output=str(forged_response), response_sha256=run.digest(forged_response),
                          completed_count=len(prompt_rows) if state == "partial" else 0,
                          not_evaluated_count=0 if state == "partial" else len(prompt_rows),
                          completion_status="complete" if state == "partial" else "partial")
            forged_path = self.root / f"claim-schema-{name}-forged-manifest.json"
            run.write(forged_path, forged)
            with self.assertRaisesRegex(ValueError, "counts do not match response rows"):
                run._route_b_execution_context(forged_path, prepared, forged_response, self.root)
            with self.assertRaisesRegex(ValueError, "counts do not match response rows"):
                _validate_route_b_execution({**forged, "artifact_path": str(forged_path),
                                             "artifact_sha256": run.digest(forged_path)}, "case900",
                                            bundle_path=prepared["bundle"]["bundle"])

    def test_route_b_execution_manifest_rejects_prompt_unit_swap(self):
        from practical_packet_run import _route_b_execution_context

        prompts = self.root / "swap-prompts.jsonl"
        prompts.write_text(
            '{"prompt_id":"p1","route_b_unit_id":"u1","case_id":"case-x"}\n'
            '{"prompt_id":"p2","route_b_unit_id":"u2","case_id":"case-x"}\n',
            encoding="utf-8",
        )
        response = self.root / "swap-response.json"
        run.write(response, {"case_id": "case-x", "candidates": [
            {"prompt_id": "p1", "route_b_unit_id": "u1", "parse_status": "json", "raw_text": '{"claims":[]}', "model_json": {"claims": []}},
            {"prompt_id": "p2", "route_b_unit_id": "u2", "parse_status": "json", "raw_text": '{"claims":[]}', "model_json": {"claims": []}},
        ]})
        raw_hash = hashlib.sha256(b'{"claims":[]}').hexdigest()
        manifest = {
            "contract": "agentic-route-b-execution-1.1.0", "case_id": "case-x",
            "prompt_count": 2, "completion_status": "complete", "completed_count": 2,
            "not_evaluated_count": 0, "prompt_sha256": run.digest(prompts),
            "response_sha256": run.digest(response),
            "partial_response_provenance": [{"attempt_index": 0, "file": str(response.resolve()), "sha256": run.digest(response), "executor": "TEST", "model": "fixture"}],
            "retry_history": [],
            "response_roster": [
                {"prompt_id": "p1", "route_b_unit_id": "u1", "parse_status": "json", "raw_text_sha256": raw_hash},
                {"prompt_id": "p2", "route_b_unit_id": "u2", "parse_status": "json", "raw_text_sha256": raw_hash},
            ],
        }
        good_manifest = self.root / "swap-manifest-good.json"
        run.write(good_manifest, manifest)
        good_out = self.root / "swap-good"
        good_out.mkdir()
        _route_b_execution_context(
            good_manifest, {"case_id": "case-x", "bundle": {"prompts": str(prompts)}}, response, good_out,
        )
        manifest["response_roster"][0]["route_b_unit_id"] = "u2"
        bad_manifest = self.root / "swap-manifest-bad.json"
        run.write(bad_manifest, manifest)
        bad_out = self.root / "swap-bad"
        bad_out.mkdir()
        with self.assertRaisesRegex(ValueError, "mapping"):
            _route_b_execution_context(
                bad_manifest, {"case_id": "case-x", "bundle": {"prompts": str(prompts)}}, response, bad_out,
            )

    def test_failed_provider_attempt_can_resume_from_same_frozen_review_stage(self):
        reviewed = self.review(self.prepare())
        with patch("candidate_qtl.blast_evidence.run_panel_blastp", side_effect=RuntimeError("local provider unavailable")):
            with self.assertRaisesRegex(RuntimeError, "local provider unavailable"):
                self.finish(reviewed, out_name="provider-failed")
        failed = self.root / "provider-failed" / "handoff" / "failure.json"
        self.assertTrue(failed.is_file())
        self.assertIn("local provider unavailable", failed.read_text(encoding="utf-8"))
        run.load_stage(Path(reviewed["manifest"]), reviewed["sha256"], "awaiting_independent_reviews")
        with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
            resumed = self.finish(reviewed, out_name="provider-resumed")
        self.assertEqual("first_stage_frozen", resumed["state"])
        self.assertTrue(Path(resumed["handoff"]["wet"]["wet_candidates"]).is_file())
        self.assertTrue(pd.read_csv(resumed["handoff"]["b"]["review_candidates"], sep="\t").empty)
        self.assertEqual({"G1"}, set(pd.read_csv(resumed["handoff"]["wet"]["wet_candidates"], sep="\t").candidate_gene_id))

    def test_no_source_b_review_is_marked_unexecuted(self):
        from scripts.agentic_no_source_review import build
        reviewed = self.review(self.prepare())
        prompt = Path(reviewed["sessions"]["b"]["prompt"])
        response = build(prompt, "b")
        self.assertEqual("not_executed", response["execution_status"])
        self.assertEqual({"not_evaluated"}, {row["decision"] for row in response["reviews"]})

    def test_empty_source_auto_finish_matches_explicit_review_and_preserves_freezes(self):
        reviewed = self.review(self.prepare())
        self.assertEqual([], reviewed["independent_review_routes"])
        self.assertEqual({"a", "b"}, set(reviewed["machine_reviews"]))
        before = {p: p.read_bytes() for folder in ("prepared", "reviewed")
                  for p in (self.root / folder).rglob("*") if p.is_file()}
        explicit = {}
        for route, spec in reviewed["machine_reviews"].items():
            payload = run.read(self.root / "reviewed" / spec["response"])
            response = self.root / (route + "-explicit.json")
            run.write(response, {"reviews": payload["reviews"]})
            explicit[route] = {"response": response.name, "reviewer_id": "TEST-review-" + route,
                               "reviewer_model": "TEST/model"}
        specs = self.root / "explicit.json"
        run.write(specs, explicit)
        raw = "".join(f"Q6UD73\t{gene}\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n" for gene in ("G1", "G2"))
        with MockBlastTools(raw_hits=raw):
            manual = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], specs, self.root / "manual")
            automatic = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "auto")
        for route in ("a", "b"):
            self.assertEqual("executed", manual["reviews"][route]["review_execution_state"])
            completed = automatic["reviews"][route]
            self.assertEqual("not_executed_no_source", completed["review_execution_state"])
            self.assertEqual(0, completed["review_counts"]["returned"])
            self.assertEqual("not_invoked", completed["identity"]["reviewer_model"])
            for result in (manual, automatic):
                review = pd.read_csv(result["handoff"][route]["review_candidates"], sep="\t")
                self.assertEqual({"G1", "G2"} if route == "a" else set(), set(review.candidate_gene_id))
        left = pd.read_csv(manual["handoff"]["wet"]["wet_candidates"], sep="\t").fillna("")
        right = pd.read_csv(automatic["handoff"]["wet"]["wet_candidates"], sep="\t").fillna("")
        # A's projected-row fingerprint binds the distinct review provenance.
        # Compare every scientific field; preserve those hashes in each output.
        projection_fingerprints = []
        for table in (left, right):
            fingerprints = {}
            def scientific_evidence(value):
                rows = json.loads(value)
                for row in rows:
                    if row.get("membership_basis") == "route_a_review_projection":
                        row.pop("input_fingerprint")
                    if row.get("homolog_projection_json"):
                        projection = json.loads(row["homolog_projection_json"])
                        fingerprint = projection["evidence"].pop("input_fingerprint")
                        self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")
                        self.assertEqual(row["candidate_gene_id"], projection["candidate_gene_id"])
                        fingerprints.setdefault(row["candidate_gene_id"], set()).add(fingerprint)
                        row["homolog_projection_json"] = json.dumps(projection, sort_keys=True)
                return json.dumps(rows, sort_keys=True)
            table["relevance_evidence_json"] = table.relevance_evidence_json.map(scientific_evidence)
            projection_fingerprints.append(fingerprints)
        for fingerprints in projection_fingerprints:
            self.assertEqual({"G1", "G2"}, set(fingerprints))
            self.assertTrue(all(len(values) == 1 for values in fingerprints.values()))
        for gene in ("G1", "G2"):
            self.assertNotEqual(projection_fingerprints[0][gene], projection_fingerprints[1][gene])
        pd.testing.assert_frame_equal(left, right)
        self.assertEqual({"G1", "G2"}, set(right.candidate_gene_id))
        self.assertFalse(right.route_a_supported.any())
        self.assertFalse(right.route_b_supported.any())
        self.assertEqual(before, {p: p.read_bytes() for p in before})

    def test_nonempty_source_route_requires_explicit_review_and_empty_route_is_automatic(self):
        prepared = self.prepare()
        self.review(prepared)  # Produce the fixed B response once.
        source = self.root / "source.json"
        run.write(source, {"id": "12345", "pmid": "12345", "title": "Synthetic source",
                           "abstractText": "In Arabidopsis thaliana, loss of Q6UD73 caused fasciation."})
        reviewed = run.review(Path(prepared["manifest"]), prepared["sha256"], self.root / "b_response.json",
                              self.root / "mixed", source_cache=source)
        self.assertEqual(["a"], reviewed["independent_review_routes"])
        self.assertEqual({"b"}, set(reviewed["machine_reviews"]))
        with self.assertRaisesRegex(ValueError, "required for routes: a"):
            run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "missing-a")
        self.assertFalse((self.root / "missing-a").exists())
        prompt = run.read(Path(reviewed["sessions"]["a"]["prompt"]))
        response = self.root / "a-only-response.json"
        run.write(response, {"reviews": [{"hypothesis_id": row["hypothesis_id"], "decision": "not_evaluated",
                                         "reason": "TEST_ONLY fixed unresolved review"} for row in prompt["hypotheses"]]})
        specs = self.root / "a-only.json"
        run.write(specs, {"a": {"response": response.name, "reviewer_id": "TEST-review-a", "reviewer_model": "TEST/model"}})
        with MockBlastTools(raw_hits="Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"):
            result = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], specs, self.root / "mixed-finished")
        self.assertEqual("executed", result["reviews"]["a"]["review_execution_state"])
        self.assertEqual("not_executed_no_source", result["reviews"]["b"]["review_execution_state"])

    def test_legacy_stage_still_requires_explicit_reviews(self):
        reviewed = self.review(self.prepare())
        anchor = Path(reviewed["manifest"])
        stage = run.read(anchor)
        stage.pop("machine_reviews")
        stage.pop("independent_review_routes")
        anchor.write_text(json.dumps(stage), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "required for routes: a, b"):
            run.finish(anchor, run.digest(anchor), None, self.root / "legacy-missing")

    def test_changed_machine_response_is_rejected_before_finish(self):
        reviewed = self.review(self.prepare())
        response = self.root / "reviewed" / reviewed["machine_reviews"]["b"]["response"]
        response.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "stage artifact mismatch"):
            run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "tampered")
        self.assertFalse((self.root / "tampered").exists())

    def test_packets_through_projection_and_reviews_to_wet(self):
        prepared = self.prepare()
        reviewed = self.review(prepared)
        raw = "".join(f"Q6UD73\t{gene}\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n" for gene in ("G1", "G2"))
        with MockBlastTools(raw_hits=raw):
            result = self.finish(reviewed)
        self.assertEqual("first_stage_frozen", result["state"])
        projected = pd.read_csv(result["handoff"]["a"]["review_candidates"], sep="\t")
        self.assertEqual({"G1", "G2"}, set(projected["candidate_gene_id"]))
        b_review = pd.read_csv(result["handoff"]["b"]["review_candidates"], sep="\t")
        self.assertTrue(b_review.empty)  # #901: missing source confirmation stays in the B ledger.
        wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t")
        self.assertEqual({"G1", "G2"}, set(wet["candidate_gene_id"]))
        self.assertEqual({"route_review_candidate"}, set(wet["membership_basis"]))
        self.assertFalse(wet["route_a_supported"].any())
        self.assertFalse(wet["route_b_supported"].any())
        report_path = Path(result["report"])
        report = report_path.read_text(encoding="utf-8")
        self.assertIn("Canonical Wet: present", report)
        self.assertIn("Direct C execution: `not_evaluated_no_input`", report)
        self.assertIn("Candidate rows: 2", report)
        self.assertIn("no causal-accuracy claim", report)
        self.assertEqual(run.digest(report_path), run.read(Path(result["manifest"]))["files"]["report.md"])
        run.load_stage(Path(prepared["manifest"]), prepared["sha256"], "awaiting_b_response")
        run.load_stage(Path(reviewed["manifest"]), reviewed["sha256"], "awaiting_independent_reviews")

    @unittest.skipUnless(sys.platform == "linux", "ordinary extraction requires Linux/WSL and pysam")
    def test_ordinary_prepare_preserves_distinct_fasta_and_roster_ids(self):
        fasta = self.root / "colon.fa"
        fasta.write_text(">chr1\nATGAAATAA" + "N" * 10 + "ATGAAATAA\n", encoding="utf-8")
        gff = self.root / "colon.gff3"
        lines = ["##gff-version 3"]
        for gene, start, end in (("gene:A", 1, 9), ("gene_A", 20, 28)):
            lines.extend((
                f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
                f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gene}.t1",
            ))
        gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
        config = run.read(self.root / "config.json")
        for key in ("packet_a", "packet_b", "annotations", "orthology", "literature"):
            config.pop(key)
        config.update({"case_id": "case900", "focal_species": "Pisum sativum",
                       "phenotype": "fasciation", "interval": "chr1:1-28",
                       "reference_fasta": "colon.fa", "annotation_gff3": "colon.gff3",
                       "reference_storage": "snapshot"})
        run.write(self.root / "colon-config.json", config)
        prepared = run.prepare(self.root / "colon-config.json", self.root / "colon-prepared")
        packet = Path(prepared["inputs"]["packet_b"])
        with (packet / "candidates.tsv").open(encoding="utf-8", newline="") as handle:
            roster = {row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
        with (packet / "candidate_proteins.fa").open(encoding="utf-8") as handle:
            proteins = {record.id for record in SeqIO.parse(handle, "fasta")}
        self.assertEqual({"gene:A", "gene_A"}, roster)
        self.assertEqual(roster, proteins)

    def test_packet_protein_roster_mismatch_is_rejected(self):
        (self.root / "b" / "candidate_proteins.fa").write_text(">G1\nMAAA\n>not-in-roster\nMAAA\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "packet protein roster must match"):
            self.prepare()

    @unittest.skipUnless(sys.platform == "linux" and shutil.which("blastp") and shutil.which("makeblastdb"),
                         "requires Linux/WSL, pysam, and local BLAST+")
    def _ordinary_prepared(self):
        """Ordinary FASTA/GFF/QTL inputs through the real prepare CLI (30 genes)."""
        from tests.test_blast_evidence import S

        codons = {
            "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
            "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
            "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
            "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
        }
        cds = "".join(codons[aa] for aa in S) + "TAA"
        gene_count = 30
        stride = len(cds) + 10
        genome = list("N" * (gene_count * stride))
        gff_lines = ["##gff-version 3"]
        for index in range(1, gene_count + 1):
            start = (index - 1) * stride + 1
            end = start + len(cds) - 1
            genome[start - 1:end] = cds
            gene = f"G{index}"
            gff_lines.extend([
                f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
                f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gene}.t1",
            ])
        fasta = self.root / "ordinary.fa"
        fasta.write_text(">1\n" + "".join(genome) + "\n", encoding="utf-8")
        gff = self.root / "ordinary.gff3"
        gff.write_text("\n".join(gff_lines) + "\n", encoding="utf-8")
        gene_ids = [f"G{index}" for index in range(1, gene_count + 1)]
        pd.DataFrame([{
            "gene_id": gene, "family_id": f"FAM{index}", "evolutionary_family_id": f"EV{index}",
            "orthogroup": f"OG{index}", "subfamily_id": f"SUB{index}",
            "orthology_reference_gene": "REF1", "domain_architecture": "kinase",
            "broad_function_bins": "signaling",
        } for index, gene in enumerate(gene_ids, 1)]).to_csv(self.root / "annotations.tsv", sep="\t", index=False)
        pd.DataFrame([{
            "gene_id": gene, "reference_gene_id": "REF1", "orthogroup": f"OG{index}",
            "relationship": "homolog", "provider": "TEST_ONLY",
        } for index, gene in enumerate(gene_ids, 1)]).to_csv(self.root / "orthology.tsv", sep="\t", index=False)
        literature = dict(fixtures.LITERATURE_ROWS[0])
        literature["interval_gene_ids"] = ";".join(gene_ids)
        pd.DataFrame([literature]).to_csv(self.root / "literature.tsv", sep="\t", index=False)
        ordinary_config = run.read(self.root / "config.json")
        ordinary_config.pop("packet_a")
        ordinary_config.pop("packet_b")
        ordinary_config.update({"case_id": "case900", "focal_species": "Pisum sativum",
                                "phenotype": "fasciation", "interval": f"chr1:1-{len(genome)}",
                                "reference_fasta": "ordinary.fa", "annotation_gff3": "ordinary.gff3",
                                "reference_storage": "snapshot"})
        pd.DataFrame([{"gene_id": "G1", "function_description": "meristem receptor",
                       "source_id": "TEST_ONLY:reference-function"}]).to_csv(
            self.root / "ordinary_reference.tsv", sep="\t", index=False)
        ordinary_config["reference_corpus"] = "ordinary_reference.tsv"
        (self.root / "config.json").write_text(json.dumps(ordinary_config), encoding="utf-8")

        rejected = subprocess.run(
            [sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
             "--config", str(self.root / "config.json"), "--out-dir", str(self.root / "ordinary_unmapped")],
            capture_output=True, text=True,
        )
        self.assertNotEqual(0, rejected.returncode, "seqid mismatch must fail without an explicit mapping")
        ordinary_config["seqid_map"] = "chr1=1"
        (self.root / "config.json").write_text(json.dumps(ordinary_config), encoding="utf-8")
        prepared_result = subprocess.run(
            [sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
             "--config", str(self.root / "config.json"), "--out-dir", str(self.root / "ordinary_prepared")],
            capture_output=True, text=True,
        )
        self.assertEqual(0, prepared_result.returncode, prepared_result.stderr)
        prepared = json.loads(prepared_result.stdout)
        self.assertEqual("case900", prepared["ordinary_case"]["case_id"])
        self.assertEqual(gene_count, json.loads(Path(prepared["inputs"]["extraction_manifest"]).read_text(encoding="utf-8"))["candidate_count"])
        self.assertIn("meristem receptor", Path(prepared["bundle"]["prompts"]).read_text(encoding="utf-8"))
        return prepared, gene_count

    def test_issue758_ordinary_reference_inputs_reach_adopted_wet_without_handwritten_candidate_tsv(self):
        """Ordinary FASTA/GFF/QTL inputs preserve the fixed v2 independent-review replay."""
        from scripts.agentic_route_b import merge as merge_route_b
        from scripts.agentic_route_b import prepare as prepare_route_b

        prepared, gene_count = self._ordinary_prepared()

        shard_root = self.root / "ordinary_b_shards"
        shard_manifest = prepare_route_b(
            prepared["bundle"]["prompts"], str(shard_root), max_prompts=8,
        )
        self.assertGreaterEqual(len(shard_manifest["shards"]), 3)
        claim = [{
            "claim_type": "mutant_phenotype_bridge", "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
            "reference_gene": "REF1", "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup", "source_ids": ["PMID:123"],
            "alternative_hypothesis": "related copy", "uncertainty_note": "transfer unresolved",
        }]
        initial_attempts = []
        retry_attempts = []
        for index, shard in enumerate(shard_manifest["shards"]):
            shard_rows = [
                json.loads(line) for line in (shard_root / shard["file"]).read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            rows_path = self.root / f"ordinary_b_rows_{index}.json"
            shard_prompts = self.root / f"ordinary_b_prompts_{index}.jsonl"
            shard_prompts.write_text(
                "".join(json.dumps(row) + "\n" for row in shard_rows), encoding="utf-8"
            )
            fixtures.PacketToCanonicalTests()._responses(shard_prompts, claim, rows_path)
            shard_candidates = run.read(rows_path)["candidates"]
            if index == 1:
                failed_attempt = self.root / f"ordinary_b_attempt_{index}_failed.json"
                run.write(failed_attempt, {
                    "case_id": "case900", "executor": f"TEST-b-shard-{index}", "model": "TEST/model",
                    "response_id": f"shard-{index}-failed",
                    "responses": [
                        {"prompt_id": row["prompt_id"], "status": "failed", "reason": "fixture shard timeout"}
                        for row in shard_candidates
                    ],
                })
                retry_attempt = self.root / f"ordinary_b_attempt_{index}_retry.json"
                run.write(retry_attempt, {
                    "case_id": "case900", "executor": f"TEST-b-shard-{index}-retry", "model": "TEST/model",
                    "response_id": f"shard-{index}-retry", "responses": shard_candidates,
                })
                initial_attempts.append(failed_attempt)
                retry_attempts.append(retry_attempt)
            else:
                attempt = self.root / f"ordinary_b_attempt_{index}.json"
                run.write(attempt, {
                    "case_id": "case900", "executor": f"TEST-b-shard-{index}", "model": "TEST/model",
                    "response_id": f"shard-{index}", "responses": shard_candidates,
                })
                initial_attempts.append(attempt)
        first_merged = self.root / "ordinary_b_first_merged.json"
        first_execution = merge_route_b(
            prepared["bundle"]["prompts"], [str(path) for path in initial_attempts], str(first_merged),
            executor="TEST-b-parent", model="TEST/model",
        )
        self.assertEqual("partial", first_execution["completion_status"])
        self.assertGreater(first_execution["not_evaluated_count"], 0)
        merged = self.root / "ordinary_b_merged.json"
        execution = merge_route_b(
            prepared["bundle"]["prompts"],
            [str(first_merged), *[str(path) for path in retry_attempts]], str(merged),
            executor="TEST-b-parent", model="TEST/model",
        )
        self.assertEqual("complete", execution["completion_status"])
        self.assertGreaterEqual(len(execution["retry_history"]), 1)
        execution_manifest = Path(merged).with_name(Path(merged).stem + "_manifest.json")
        review_result = subprocess.run(
            [sys.executable, str(HERE / "practical_packet_run.py"), "review",
             "--prepared", prepared["manifest"], "--sha256", prepared["sha256"],
             "--responses", str(merged), "--b-execution-manifest", str(execution_manifest),
             "--out-dir", str(self.root / "ordinary_reviewed")],
            capture_output=True, text=True,
        )
        self.assertEqual(0, review_result.returncode, review_result.stderr)
        reviewed = json.loads(review_result.stdout)
        self.assertEqual("complete", reviewed["route_b_execution"]["completion_status"])

        c_status = self.root / "ordinary_c_status.json"
        run.write(c_status, {"case_id": "case900", "state": "provider_failure", "reason": "fixture C provider unavailable"})
        # The ordinary CLI consumes sealed empty-source records without starting
        # reviewers; later Wet boundary comparison still has independent answers.
        self.assertEqual([], reviewed["independent_review_routes"])
        # Preserve this fixed-answer independent-review replay as an explicit v2
        # generation. The new v3 original-assessment CLI has its own scope fixture.
        finish_v2 = ("import runpy,sys; from pathlib import Path; "
                     "sys.path.insert(0,str(Path(sys.argv[1]).resolve().parent)); "
                     "import candidate_qtl.boundary_tier_review as policy; "
                     "policy.DEFAULT_SELECTION_POLICY=policy.EVIDENCE_SCOPED_POLICY; "
                     "runpy.run_path(sys.argv.pop(1),run_name='__main__')")
        finished = subprocess.run([
            sys.executable, "-c", finish_v2, str(HERE / "practical_packet_run.py"), "finish",
            "--reviewed", reviewed["manifest"], "--sha256", reviewed["sha256"],
            "--route-c-status", str(c_status), "--max-adaptive-rounds", "2",
            "--out-dir", str(self.root / "ordinary_finished")], capture_output=True, text=True)
        self.assertEqual(0, finished.returncode, finished.stderr)
        result = json.loads(finished.stdout)
        for route in ("a", "b"):
            self.assertEqual("not_executed_no_source", result["reviews"][route]["review_execution_state"])
        self.assertEqual("complete", result["route_b_execution"]["completion_status"])
        self.assertEqual("provider_failure", result["route_c_execution"]["state"])
        first_report = Path(result["report"]).read_text(encoding="utf-8")
        self.assertIn("Route B execution manifest", first_report)
        self.assertIn("Route C execution status", first_report)
        wet_path = Path(result["handoff"]["wet"]["wet_candidates"])
        wet = pd.read_csv(wet_path, sep="\t", dtype=str).fillna("")
        self.assertEqual({f"G{index}" for index in range(1, gene_count + 1)}, set(wet["candidate_gene_id"]))
        self.assertTrue(result["handoff"]["wet"]["adaptive_synthesis_request"])

        first_stage = run.read(Path(result["manifest"]))
        self.assertTrue(any(name.endswith("route_b_execution_manifest.json") for name in first_stage["files"]))
        self.assertTrue(any(name.endswith("route_c_status.json") for name in first_stage["files"]))
        source_wet_bytes = wet_path.read_bytes()
        scope = run.read(Path(result["handoff"]["wet"]["boundary_scope_manifest"]))
        self.assertEqual("evidence_scoped_v2", scope["policy_version"])
        self.assertEqual(set(wet["candidate_gene_id"]), set(scope["eligible_gene_ids"]))
        proposal, independent = self.agentic_boundary_round(
            result["handoff"]["wet"]["adaptive_synthesis_request"],
            self.root / "ordinary_boundary", adopt=8,
        )
        handoff_root = Path(result["manifest"]).parent / "handoff"
        adopted_dir = self.root / "ordinary_adopted"
        command = [
            sys.executable, str(HERE.parents[1] / "scripts" / "claim_review_session.py"),
            "finish-wet-review", "--handoff", str(handoff_root),
            "--completion-sha256", run.digest(handoff_root / "completion.json"),
            "--out", str(adopted_dir), "--proposals", str(proposal),
            "--independent-review", str(independent),
            "--executor-id", "TEST-boundary-executor", "--executor-model", "TEST/model",
            "--reviewer-id", "TEST-boundary-reviewer", "--reviewer-model", "TEST/model",
        ]
        b_manifest = Path(result["route_b_execution"]["path"])
        b_manifest_bytes = b_manifest.read_bytes()
        b_manifest.write_text("{}\n", encoding="utf-8")
        tampered_result = subprocess.run(command, cwd=HERE.parents[1], capture_output=True, text=True)
        self.assertNotEqual(0, tampered_result.returncode)
        b_manifest.write_bytes(b_manifest_bytes)
        adopted_result = subprocess.run(command, cwd=HERE.parents[1], capture_output=True, text=True)
        self.assertEqual(0, adopted_result.returncode, adopted_result.stderr)
        adopted = json.loads(adopted_result.stdout)
        self.assertEqual("completed_evidence_scoped_selection", adopted["state"])
        self.assertEqual("evidence_scoped_v2", adopted["boundary_selection"]["policy_version"])
        self.assertEqual(8, adopted["adopted_candidate_count"])
        self.assertIsNone(adopted["next_boundary_scope_manifest"])
        self.assertIsNone(adopted["wet"]["adaptive_synthesis_request"])
        adopted_path = Path(adopted["adopted_wet_candidates"])
        self.assertEqual(adopted["adopted_wet_sha256"], run.digest(adopted_path))
        adopted_ids = set(pd.read_csv(adopted_path, sep="\t", dtype=str)["candidate_gene_id"])
        self.assertEqual(8, len(adopted_ids))
        self.assertTrue(adopted_ids <= set(wet["candidate_gene_id"]))
        self.assertEqual(source_wet_bytes, wet_path.read_bytes())
        deferred = pd.read_csv(adopted["boundary_selection"]["deferral_ledger"]["path"], sep="\t", dtype=str).fillna("")
        self.assertEqual(set(wet["candidate_gene_id"]) - adopted_ids, set(deferred["candidate_gene_id"]))
        final_report = Path(adopted["report"]["path"]).read_text(encoding="utf-8")
        self.assertIn("route_b_execution", final_report)
        self.assertIn("provider_failure", final_report)
        self.assertIn("Wet SHA-256", final_report)
        self.assertIn(adopted["adopted_wet_sha256"], final_report)
        self.assertEqual(adopted["report"]["sha256"], run.digest(Path(adopted["report"]["path"])))
        self.assertEqual(execution["response_sha256"], result["route_b_execution"]["response_sha256"])

    def shared_reference_config(self):
        import pysam
        fasta = self.root / "shared.fa"
        fasta.write_text(">1\nATGAAATAA" + "N" * 10 + "ATGAAATAA" + "N" * 10 + "ATGAAATAA\n", encoding="utf-8")
        gff = self.root / "shared.gff3"
        lines = ["##gff-version 3"]
        for gene, start, end in (("G1", 1, 9), ("G2", 20, 28), ("outside", 39, 47)):
            lines.extend((
                f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
                f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gene}.t1",
            ))
        gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
        pysam.faidx(str(fasta))
        config = run.read(self.root / "config.json")
        config.pop("packet_a")
        config.pop("packet_b")
        config.update(case_id="case900", focal_species="Pisum sativum", phenotype="fasciation",
                      interval="chr1:5-28", seqid_map="chr1=1", reference_fasta="shared.fa",
                      annotation_gff3="shared.gff3", reference_storage="shared_readonly",
                      reference_version="synthetic-v2", reference_corpus="b/reference_corpus.tsv")
        path = self.root / "shared-config.json"
        run.write(path, config)
        return path, config, fasta, gff

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_ac_first_ordinary_finish_records_optional_route_b_not_run(self):
        """A/C-first (2026-11-08): optional B deliberately not run reaches Wet honestly."""
        from scripts.agentic_route_b import merge as merge_route_b

        prepared, _ = self._ordinary_prepared()
        not_run = self.root / "b_not_run.json"
        execution = merge_route_b(prepared["bundle"]["prompts"], [], str(not_run),
                                  executor="not_run", model="not_run")
        self.assertEqual(0, execution["completed_count"])
        manifest = Path(not_run).with_name(Path(not_run).stem + "_manifest.json")
        anchor = Path(prepared["manifest"])
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            run.review(anchor, prepared["sha256"], not_run, self.root / "both",
                       b_execution_manifest=manifest, route_b_not_run="A/C-first")
        with self.assertRaisesRegex(ValueError, "non-empty reason"):
            run.review(anchor, prepared["sha256"], not_run, self.root / "blank", route_b_not_run=" ")
        rows = run.read(not_run)["candidates"]
        mixed = self.root / "mixed_not_run.json"
        fixtures.PacketToCanonicalTests()._responses(Path(prepared["bundle"]["prompts"]),
            [{"claim_type": "mechanistic_review_seed"}, {}], mixed)
        with self.assertRaisesRegex(ValueError, "all be not evaluated"):
            run._route_b_not_run_context("optional B omitted", mixed)
        answered = self.root / "b_answered.json"
        fixtures.PacketToCanonicalTests()._responses(Path(prepared["bundle"]["prompts"]), [], answered)
        with self.assertRaisesRegex(ValueError, "all be not evaluated"):
            run.review(anchor, prepared["sha256"], answered, self.root / "relabelled",
                       route_b_not_run="A/C-first")
        self.assertTrue(rows)

        result = subprocess.run(
            [sys.executable, str(HERE / "practical_packet_run.py"), "review",
             "--prepared", prepared["manifest"], "--sha256", prepared["sha256"],
             "--responses", str(not_run), "--route-b-not-run", "A/C-first: optional B not requested",
             "--out-dir", str(self.root / "ac_reviewed")], capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
        reviewed = json.loads(result.stdout)
        self.assertEqual("not_run_optional", reviewed["route_b_execution"]["state"])
        self.assertEqual(len(rows), reviewed["route_b_execution"]["not_evaluated_count"])
        c_status = self.root / "ac_c_status.json"
        run.write(c_status, {"case_id": "case900", "state": "not_evaluated", "reason": "fixture C not run"})
        finished = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "ac_finished",
                              route_c_status=c_status, max_adaptive_rounds=2)
        self.assertEqual("not_run_optional", finished["route_b_execution"]["state"])
        self.assertEqual(0, finished["route_b_execution"]["completed_count"])
        self.assertNotEqual("route_b_declared_but_unverified", finished["handoff"]["wet"]["input_coverage_state"])
        report = Path(finished["report"]).read_text(encoding="utf-8")
        self.assertIn("Route B execution: `not_run_optional`", report)
        self.assertIn("A/C-first: optional B not requested", report)
        wet = pd.read_csv(finished["handoff"]["wet"]["wet_candidates"], sep="\t", dtype=str).fillna("")
        self.assertNotIn("interval_review_coverage", set(wet.get("membership_basis", pd.Series(dtype=str))))
        audit = run.audit(Path(finished["manifest"]), finished["sha256"])
        self.assertEqual("not_run_optional", audit["route_b_execution"]["state"])

        from tests._issue664_helpers import manifest_for
        c_members = self.root / "empty_c.tsv"
        pd.DataFrame(columns=["candidate_gene_id", "support_status"]).to_csv(c_members, sep="\t", index=False)
        c_manifest = manifest_for(c_members, route="c", case_id="case900")
        query_result = self.root / "native_query.json"
        run.write(query_result, {"query": "fixture phenotype reported gene", "results": []})
        receipt = self.root / "native_receipt.json"
        run.write(receipt, {"contract": "route-c-native-research-receipt-1.0.0", "case_id": "case900",
            "worker_id": "fixture-c", "model": "fixture", "interval_blind": True, "state": "completed",
            "finding_count": 0, "no_findings_rationale": "No relevant direct report in the bounded search",
            "queries": [{"query": "fixture phenotype reported gene", "status": "no_hit",
                         "result_path": str(query_result), "result_sha256": run.digest(query_result)}], "read_sources": []})
        empty_status = self.root / "native_empty_status.json"
        run.write(empty_status, {"case_id": "case900", "state": "completed_empty",
            "research_receipt_path": str(receipt), "research_receipt_sha256": run.digest(receipt)})
        empty_finished = run.finish(Path(reviewed["manifest"]), reviewed["sha256"], None, self.root / "empty_c_finished",
            route_c_members=c_members, route_c_manifest=c_manifest, route_c_status=empty_status)
        sealed = empty_finished["handoff"]["execution_context"]["route_c_execution"]
        self.assertEqual("completed_empty", sealed["state"])
        sealed_receipt = run.read(Path(sealed["research_receipt_path"]))
        self.assertTrue(Path(sealed_receipt["queries"][0]["result_path"]).is_relative_to(self.root / "empty_c_finished"))
        self.assertEqual(run.digest(receipt), sealed["original_research_receipt_sha256"])

    def test_shared_reference_cli_matches_snapshot_and_resumes_without_source(self):
        path, config, fasta, gff = self.shared_reference_config()
        index = Path(str(fasta) + ".fai")
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (fasta, gff, index)}
        command = [sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
                   "--config", str(path), "--out-dir", str(self.root / "prepared")]
        executed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(0, executed.returncode, executed.stderr)
        prepared = json.loads(executed.stdout)
        asset = run.read(Path(prepared["inputs"]["reference_asset_manifest"]))
        self.assertEqual("synthetic-v2", asset["version"])
        self.assertEqual("chr1", asset["gff_seqid"])
        self.assertEqual("1", asset["fasta_seqid"])
        for key, source in (("fasta", fasta), ("gff3", gff)):
            self.assertEqual(run.digest(source), asset[key]["sha256"])
            self.assertEqual("strong", asset[key]["strength"])
        self.assertFalse(any(name.startswith(("reference_fasta", "annotation_gff3"))
                             for name in prepared["files"]))
        config["reference_storage"] = "snapshot"
        path.write_text(json.dumps(config), encoding="utf-8")
        copied = run.prepare(path, self.root / "snapshot")
        b = Path(prepared["inputs"]["packet_b"])
        for name in ("candidates.tsv", "transcripts.tsv", "id_map.tsv", "candidate_genomic.fa",
                     "candidate_cds.fa", "candidate_transcript_cds.fa", "candidate_proteins.fa"):
            self.assertEqual((b / name).read_bytes(), (Path(copied["inputs"]["packet_b"]) / name).read_bytes())
        self.assertEqual({"G1", "G2"}, set(pd.read_csv(b / "candidates.tsv", sep="\t")["gene_id"]))
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before})
        for p in before:
            p.unlink()
        from scripts.agentic_route_b import merge
        responses = self.root / "raw-b.json"
        fixtures.PacketToCanonicalTests()._responses(prepared["bundle"]["prompts"], [{
            "claim_type": "mechanistic_review_seed", "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
            "reference_gene": "REF1", "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup", "source_ids": [],
            "alternative_hypothesis": "related copy", "uncertainty_note": "transfer unresolved",
        }], responses)
        payload = run.read(responses)
        payload.update(case_id="case900", executor="TEST-b", model="TEST/model")
        payload["responses"] = payload.pop("candidates")
        responses.write_text(json.dumps(payload), encoding="utf-8")
        merged = self.root / "merged-b.json"
        merge(prepared["bundle"]["prompts"], [str(responses)], str(merged),
                          executor="TEST-b", model="TEST/model")
        reviewed = run.review(Path(prepared["manifest"]), prepared["sha256"], merged,
                              self.root / "reviewed", merged.with_name("merged-b_manifest.json"))
        finished = self.finish(reviewed, cli=True)
        self.assertTrue(Path(finished["handoff"]["wet"]["wet_candidates"]).is_file())
        (b / "candidate_genomic.fa").write_text(">changed\nAAA\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            run.load_stage(Path(prepared["manifest"]), prepared["sha256"], "awaiting_b_response")

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_ordinary_default_is_readonly_and_matches_snapshot(self):
        path, config, fasta, gff = self.shared_reference_config()
        config.pop("reference_storage")
        path.write_text(json.dumps(config), encoding="utf-8")
        sources = (fasta, gff, Path(str(fasta) + ".fai"))
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in sources}
        prepared = run.prepare(path, self.root / "default-readonly")
        self.assertEqual("shared_readonly", prepared["ordinary_case"]["reference_storage"])
        self.assertFalse(any(name.startswith(("inputs/reference_fasta", "inputs/annotation_gff3"))
                             for name in prepared["files"]))
        self.assertNotIn("managed_asset", run.read(Path(prepared["inputs"]["reference_asset_manifest"])))
        config["reference_storage"] = "snapshot"
        path.write_text(json.dumps(config), encoding="utf-8")
        snapshot = run.prepare(path, self.root / "explicit-snapshot")
        for name in ("candidates.tsv", "transcripts.tsv", "id_map.tsv", "candidate_genomic.fa",
                     "candidate_cds.fa", "candidate_transcript_cds.fa", "candidate_proteins.fa"):
            self.assertEqual((Path(prepared["inputs"]["packet_b"]) / name).read_bytes(),
                             (Path(snapshot["inputs"]["packet_b"]) / name).read_bytes())
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in sources})

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_default_readonly_requires_version_and_does_not_create_missing_index(self):
        path, config, fasta, _ = self.shared_reference_config()
        config.pop("reference_storage")
        config.pop("reference_version")
        path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "reference_version"):
            run.prepare(path, self.root / "default-no-version")
        self.assertFalse((self.root / "default-no-version").exists())
        config["reference_version"] = "synthetic-v2"
        path.write_text(json.dumps(config), encoding="utf-8")
        index = Path(str(fasta) + ".fai")
        index.unlink()
        with patch("pysam.faidx", side_effect=AssertionError("must not write source index")):
            with self.assertRaisesRegex(ValueError, "existing .fai"):
                run.prepare(path, self.root / "default-no-index")
        self.assertFalse(index.exists())

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_ordinary_default_rejects_stale_index_without_sealing_stage(self):
        path, config, fasta, gff = self.shared_reference_config()
        config.pop("reference_storage")
        path.write_text(json.dumps(config), encoding="utf-8")
        # Preserve bases and lengths but change the offset by lengthening a header.
        content = fasta.read_text(encoding="utf-8")
        fasta.write_text(content.replace(">1", ">1 longer description", 1), encoding="utf-8")
        sources = (fasta, gff, Path(str(fasta) + ".fai"))
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in sources}
        out = self.root / "stale-index"
        executed = subprocess.run([sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
                                   "--config", str(path), "--out-dir", str(out)], capture_output=True, text=True)
        self.assertNotEqual(0, executed.returncode)
        self.assertIn("FASTA/.fai layout mismatch", executed.stderr)
        self.assertFalse((out / "stage.json").exists())
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in sources})

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_shared_reference_missing_index_never_builds_source_index(self):
        path, config, fasta, gff = self.shared_reference_config()
        index = Path(str(fasta) + ".fai")
        index.unlink()
        before = {p: p.read_bytes() for p in (fasta, gff)}
        with patch("pysam.faidx", side_effect=AssertionError("must not index shared input")):
            with self.assertRaisesRegex(ValueError, "existing .fai"):
                run.prepare(path, self.root / "missing-index")
        self.assertFalse(index.exists())
        self.assertEqual(before, {p: p.read_bytes() for p in before})

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_shared_reference_requires_explicit_version_and_rejects_symlink(self):
        path, config, fasta, _ = self.shared_reference_config()
        config.pop("reference_version")
        path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "reference_version"):
            run.prepare(path, self.root / "no-version")
        self.assertFalse((self.root / "no-version").exists())
        link = self.root / "link.fa"
        link.symlink_to(fasta)
        config.update(reference_version="v2", reference_fasta=str(link))
        path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "non-symlink"):
            run.prepare(path, self.root / "linked-source")

    def managed_reference_config(self):
        path, config, fasta, gff = self.shared_reference_config()
        Path(str(fasta) + ".fai").unlink()
        config.update(reference_storage="managed_cache", reference_cache_dir="reference-cache")
        path.write_text(json.dumps(config), encoding="utf-8")
        return path, config, fasta, gff

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_managed_reference_cli_reuses_extraction_and_seals_independent_case_files(self):
        path, config, fasta, gff = self.managed_reference_config()
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (fasta, gff)}
        command = [sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
                   "--config", str(path), "--out-dir", str(self.root / "prepared")]
        executed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(0, executed.returncode, executed.stderr)
        first = json.loads(executed.stdout)
        lineage = run.read(Path(first["inputs"]["reference_asset_manifest"]))
        self.assertFalse(lineage["managed_asset"]["cache_hit"])
        self.assertFalse(lineage["interval_cache"]["hit"])
        self.assertFalse(Path(str(fasta) + ".fai").exists())
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before})
        config.update(case_id="case901", phenotype="different observed phenotype")
        path.write_text(json.dumps(config), encoding="utf-8")
        with patch("pysam.faidx", side_effect=AssertionError("shared index must be reused")), \
             patch("candidate_qtl.reference_extract.iter_gff_rows",
                   side_effect=AssertionError("same interval must not be parsed again")):
            second = run.prepare(path, self.root / "reused")
        reused = run.read(Path(second["inputs"]["reference_asset_manifest"]))
        self.assertTrue(reused["managed_asset"]["cache_hit"])
        self.assertTrue(reused["interval_cache"]["hit"])
        self.assertEqual(lineage["interval_cache"]["key"], reused["interval_cache"]["key"])
        first_b, second_b = Path(first["inputs"]["packet_b"]), Path(second["inputs"]["packet_b"])
        self.assertEqual("case901", run.read(second_b / "case.json")["case_id"])
        self.assertEqual("different observed phenotype", run.read(second_b / "phenotype.json")["phenotype"])
        names = ("candidates.tsv", "transcripts.tsv", "id_map.tsv", "candidate_genomic.fa",
                 "candidate_cds.fa", "candidate_transcript_cds.fa", "candidate_proteins.fa")
        for name in names:
            self.assertEqual((first_b / name).read_bytes(), (second_b / name).read_bytes())
            self.assertNotEqual((first_b / name).stat().st_ino, (second_b / name).stat().st_ino)
        extraction = run.read(second_b / "extract_manifest.json")
        self.assertTrue(all(Path(value).parent == second_b for value in extraction["outputs"].values()))
        config.pop("reference_cache_dir")
        config.update(reference_storage="snapshot", case_id="case900", phenotype="fasciation")
        path.write_text(json.dumps(config), encoding="utf-8")
        snapshot = run.prepare(path, self.root / "snapshot")
        for name in names:
            self.assertEqual((first_b / name).read_bytes(),
                             (Path(snapshot["inputs"]["packet_b"]) / name).read_bytes())
        # Frozen stages consume independent case artifacts, even when every
        # external reference and its reusable cache are subsequently unavailable.
        for source in before:
            source.unlink()
        shutil.rmtree(self.root / "reference-cache")
        from scripts.agentic_route_b import merge
        responses = self.root / "managed-b.json"
        fixtures.PacketToCanonicalTests()._responses(first["bundle"]["prompts"], [{
            "claim_type": "mechanistic_review_seed", "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
            "reference_gene": "REF1", "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "orthogroup", "source_ids": [],
            "alternative_hypothesis": "related copy", "uncertainty_note": "transfer unresolved",
        }], responses)
        payload = run.read(responses)
        payload.update(case_id="case900", executor="TEST-b", model="TEST/model")
        payload["responses"] = payload.pop("candidates")
        responses.write_text(json.dumps(payload), encoding="utf-8")
        merged = self.root / "managed-b-merged.json"
        merge(first["bundle"]["prompts"], [str(responses)], str(merged), executor="TEST-b", model="TEST/model")
        reviewed = run.review(Path(first["manifest"]), first["sha256"], merged, self.root / "reviewed",
                              merged.with_name("managed-b-merged_manifest.json"))
        with MockBlastTools(raw_hits=""):
            finished = self.finish(reviewed)
        self.assertTrue(Path(finished["handoff"]["wet"]["wet_candidates"]).is_file())
        (second_b / "candidate_genomic.fa").write_text(">changed\nAAA\n", encoding="utf-8")
        run.load_stage(Path(first["manifest"]), first["sha256"], "awaiting_b_response")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            run.load_stage(Path(second["manifest"]), second["sha256"], "awaiting_b_response")

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_managed_config_rejects_missing_cache_or_overlapping_output(self):
        path, config, _, _ = self.managed_reference_config()
        config.pop("reference_cache_dir")
        path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "reference_cache_dir"):
            run.prepare(path, self.root / "missing-cache")
        self.assertFalse((self.root / "missing-cache").exists())
        config["reference_cache_dir"] = "overlap"
        path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "separate"):
            run.prepare(path, self.root / "overlap" / "case")
        self.assertFalse((self.root / "overlap").exists())

    @unittest.skipUnless(sys.platform == "linux", "requires pysam/WSL")
    def test_shared_protein_db_keeps_external_search_target_and_rejects_changes(self):
        path, config, _, _ = self.managed_reference_config()
        db = self.root / "proteins" / "multi-species"
        db.parent.mkdir()
        for suffix in ("pin", "phr", "psq"):
            Path(str(db) + "." + suffix).write_bytes(b"synthetic " + suffix.encode())
        config.update(protein_db=str(db), protein_db_storage="shared_readonly", protein_db_version="fixture-v1")
        config.pop("annotations")
        config.pop("orthology")
        path.write_text(json.dumps(config), encoding="utf-8")
        calls = []
        def unavailable(query, database, *args):
            calls.append((query, database))
            self.assertEqual(db, database)
            self.assertEqual({"G1", "G2"}, {record.id for record in SeqIO.parse(query, "fasta")})
            return [], {"available": False, "evaluation_state": "not_evaluated_provider_unavailable",
                        "reason": "synthetic unavailable backend"}
        with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=unavailable):
            prepared = run.prepare(path, self.root / "shared-db")
        self.assertEqual(1, len(calls))
        self.assertEqual(str(db), prepared["inputs"]["protein_db"])
        self.assertFalse((self.root / "shared-db/inputs/protein_db").exists())
        saved = run.read(Path(prepared["inputs"]["protein_db_asset_manifest"]))
        self.assertEqual("fixture-v1", saved["version"])
        self.assertEqual({db.name + "." + suffix for suffix in ("pin", "phr", "psq")}, set(saved["file_roster"]))
        b0_manifest = run.read(Path(prepared["b0"]["manifest"]))
        self.assertEqual("not_evaluated_provider_unavailable", b0_manifest["providers"]["homology_search"]["evaluation_state"])
        def changed(*args):
            result = unavailable(*args)
            Path(str(db) + ".psq").write_bytes(b"changed database")
            return result
        with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "changed"):
                run.prepare(path, self.root / "changed-db")
        self.assertFalse((self.root / "changed-db/stage.json").exists())

    def test_ordinary_config_rejects_missing_or_unknown_inputs(self):
        config_path = self.root / "config.json"
        config = run.read(config_path)
        config.pop("packet_a")
        config.pop("packet_b")
        config.update({"case_id": "case900", "focal_species": "Pisum sativum",
                       "phenotype": "fasciation", "interval": "chr1:1-9999",
                       "reference_fasta": "ordinary.fa"})
        config_path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "missing or unknown"):
            run.prepare(config_path, self.root / "missing-input-output")
        config["annotation_gff3"] = "ordinary.gff3"
        config["unexpected"] = "ignored.tsv"
        config_path.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "missing or unknown"):
            run.prepare(config_path, self.root / "unknown-input-output")

    @unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "local BLAST+ required")
    def test_real_blast_packet_handoff(self):
        result = self.finish(self.review(self.prepare()), cli=True)
        projected = pd.read_csv(result["handoff"]["a"]["review_candidates"], sep="\t")
        self.assertEqual({"G1", "G2"}, set(projected["candidate_gene_id"]))

    def test_missing_blast_does_not_create_a_negative(self):
        reviewed = self.review(self.prepare())
        with patch("candidate_qtl.blast_evidence.shutil.which", return_value=None):
            result = self.finish(reviewed)
        wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t")
        self.assertTrue(wet.empty)
        self.assertEqual("review_required_no_supported_candidates", result["handoff"]["wet"]["handoff_state"])
        inventory = pd.read_csv(self.root / "finished/handoff/wet/wet_relevance_inventory.tsv", sep="\t")
        self.assertEqual({"G1", "G2"}, set(inventory.candidate_gene_id))
        evidence = pd.read_csv(self.root / "finished/handoff/blast/blast_evidence.tsv", sep="\t")
        self.assertTrue(evidence["provider_status"].str.startswith("not_evaluated").all())

    def test_report_distinguishes_unavailable_wet_from_zero_candidates(self):
        out = self.root / "report-only"
        out.mkdir()
        path = Path(run.finish_report(out, "case900", "wet_unavailable", {"wet": {}}))
        report = path.read_text(encoding="utf-8")
        self.assertIn("unavailable (not evaluated)", report)
        self.assertIn("no candidate count is inferred", report)
        self.assertNotIn("Candidate rows: 0", report)

    def test_report_uses_requested_boundary_state_for_next_action(self):
        out = self.root / "boundary-report"
        out.mkdir()
        path = Path(run.finish_report(out, "case`900", "first_stage_frozen", {
            "wet": {"handoff_state": "boundary_review_requested",
                    "second_stage_state": "requested",
                    "boundary_scope_manifest": {"path": "scope`review.json"}},
        }))
        report = path.read_text(encoding="utf-8")
        self.assertIn("Execution state: `first_stage_frozen`", report)
        self.assertIn("Wet handoff state: `boundary_review_requested`", report)
        self.assertIn("Complete the requested boundary review", report)
        self.assertIn("scope\\`review.json", report)
        self.assertIn("case\\`900", report)

    def test_report_explains_wet_membership_and_unreached_trace_without_negative_inference(self):
        out = self.root / "readable-report"
        out.mkdir()
        wet_path = out / "wet_candidates.tsv"
        wet_path.write_text(
            "candidate_gene_id\tmembership_basis\trelevance_tier\trelevance_basis\t"
            "mapping_tier\tcopy_state\tsupport_state\troute_provenance\troute_support\tsource_ids\t"
            "tier_reason\treview_reasons\n"
            "G1\troute_review_candidate\tT2_trait_or_mechanism_linked\texplicit_mechanism_bridge\t"
            "tier2_reference_phenotype\tparalog_ambiguous\treview\tA;B\treview_only\tPMID:123\t"
            "generic tier rationale\tclaim remains under review\n",
            encoding="utf-8")
        with wet_path.open("a", encoding="utf-8") as stream:
            stream.write("G3\tinterval_review_coverage\tT3_functionally_plausible\t\t\t\tcontext\t\t\t\tfill tier rationale\t\n")
            for index in range(4, 33):
                stream.write(f"G{index}\tinterval_review_coverage\tT4_identity_only\t\t\t\tcontext\t\t\t\tidentity tier rationale\t\n")
        trace_path = out / "wet_candidate_trace.tsv"
        trace_path.write_text(
            "candidate_gene_id\tselection_disposition\trelevance_tier\tsupport_state\tbiological_exclusion\n"
            "G1\tretained\tT2_trait_or_mechanism_linked\treview\tnot_established\n"
            "G2\tnot_reached_lower_tier\tT4_identity_only\tcontext\tnot_established\n"
            "G3\tretained_pending_boundary_review\tT3_functionally_plausible\tcontext\tnot_established\n",
            encoding="utf-8")
        membership_manifest_path = out / "wet_candidates_manifest.json"
        run.write(membership_manifest_path, {"candidate_trace": str(trace_path),
            "candidate_trace_sha256": run.digest(trace_path), "relevance_inventory": "",
            "relevance_selection_applied": True})
        handoff_manifest_path = out / "wet_candidate_handoff_manifest.json"
        run.write(handoff_manifest_path, {"candidate_trace": "not-the-canonical-trace-location"})
        path = Path(run.finish_report(out, "case900", "first_stage_frozen", {
            "adopted_stage": "first_stage",
            "wet": {"handoff_state": "boundary_review_requested", "second_stage_state": "requested",
                    "wet_candidates": str(wet_path), "wet_candidates_sha256": run.digest(wet_path),
                    "manifest": str(handoff_manifest_path),
                    "wet_candidates_manifest": str(membership_manifest_path)}}))
        report = path.read_text(encoding="utf-8")
        self.assertIn("## Wet candidates", report)
        self.assertIn("G1", report)
        self.assertIn("pending boundary review", report)
        self.assertIn("First-stage lower-tier genes not reached: `G2`", report)
        self.assertIn("Showing 30 of 31 Wet genes", report)
        self.assertNotIn("| G32 |", report)
        self.assertIn("route_review_candidate", report)
        self.assertIn("fill tier rationale", report)
        g1_line = next(line for line in report.splitlines() if line.startswith("| G1 |"))
        self.assertIn("claim remains under review", g1_line)
        self.assertNotIn("generic tier rationale", g1_line)
        self.assertIn("T2_trait_or_mechanism_linked", report)
        self.assertIn("tier2_reference_phenotype", report)
        self.assertIn("paralog_ambiguous", report)
        self.assertIn("PMID:123", report)
        self.assertIn("not_reached_lower_tier=1", report)
        self.assertIn("retained_pending_boundary_review=1", report)
        self.assertIn("Neither state is a biological exclusion or a negative finding", report)
        self.assertIn("Candidate trace", report)
        self.assertIn(run.digest(trace_path), report)
        self.assertIn("Wet candidates manifest", report)
        self.assertIn("## Interval coverage", report)
        adopted_out = self.root / "adopted-report"
        adopted_out.mkdir()
        adopted_path = Path(run.finish_report(adopted_out, "case900", "completed_boundary_selection", {
            "adopted_stage": "boundary_tier_review",
            "wet": {"handoff_state": "completed_boundary_selection", "second_stage_state": "not_required",
                    "wet_candidates": str(wet_path), "wet_candidates_sha256": run.digest(wet_path),
                    "manifest": str(handoff_manifest_path),
                    "wet_candidates_manifest": str(membership_manifest_path)}}))
        adopted_report = adopted_path.read_text(encoding="utf-8")
        self.assertIn("adopted (boundary_tier_review)", adopted_report)
        g3_adopted_line = next(line for line in adopted_report.splitlines() if line.startswith("| G3 |"))
        self.assertIn("adopted (boundary_tier_review)", g3_adopted_line)
        self.assertNotIn("pending boundary review", g3_adopted_line)
        self.assertIn("Trace is the first-stage coverage ledger", adopted_report)
        self.assertIn("Final decision reasons: not provided", adopted_report)

    def test_supported_a_reaches_family_copy_while_b_without_source_stays_in_ledger(self):
        from candidate_qtl.agent_claim_review import AXES
        text = "In Arabidopsis thaliana, loss of Q6UD73 caused fasciation."
        run.write(self.root / "sources.json", {"id": "12345", "pmid": "12345",
            "title": "Q6UD73 meristem", "abstractText": text})
        config_path = self.root / "config.json"
        config = run.read(config_path)
        config["source_cache"] = "sources.json"
        config_path.write_text(json.dumps(config))
        reviewed = self.review(self.prepare())
        prompt = run.read(Path(reviewed["sessions"]["a"]["prompt"]))
        source = prompt["source_review_records"]["records"]["PMID:12345"]
        evidence = {"hypothesis_id": "H1", "citation": "PMID:12345", "source_record_id": "12345",
            "sections_sha256": source["sections_sha256"], "subject_gene_identity": "Q6UD73",
            "reported_gene_identity": "Q6UD73", "organism": "Arabidopsis thaliana",
            "reported_organism": "Arabidopsis thaliana", "organism_quote": "In Arabidopsis thaliana",
            "section": "abstractText", "quoted_support": text, "quote_start": 0, "quote_end": len(text),
            "experiment_perturbation": "loss of Q6UD73", "reported_phenotype": "fasciation",
            "predicate_relation": "causal_perturbation_supports_target_trait", "decision": "supported",
            "reason": "TEST_ONLY supported source fixture", "axes": dict.fromkeys(AXES, True)}
        raw = "Q6UD73\tG1\t100\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200\n"
        with MockBlastTools(raw_hits=raw):
            result = self.finish(reviewed, evidence)
        self.assertGreater(result["handoff"]["family"]["family_count"], 0)
        self.assertTrue((self.root / "finished/handoff/copy/family_copy_discrimination.tsv").is_file())
        wet = pd.read_csv(result["handoff"]["wet"]["wet_candidates"], sep="\t").set_index("candidate_gene_id")
        self.assertTrue(wet.loc["G1", "route_a_supported"])
        self.assertEqual({"G1"}, set(wet.index))
        b_review = pd.read_csv(result["handoff"]["b"]["review_candidates"], sep="\t")
        self.assertTrue(b_review.empty)
        self.assertFalse(wet.route_b_supported.any())

    def test_missing_b_response_roster_rejected(self):
        prepared = self.prepare()
        run.write(self.root / "empty-response.json", {"candidates": []})
        with self.assertRaisesRegex(ValueError, "roster|missing"):
            run.review(Path(prepared["manifest"]), prepared["sha256"],
                       self.root / "empty-response.json", self.root / "reviewed")

    def test_changed_preparation_refused_before_review(self):
        prepared = self.prepare()
        Path(prepared["bundle"]["prompts"]).write_text("changed")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            run.review(Path(prepared["manifest"]), prepared["sha256"], self.root / "absent", self.root / "reviewed")
        self.assertFalse((self.root / "reviewed").exists())

    def test_cross_case_packets_rejected(self):
        path = self.root / "a/phenotype.json"
        data = run.read(path)
        data["case_id"] = "different"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "packet mismatch"):
            self.prepare()

    def test_cli_prepare_and_hash_anchor(self):
        result = subprocess.run([sys.executable, str(HERE / "practical_packet_run.py"), "prepare",
            "--config", str(self.root / "config.json"), "--out-dir", str(self.root / "cli")],
            capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
        state = json.loads(result.stdout)
        self.assertEqual(state["sha256"], run.digest(Path(state["manifest"])))


    def test_readonly_audit_reports_observed_counts_not_scientific_success(self):
        root = self.root / "ordinary-audit"
        prepared_root = root / "prepared"
        prepared_root.mkdir(parents=True)
        prepared = run.seal(prepared_root, {
            "state": "awaiting_b_response", "case_id": "audit-case",
        })
        reviewed_root = root / "reviewed"
        for route, sources in (("a", {}), ("b", {"PMID:123": {"title": "test"}})):
            directory = reviewed_root / route
            directory.mkdir(parents=True)
            run.write(directory / "prompt.json", {
                "source_review_records": {"records": sources},
            })
        reviewed = run.seal(reviewed_root, {
            "state": "awaiting_independent_reviews", "case_id": "audit-case",
            "prepared": prepared["manifest"], "prepared_sha256": prepared["sha256"],
            "sessions": {route: {"prompt": str(reviewed_root / route / "prompt.json")}
                         for route in ("a", "b")},
            "independent_review_routes": ["b"], "machine_reviews": {"a": {"reviewer_model": "not_invoked"}},
        })
        final_root = root / "finished"
        final_root.mkdir()
        b_path = final_root / "route_b_execution_manifest.json"
        run.write(b_path, {
            "case_id": "audit-case", "prompt_count": 2, "completed_count": 1,
            "not_evaluated_count": 1,
            "partial_response_provenance": [{"attempt_index": 0}], "retry_history": [],
        })
        wet_path = final_root / "wet_candidates.tsv"
        wet_path.write_text("candidate_gene_id\tmembership_basis\nG1\troute_review_candidate\n"
                            "G2\tinterval_review_coverage\n", encoding="utf-8")
        finished = run.seal(final_root, {
            "state": "first_stage_frozen", "case_id": "audit-case",
            "review_stage": {"path": reviewed["manifest"], "sha256": reviewed["sha256"]},
            "route_b_execution": {"state": "partial", "path": str(b_path), "sha256": run.digest(b_path)},
            "handoff": {"adopted_stage": "first_stage", "wet": {
                "handoff_state": "boundary_review_requested",
                "second_stage_state": "requested",
                "wet_candidates": str(wet_path), "wet_candidates_sha256": run.digest(wet_path),
            }},
        })
        seal_before = Path(finished["manifest"]).read_bytes()
        summary = run.audit(Path(finished["manifest"]), finished["sha256"])
        self.assertEqual(summary["case_id"], "audit-case")
        self.assertEqual(summary["source_records_presented_to_review"], {"a": 0, "b": 1})
        self.assertEqual(summary["route_b_execution"]["not_evaluated_count"], 1)
        self.assertEqual(summary["route_b_execution"]["completed_count"], 1)
        self.assertEqual(summary["no_source_machine_not_executed_routes"], ["a"])
        self.assertEqual(summary["independent_review_required_routes"], ["b"])
        self.assertEqual(summary["wet"]["candidate_count"], 2)
        self.assertEqual(summary["wet"]["handoff_state"], "boundary_review_requested")
        self.assertIn("truth_retention", summary["not_inferred"])
        self.assertIn("walltime", summary["not_inferred"])
        cli = subprocess.run([sys.executable, str(HERE / "practical_packet_run.py"), "audit",
            "--finished", finished["manifest"], "--sha256", finished["sha256"]],
            capture_output=True, text=True)
        self.assertEqual(cli.returncode, 0, cli.stderr)
        self.assertEqual(json.loads(cli.stdout), summary)
        self.assertEqual(Path(finished["manifest"]).read_bytes(), seal_before)
        self.assertFalse((root / "audit-output").exists())
        wet_path.write_text(wet_path.read_text(encoding="utf-8") + "G3\tunknown\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            run.audit(Path(finished["manifest"]), finished["sha256"])

    def test_readonly_audit_does_not_infer_zero_wet_or_complete_b(self):
        root = self.root / "audit-missing"
        root.mkdir()
        prepared_root, reviewed_root, finished_root = (root / x for x in ("prepared", "reviewed", "finished"))
        prepared_root.mkdir()
        reviewed_root.mkdir()
        finished_root.mkdir()
        prepared = run.seal(prepared_root, {"case_id": "c1", "state": "awaiting_b_response"})
        reviewed = run.seal(reviewed_root, {
            "case_id": "c1", "state": "awaiting_independent_reviews",
            "prepared": prepared["manifest"], "prepared_sha256": prepared["sha256"],
        })
        finished = run.seal(finished_root, {
            "case_id": "c1", "state": "first_stage_frozen",
            "review_stage": {"path": reviewed["manifest"], "sha256": reviewed["sha256"]},
            "route_b_execution": {"state": "not_recorded"},
            "handoff": {"wet": {"handoff_state": "boundary_review_requested"}},
        })
        summary = run.audit(Path(finished["manifest"]), finished["sha256"])
        self.assertIsNone(summary["wet"]["candidate_count"])
        self.assertIsNone(summary["route_b_execution"]["completed_count"])
        self.assertEqual(summary["route_b_execution"]["state"], "not_recorded")


if __name__ == "__main__":
    unittest.main()
