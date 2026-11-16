"""Synthetic scope and isolated existing-selection consumer regressions."""
import copy
import csv
import hashlib
import html
import json
import random
import re
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from candidate_qtl.boundary_tier_review import (
    partition_relevance_tiers, freeze_boundary_scope, consume_boundary_selection,
    write_wet_review_report, _load_scope, EVIDENCE_SCOPED_POLICY,
)
from candidate_qtl.relevance_tiers import T1,T2,T3,T4,UNRESOLVED,NEW_TIERS


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(*counts):
    return [dict(candidate_gene_id=f"t{t}_{i:03}",relevance_tier=tier,
                 support_state="review" if t>1 else "verified",copy_state="co_preferred",
                 source_ids=f"PMID:{t+1}",rank_contribution="none")
            for t,(tier,count) in enumerate(zip((T1,T2,T3,T4,UNRESOLVED),counts)) for i in range(count)]


def write(path,rows,columns=None):
    columns=columns or list(rows[0])
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    return str(path)


def candidate_table_rows(report):
    section = report.split("## 候補一覧", 1)[1].split("## 判断・来歴の詳細", 1)[0]
    table = [line for line in section.splitlines() if line.startswith("| ")]
    return [[html.unescape(re.sub(r"\\([\\|`*_\[\]])", r"\1", cell.strip()))
             for cell in re.split(r"(?<!\\)\|", line[1:-1])]
            for line in table[2:]]


class PartitionTests(unittest.TestCase):
    def test_evidence_scoped_pool_keeps_every_lower_seed_and_excludes_identity_only(self):
        rows = inventory(0, 26, 3, 300, 300)
        eligible = {row["candidate_gene_id"] for row in rows if row["relevance_tier"] in {T2, T3}}
        result = partition_relevance_tiers(rows, policy_version=EVIDENCE_SCOPED_POLICY,
                                           eligible_gene_ids=eligible)
        self.assertEqual(len(result["frozen"]), 0)
        self.assertEqual(len(result["boundary"]), 29)
        self.assertEqual(len(result["lower"]), 600)
        self.assertEqual({T2, T3}, {row["relevance_tier"] for row in result["boundary"]})
        self.assertEqual(eligible, set(result["eligible_gene_ids"]))

    def test_evidence_scoped_short_pool_does_not_pad_from_identity_only(self):
        rows = inventory(0, 8, 4, 200, 200)
        eligible = {row["candidate_gene_id"] for row in rows if row["relevance_tier"] in {T2, T3}}
        result = partition_relevance_tiers(rows, policy_version=EVIDENCE_SCOPED_POLICY,
                                           eligible_gene_ids=eligible)
        self.assertEqual(12, len(result["frozen"]))
        self.assertFalse(result["boundary"])
        self.assertEqual(400, len(result["lower"]))

    def test_3_8_35_only_t3_boundary(self):
        result=partition_relevance_tiers(inventory(3,8,35,10))
        self.assertEqual([len(result[k]) for k in ("frozen","boundary","lower")],[11,35,10])
        self.assertEqual(result["boundary_tier"],T3)

    def test_21_to_25_never_cut_to_20(self):
        for count in range(21,26):
            result=partition_relevance_tiers(inventory(3,count-3,30))
            self.assertEqual(len(result["frozen"]),count)
            self.assertFalse(result["boundary"])

    def test_t1_only_overflow_is_not_reviewed(self):
        result=partition_relevance_tiers(inventory(31,50))
        self.assertEqual(len(result["frozen"]),31)
        self.assertEqual(result["state"],"known_rescue_over_upper_retained")
        self.assertFalse(result["boundary"])

    def test_unknown_trait_advances_to_identity_and_unresolved(self):
        result=partition_relevance_tiers(inventory(0,0,2,3,30))
        self.assertEqual(len(result["frozen"]),5)
        self.assertEqual(result["boundary_tier"],UNRESOLVED)
        self.assertEqual(len(result["boundary"]),30)

    def test_large_families_do_not_collapse_copies(self):
        result=partition_relevance_tiers(inventory(0,35))
        self.assertEqual(len(result["boundary"]),35)
        self.assertTrue(all(r["copy_state"]=="co_preferred" for r in result["boundary"]))

    def test_permutation_id_rename_and_input_immutable(self):
        rows=inventory(3,8,35,12);original=copy.deepcopy(rows)
        expected=partition_relevance_tiers(rows)
        random.Random(1).shuffle(rows)
        self.assertEqual(partition_relevance_tiers(rows),expected)
        self.assertEqual(sorted(rows,key=lambda r:r["candidate_gene_id"]),original)
        for row in rows: row["candidate_gene_id"]="reversed_"+row["candidate_gene_id"][::-1]
        changed=partition_relevance_tiers(rows)
        self.assertEqual([len(changed[k]) for k in ("frozen","boundary","lower")],[11,35,12])

    def test_no_arbitrary_frozen_lower_tier_prefix(self):
        with self.assertRaisesRegex(ValueError,"skip unreviewed"):
            partition_relevance_tiers(inventory(3,8,35),frozen_gene_ids=["t2_000"])

    def test_invalid_inputs_fail_closed(self):
        for target,upper in ((0,25),(26,25),(True,25),(20,25.5)):
            with self.assertRaises(ValueError): partition_relevance_tiers([],target=target,upper=upper)
        for rows in ([{"candidate_gene_id":"g","support_tier":"tier1_verified_support"}],
                     [dict(candidate_gene_id="g",relevance_tier=T1)]*2):
            with self.assertRaises(ValueError): partition_relevance_tiers(rows)


class ReportBindingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.wet = self.root / "wet.tsv"
        write(self.wet, [{
            "candidate_gene_id": "g1", "relevance_tier": "T1",
            "support_state": "review", "membership_basis": "route_review_candidate",
        }])
        self.decisions = self.root / "decisions.tsv"
        write(self.decisions, [{
            "candidate_gene_id": "g1", "applied_decision": "adopt",
            "independent_reason": "retained for review",
        }])
        self.scope = self.root / "next-scope.json"
        self.scope.write_text(json.dumps({"state": "next"}), encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def test_candidate_table_preserves_roles_values_and_missing_metadata(self):
        # Synthetic namespaces deliberately distinguish the focal locus,
        # original reported gene and sequence proxy; no biological assertion.
        paper = {
            "source_id": "PMID:123456", "source_title": "MtA phenotype report",
            "reported_gene": "MtA", "reported_gene_id": "MtA.1",
            "reported_species": "Medicago truncatula", "reference_gene_id": "AtB",
            "gene_correspondence": "homolog", "initial_tier": "T2",
            "initial_tier_reason": "The original report motivates the focal hypothesis.",
        }
        claims = [{
            "candidate_gene_id": "PsC_z", "reference_gene": "AtB;AtC",
            "reference_species": "Arabidopsis thaliana", "source_id": "PMID:123456",
            "paper_reports_json": json.dumps([paper]),
        }, {
            "candidate_gene_id": "PsC_z", "reference_gene": "OsD",
            "reference_species": "Oryza sativa", "source_id": "PMID:123456",
            "paper_reports_json": json.dumps([{**paper, "reported_gene": "MtD", "reported_gene_id": "MtD.1",
                "reported_species": "Lotus japonicus", "reference_gene_id": "OsD"}]),
        }, {
            "candidate_gene_id": "PsC_z", "reference_gene": "LjE",
            "reference_species": "Lotus japonicus", "source_id": "PMID:123456",
            "paper_reports_json": json.dumps([{**paper, "reported_gene": "MtF", "reported_gene_id": "MtF.1",
                "reference_gene_id": "UnboundX"}]),
        }]
        rows = [{
            "candidate_gene_id": "PsC_z", "gene_name": "PsC_name",
            "chrom": "chr05", "start": "101", "end": "300",
            "relevance_tier": NEW_TIERS[1], "support_state": "review",
            "membership_basis": "route_review_candidate", "copy_state": "paralog_ambiguous",
            # These independently aggregated fields cannot be zipped into pairs.
            "reference_gene": "AtB;AtC;OsD;LjE", "reference_species": "Oryza sativa;Arabidopsis thaliana;Lotus japonicus",
            "source_ids": "PMID:123456;DOI:10.1234/qtl.test;PMCID:PMC1234567",
            "tier_source_id": "DOI:10.1234/qtl.part;2",
            "paper_evidence_summary_json": json.dumps([
                {"source_id": source} for source in ("PMID:123456", "DOI:10.1234/qtl.test", "PMCID:PMC1234567",
                                                     "DOI:10.1234/qtl.part;2")]),
            "tier_reason": "候補機構 | copy\nの解決待ち",
            "relevance_evidence_json": json.dumps(claims, ensure_ascii=False),
        }]
        for gene, tier in (("PsC_a", UNRESOLVED), ("PsC_m", NEW_TIERS[5])):
            rows.append({**{key: "" for key in rows[0]}, "candidate_gene_id": gene,
                         "relevance_tier": tier})
        write(self.wet, rows)
        original = self.wet.read_bytes()
        result = write_wet_review_report(
            str(self.wet), str(self.root / "report.md"), "case", "boundary_review_requested", "first_stage",
            adopted_gene_ids=["PsC_z"], pending_gene_ids=["PsC_a"],
        )
        report = Path(result["path"]).read_text(encoding="utf-8")
        self.assertIn("| 遺伝子名・ID | 位置 | 相同遺伝子（種・ID） | Tier | 論文 | 一言メモ |", report)
        displayed = candidate_table_rows(report)
        self.assertEqual(3, len(displayed))
        self.assertTrue(all(len(row) == 6 for row in displayed))
        for actual, supplied in zip(displayed, rows):
            self.assertIn(supplied["candidate_gene_id"], actual[0])
            self.assertEqual(supplied["relevance_tier"], actual[3])
        self.assertIn("PsC_name", displayed[0][0])
        self.assertNotIn("MtA", displayed[0][0])
        self.assertNotIn("AtB", displayed[0][0])
        self.assertEqual("chr05:101-300", displayed[0][1])
        homologs = displayed[0][2]
        original_report, proxy = homologs.split("報告:", 1)[1].split("→", 1)
        self.assertIn("MtA", original_report)
        self.assertIn("Medicago truncatula", original_report)
        self.assertNotIn("Arabidopsis thaliana", original_report)
        self.assertRegex(proxy, r"参照:\s*AtB[^;→]*Arabidopsis thaliana")
        self.assertNotIn("AtC", proxy.split(";", 1)[0])
        self.assertNotRegex(proxy, r"参照:\s*AtB[^;→]*Medicago truncatula")
        self.assertNotRegex(homologs, r"AtB[^;→]*Oryza sativa")
        self.assertNotRegex(homologs, r"OsD[^;→]*Arabidopsis thaliana")
        self.assertIn("MtD", homologs)
        self.assertIn("Lotus japonicus", homologs)
        self.assertRegex(homologs, r"参照:\s*OsD[^;→]*Oryza sativa")
        unbound = homologs.split("原報の対応先ID:", 1)[1].split(";", 1)[0].strip()
        self.assertEqual("UnboundX（参照との対応は未記録）", unbound)
        self.assertNotIn("Lotus japonicus", unbound)
        self.assertNotIn("→ 参照: UnboundX", homologs)
        self.assertIn("参照: LjE（Lotus japonicus）", homologs)
        for url in ("https://pubmed.ncbi.nlm.nih.gov/123456/", "https://doi.org/10.1234/qtl.test",
                    "https://pmc.ncbi.nlm.nih.gov/articles/PMC1234567/", "https://doi.org/10.1234/qtl.part%3B2"):
            self.assertIn(url, displayed[0][4])
        self.assertNotIn("(https://doi.org/10.1234/qtl.part)", displayed[0][4])
        self.assertIn("候補機構 | copy の解決待ち", displayed[0][5])
        self.assertIn("paralog_ambiguous", displayed[0][5])
        self.assertEqual(["未記録"] * 2, [displayed[1][1], displayed[1][4]])
        self.assertEqual("未記録", displayed[1][2])
        self.assertTrue(displayed[1][5].startswith("保留中"))
        self.assertIn("未記録", displayed[1][5])
        self.assertIn("保持（採用状態未指定）", displayed[2][5])
        self.assertEqual(original, self.wet.read_bytes())
        self.assertEqual(sha(self.wet), result["wet_candidates_sha256"])
        self.assertEqual((3, 1, 1), (result["candidate_count"], result["adopted_candidate_count"], result["pending_candidate_count"]))

    def test_report_rejects_decision_or_next_scope_generation_mismatch(self):
        with self.assertRaisesRegex(ValueError, "sealed selection hash"):
            write_wet_review_report(
                str(self.wet), str(self.root / "bad-decisions.md"), "case", "state", "stage",
                review_decisions_path=self.decisions,
                review_decisions_sha256="0" * 64,
            )
        with self.assertRaisesRegex(ValueError, "sealed SHA-256"):
            write_wet_review_report(
                str(self.wet), str(self.root / "missing-decision-sha.md"), "case", "state", "stage",
                review_decisions_path=self.decisions,
            )
        with self.assertRaisesRegex(ValueError, "changed after selection"):
            write_wet_review_report(
                str(self.wet), str(self.root / "bad-scope.md"), "case", "state", "stage",
                next_boundary_scope_manifest={"path": str(self.scope), "sha256": "0" * 64},
            )
        with self.assertRaisesRegex(ValueError, "sealed SHA-256"):
            write_wet_review_report(
                str(self.wet), str(self.root / "missing-scope-sha.md"), "case", "state", "stage",
                next_boundary_scope_manifest={"path": str(self.scope)},
            )

    def test_report_prefers_recorded_boundary_reasons_without_changing_inputs(self):
        write(self.wet, [{
            "candidate_gene_id": "g1", "relevance_tier": T2,
            "tier_reason": "initial tier rationale", "reason": "initial proposal rationale",
        }])
        write(self.decisions, [{
            "candidate_gene_id": "g1", "applied_decision": "adopt",
            "proposal_reason": "retain this copy | boundary comparison",
            "independent_reason": "independent recorded reason",
        }])
        wet_before, decisions_before = self.wet.read_bytes(), self.decisions.read_bytes()
        report = self.root / "adopted.md"
        result = write_wet_review_report(
            str(self.wet), str(report), "case", "completed_evidence_scoped_selection",
            "boundary_tier_review", adopted_gene_ids={"g1"}, review_history_state="complete",
            review_decisions_path=self.decisions, review_decisions_sha256=sha(self.decisions),
        )
        text = report.read_text(encoding="utf-8")
        self.assertIn("retain this copy \\| boundary comparison; independent recorded reason", text)
        # The six-column view retains the original Tier rationale.
        # The adoption detail must prefer the recorded boundary decision.
        self.assertIn("initial tier rationale", text.split("## 判断・来歴の詳細")[0])
        detail = text.split("## 判断・来歴の詳細")[1]
        self.assertNotIn("initial tier rationale", detail)
        self.assertNotIn("initial proposal rationale", detail)
        self.assertEqual(wet_before, self.wet.read_bytes())
        self.assertEqual(decisions_before, self.decisions.read_bytes())
        self.assertEqual(sha(self.wet), result["wet_candidates_sha256"])
        self.assertEqual(1, result["adopted_candidate_count"])
        self.assertEqual(0, result["pending_candidate_count"])

    def test_report_keeps_initial_reason_when_decision_reason_is_unavailable(self):
        write(self.wet, [{"candidate_gene_id": "g1", "tier_reason": "initial tier rationale"}])
        write(self.decisions, [{"candidate_gene_id": "g1", "applied_decision": "adopt"}])
        for with_decision in (False, True):
            with self.subTest(with_decision=with_decision):
                report = self.root / f"fallback-{with_decision}.md"
                write_wet_review_report(
                    str(self.wet), str(report), "case", "boundary_review_requested", "first_stage",
                    **({"review_decisions_path": self.decisions,
                        "review_decisions_sha256": sha(self.decisions)} if with_decision else {}),
                )
                self.assertIn("initial tier rationale", report.read_text(encoding="utf-8"))

    def test_completed_evidence_scoped_report_points_to_experimental_confirmation(self):
        report = self.root / "completed.md"
        write_wet_review_report(
            str(self.wet), str(report), "case", "completed_evidence_scoped_selection",
            "boundary_tier_review", adopted_gene_ids={"g1"}, review_history_state="complete",
            review_decisions_path=self.decisions, review_decisions_sha256=sha(self.decisions),
        )
        text = report.read_text(encoding="utf-8")
        self.assertIn("- next action: Inspect the adopted Wet set and plan experimental confirmation.", text)
        self.assertNotIn("before Wet adoption", text)


class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.root=Path(self.temp.name)
        self.source=self.root/"inventory.tsv"
        write(self.source,inventory(3,8,35,12,5))
        self.scope=freeze_boundary_scope(str(self.source),str(self.root/"scope"),
                    inventory_sha256=sha(self.source),case_id="test")
    def tearDown(self): self.temp.cleanup()

    def review(self,keep=10,*,disagreement=False,scope=None,name="review"):
        scope=scope or self.scope
        boundary=pd.read_csv(scope["artifacts"]["boundary"]["path"],sep="\t").to_dict("records")
        proposal=[]
        for index,row in enumerate(boundary):
            proposal.append(dict(candidate_gene_id=row["candidate_gene_id"],decision="adopt" if index<keep else "defer",
                reason="synthetic mechanistic comparison",executor_id="executor",executor_model="model-a",
                source_wet_sha256=scope["artifacts"]["boundary"]["sha256"]))
        p=self.root/f"{name}.tsv";write(p,proposal)
        independent=[dict(candidate_gene_id=r["candidate_gene_id"],decision=("adopt" if disagreement else r["decision"]),
                     reason="independent synthetic comparison",reviewer_id="reviewer",reviewer_model="model-b",
                     source_wet_sha256=r["source_wet_sha256"],proposal_sha256=sha(p)) for r in proposal]
        i=self.root/f"{name}-independent.tsv";write(i,independent)
        return str(p),str(i)

    def consume(self,proposal,independent,*,scope=None,name="selected",**kwargs):
        scope=scope or self.scope
        return consume_boundary_selection(scope["manifest"],proposal,str(self.root/name),
            scope_sha256=scope["manifest_sha256"],case_id="test",reviewer_id="reviewer",reviewer_model="model-b",
            executor_id="executor",executor_model="model-a",independent_review_arg=independent,**kwargs)

    def test_post_review_provisional_preserves_frozen_and_protected_c_rows(self):
        rows = inventory(3,30,4,35,5)
        source = self.root / "protected-source.tsv"
        write(source, rows)
        protected = "t3_000"
        scope = freeze_boundary_scope(str(source), str(self.root / "protected-scope"),
            inventory_sha256=sha(source), case_id="test", protected_gene_ids=[protected])
        screen = self.root / "screen.tsv"
        write(screen, [{"case_id":"test", "unit_id":"full", "member_interval_genes":";".join(r["candidate_gene_id"] for r in rows),
                       "screen_state":"screened", "answer_sha256":"a"*64}])
        p,i = self.review(keep=2, scope=scope, name="protected-review")
        result = self.consume(p,i,scope=scope,allow_provisional_wet=True,
            provisional_screen_ledger_arg=str(screen),provisional_screen_ledger_sha256=sha(screen))
        self.assertEqual("provisional_adopted_wet",result["state"])
        final = pd.read_csv(result["wet_candidates"]["path"],sep="\t",dtype=str).to_dict("records")
        self.assertEqual(10,len(final))
        self.assertTrue({protected,*[f"t0_{i:03}" for i in range(3)]} <= {r["candidate_gene_id"] for r in final})
        self.assertTrue(all(row in rows for row in final))

    def test_identity_comparison_evidence_blocks_post_review_provisional(self):
        rows = inventory(0,30,4,35,5)
        for row in rows:
            row["review_claim_ids"] = "identity-comparison" if row["candidate_gene_id"] == "t3_001" else ""
        source = self.root / "evidence-source.tsv"
        write(source, rows)
        scope = freeze_boundary_scope(str(source),str(self.root / "evidence-scope"),inventory_sha256=sha(source),case_id="test")
        screen = self.root / "screen.tsv"
        write(screen,[{"case_id":"test", "unit_id":"full", "member_interval_genes":";".join(r["candidate_gene_id"] for r in rows),
                       "screen_state":"screened", "answer_sha256":"a"*64}])
        p,i = self.review(keep=2,scope=scope,name="evidence-review")
        result = self.consume(p,i,scope=scope,allow_provisional_wet=True,
            provisional_screen_ledger_arg=str(screen),provisional_screen_ledger_sha256=sha(screen))
        self.assertEqual("next_boundary_review_required",result["state"])
        self.assertEqual("identity_tier_has_comparison_evidence",result["provisional_wet"]["reason"])
        self.assertEqual(41,result["output_candidate_count"])

    def test_scope_hashes_and_no_overwrite(self):
        loaded,_,plan=_load_scope(self.scope["manifest"],self.scope["manifest_sha256"],"test")
        self.assertEqual(len(plan["frozen"]),11)
        self.assertEqual(loaded["inventory"]["sha256"],sha(self.source))
        with self.assertRaises(FileExistsError):
            freeze_boundary_scope(str(self.source),str(self.root/"scope"),inventory_sha256=sha(self.source),case_id="test")

    def test_unversioned_scope_replays_legacy_lower_tier_continuation(self):
        original_bytes = Path(self.scope["manifest"]).read_bytes()
        legacy = json.loads(original_bytes)
        legacy.pop("policy_version")
        legacy.pop("eligible_gene_ids")
        path = self.root / "unversioned-scope.json"
        path.write_text(json.dumps(legacy), encoding="utf-8")
        legacy.update(manifest=str(path), manifest_sha256=sha(path))
        proposal, independent = self.review(keep=8, scope=legacy, name="unversioned-review")
        result = self.consume(proposal, independent, scope=legacy, name="unversioned-selected")
        self.assertEqual("next_boundary_review_required", result["state"])
        next_scope = json.loads(Path(result["next_scope_manifest"]["path"]).read_text())
        self.assertEqual({"frozen": 19, "boundary": 12, "lower": 5}, next_scope["counts"])
        self.assertEqual("legacy_whole_tier_v1", next_scope["policy_version"])
        self.assertEqual(original_bytes, Path(self.scope["manifest"]).read_bytes())

    def test_source_tamper_rejected_before_consumer(self):
        p,i=self.review();self.source.write_text(self.source.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"hash mismatch"):
            self.consume(p,i)
        self.assertFalse((self.root/"selected").exists())

    def test_higher_row_change_rejected_even_with_updated_artifact_hash(self):
        path=Path(self.scope["artifacts"]["frozen"]["path"])
        path.write_text(path.read_text().replace("co_preferred","preferred"))
        manifest=Path(self.scope["manifest"]);payload=json.loads(manifest.read_text())
        payload["artifacts"]["frozen"]["sha256"]=sha(path);manifest.write_text(json.dumps(payload))
        with self.assertRaisesRegex(ValueError,"rows differ"):
            _load_scope(str(manifest),sha(manifest),"test")

    def test_higher_tier_review_row_rejected(self):
        p,i=self.review()
        rows=pd.read_csv(p,sep="\t").to_dict("records");rows[0]["candidate_gene_id"]="t0_000";write(Path(p),rows)
        with self.assertRaisesRegex(ValueError,"higher/lower"):
            self.consume(p,i)
        self.assertFalse((self.root/"selected").exists())

    def test_existing_consumer_preserves_higher_and_original_bytes(self):
        source=self.source.read_bytes();frozen=Path(self.scope["artifacts"]["frozen"]["path"]).read_bytes()
        p,i=self.review(keep=10);result=self.consume(p,i)
        self.assertEqual(result["state"],"within_target_range")
        self.assertEqual(result["output_candidate_count"],21)
        self.assertEqual(self.source.read_bytes(),source)
        self.assertEqual(Path(self.scope["artifacts"]["frozen"]["path"]).read_bytes(),frozen)
        final=pd.read_csv(result["wet_candidates"]["path"],sep="\t",dtype=str).to_dict("records")
        self.assertEqual(final[:11],pd.read_csv(self.scope["artifacts"]["frozen"]["path"],sep="\t",dtype=str).to_dict("records"))
        self.assertEqual(len(pd.read_csv(result["deferral_ledger"]["path"],sep="\t")),25)

    def test_request_policy_separates_adoption_from_support_and_copy_verification(self):
        from candidate_qtl.boundary_tier_review import build_boundary_review_request

        request = build_boundary_review_request(
            self.scope["manifest"], str(self.root / "request-policy"),
            scope_sha256=self.scope["manifest_sha256"], case_id="test",
        )
        payload = json.loads(Path(request["request"]).read_text())
        policy = payload["policy"]
        self.assertTrue(any("adopt means retain" in item and "does not promote support" in item for item in policy))
        self.assertTrue(any("unresolved copy alone" in item and "not_evaluated or defer" in item for item in policy))
        self.assertTrue(any("explicit agreed experiment-priority" in item for item in policy))

    def test_agreed_adoption_keeps_ambiguous_review_copies_without_support_promotion(self):
        boundary = pd.DataFrame([
            {"candidate_gene_id": "G1", "relevance_tier": "T3_functionally_plausible",
             "support_state": "review", "copy_state": "paralog_ambiguous", "route": "a"},
            {"candidate_gene_id": "G2", "relevance_tier": "T3_functionally_plausible",
             "support_state": "review", "copy_state": "paralog_ambiguous", "route": "a"},
        ])
        source = self.root / "ambiguous-copy-input.tsv"
        boundary.to_csv(source, sep="\t", index=False)
        scope = freeze_boundary_scope(str(source), str(self.root / "ambiguous-copy-scope"),
                                      inventory_sha256=sha(source), case_id="ambiguous-copy",
                                      target=1, upper=1)
        boundary_path = Path(scope["artifacts"]["boundary"]["path"])
        self.assertEqual(2, len(pd.read_csv(boundary_path, sep="\t")))
        wet_hash = sha(boundary_path)
        proposal = self.root / "ambiguous-copy-proposal.tsv"
        pd.DataFrame([
            {"candidate_gene_id": gene, "decision": "adopt", "reason": "retain plausible copy pending resolution",
             "source_wet_sha256": wet_hash, "executor_id": "executor", "executor_model": "model-a"}
            for gene in ("G1", "G2")
        ]).to_csv(proposal, sep="\t", index=False)
        proposal_hash = sha(proposal)
        independent = self.root / "ambiguous-copy-independent.tsv"
        pd.DataFrame([
            {"candidate_gene_id": gene, "decision": "adopt", "reason": "retain both copies for Wet review",
             "source_wet_sha256": wet_hash, "proposal_sha256": proposal_hash,
             "reviewer_id": "reviewer", "reviewer_model": "model-b"}
            for gene in ("G1", "G2")
        ]).to_csv(independent, sep="\t", index=False)
        result = consume_boundary_selection(
            scope["manifest"], str(proposal), str(self.root / "ambiguous-copy-out"),
            scope_sha256=scope["manifest_sha256"], case_id="ambiguous-copy",
            reviewer_id="reviewer", reviewer_model="model-b", executor_id="executor",
            executor_model="model-a", independent_review_arg=str(independent),
        )
        selected = pd.read_csv(result["wet_candidates"]["path"], sep="\t", dtype=str)
        self.assertEqual({"G1", "G2"}, set(selected["candidate_gene_id"]))
        self.assertTrue((selected["support_state"] == "review").all())
        self.assertTrue((selected["copy_state"] == "paralog_ambiguous").all())
        ledger = pd.read_csv(result["deferral_ledger"]["path"], sep="\t", dtype=str)
        self.assertEqual(0, len(ledger))

    def test_evidence_scoped_consumer_accepts_mixed_t3_t4_seeds_with_ineligible_same_tier_rows(self):
        rows = inventory(1, 26, 2, 2, 4)
        eligible = {row["candidate_gene_id"] for row in rows if row["relevance_tier"] == T2}
        eligible.add("t4_000")
        eligible.add("t3_000")  # T3 seed below the T2 overflow.
        source = self.root / "evidence-scoped-source.tsv"
        write(source, rows)
        scope = freeze_boundary_scope(
            str(source), str(self.root / "evidence-scoped-scope"),
            inventory_sha256=sha(source), case_id="test",
            policy_version=EVIDENCE_SCOPED_POLICY, eligible_gene_ids=eligible,
        )
        boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
        lower = pd.read_csv(scope["artifacts"]["lower"]["path"], sep="\t")
        self.assertIn(T4, set(boundary.relevance_tier))
        self.assertTrue({"t3_001", "t4_001"} <= set(lower.candidate_gene_id))
        proposal, independent = self.review(keep=len(boundary), scope=scope, name="mixed-evidence")
        result = self.consume(proposal, independent, scope=scope, name="mixed-evidence-selected")
        selected = pd.read_csv(result["wet_candidates"]["path"], sep="\t")
        self.assertEqual(len(boundary) + 1, len(selected))  # Frozen T1 plus all eligible review seeds.
        self.assertIsNone(result["next_scope_manifest"])
        self.assertEqual("completed_evidence_scoped_selection", result["state"])

    def test_under_target_advances_and_next_boundary_is_t4(self):
        p,i=self.review(keep=8);result=self.consume(p,i)
        self.assertEqual(result["state"],"next_boundary_review_required")
        self.assertEqual(result["next_boundary_tier"],T4)
        next_path=Path(result["next_scope_manifest"]["path"])
        next_scope=json.loads(next_path.read_text());next_scope.update(manifest=str(next_path),manifest_sha256=sha(next_path))
        self.assertEqual(next_scope["counts"],{"frozen":19,"boundary":12,"lower":5})
        p2,i2=self.review(keep=1,scope=next_scope,name="next-review")
        second=self.consume(p2,i2,scope=next_scope,name="second-selection")
        self.assertEqual(second["output_candidate_count"],20)
        self.assertEqual(second["state"],"within_target_range")

    def test_disagreement_retains_overflow_instead_of_count_pruning(self):
        p,i=self.review(keep=10,disagreement=True);result=self.consume(p,i)
        self.assertEqual(result["output_candidate_count"],46)
        self.assertEqual(result["state"],"boundary_review_incomplete_retained")

    def test_deferral_without_independent_review_still_rejected(self):
        p,_=self.review()
        with self.assertRaisesRegex(ValueError,"independent review"):
            self.consume(p,None)
        self.assertFalse((self.root/"selected/wet_candidates.tsv").exists())

    def test_case_and_scope_hash_mismatch(self):
        with self.assertRaisesRegex(ValueError,"manifest hash"):
            _load_scope(self.scope["manifest"],"0"*64,"test")
        with self.assertRaisesRegex(ValueError,"case mismatch"):
            _load_scope(self.scope["manifest"],self.scope["manifest_sha256"],"another")

if __name__=="__main__": unittest.main()
