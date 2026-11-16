"""#1005: light T4/T5 follow-up planning, without changing canonical Tier/Wet."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.agentic_route_b_initial import summarize_initial
from scripts.agentic_route_b_triage import CONTRACT, plan_light_followup
from tests.test_agentic_route_b_research import fixture

ROOT = Path(__file__).resolve().parents[1]


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


class LightTierFollowupTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.bundle, _ = fixture()
        self.bundle["units"][0]["genes"][0]["reference_species"] = "Arabidopsis thaliana"
        self.bundle["units"][0]["genes"][1]["reference_species"] = "Arabidopsis thaliana"
        self.bundle_path = self.root / "bundle.json"
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())

    def judgment(self, gid="g1", state="assessed_claim", pointer="/units/0/genes/0",
                 fields=None, uid="u1"):
        return {
            "gene_id": gid, "unit_id": uid, "state": state,
            "reason": "Original researcher identified a specific functional bridge.",
            "used_material_pointers": [pointer] if pointer else [],
            "assessed_fields": fields if fields is not None else ["reference_function"],
        }

    def initial(self, assessments=None):
        path = self.root / "assessments.json"
        path.write_text(json.dumps({"bundle_sha256": self.bundle_sha,
                                    "assessments": assessments if assessments is not None else [self.judgment()]}),
                        encoding="utf-8")
        summary_dir = self.root / "initial"
        summarize_initial(str(self.bundle_path), self.bundle_sha, str(path), str(summary_dir))
        summary = summary_dir / "initial_assessment_summary.json"
        return summary

    @staticmethod
    def screen(members=("g1",), ref="ref1", species="Arabidopsis thaliana",
               state="no_concrete_upper_lead", **changes):
        row = {
            "screen_id": "S1", "reference_gene": ref, "reference_species": species,
            "member_gene_ids": list(members), "state": state,
            "retrieval_state": "completed", "queries": ["REF1 phenotype mutant plant development"],
            "checked_species_scopes": ["focal_species", "other_plants"],
            "screened_source_ids": [], "reason": "Search had no specific upper-Tier lead in the inspected results.",
        }
        row.update(changes)
        return row

    @staticmethod
    def assessed(gid="g1", tier="T4", ids=("S1",), uid="u1"):
        return {"gene_id": gid, "unit_id": uid, "provisional_tier": tier,
                "tier_reason": "Source-bound function suggests a specific conditional mechanism.",
                "screen_ids": list(ids)}

    def screening(self, summary, screens=None, decisions=None):
        data = {
            "contract": CONTRACT, "case_id": "fixture",
            "initial_summary_sha256": sha(summary.read_bytes()),
            "reference_screens": screens if screens is not None else [self.screen()],
            "gene_assessments": decisions if decisions is not None else [self.assessed()],
        }
        path = self.root / f"screening-{len(list(self.root.glob('screening-*.json'))):03d}.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    def run_plan(self, summary, path):
        out = self.root / "followup"
        return plan_light_followup(str(summary), str(path), str(out))

    def output(self):
        return json.loads((self.root / "followup/followup_plan.json").read_text(encoding="utf-8"))

    def test_grounded_t4_no_lead_allows_only_optional_deep_read_stop(self):
        summary = self.initial()
        screen = self.screening(summary)
        before_bundle = self.bundle_path.read_bytes()
        result = self.run_plan(summary, screen)
        rows = {g["gene_id"]: g for g in self.output()["genes"]}
        self.assertEqual("light_stop_eligible", rows["g1"]["research_action"])
        self.assertEqual("not_evaluated", rows["g2"]["research_action"])
        self.assertEqual("not_evaluated", rows["g3"]["research_action"])
        self.assertEqual("not_evaluated", rows["g4"]["research_action"])
        self.assertEqual("requires_assessment", self.output()["units"][0]["research_action"])
        self.assertIn("No gene is excluded", self.output()["limitations"])
        self.assertEqual(self.bundle_sha, sha(self.bundle_path.read_bytes()))
        self.assertEqual(before_bundle, self.bundle_path.read_bytes())
        self.assertEqual(sha((self.root / "followup/followup_plan.json").read_bytes()),
                         result["plan"]["sha256"])
        self.assertEqual(4, result["gene_count"])

    def test_one_shared_reference_screen_covers_two_interpreted_copies(self):
        self.bundle["units"][0]["genes"][1]["reference_gene"] = "ref1"
        self.bundle["units"][0]["reference_material_cards"][0]["reference_species"] = ""
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        assessments = [
            self.judgment(),
            self.judgment("g2", pointer="/units/0/genes/1"),
            self.judgment("g2", pointer="/units/0/reference_material_cards/0",
                          fields=["function_description"]),
        ]
        summary = self.initial(assessments)
        screen = self.screening(summary, screens=[self.screen(("g1", "g2"))],
                                decisions=[self.assessed("g1"), self.assessed("g2")])
        self.run_plan(summary, screen)
        output = self.output()
        self.assertEqual(1, len(output["reference_screens"]))
        self.assertEqual(["g1", "g2"], output["reference_screens"][0]["member_gene_ids"])
        self.assertEqual("light_stop_eligible", output["units"][0]["research_action"])
        self.assertEqual(2, output["gene_actions"]["light_stop_eligible"])
        self.assertEqual({"g1", "g2"}, set(output["units"][0]["member_gene_ids"]))

    def test_unread_go_and_specific_lead_force_targeted_followup(self):
        self.bundle["units"][0]["genes"][0]["go_terms"] = "meristem specification"
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial()
        screen = self.screening(summary)
        self.run_plan(summary, screen)
        row = self.output()["genes"][0]
        self.assertEqual("targeted_followup", row["research_action"])
        self.assertIn("go_terms", row["unread_relevant_material"]["/units/0/genes/0"])
        # If a gene already carries a concrete upper lead, the no-lead screen
        # cannot suppress additional targeted evaluation.
        other = self.root / "second"
        second_rows = [self.judgment(), self.judgment(state="assessed_concrete_lead")]
        other.mkdir()
        # One documented lead is a distinct original assessment (no duplicate).
        second_rows[1]["reason"] = "A different cited developmental observation."
        alt_initial = self.root / "initial2"
        answers = self.root / "answers2.json"
        answers.write_text(json.dumps({"bundle_sha256": self.bundle_sha, "assessments": second_rows}))
        summarize_initial(str(self.bundle_path), self.bundle_sha, str(answers), str(alt_initial))
        summary2 = alt_initial / "initial_assessment_summary.json"
        screen2 = self.screening(summary2)
        plan_light_followup(str(summary2), str(screen2), str(other / "out"))
        alt = json.loads((other / "out/followup_plan.json").read_text())
        self.assertEqual("targeted_followup", alt["genes"][0]["research_action"])

    def test_unassessed_and_failed_screen_are_not_negative(self):
        summary = self.initial()
        for state, expected in (("not_assessed", "needs_light_screen"),
                                ("retrieval_failed", "not_evaluated")):
            with self.subTest(state=state):
                base = self.root / state
                base.mkdir()
                screen = self.screening(summary, screens=[self.screen(state=state)])
                output_dir = base / "plan"
                plan_light_followup(str(summary), str(screen), str(output_dir))
                rows = json.loads((output_dir / "followup_plan.json").read_text())["genes"]
                self.assertEqual(expected, rows[0]["research_action"])

    def test_other_unscreened_reference_prevents_early_stop(self):
        self.bundle["units"][0]["reference_material_cards"][0].update(
            {"interval_gene_id": "g1", "reference_gene": "refX", "reference_species": "Medicago truncatula"})
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial([self.judgment(),
                                self.judgment(pointer="/units/0/reference_material_cards/0",
                                              fields=["function_description"])])
        screen = self.screening(summary)
        self.run_plan(summary, screen)
        row = self.output()["genes"][0]
        self.assertEqual("needs_light_screen", row["research_action"])
        self.assertIn(["refX", "medicago truncatula"], row["unscreened_reference_pairs"])

    def test_t5_requires_annotation_and_t6_never_light_stops(self):
        summary = self.initial()
        screen = self.screening(summary, decisions=[self.assessed(tier="T5")])
        self.run_plan(summary, screen)
        self.assertEqual("needs_grounded_assessment", self.output()["genes"][0]["research_action"])
        screen2 = self.screening(summary, decisions=[self.assessed(tier="T6")])
        out2 = self.root / "plan2"
        plan_light_followup(str(summary), str(screen2), str(out2))
        self.assertEqual("not_evaluated", json.loads((out2 / "followup_plan.json").read_text())["genes"][0]["research_action"])

    def test_uninspected_frozen_literature_link_blocks_stop(self):
        self.bundle["units"][0]["literature_link_candidates"] = [
            {"reference_gene": "ref1", "reference_species": "Arabidopsis thaliana",
             "source_id": "PMID:987654", "title": "Potential developmental phenotype study"}
        ]
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial()
        screen = self.screening(summary)
        self.run_plan(summary, screen)
        gene = self.output()["genes"][0]
        self.assertEqual("needs_light_screen", gene["research_action"])
        self.assertEqual(["PMID:987654"], gene["unscreened_literature_handles"])
        # The originating screen may report that it actually checked the
        # known source; an uninspected link is not a scientific negative.
        checked = self.screening(summary, screens=[self.screen(
            screened_source_ids=["PMID:987654"])])
        out2 = self.root / "checked-paper"
        plan_light_followup(str(summary), str(checked), str(out2))
        gene2 = json.loads((out2 / "followup_plan.json").read_text())["genes"][0]
        self.assertEqual("light_stop_eligible", gene2["research_action"])

    def test_unresolved_frozen_query_status_blocks_stop(self):
        self.bundle["units"][0]["literature_query_status"] = [
            {"reference_gene": "ref1", "retrieval_status": "provider_failure"}
        ]
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial()
        screen = self.screening(summary)
        self.run_plan(summary, screen)
        gene = self.output()["genes"][0]
        self.assertEqual("needs_light_screen", gene["research_action"])
        self.assertEqual(["provider_failure"], gene["unresolved_literature_query_states"])

    def test_explicit_t5_annotation_could_light_stop_after_reading_function_too(self):
        self.bundle["units"][0]["genes"][0]["go_terms"] = "documented meristem initiation"
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        assessment = self.judgment(fields=["reference_function", "go_terms"])
        summary = self.initial([assessment])
        screen = self.screening(summary, decisions=[self.assessed(tier="T5")])
        self.run_plan(summary, screen)
        self.assertEqual("light_stop_eligible", self.output()["genes"][0]["research_action"])

    def test_upper_tier_lead_and_previously_known_upper_tier_never_light_stop(self):
        summary = self.initial()
        screen = self.screening(summary, screens=[self.screen(
            state="concrete_upper_lead",
            upper_lead="Specific reported developmental mutant phenotype in another species.")])
        self.run_plan(summary, screen)
        self.assertEqual("targeted_followup", self.output()["genes"][0]["research_action"])
        # A separate frozen planning output with a T1 lead is not lower-Tier
        # research, even though the paper-to-member inference is still provisional.
        other = self.root / "existing-upper"
        screen2 = self.screening(summary, decisions=[self.assessed(tier="T1")])
        plan_light_followup(str(summary), str(screen2), str(other))
        output = json.loads((other / "followup_plan.json").read_text())
        self.assertEqual("retain_upper_evidence", output["genes"][0]["research_action"])

    def test_bad_hash_species_member_and_unperformed_search_rejected(self):
        summary = self.initial()
        cases = [
            ([self.screen(species="Medicago truncatula")], [self.assessed()]),
            ([self.screen(members=("g4",))], [self.assessed()]),
            ([self.screen(queries=[])], [self.assessed()]),
            ([self.screen(retrieval_state="failed")], [self.assessed()]),
            ([self.screen()], [self.assessed(uid="u2")]),
            ([self.screen()], [self.assessed(ids=("unknown",))]),
        ]
        for screens, decisions in cases:
            with self.subTest(screens=screens, decisions=decisions):
                candidate = self.screening(summary, screens=screens, decisions=decisions)
                with self.assertRaises(ValueError):
                    self.run_plan(summary, candidate)
                self.assertFalse((self.root / "followup").exists())
        good = self.screening(summary)
        tampered = json.loads(good.read_text())
        tampered["initial_summary_sha256"] = "0000"
        good.write_text(json.dumps(tampered))
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            self.run_plan(summary, good)

    def test_canonical_cli_produces_auditable_plan(self):
        summary = self.initial()
        screen = self.screening(summary)
        p = subprocess.run([sys.executable, str(ROOT / "scripts/agentic_route_b.py"), "plan-light-followup",
                            "--summary", str(summary), "--screening", str(screen),
                            "--out", str(self.root / "followup")],
                           text=True, capture_output=True)
        self.assertEqual(0, p.returncode, p.stderr)
        self.assertEqual(1, json.loads(p.stdout)["gene_actions"]["light_stop_eligible"])
        self.assertTrue((self.root / "followup/followup_plan.json").is_file())

    def test_omitted_bound_upper_lead_or_failed_screen_cannot_authorize_stop(self):
        summary = self.initial()
        for state, expected in (("concrete_upper_lead", "targeted_followup"),
                                ("retrieval_failed", "not_evaluated"),
                                ("not_assessed", "needs_light_screen")):
            with self.subTest(state=state):
                screens = [self.screen(), self.screen(screen_id="S2", state=state,
                           upper_lead="A reported concrete mutant phenotype." if state == "concrete_upper_lead" else "")]
                screen = self.screening(summary, screens=screens)
                out = self.root / state
                plan_light_followup(str(summary), str(screen), str(out))
                row = json.loads((out / "followup_plan.json").read_text())["genes"][0]
                self.assertEqual(expected, row["research_action"])
                self.assertIn("S2", row["unlinked_screen_ids"])

    def test_focal_and_other_species_screen_scope_must_be_recorded(self):
        summary = self.initial()
        for scopes in ([], ["focal_species"], ["other_plants"]):
            with self.subTest(scopes=scopes):
                screen = self.screening(summary, screens=[self.screen(checked_species_scopes=scopes)])
                out = self.root / ("scope-" + str(len(list(self.root.glob('scope-*')))))
                plan_light_followup(str(summary), str(screen), str(out))
                row = json.loads((out / "followup_plan.json").read_text())["genes"][0]
                self.assertEqual("needs_light_screen", row["research_action"])

    def test_same_paper_screened_for_one_reference_cannot_close_another(self):
        self.bundle["units"][0]["reference_material_cards"].append({
            "function_description": "A distinct reference context.", "interval_gene_id": "g1",
            "reference_gene": "refX", "reference_species": "Medicago truncatula"})
        self.bundle["units"][0]["literature_link_candidates"] = [{
            "reference_gene": "refX", "reference_species": "Medicago truncatula", "source_id": "PMID:42"}]
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial([self.judgment(), self.judgment(
            pointer="/units/0/reference_material_cards/1", fields=["function_description"])])
        screens = [self.screen(screened_source_ids=["PMID:42"]),
                   self.screen(screen_id="S2", ref="refX", species="Medicago truncatula")]
        screen = self.screening(summary, screens=screens,
                                decisions=[self.assessed(ids=("S1", "S2"))])
        self.run_plan(summary, screen)
        row = self.output()["genes"][0]
        self.assertEqual("needs_light_screen", row["research_action"])
        self.assertEqual(["PMID:42"], row["unscreened_literature_handles"])

    def test_complementary_scopes_for_same_reference_can_share_search(self):
        summary = self.initial()
        screen = self.screening(summary, screens=[
            self.screen(checked_species_scopes=["focal_species"]),
            self.screen(screen_id="S2", checked_species_scopes=["other_plants"])],
            decisions=[self.assessed(ids=("S1", "S2"))])
        self.run_plan(summary, screen)
        self.assertEqual("light_stop_eligible", self.output()["genes"][0]["research_action"])

    def test_missing_frozen_query_status_is_not_completed_retrieval(self):
        self.bundle["units"][0]["literature_query_status"] = [{"reference_gene": "ref1"}]
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf8")
        self.bundle_sha = sha(self.bundle_path.read_bytes())
        summary = self.initial()
        self.run_plan(summary, self.screening(summary))
        row = self.output()["genes"][0]
        self.assertEqual("needs_light_screen", row["research_action"])
        self.assertEqual(["not_assessed"], row["unresolved_literature_query_states"])

    def test_contradictory_no_lead_record_preserves_the_explicit_lead(self):
        summary = self.initial()
        screen = self.screening(summary, screens=[self.screen(upper_lead="A specific mutant phenotype.")])
        self.run_plan(summary, screen)
        self.assertEqual("targeted_followup", self.output()["genes"][0]["research_action"])


if __name__ == "__main__":
    unittest.main()
