"""Documentation guards for recall-oriented Wet budget use (#732).

These are wording/ownership checks, not runtime or biological performance tests.
The production fill algorithm is intentionally unchanged by this contract PR.
"""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
POLICY = "docs/CANONICAL-WET-MEMBERSHIP.ja.md"


class WetRecallBudgetContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = (ROOT / POLICY).read_text(encoding="utf-8")
        cls.core = (ROOT / "CORE-README.ja.md").read_text(encoding="utf-8")
        cls.agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")

    def test_count_is_neither_quota_nor_truncation(self):
        sentence = "20件は目安であり、最低件数・切断上限・調査完了の条件ではありません。20件を超える既存候補は数字だけで切らない。"
        for content in (self.policy, self.core):
            self.assertIn(sentence, content)
        self.assertIn("Twenty is a guide, not a minimum count, truncation cap, or completion condition", self.agents)
        self.assertIn("never truncate retained candidates merely to reach 20", self.agents)

    def test_weak_fill_does_not_require_a_strong_mechanistic_claim(self):
        self.assertIn("補充に本命と同じ直接証拠・claim確認・具体的な機構説明を要求しません", self.policy)
        self.assertIn("rank_contribution=none", self.policy)
        self.assertIn("支持フラグや架空の引用は付けません", self.policy)
        self.assertIn("supported_only", self.policy)
        self.assertIn("診断・比較用", self.policy)

    def test_policy_owner_is_shared_instead_of_a_new_workflow(self):
        self.assertIn(POLICY, self.core)
        self.assertIn(POLICY, self.agents)
        self.assertIn("固定比率やfamilyごとの枠", self.policy)
        self.assertIn("新スコアを設ける指示ではありません", self.policy)

    def test_eligibility_and_unfinished_work_remain_visible(self):
        for token in ("tier4_dark", "unfilled_budget_slots", "gene ID順", "#722", "#732"):
            self.assertIn(token, self.policy)
        self.assertIn("全区間に20遺伝子あるだけでは常に20件になるとは限りません", self.policy)
        self.assertIn("20件を埋めただけで未知の正解保持を保証できるわけではありません", self.policy)
        self.assertIn("#732のreview込みconsumer", self.policy)

    def test_machine_checks_are_not_biological_proof(self):
        for phrase in ("研究者役LLM / agent", "生物学的な証明と呼ばない", "支持状態", "採用/保留"):
            self.assertIn(phrase, self.policy)
        self.assertIn("明示的な保留を通常fillで黙って復活させません", self.policy)

    def test_second_stage_trigger_and_benchmark_pause_are_preserved(self):
        self.assertIn("`<=25` can complete without a boundary request, including 8 or 12 eligible genes", self.agents)
        self.assertIn("not a truncation cap", self.agents)
        self.assertIn("boundary review", self.policy)
        self.assertIn("round budget exhausted", self.policy)


if __name__ == "__main__":
    unittest.main()
