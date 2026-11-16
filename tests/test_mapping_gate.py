import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.mapping_gate import mapping_adequacy_gate


class MappingGateTests(unittest.TestCase):
    def _run(self, root: Path, values: dict[str, object]):
        pd.DataFrame([
            {"metric": metric, "value": value, "source": "PMID:1"}
            for metric, value in values.items()
        ]).to_csv(root / "metrics.tsv", sep="\t", index=False)
        return mapping_adequacy_gate(str(root / "metrics.tsv"), str(root / "out"))

    def test_ready_and_probability_guard(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = self._run(Path(tmp), {
                "interval_size_bp": 5_000_000, "interval_size_cm": 5,
                "marker_count": 20, "informative_recombinant_count": 8,
                "lod_peak_count": 1, "peak_stability": "true",
                "phenotype_direction_consistent": "true", "candidate_gene_count": 80,
            })
            self.assertTrue(result["candidate_triage_allowed"])
            gate = json.loads(Path(result["gate"]).read_text(encoding="utf-8"))
            self.assertFalse(gate["causal_probability"])
            self.assertFalse(gate["ranking_score_modified"])

    def test_linked_low_density_and_recombinant_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = self._run(Path(tmp), {
                "interval_size_bp": 30_000_000, "interval_size_cm": 2,
                "marker_count": 2, "informative_recombinant_count": 1,
                "lod_peak_count": 2, "peak_stability": "true",
                "phenotype_direction_consistent": "true", "candidate_gene_count": 300,
            })
            gate = json.loads(Path(result["gate"]).read_text(encoding="utf-8"))
            self.assertEqual("linked_qtl_suspected", gate["primary_status"])
            self.assertIn("marker_density_insufficient", gate["statuses"])
            self.assertIn("insufficient_informative_recombinants", gate["statuses"])
            self.assertIn("low_recombination_region", gate["statuses"])
            self.assertFalse(gate["candidate_triage_allowed"])

    def test_missing_input_is_not_scientific_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = self._run(Path(tmp), {"interval_size_bp": 1_000_000})
            gate = json.loads(Path(result["gate"]).read_text(encoding="utf-8"))
            self.assertEqual("mapping_model_review_required", gate["primary_status"])
            self.assertIn("marker_count", gate["missing_required_inputs"])


if __name__ == "__main__":
    unittest.main()
