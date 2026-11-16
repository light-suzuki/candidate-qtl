from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class LlmUsageCliTests(unittest.TestCase):
    def test_cli_aggregates_provider_reported_usage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "llm_raw_responses.jsonl"
            raw.write_text(
                json.dumps({
                    "prompt_id": "p1",
                    "gene_id": "g1",
                    "response": {"usage": {"input_tokens": 25, "output_tokens": 5, "total_tokens": 30}},
                }) + "\n",
                encoding="utf-8",
            )
            llm_manifest = root / "llm_run_manifest.json"
            llm_manifest.write_text(json.dumps({
                "provider": "openai-responses",
                "model": "fixture-model",
                "prompt_count": 1,
            }) + "\n", encoding="utf-8")

            self.assertEqual(0, main([
                "llm-usage",
                "--raw-responses", str(raw),
                "--llm-manifest", str(llm_manifest),
                "--out-dir", str(root / "usage"),
                "--case-id", "cli-usage",
            ]))
            manifest = json.loads((root / "usage" / "llm_usage_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("cli-usage", manifest["case_id"])
            self.assertEqual(1, manifest["call_count"])
            self.assertEqual(25, manifest["input_tokens"])
            self.assertEqual(5, manifest["output_tokens"])
            self.assertEqual(30, manifest["total_tokens"])
            self.assertEqual("complete", manifest["usage_status"])
            self.assertFalse(manifest["token_estimation_performed"])
            self.assertFalse(manifest["price_estimation_performed"])


if __name__ == "__main__":
    unittest.main()
