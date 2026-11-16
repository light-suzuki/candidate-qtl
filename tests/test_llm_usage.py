from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.llm_usage import summarize_llm_usage


class LlmUsageTests(unittest.TestCase):
    def _write_raw(self, path: Path, responses: list[dict[str, object]]) -> None:
        with path.open("w", encoding="utf-8") as handle:
            for index, response in enumerate(responses, 1):
                handle.write(json.dumps({
                    "prompt_id": f"p{index}",
                    "gene_id": f"g{index}",
                    "response": response,
                }) + "\n")

    def test_openai_responses_usage_is_aggregated_without_estimation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.jsonl"
            self._write_raw(raw, [
                {"usage": {"input_tokens": 100, "output_tokens": 20, "total_tokens": 120}},
                {"usage": {"input_tokens": 80, "output_tokens": 10, "total_tokens": 90}},
            ])
            manifest = root / "llm.json"
            manifest.write_text(json.dumps({"provider": "openai-responses", "model": "m", "prompt_count": 2}), encoding="utf-8")
            result = summarize_llm_usage(str(raw), str(root / "out"), llm_manifest_arg=str(manifest), case_id="usage")
            self.assertEqual(2, result["call_count"])
            self.assertEqual(180, result["input_tokens"])
            self.assertEqual(30, result["output_tokens"])
            self.assertEqual(210, result["total_tokens"])
            self.assertEqual("complete", result["usage_status"])
            data = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertTrue(data["token_counts_are_provider_reported"])
            self.assertFalse(data["token_estimation_performed"])
            self.assertFalse(data["price_estimation_performed"])

    def test_chat_completions_usage_shape_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.jsonl"
            self._write_raw(raw, [
                {"usage": {"prompt_tokens": 70, "completion_tokens": 9, "total_tokens": 79}},
            ])
            result = summarize_llm_usage(str(raw), str(root / "out"))
            frame = pd.read_csv(result["usage"], sep="\t", dtype=str).fillna("")
            self.assertEqual("chat_completions", frame.iloc[0]["usage_shape"])
            self.assertEqual("70", frame.iloc[0]["input_tokens"])
            self.assertEqual("9", frame.iloc[0]["output_tokens"])

    def test_missing_usage_remains_not_available_instead_of_being_estimated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.jsonl"
            self._write_raw(raw, [{"provider": "lmstudio-cli", "choices": [{"message": {"content": "{}"}}]}])
            result = summarize_llm_usage(str(raw), str(root / "out"))
            self.assertEqual("not_available", result["usage_status"])
            self.assertEqual(0, result["total_tokens"])
            frame = pd.read_csv(result["usage"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_available", frame.iloc[0]["usage_status"])
            self.assertEqual("", frame.iloc[0]["total_tokens"])

    def test_manifest_prompt_count_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.jsonl"
            self._write_raw(raw, [{"usage": {"input_tokens": 10, "output_tokens": 1}}])
            manifest = root / "llm.json"
            manifest.write_text(json.dumps({"provider": "openai-responses", "model": "m", "prompt_count": 2}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "roster mismatch"):
                summarize_llm_usage(str(raw), str(root / "out"), llm_manifest_arg=str(manifest))

    def test_duplicate_prompt_ids_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.jsonl"
            raw.write_text(
                json.dumps({"prompt_id": "p", "gene_id": "g1", "response": {}}) + "\n" +
                json.dumps({"prompt_id": "p", "gene_id": "g2", "response": {}}) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                summarize_llm_usage(str(raw), str(root / "out"))


if __name__ == "__main__":
    unittest.main()
