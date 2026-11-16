"""#688: canonical source text remains reviewable without changing support."""

import hashlib
import json
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.source_review_capture import capture_source, content_sha256
from candidate_qtl.source_verify import source_verify


def _preserves_exact_sections_and_separates_available_from_classified_scope():
    source = {
        "id": "12345", "title": "  Source title\n", "abstractText": "<p>A  claim.</p>",
        "fullText": "Body\n  with spacing.\n", "keywordList": {"keyword": ["test"]},
    }
    result = capture_source(source, citation="PMID:12345", citation_type="pmid",
                            resolved_id="12345", classification_scope="title_abstract")
    assert result["sections"]["fullText"] == source["fullText"]
    assert result["sections"]["abstractText"] == source["abstractText"]
    assert result["sections"]["keywordList"] == source["keywordList"]
    assert result["available_text_scope"] == "full_text_cached"
    assert result["classification_scope"] == "title_abstract"
    assert result["sections_sha256"] == content_sha256(result["sections"])
    assert result["semantic_binding_state"] == "not_evaluated"


def _inputs(tmp_path, reverse=False):
    records = [
        {"hypothesis_id": "h1", "reference_gene_name": "ABC1", "source_ids": "PMID:12345"},
        {"hypothesis_id": "h2", "reference_gene_name": "ZZZ2", "source_ids": "PMID:12345"},
        {"hypothesis_id": "h3", "reference_gene_name": "ABC1", "source_ids": "PMID:98765"},
    ]
    path = tmp_path / ("reversed.tsv" if reverse else "mechanisms.tsv")
    pd.DataFrame(records[::-1] if reverse else records).to_csv(path, sep="\t", index=False)
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({
        "id": "12345", "pmid": "12345", "title": "ABC1 function",
        "abstractText": "An ABC1 mutant has altered organ development.\n",
    }), encoding="utf-8")
    return path, cache


def _shared_paper_keeps_supported_and_unmatched_hypotheses(tmp_path):
    mechanism, cache = _inputs(tmp_path)
    result = source_verify(str(mechanism), str(tmp_path / "out"), cache_arg=str(cache))
    capture = json.loads(Path(result["source_review_records"]).read_text(encoding="utf-8"))
    assert capture["membership_effect"] == "none"
    assert list(capture["records"]) == ["PMID:12345"]
    row = capture["records"]["PMID:12345"]
    assert row["hypothesis_ids"] == ["h1", "h2"]
    assert row["source_record_id"] == "12345"
    validation = pd.read_csv(result["source_validation"], sep="\t").set_index("hypothesis_id")
    assert validation.loc["h1", "verification_status"] == "verified"
    assert validation.loc["h2", "verification_status"] == "review_required"
    assert validation.loc["h3", "review_reason"] == "verification_not_performed_offline"
    mechanism2, _ = _inputs(tmp_path, reverse=True)
    result2 = source_verify(str(mechanism2), str(tmp_path / "out2"), cache_arg=str(cache))
    assert Path(result["source_review_records"]).read_bytes() == Path(result2["source_review_records"]).read_bytes()


def _cli_exposes_fingerprinted_source_capture(tmp_path):
    mechanism, cache = _inputs(tmp_path)
    out = tmp_path / "cli"
    assert main(["source-verify", "--verification-policy", "content_v1", "--mechanisms", str(mechanism), "--cache", str(cache),
                 "--out-dir", str(out), "--json"]) == 0
    manifest = json.loads((out / "source_verification_manifest.json").read_text())
    artifact = out / "source_review_records.json"
    assert manifest["outputs"]["source_review_records"] == str(artifact)
    assert manifest["source_review_records_fingerprint"]["sha256"] == hashlib.sha256(artifact.read_bytes()).hexdigest()


def _missing_title_only_and_nontext_values_never_gain_scope():
    for source, expected in [({}, "none"), ({"title": "Only a title"}, "title_only"),
                             ({"title": "Title", "abstractText": None, "fullText": []}, "title_only")]:
        row = capture_source(source, citation="PMID:12345", citation_type="pmid",
                             resolved_id="12345", classification_scope=expected)
        assert row["available_text_scope"] == expected
        assert row["semantic_binding_state"] == "not_evaluated"


def _offline_without_cache_emits_empty_explicit_review_artifact(tmp_path):
    mechanism, _ = _inputs(tmp_path)
    result = source_verify(str(mechanism), str(tmp_path / "offline"))
    capture = json.loads(Path(result["source_review_records"]).read_text(encoding="utf-8"))
    assert capture["records"] == {}
    assert capture["membership_effect"] == "none"


class SourceReviewCaptureTests(unittest.TestCase):
    def test_preserves_exact_sections(self):
        _preserves_exact_sections_and_separates_available_from_classified_scope()

    def test_shared_paper_associations_and_order_invariance(self):
        with tempfile.TemporaryDirectory() as directory:
            _shared_paper_keeps_supported_and_unmatched_hypotheses(Path(directory))

    def test_cli_capture_reachability(self):
        with tempfile.TemporaryDirectory() as directory:
            _cli_exposes_fingerprinted_source_capture(Path(directory))

    def test_missing_scope(self):
        _missing_title_only_and_nontext_values_never_gain_scope()

    def test_offline_empty_capture(self):
        with tempfile.TemporaryDirectory() as directory:
            _offline_without_cache_emits_empty_explicit_review_artifact(Path(directory))


class LazySourceCaptureTests(unittest.TestCase):
    def test_shared_classified_paper_captured_once_without_losing_bindings(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mechanisms, cache = _inputs(root)
            with mock.patch("candidate_qtl.source_verify.capture_source", wraps=capture_source) as capture:
                result = source_verify(str(mechanisms), str(root / "out"), cache_arg=str(cache))
            self.assertEqual(capture.call_count, 1)
            records = json.loads(Path(result["source_review_records"]).read_text(encoding="utf-8"))["records"]
            source = json.loads(cache.read_text(encoding="utf-8"))
            expected = capture_source(source, citation="PMID:12345", citation_type="pmid",
                                      resolved_id="12345", classification_scope=records["PMID:12345"]["classification_scope"])
            expected["hypothesis_ids"] = ["h1", "h2"]
            self.assertEqual(records, {"PMID:12345": expected})
            validation = pd.read_csv(result["source_validation"], sep="\t").set_index("hypothesis_id")
            self.assertEqual(validation.loc["h1", "verification_status"], "verified")
            self.assertEqual(validation.loc["h2", "verification_status"], "review_required")
            self.assertEqual(validation.loc["h3", "review_reason"], "verification_not_performed_offline")

    def test_shared_review_only_paper_captured_once_keeps_lineage_and_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mechanisms, cache = _inputs(root)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame["review_only_source_ids"] = ["PMID:12345", "PMID:12345", ""]
            frame.loc[:1, "source_ids"] = ""
            frame.to_csv(mechanisms, sep="\t", index=False)
            with mock.patch("candidate_qtl.source_verify.capture_source", wraps=capture_source) as capture:
                result = source_verify(str(mechanisms), str(root / "out"), cache_arg=str(cache))
            self.assertEqual(capture.call_count, 1)
            payload = json.loads(Path(result["source_review_records"]).read_text(encoding="utf-8"))
            record = payload["records"]["PMID:12345"]
            self.assertEqual(record["hypothesis_ids"], ["h1", "h2"])
            self.assertEqual(record["review_only_hypothesis_ids"], ["h1", "h2"])
            self.assertEqual(record["review_only_lineage"], [
                {"hypothesis_id": "h1", "parent_hypothesis_id": "h1"},
                {"hypothesis_id": "h2", "parent_hypothesis_id": "h2"}])
            self.assertEqual(record["classification_scope"], "review_only_not_classified")
            source = json.loads(cache.read_text(encoding="utf-8"))
            expected = capture_source(source, citation="PMID:12345", citation_type="pmid",
                                      resolved_id="12345", classification_scope="review_only_not_classified")
            for key in ("sections", "sections_sha256", "provider_record_sha256", "source_identifiers"):
                self.assertEqual(record[key], expected[key])
            self.assertEqual(payload["membership_effect"], "none")

    def test_review_only_reuse_does_not_replace_classified_first_capture(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mechanisms, cache = _inputs(root)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame["review_only_source_ids"] = ["", "PMID:12345", ""]
            frame.to_csv(mechanisms, sep="\t", index=False)
            with mock.patch("candidate_qtl.source_verify.capture_source", wraps=capture_source) as capture:
                result = source_verify(str(mechanisms), str(root / "out"), cache_arg=str(cache))
            self.assertEqual(capture.call_count, 1)
            record = json.loads(Path(result["source_review_records"]).read_text(encoding="utf-8"))["records"]["PMID:12345"]
            self.assertNotEqual(record["classification_scope"], "review_only_not_classified")
            self.assertEqual(record["hypothesis_ids"], ["h1", "h2"])
            self.assertEqual(record["review_only_hypothesis_ids"], ["h2"])
