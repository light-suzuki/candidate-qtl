import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.claim_review_session import _validate_empty_c_research, _validate_route_c_execution, sha256


class EmptyCResearchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.result = self.root / "query.json"
        self.result.write_text('{"query":"phenotype gene","results":[]}', encoding="utf-8")
        self.receipt = {"contract": "route-c-native-research-receipt-1.0.0", "case_id": "case",
            "worker_id": "c-worker", "model": "recorded-model", "interval_blind": True,
            "state": "completed", "finding_count": 0, "no_findings_rationale": "No relevant direct report found",
            "queries": [{"query": "phenotype gene", "status": "no_hit", "result_path": str(self.result),
                         "result_sha256": sha256(self.result)}], "read_sources": []}

    def context(self):
        path = self.root / "receipt.json"
        path.write_text(json.dumps(self.receipt), encoding="utf-8")
        return {"research_receipt_path": str(path), "research_receipt_sha256": sha256(path)}

    def test_native_no_hit_is_allowed_with_actual_sealed_payload(self):
        _validate_empty_c_research(self.context(), "case")

    def test_existing_search_generation_requires_actual_cassette_payload(self):
        import pandas as pd
        from io import BytesIO
        from unittest.mock import patch
        from candidate_qtl.route_c_direct_search_cassette import direct_literature_search

        class Response(BytesIO):
            status = 200

        scope = self.root / "scope.tsv"
        pd.DataFrame([{"query_id": "q1", "reference_species": "Pisum sativum",
                       "species_distance": "S1_same_species", "phenotype_query": "test phenotype"}]).to_csv(
                           scope, sep="\t", index=False)
        with patch("candidate_qtl.literature_cassette.urllib.request.urlopen",
                   side_effect=lambda *_a, **_k: Response(b'{"resultList":{"result":[]}}')):
            search = direct_literature_search("test phenotype", "Pisum sativum", str(scope),
                str(self.root / "search"), allow_external_api=True, case_id="case")
        context = {"search_manifest_path": search["manifest"],
                   "search_manifest_sha256": sha256(search["manifest"]), "reason": "No relevant reported fact found"}
        _validate_empty_c_research(context, "case")
        manifest = json.loads(Path(search["manifest"]).read_text())
        cassette = next(Path(manifest["outputs"]["cassette_directory"]).glob("*.json"))
        cassette.write_text("{}", encoding="utf-8")
        with self.assertRaises(ValueError):
            _validate_empty_c_research(context, "case")

    def test_completed_empty_requires_sealed_status_and_empty_members(self):
        import pandas as pd
        from tests._issue664_helpers import manifest_for
        members = self.root / "members.tsv"
        pd.DataFrame(columns=["candidate_gene_id", "support_status"]).to_csv(members, sep="\t", index=False)
        producer = manifest_for(members, route="c", case_id="case")
        status = {"state": "completed_empty", "case_id": "case", **self.context()}
        path = self.root / "status.json"
        path.write_text(json.dumps(status), encoding="utf-8")
        status.update(status_path=str(path), status_sha256=sha256(path))
        actual = _validate_route_c_execution(status, "case", str(members), str(producer))
        self.assertEqual("completed_empty", actual["state"])
        pd.DataFrame([{"candidate_gene_id": "g1", "support_status": "review"}]).to_csv(members, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "member row count"):
            _validate_route_c_execution(status, "case", str(members), str(producer))

    def test_retrieved_requires_saved_read_source_and_can_be_judged_irrelevant(self):
        self.receipt["queries"][0]["status"] = "retrieved"
        with self.assertRaisesRegex(ValueError, "read-source"):
            _validate_empty_c_research(self.context(), "case")
        self.receipt["read_sources"] = [{"path": str(self.result), "sha256": sha256(self.result)}]
        _validate_empty_c_research(self.context(), "case")

    def test_failed_tampered_missing_and_unattested_search_are_rejected(self):
        for field, replacement in (("queries", []), ("model", ""), ("interval_blind", False)):
            original = self.receipt[field]
            self.receipt[field] = replacement
            with self.assertRaises(ValueError):
                _validate_empty_c_research(self.context(), "case")
            self.receipt[field] = original
        self.receipt["queries"][0]["status"] = "provider_unavailable"
        with self.assertRaises(ValueError):
            _validate_empty_c_research(self.context(), "case")
        self.receipt["queries"][0]["status"] = "no_hit"
        context = self.context()
        self.result.write_text("tampered", encoding="utf-8")
        with self.assertRaises(ValueError):
            _validate_empty_c_research(context, "case")
        with self.assertRaises(ValueError):
            _validate_empty_c_research({}, "case")
        self.result.write_text("{}", encoding="utf-8")
        self.receipt["queries"][0]["result_sha256"] = sha256(self.result)
        with self.assertRaisesRegex(ValueError, "placeholder"):
            _validate_empty_c_research(self.context(), "case")
