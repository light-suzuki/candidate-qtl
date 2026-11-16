from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import tempfile
import unittest
from unittest import mock

from candidate_qtl import llm_pathway_legacy as pathway
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.reference_gene_resolution_legacy import _cached_gene_symbol_candidates_exact


def _record(accession: str, pmid: str | None, sequence: str = "M" * 12, gene: str = "DUP1") -> dict:
    value = {
        "primaryAccession": accession,
        "organism": {"scientificName": "Arabidopsis thaliana"},
        "genes": [{"geneName": {"value": gene}}],
        "sequence": {"value": sequence},
    }
    if pmid:
        value["references"] = [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": pmid}]}}]
    return value


def _mechanism(path: Path) -> None:
    row = {column: "" for column in MECHANISM_COLUMNS}
    row.update({"hypothesis_id": "mech_001", "case_id": "fanout", "reference_species": "Arabidopsis thaliana",
                "reference_gene_name": "DUP1", "gene_family": "test", "source_ids": ""})
    pd.DataFrame([row])[MECHANISM_COLUMNS].to_csv(path, sep="\t", index=False)


class ReferenceGeneFanoutTests(unittest.TestCase):
    def setUp(self) -> None:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp_path = Path(tmp.name)

    def test_ambiguous_gene_fanout_preserves_lineage_and_per_accession_sources(self) -> None:
        tmp_path = self.tmp_path
        mechanisms = tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        cache = tmp_path / "records.json"
        cache.write_text(json.dumps({"results": [_record("Q9XYZ1", "111"), _record("Q9XYZ2", None)]}), encoding="utf-8")
        out = tmp_path / "out"

        result = reference_gene_resolve(str(mechanisms), str(out), record_cache_arg=str(cache))
        resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
        children = resolved[resolved["parent_hypothesis_id"] == "mech_001"].sort_values("hypothesis_id")
        assert list(children["hypothesis_id"]) == ["mech_001__ref_q9xyz1", "mech_001__ref_q9xyz2"]
        assert list(children["proposed_stable_id"]) == ["Q9XYZ1", "Q9XYZ2"]
        assert set(children["source_ids"]) == {"", "PMID:111"}
        assert all("parent_hypothesis_id" in value for value in children["retrieval_provenance"])
        assert result["provider_resolved_count"] == 1
        ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
        assert set(ledger["hypothesis_id"]) == {"mech_001__ref_q9xyz1", "mech_001__ref_q9xyz2"}


    def test_candidate_matching_is_casefolded_but_species_remains_exact(self) -> None:
        payload = {"results": [_record("Q9XYZ1", "111", gene="FAN1"), _record("Q9XYZ2", "111", gene="FAN1")]}
        wrong_species = _record("Q9XYZ3", "111", gene="FAN1")
        wrong_species["organism"]["scientificName"] = "Arabidopsis lyrata"
        payload["results"].append(wrong_species)
        response = mock.MagicMock()
        response.read.return_value = json.dumps(payload).encode()
        response.headers = {}
        response.__enter__.return_value = response
        with mock.patch.object(pathway.urllib.request, "urlopen", return_value=response):
            candidates, _, _, _ = pathway._resolve_accession_candidates("fan1", "arabidopsis  thaliana")
        assert candidates == ["Q9XYZ1", "Q9XYZ2"]
        cache = {record["primaryAccession"]: record for record in reversed(payload["results"])}
        assert _cached_gene_symbol_candidates_exact(cache, "fan1", "arabidopsis  thaliana") == candidates


    def test_empty_candidate_query_has_status_without_network(self) -> None:
        with mock.patch.object(pathway.urllib.request, "urlopen", side_effect=AssertionError("unexpected network")):
            assert pathway._resolve_accession_candidates("", "Arabidopsis thaliana") == ([], "no_gene_name_to_resolve", "", 0)
