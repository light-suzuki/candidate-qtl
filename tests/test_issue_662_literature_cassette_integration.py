from __future__ import annotations

import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.literature_cassette import LiteratureRetrieval
from candidate_qtl.reverse_literature_cassette import reverse_literature
from candidate_qtl.route_c_direct_search_cassette import direct_literature_search


class _Response(BytesIO):
    status = 200

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()
        return False


def _response(*_args, **_kwargs):
    return _Response(
        json.dumps(
            {
                "resultList": {
                    "result": [
                        {
                            "pmid": "123",
                            "title": "A source-backed phenotype report",
                            "abstractText": "A mutant allele produced a fasciated stem.",
                            "journalTitle": "Plant Journal",
                            "pubYear": "2026",
                            "authorString": "Example A",
                        }
                    ]
                }
            }
        ).encode("utf-8")
    )


def _retrieval(
    query: str,
    *,
    status: str = "retrieved",
    origin: str = "live_europe_pmc",
    execution: str = "immutable_cassette",
) -> LiteratureRetrieval:
    articles = (
        {
            "pmid": "1",
            "title": "Paper",
            "abstractText": "fasciated stem",
            "journalTitle": "J",
            "pubYear": "2026",
            "authorString": "A",
        },
    ) if status == "retrieved" else ()
    return LiteratureRetrieval(
        query=query,
        request_key_sha256="a" * 64,
        retrieval_status=status,
        articles=articles,
        provider="europe_pmc_rest",
        retrieval_source=origin,
        provider_retrieved_at="2026-09-03T00:00:00+00:00",
        raw_response_sha256="b" * 64,
        cassette_path="/cache/a.json",
        execution_source=execution,
        http_status=200,
        attempt_count=1,
        error="" if status == "retrieved" else "provider down",
    )


class LiteratureCassetteReplayIntegrationTests(unittest.TestCase):
    def test_route_b_live_then_offline_replay_is_byte_stable(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            lookup = root / "lookup.tsv"
            pd.DataFrame(
                [
                    {
                        "gene_id": "Psat.G1",
                        "reference_gene": "CLV1",
                        "reference_function": "meristem receptor",
                    }
                ]
            ).to_csv(lookup, sep="\t", index=False, lineterminator="\n")
            out = root / "route-b"
            with patch(
                "candidate_qtl.literature_cassette.urllib.request.urlopen",
                _response,
            ):
                first = reverse_literature(
                    str(lookup),
                    str(out),
                    phenotype="fasciated stem",
                    allow_external_api=True,
                    case_id="case1",
                )
            scientific = {
                name: Path(first[name]).read_bytes()
                for name in ("literature", "status", "query_provenance")
            }
            second = reverse_literature(
                str(lookup),
                str(out),
                phenotype="fasciated stem",
                allow_external_api=False,
                case_id="case1",
            )
            for name, expected in scientific.items():
                self.assertEqual(expected, Path(second[name]).read_bytes())
            run = json.loads(Path(second["run_envelope"]).read_text())
            self.assertEqual(
                "immutable_cassette",
                run["query_execution_sources"][0]["execution_source"],
            )

    def test_direct_c_live_then_offline_replay_is_byte_stable(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scope = root / "scope.tsv"
            pd.DataFrame(
                [
                    {
                        "query_id": "pea",
                        "reference_species": "Pisum sativum",
                        "species_distance": "S1_same_species",
                        "phenotype_query": "fasciated stem",
                    }
                ]
            ).to_csv(scope, sep="\t", index=False, lineterminator="\n")
            out = root / "route-c"
            with patch(
                "candidate_qtl.literature_cassette.urllib.request.urlopen",
                _response,
            ):
                first = direct_literature_search(
                    "fasciated stem",
                    "Pisum sativum",
                    str(scope),
                    str(out),
                    allow_external_api=True,
                    case_id="case1",
                )
            scientific = {
                name: Path(first[name]).read_bytes()
                for name in ("candidates", "status", "query_provenance")
            }
            second = direct_literature_search(
                "fasciated stem",
                "Pisum sativum",
                str(scope),
                str(out),
                allow_external_api=False,
                case_id="case1",
            )
            for name, expected in scientific.items():
                self.assertEqual(expected, Path(second[name]).read_bytes())
            run = json.loads(Path(second["run_envelope"]).read_text())
            self.assertEqual(
                "immutable_cassette",
                run["query_execution_sources"][0]["execution_source"],
            )

    def test_route_b_provider_failure_is_per_query_and_does_not_abort_batch(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            lookup = root / "lookup.tsv"
            pd.DataFrame(
                [
                    {
                        "gene_id": "G1",
                        "reference_gene": "CLV1",
                        "reference_function": "receptor",
                    },
                    {
                        "gene_id": "G2",
                        "reference_gene": "WUS",
                        "reference_function": "transcription factor",
                    },
                ]
            ).to_csv(lookup, sep="\t", index=False, lineterminator="\n")

            def retrieve(query, *_args, **_kwargs):
                return (
                    _retrieval(query)
                    if "CLV1" in query
                    else _retrieval(query, status="provider_unavailable")
                )

            with patch(
                "candidate_qtl.reverse_literature_cassette.retrieve_europe_pmc",
                retrieve,
            ):
                result = reverse_literature(
                    str(lookup),
                    str(root / "out"),
                    phenotype="fasciated stem",
                    allow_external_api=True,
                )
            status = pd.read_csv(
                result["status"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(
                {"retrieved", "not_evaluated_provider_unavailable"},
                set(status["query_status"]),
            )
            literature = pd.read_csv(
                result["literature"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(1, len(literature))
            self.assertEqual(
                "europe_pmc_rest:live_europe_pmc",
                literature.loc[0, "provider"],
            )


if __name__ == "__main__":
    unittest.main()
