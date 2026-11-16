"""Offline regressions for UniProt input/provider error taxonomy (#655)."""

from __future__ import annotations

import json
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.llm_pathway import hypothesis_verify
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.uniprot_error_policy import (
    capture_uniprot_entry_events,
    fetch_uniprot_entry,
)


class UniProtEntryPolicyTests(unittest.TestCase):
    def _http_error(self, code: int) -> urllib.error.HTTPError:
        return urllib.error.HTTPError(
            "https://rest.uniprot.org/uniprotkb/P12345.json",
            code,
            f"HTTP {code}",
            None,
            None,
        )

    def test_invalid_syntax_is_rejected_before_network(self) -> None:
        with patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen") as urlopen:
            with capture_uniprot_entry_events() as events:
                self.assertIsNone(fetch_uniprot_entry("Pisum_sativum__NIN"))
        urlopen.assert_not_called()
        self.assertEqual("invalid_syntax", events[0].state)
        self.assertFalse(events[0].network_attempted)

    def test_http_400_is_input_error_and_404_is_definitive_not_found(self) -> None:
        for code, expected in ((400, "invalid_http_400"), (404, "not_found_http_404")):
            with self.subTest(code=code):
                with patch(
                    "candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                    side_effect=self._http_error(code),
                ):
                    with capture_uniprot_entry_events() as events:
                        self.assertIsNone(fetch_uniprot_entry("P12345"))
                self.assertEqual(expected, events[0].state)
                self.assertTrue(events[0].network_attempted)

    def test_transient_and_provider_errors_raise_not_evaluated_outcomes(self) -> None:
        cases = [
            (self._http_error(429), "rate_limited"),
            (self._http_error(500), "provider_unavailable"),
            (TimeoutError("timed out"), "provider_unavailable"),
        ]
        for error, expected in cases:
            with self.subTest(expected=expected, error=type(error).__name__):
                with patch(
                    "candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                    side_effect=error,
                ):
                    with capture_uniprot_entry_events() as events:
                        with self.assertRaises(RuntimeError):
                            fetch_uniprot_entry("P12345")
                self.assertEqual(expected, events[0].state)


class HypothesisVerifyErrorTaxonomyTests(unittest.TestCase):
    def test_http_400_is_not_reported_as_provider_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame(
                [
                    {
                        "hypothesis_id": "hyp_001",
                        "case_id": "case-655",
                        "gene_family": "NIN",
                        "reference_gene_name": "NIN",
                        "reference_species": "Medicago truncatula",
                        "proposed_stable_id": "P12345",
                    }
                ]
            ).to_csv(hypotheses, sep="\t", index=False, lineterminator="\n")
            error = urllib.error.HTTPError(
                "https://rest.uniprot.org/uniprotkb/P12345.json",
                400,
                "bad request",
                None,
                None,
            )
            with patch(
                "candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                side_effect=error,
            ):
                result = hypothesis_verify(
                    str(hypotheses),
                    str(root / "out"),
                    "Pisum sativum",
                    allow_external_api=True,
                    min_panel_size=1,
                )

            validation = pd.read_csv(result["validation"], sep="\t", dtype=str).fillna("")
            row = validation.iloc[0]
            self.assertEqual("rejected", row["verification_status"])
            self.assertEqual("invalid_submitted_accession", row["rejection_reason"])
            self.assertEqual("invalid_http_400", row["accession_resolution_status"])
            self.assertNotEqual("accession_provider_unavailable", row["review_reason"])

            reviews = pd.read_csv(
                result["review_required"], sep="\t", dtype=str
            ).fillna("")
            target_reviews = reviews[reviews["target"].eq("hyp_001")]
            self.assertEqual(1, len(target_reviews))
            self.assertEqual(
                "proposal_rejected:invalid_submitted_accession",
                target_reviews.iloc[0]["reason"],
            )
            self.assertFalse(
                target_reviews.astype(str)
                .apply(lambda column: column.str.contains("provider_unavailable"))
                .any()
                .any()
            )

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["uniprot_entry_outcome_counts"]["invalid_http_400"])
            self.assertEqual(["hyp_001"], manifest["invalid_submitted_accession_hypothesis_ids"])


class ReferenceGeneResolveErrorTaxonomyTests(unittest.TestCase):
    def test_http_400_is_evaluated_input_error_and_invalid_id_is_cleared(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            row = {column: "" for column in MECHANISM_COLUMNS}
            row.update(
                {
                    "hypothesis_id": "hyp_001",
                    "case_id": "case-655",
                    "reference_gene_name": "",
                    "reference_species": "Medicago truncatula",
                    "proposed_stable_id": "P12345",
                    "source_ids": "",
                }
            )
            pd.DataFrame([row], columns=MECHANISM_COLUMNS).to_csv(
                mechanisms, sep="\t", index=False, lineterminator="\n"
            )
            error = urllib.error.HTTPError(
                "https://rest.uniprot.org/uniprotkb/P12345.json",
                400,
                "bad request",
                None,
                None,
            )
            with patch(
                "candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                side_effect=error,
            ):
                result = reference_gene_resolve(
                    str(mechanisms),
                    str(root / "out"),
                    allow_external_api=True,
                    focal_species="Pisum sativum",
                )

            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            ledger_row = ledger.iloc[0]
            self.assertEqual("evaluated", ledger_row["evaluation_status"])
            self.assertEqual("invalid_submitted_accession", ledger_row["review_reason"])
            self.assertEqual("invalid_http_400", ledger_row["accession_resolution_status"])
            self.assertEqual("", ledger_row["provider_error"])
            resolved = pd.read_csv(
                result["mechanisms_resolved"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual("", resolved.iloc[0]["proposed_stable_id"])

            reviews = pd.read_csv(
                result["review_required"], sep="\t", dtype=str
            ).fillna("")
            target_reviews = reviews[reviews["target"].eq("hyp_001")]
            self.assertEqual(1, len(target_reviews))
            self.assertEqual(
                "reference_gene_resolution_required:invalid_submitted_accession",
                target_reviews.iloc[0]["reason"],
            )
            self.assertFalse(
                target_reviews.astype(str)
                .apply(lambda column: column.str.contains("provider_unavailable"))
                .any()
                .any()
            )

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(["hyp_001"], manifest["invalid_submitted_accession_hypothesis_ids"])
            self.assertEqual(0, manifest["evaluation_status_counts"]["not_evaluated"])


if __name__ == "__main__":
    unittest.main()
