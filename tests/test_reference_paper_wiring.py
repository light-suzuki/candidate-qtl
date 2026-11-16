"""Real A resolver/source/sequence/projection/Wet seam, with fixed provider data."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.blast_evidence import run_panel_blastp
from candidate_qtl.blast_panel import freeze_blast_panel, panel_from_verified_bait_panel
from candidate_qtl.family_candidate_bridge import family_candidate_bridge
from candidate_qtl.llm_pathway import hypothesis_verify, bait_panel_freeze
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.dual_route import phenotype_first_status
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.relevance_tiers import NEW_TIERS
from candidate_qtl.source_verify import source_verify
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_blast_evidence import FIXTURE, MockBlastTools, S
from tests.test_reference_gene_resolution import _mechanism_row, _uniprot_record, _write_cache, _write_mechanisms


class ReferencePaperWiringTests(unittest.TestCase):
    def test_external_id_without_a_name_resolves_from_explicit_record_xref(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_mechanisms(root / "input.tsv", [_mechanism_row("h1", reference_gene_name="",
                proposed_stable_id="genbank:ab123456", reference_species="Plant alpha", source_ids="PMID:1234567")])
            record = _uniprot_record("Q00001", "Plant alpha", S, gene_name="A different display name")
            record["uniProtKBCrossReferences"] = [{"database": "EMBL", "id": "AB123456"}]
            _write_cache(root / "records.json", [record])
            result = reference_gene_resolve(str(root / "input.tsv"), str(root / "resolve"), record_cache_arg=str(root / "records.json"))
            rows = pd.read_csv(result["mechanisms_resolved"], sep="\t").fillna("")
            self.assertEqual({"Q00001"}, set(rows.proposed_stable_id))
            self.assertIn("PMID:1234567", rows.iloc[0].source_ids)

    def test_genbank_reference_bibliography_sequence_and_both_interval_copies_reach_wet_t1(self):
        self._reference_bibliography_to_wet("GenBank:AB123456")

    def test_paper_specific_external_id_reaches_wet_without_parent_stable_id(self):
        self._reference_bibliography_to_wet("")

    def _reference_bibliography_to_wet(self, proposed_stable_id):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = {"source_id": "PMID:1234567", "source_title": "A synthetic phenotype observation",
                      "reported_gene": "Old locus", "reported_gene_id": "GenBank:AB123456",
                      "reported_species": "Plant alpha", "reference_gene_id": "GenBank:AB123456",
                      "gene_correspondence": "same_gene", "initial_tier_reason": "Same focal locus under an older reference ID; equivalent observed phenotype.",
                      "paper_checklist": {"focal_species": "yes", "same_gene_as_focal": "yes",
                          "phenotype_reported": "yes", "equivalent_target_phenotype": "yes", "identity_resolved": "yes"}}
            _write_mechanisms(root / "input.tsv", [_mechanism_row("h1", reference_gene_name="NODE1",
                reference_species="Plant alpha", proposed_stable_id=proposed_stable_id, source_ids="PMID:1234567",
                **{PAPER_FIELD: normalize_paper_reports([report])})])
            records = [_uniprot_record(accession, "Plant alpha", S, pmid="9999999", gene_name="NODE1")
                       for accession in ("Q00001", "Q00002")]
            for record, external in zip(records, ("AB123456", "AB654321")):
                record["uniProtKBCrossReferences"] = [{"database": "EMBL", "id": external}]
            _write_cache(root / "references.json", records)
            resolved = reference_gene_resolve(str(root / "input.tsv"), str(root / "resolve"),
                                              record_cache_arg=str(root / "references.json"))
            frame = pd.read_csv(resolved["mechanisms_resolved"], sep="\t", dtype=str).fillna("").set_index("proposed_stable_id")
            if proposed_stable_id:
                self.assertIn("PMID:1234567", frame.loc["Q00001", "source_ids"])
            self.assertNotIn("PMID:1234567", frame.loc["Q00002", "source_ids"])
            original = json.loads(frame.loc["Q00001", PAPER_FIELD])[0]
            self.assertEqual("Q00001", original["reference_gene_id"])
            self.assertEqual("GenBank:AB123456", original["reported_gene_id"])
            bibliography = root / "bibliography.json"
            bibliography.write_text(json.dumps({"records": [{"id": "1234567", "title": report["source_title"]}],
                "gene_records": {"embl:AB123456": {"ids": ["GenBank:AB123456"], "names": ["Old locus"], "species": ["Plant alpha"]}}}))
            sources = source_verify(resolved["mechanisms_resolved"], str(root / "sources"),
                                    cache_arg=str(bibliography), verification_policy=INITIAL_POLICY)
            sourced = pd.read_csv(sources["mechanisms_with_verified_sources"], sep="\t", dtype=str).fillna("").set_index("proposed_stable_id")
            self.assertEqual("resolved", sourced.loc["Q00001", "bibliography_state"])
            self.assertEqual("reference_endpoint_changed", sourced.loc["Q00002", "bibliography_state"])
            self.assertEqual("matched", json.loads(sourced.loc["Q00001", "gene_identity_checks_json"])[0]["state"])
            verified = hypothesis_verify(sources["mechanisms_with_verified_sources"], str(root / "verify"),
                                         focal_species="Plant alpha", record_cache_arg=str(root / "references.json"), min_panel_size=1)
            panel = bait_panel_freeze(sources["mechanisms_with_verified_sources"], verified["validation"],
                                      verified["verified_sequences"], str(root / "panel"), min_panel_size=1,
                                      reviewable_sequences_arg=verified["reviewable_sequences"])
            projection = pd.read_csv(panel["projection_panel"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"Q00001"}, set(projection.bait_id))
            skipped = sourced.loc["Q00002", "hypothesis_id"]
            reviews = pd.read_csv(panel["review_required"], sep="\t", dtype=str).fillna("")
            self.assertIn(("projection_blocked:bibliography_reference_endpoint_changed", skipped),
                          set(zip(reviews.reason, reviews.target)))
            blast_input = root / "blast-input.tsv"
            panel_from_verified_bait_panel(panel["projection_panel"], verified["reviewable_sequences"]).to_csv(blast_input, sep="\t", index=False)
            frozen = freeze_blast_panel(str(blast_input), str(root / "freeze"), min_panel_size=1)
            raw = "".join(f"Q00001\t{gene}\t99\t120\t120\t120\t1\t120\t1\t120\t1e-40\t240\n" for gene in ("g1", "g2"))
            with MockBlastTools(raw_hits=raw):
                search = run_panel_blastp(frozen["manifest"], str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            status = phenotype_first_status(str(root / "status"), hypotheses_arg=sources["mechanisms_with_verified_sources"],
                                             panel_arg=panel["projection_panel"], case_id="reference-paper")
            atlas = root / "atlas.tsv"
            pd.DataFrame({"gene_id": ["g1", "g2", *[f"unused_{i}" for i in range(998)]]}).to_csv(atlas, sep="\t", index=False)
            bridge = family_candidate_bridge(str(root / "bridge"), route_a_hypotheses_arg=status["hypotheses"],
                route_b_atlas_arg=str(atlas), genes_arg=str(atlas), route_a_projection_panel_arg=panel["projection_panel"],
                route_a_blast_arg=search["evidence"], case_id="reference-paper")
            wet = wet_candidate_handoff(str(root / "wet"), route_a_bridge_arg=bridge["bridge_tsv"], route_a_manifest_arg=bridge["manifest"],
                review_candidates_arg=bridge["review_candidates"], review_candidates_manifest_arg=bridge["review_manifest"],
                interval_atlas_arg=str(atlas), wet_budget=20, case_id="reference-paper")
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
            self.assertEqual({"review"}, set(rows.support_state))
