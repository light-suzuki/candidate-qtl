"""Provider-curated Direct-C bait edges, offline synthetic fixtures.

These pin the contract that lets a symbol-only Direct-C record reach the Wet
review set: the BLAST query must be the one provider entry of the record's
species that cites the record's paper under that name. They are wiring tests,
not scientific accuracy measurements.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.blast_panel import freeze_blast_panel
from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding_producer import (
    prepare_claim_binding_review,
    seal_claim_binding_review,
)
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_pipeline import (
    DIRECT_PROPOSAL_COLUMNS,
    direct_rescue_freeze,
    direct_rescue_project,
)
from candidate_qtl.route_c_provider_bait import (
    EDGE_STATE,
    allelic_representative,
    build_provider_bait,
    load_provider_edges,
    provider_edge_links,
)
from candidate_qtl.route_c_review_candidates import read_direct_c_review_candidates

PMID = "12345678"
DONOR = "Synthetic donor"
SEQUENCE = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEG"
QUOTE = "The SdXYZ1 knockout abolished organ initiation in Synthetic donor."


def _entry(accession: str, *, name: str = "XYZ1", organism: str = DONOR, pmid: str = PMID,
           sequence: str = SEQUENCE, molecule: str = "mRNA", entry_type: str = "UniProtKB unreviewed (TrEMBL)") -> dict:
    return {
        "primaryAccession": accession, "entryType": entry_type,
        "organism": {"scientificName": organism},
        "genes": [{"geneName": {"value": name}}],
        "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": pmid}]},
                        "referencePositions": ["NUCLEOTIDE SEQUENCE"]}],
        "uniProtKBCrossReferences": [{"database": "EMBL", "id": f"X{accession}",
                                      "properties": [{"key": "MoleculeType", "value": molecule}]}],
        "sequence": {"value": sequence, "length": len(sequence)},
    }


class ProviderBaitTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def freeze(self, symbol: str = "SdXYZ1", locus: str = "XYZ1", cultivar: str = "") -> dict:
        proposal = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
        proposal.update(
            record_id="c1", query_phenotype="loss of organ initiation",
            reported_phenotype="loss of organ initiation", focal_species="Synthetic species",
            reference_species=DONOR, reported_gene_symbol=symbol, reported_locus=locus,
            reference_cultivar_or_accession=cultivar,
            species_distance="S5_model_or_other_plant",
            evidence_directness="D0_causal_perturbation_or_rescue",
            evidence_type="knockout", source_id=f"PMID:{PMID}",
        )
        proposals, cache, extraction = (self.root / "proposals.tsv", self.root / "cache.json",
                                         self.root / "extraction.json")
        pd.DataFrame([proposal], columns=DIRECT_PROPOSAL_COLUMNS).to_csv(proposals, sep="\t", index=False)
        cache.write_text(json.dumps([{"id": PMID, "pmid": PMID, "title": "Synthetic source",
                                      "abstractText": QUOTE, "pubYear": "2020"}]), encoding="utf-8")
        extraction.write_text(json.dumps([{"citation": f"PMID:{PMID}", "evidence_relation": "direct_genetic_test",
                                           "quoted_support": QUOTE}]), encoding="utf-8")
        return direct_rescue_freeze(str(proposals), str(self.root / "freeze"), str(cache),
                                    claim_extraction_arg=str(extraction))

    def provider(self, frozen: dict, entries: list[dict], name: str = "bait") -> dict:
        cache = self.root / f"{name}-cache"
        cache.mkdir()
        (cache / f"uniprot_lit_{PMID}.json").write_text(json.dumps({"results": entries}), encoding="utf-8")
        return build_provider_bait(frozen["frozen"], frozen["manifest"], str(self.root / name),
                                   cache_arg=str(cache))

    def search_chain(self, panel: str, target: str = "G1") -> Path:
        """Stand-in for blast-panel-search: the real freeze plus a sealed hits file."""
        frozen_panel = freeze_blast_panel(panel, str(self.root / "panel-freeze"), case_id="case",
                                          min_panel_size=1)
        search = self.root / "search"
        search.mkdir()
        hits = search / "eligible_interval_hits.tsv"
        pd.DataFrame([{"panel_id": "c1", "target_id": target, "eligibility": "eligible", "pident": "80.0",
                       "query_coverage": "0.95", "target_coverage": "0.95"}]).to_csv(hits, sep="\t", index=False)
        (search / "blast_panel_search_manifest.json").write_text(json.dumps({
            "command": "blast-panel-search", "freeze": frozen_panel["manifest"],
            "eligible_interval_hits_sha256": sha256(hits)}), encoding="utf-8")
        return hits

    def test_selects_the_entry_citing_the_paper_under_the_species_prefixed_symbol(self) -> None:
        frozen = self.freeze()
        result = self.provider(frozen, [
            _entry("A00001"),
            _entry("A00002", organism="Other plant"),
            _entry("A00003", pmid="999"),
            _entry("A00004", name="ABC9"),
        ])
        self.assertEqual(1, result["selected_count"])
        edges = pd.read_csv(result["edges"], sep="\t", dtype=str)
        self.assertEqual(["A00001"], list(edges.provider_accession))
        self.assertEqual("SdXYZ1", edges.iloc[0].source_identity_value)
        self.assertEqual("XYZ1", edges.iloc[0].matched_gene_name)
        panel = pd.read_csv(result["panel"], sep="\t", dtype=str)
        self.assertEqual(SEQUENCE, panel.iloc[0].query_sequence)
        self.assertEqual("none", panel.iloc[0].rank_contribution)

    def test_no_qualifying_or_divergent_entries_fail_closed(self) -> None:
        frozen = self.freeze()
        none = self.provider(frozen, [_entry("A00002", organism="Other plant")], name="none")
        divergent = self.provider(frozen, [_entry("A00001"), _entry("A00005", sequence="MKV" * 30)],
                                  name="divergent")
        for result, state in ((none, "no_provider_entry"), (divergent, "ambiguous_provider_entries")):
            self.assertEqual(0, result["selected_count"])
            status = pd.read_csv(result["status"], sep="\t", dtype=str).fillna("")
            self.assertEqual(state, status.iloc[0].state)
            self.assertTrue(pd.read_csv(result["panel"], sep="\t").empty)

    def test_allelic_entries_choose_reviewed_then_mrna_representative(self) -> None:
        mutant = SEQUENCE[:40] + "A" + SEQUENCE[41:]  # one substitution: an allele, not a paralog
        matches = {
            "A00009": (_entry("A00009", molecule="Genomic_DNA"), "XYZ1", PMID),
            "A00010": (_entry("A00010", sequence=mutant), "XYZ1", PMID),
        }
        self.assertEqual("A00010", allelic_representative(matches))
        matches["A00011"] = (_entry("A00011", molecule="Genomic_DNA",
                                    entry_type="UniProtKB reviewed (Swiss-Prot)"), "XYZ1", PMID)
        self.assertEqual("A00011", allelic_representative(matches))
        matches["A00012"] = (_entry("A00012", name="XYZ2"), "XYZ2", PMID)
        self.assertIsNone(allelic_representative(matches))

    def test_missing_cache_without_network_permission_raises(self) -> None:
        frozen = self.freeze()
        with self.assertRaisesRegex(FileNotFoundError, "provider search cache missing"):
            build_provider_bait(frozen["frozen"], frozen["manifest"], str(self.root / "offline"),
                                cache_arg=str(self.root / "empty-cache"))

    def test_edges_are_rechecked_against_provider_bytes_and_frozen_bait(self) -> None:
        frozen = self.freeze()
        result = self.provider(frozen, [_entry("A00001")])
        hits = self.search_chain(result["panel"])
        records = {"c1": pd.read_csv(frozen["frozen"], sep="\t", dtype=str, keep_default_na=False).iloc[0].to_dict()}
        edges_path = Path(result["edges"])
        frame = pd.read_csv(edges_path, sep="\t", dtype=str, keep_default_na=False)
        self.assertEqual(1, len(load_provider_edges(frame, edges_path, records, hits)))
        with self.assertRaisesRegex(ValueError, "require the sealed eligible hits"):
            load_provider_edges(frame, edges_path, records, None)
        other = frame.copy()
        other.loc[0, "query_sequence_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "sequence hash mismatch"):
            load_provider_edges(other, edges_path, records, hits)
        renamed = {"c1": {**records["c1"], "reported_gene_symbol": "SdABC9"}}
        with self.assertRaisesRegex(ValueError, "does not match the frozen record"):
            load_provider_edges(frame, edges_path, renamed, hits)
        hits.write_text(hits.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not the sealed search output"):
            load_provider_edges(frame, edges_path, records, hits)

    def test_link_requires_the_edge_name_for_symbol_and_locus(self) -> None:
        frame = pd.DataFrame([{"record_id": "c1", "source_identity_value": "SdXYZ1", "reference_species": DONOR,
                               "matched_gene_name": "XYZ1", "provider_accession": "A00001"}])
        record = {"record_id": "c1", "reported_gene_symbol": "SdXYZ1", "reported_locus": "XYZ1"}
        hit = {"mapping_method": "blastp", "source_identity_verification_state": EDGE_STATE,
               "source_identity_field": "reported_gene_symbol", "source_identity_value": "SdXYZ1"}
        self.assertTrue(provider_edge_links(frame, record, "gene_symbol", "SdXYZ1", hit))
        self.assertTrue(provider_edge_links(frame, record, "locus", "XYZ1", hit))
        self.assertFalse(provider_edge_links(frame, {**record, "reported_locus": "ABC9"}, "locus", "ABC9", hit))
        self.assertFalse(provider_edge_links(frame, record, "stable_gene_id", "SdXYZ1", hit))
        self.assertFalse(provider_edge_links(
            frame, record, "gene_symbol", "SdXYZ1", {**hit, "source_identity_verification_state": "verified_gene"}))

    def test_bound_symbol_record_reaches_review_only_through_the_provider_edge(self) -> None:
        self.assertEqual((0, 1), self.bound_review_counts("not_applicable"))

    def test_unverified_cultivar_stays_review_but_a_contradiction_excludes(self) -> None:
        self.assertEqual((0, 1), self.bound_review_counts("not_evaluated", cultivar="line-7"))
        self.assertEqual((0, 0), self.bound_review_counts("contradicted", cultivar="line-7"))

    def bound_review_counts(self, cultivar_state: str, cultivar: str = "") -> tuple[int, int]:
        self.root = Path(tempfile.mkdtemp(dir=self.root))
        frozen = self.freeze(cultivar=cultivar)
        self.assertEqual("eligible_direct_record",
                         pd.read_csv(frozen["frozen"], sep="\t", dtype=str).iloc[0].direct_rescue_eligibility)
        result = self.provider(frozen, [_entry("A00001")])
        hits = self.search_chain(result["panel"])
        prepared = prepare_claim_binding_review(frozen["frozen"], frozen["manifest"], str(self.root / "prepare"))
        review = pd.read_csv(prepared["review"], sep="\t", dtype=str, keep_default_na=False)
        review.loc[0, ["reported_identity_type", "reported_identity_value"]] = ["gene_symbol", "SdXYZ1"]
        review.loc[0, "source_organism_state"] = "verified_match"
        review.loc[0, "source_accession_or_cultivar_state"] = cultivar_state
        review.loc[0, "reported_identity_verification_state"] = "verified_exact"
        review.loc[0, "identity_mapping_basis"] = "exact_source_text"
        review.loc[0, "curation_method"] = "independent_agent_review"
        review.loc[0, "reviewed_by"] = "fixture-reviewer"
        review.loc[0, "binding_reason"] = "symbol named in the abstract"
        reviewed = self.root / "reviewed.tsv"
        review.to_csv(reviewed, sep="\t", index=False)
        sealed = seal_claim_binding_review(frozen["frozen"], frozen["manifest"], str(reviewed), str(self.root / "sealed"))
        atlas = self.root / "atlas.tsv"
        atlas_frame = pd.DataFrame([{"gene_id": "G1", "chrom": "chr1", "start": "1", "end": "100"}])
        atlas_frame.to_csv(atlas, sep="\t", index=False)

        def review_count(mapping: str | None, name: str) -> tuple[int, pd.DataFrame]:
            projected = direct_rescue_project(frozen["frozen"], sealed["freeze_manifest"], str(atlas), "q",
                                              str(self.root / f"{name}-projection"), str(hits),
                                              identity_mapping_arg=mapping)
            members = build_route_c_direct_members(frozen["frozen"], sealed["freeze_manifest"],
                                                   projected["projection"], projected["manifest"], str(atlas),
                                                   str(self.root / f"{name}-members"), case_id="q")
            self.assertEqual(0, members["supported_count"])
            candidates, receipt = read_direct_c_review_candidates(
                members["members"], members["manifest"], case_id="q", atlas_frame=atlas_frame, supported_ids=set())
            return receipt["review_record_count"], candidates

        # Without the provider edge the hit's query sequence is unbound to the symbol.
        plain = review_count(None, "plain")[0]
        count, candidates = review_count(result["edges"], "edge")
        if count:
            self.assertEqual("G1", candidates.iloc[0].candidate_gene_id)
            self.assertEqual("unverified_route_review", candidates.iloc[0].support_status)
            self.assertEqual("none", candidates.iloc[0].rank_contribution)
        return plain, count


if __name__ == "__main__":
    unittest.main()
