import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.citation_tokens import split_citations
from candidate_qtl.paper_reports import source_tokens
from candidate_qtl.route_c_claim_binding_producer import prepare_claim_binding_review, seal_claim_binding_review
from candidate_qtl.route_c_direct_pipeline import DIRECT_PROPOSAL_COLUMNS, direct_rescue_freeze


class CitationTokenTests(unittest.TestCase):
    def test_internal_doi_punctuation_is_not_a_list_separator(self):
        for doi in ("10.5555/neutral;a,b", "DOI:10.5555/(neutral)7<8::part>;2-3",
                    "https://doi.org/10.5555/neutral;x,y"):
            with self.subTest(doi=doi):
                self.assertEqual([doi], split_citations(doi))
                self.assertEqual(1, len(source_tokens(doi)))

    def test_explicit_multiple_citations_and_legacy_non_doi_lists(self):
        expected = ["PMID:12345678", "10.5555/neutral;a,b", "PMCID:PMC1234567", "10.5555/second;x"]
        self.assertEqual(expected, split_citations("PMID:12345678;10.5555/neutral;a,b, PMCID:PMC1234567;10.5555/second;x"))
        self.assertEqual(["12345678", "23456789", "broken"], split_citations("12345678,23456789|broken"))
        self.assertEqual(expected, split_citations([expected[:2], expected[2:]]))
        self.assertEqual(["10.5555/neutral;12345678"], split_citations("10.5555/neutral;12345678"))

    def test_semicolon_doi_survives_c_freeze_prepare_seal_and_tamper_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            doi = "10.5555/(neutral)7<8::part>;2-3"
            quote = "NeutralGene mutation caused altered growth in Synthetic species."
            proposal = {key: "" for key in DIRECT_PROPOSAL_COLUMNS}
            proposal.update(record_id="neutral-record", query_phenotype="altered growth", reported_phenotype="altered growth",
                focal_species="Synthetic species", reference_species="Synthetic species", reported_gene_id="NeutralGene",
                species_distance="S1_same_species", evidence_directness="D0_causal_perturbation_or_rescue",
                evidence_type="synthetic test", source_id=doi)
            proposals, cache, extraction = root/"proposals.tsv", root/"cache.json", root/"extraction.json"
            pd.DataFrame([proposal]).to_csv(proposals, sep="\t", index=False)
            cache.write_text(json.dumps([{"doi":doi,"title":"Neutral source","abstractText":quote}]))
            extraction.write_text(json.dumps([{"hypothesis_id":"neutral-record","citation":doi,
                "evidence_relation":"direct_genetic_test","quoted_support":quote}]))
            frozen = direct_rescue_freeze(str(proposals),str(root/"freeze"),str(cache),claim_extraction_arg=str(extraction))
            self.assertEqual([doi], pd.read_csv(frozen["frozen"],sep="\t",keep_default_na=False)["source_id"].tolist())
            checked = pd.read_csv(root/"freeze/source_verification/source_validation.tsv",sep="\t",keep_default_na=False)
            self.assertEqual([doi],checked["citation"].tolist())
            self.assertTrue(checked.iloc[0]["citation_verified"])
            prepared = prepare_claim_binding_review(frozen["frozen"],frozen["manifest"],str(root/"prepare"))
            review = pd.read_csv(prepared["review"],sep="\t",dtype=str,keep_default_na=False)
            self.assertEqual([doi],review["selected_citation"].tolist())
            self.assertEqual([doi],review["available_citations"].tolist())
            review["reviewed_by"]="synthetic-reviewer"
            reviewed = root/"reviewed.tsv"
            review.to_csv(reviewed,sep="\t",index=False)
            sealed=seal_claim_binding_review(frozen["frozen"],frozen["manifest"],str(reviewed),str(root/"sealed"))
            self.assertEqual([doi],pd.read_csv(sealed["binding"],sep="\t",keep_default_na=False)["citation"].tolist())
            for column,value in (("available_citations",doi+";corrupt"),("quoted_support",quote+" invented statement")):
                altered=review.copy();altered[column]=value
                altered.to_csv(root/"tampered.tsv",sep="\t",index=False)
                with self.assertRaises(ValueError):
                    seal_claim_binding_review(frozen["frozen"],frozen["manifest"],str(root/"tampered.tsv"),str(root/("bad-"+column)))
