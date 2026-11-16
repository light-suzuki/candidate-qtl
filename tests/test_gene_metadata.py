"""Metadata-only identity checks tolerate display differences and never reread papers."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from candidate_qtl.gene_metadata import GeneMetadataChecks, fetch_gene_metadata
from candidate_qtl.reference_name_pairs import reference_handle


class GeneMetadataTests(unittest.TestCase):
    def test_case_unicode_aliases_and_different_name_lengths_do_not_block(self):
        metadata = {"ids": ["UniProt:P12345"], "names": ["NODE1", "LONG FUNCTIONAL GENE NAME"],
                    "species": ["Plant alpha", "plant common name"]}
        with patch("candidate_qtl.gene_metadata.fetch_gene_metadata", return_value=metadata) as fetch:
            checks = GeneMetadataChecks(None, True)
            for name in ("node1", "ＮＯＤＥ１", "long functional gene name"):
                value = checks.check({"reported_gene": name, "reported_gene_id": "uniprot:p12345", "reported_species": "Plant alpha (common)"})
                self.assertEqual("matched", value["state"])
                self.assertEqual("none", value["membership_effect"])
            value = checks.check({"reported_gene": "another descriptive name", "reported_gene_id": "UniProt:P12345", "reported_species": "Plant alpha"})
            self.assertEqual("id_resolved", value["state"])
            self.assertTrue(value["id_exists"])
            fetch.assert_called_once()

    def test_species_label_is_context_and_existing_reference_metadata_is_reused(self):
        checks = GeneMetadataChecks(None, True)
        checks.seed_record({"ids": ["UniProt:P12345", "GenBank:AB123456"], "names": ["NODE1"], "species": ["renamed species"]})
        with patch("candidate_qtl.gene_metadata.fetch_gene_metadata", side_effect=AssertionError("already retrieved")):
            for identifier in ("UniProt:P12345", "GenBank:AB123456"):
                result = checks.check({"reported_gene_id": identifier, "reported_gene": "node1", "reported_species": "old species name"})
                self.assertEqual("matched", result["state"])
        self.assertEqual(0, checks.lookup_count)

    def test_wrong_returned_id_missing_id_and_provider_failure_remain_distinct(self):
        cases = [({"ids": ["UniProt:Q12345"]}, "mismatch", "returned_identifier_mismatch"),
                 ({"status": "not_found"}, "mismatch", "id_not_found"),
                 (OSError("429"), "unconfirmed", "provider_failure")]
        for metadata, state, reason in cases:
            with self.subTest(reason=reason), patch("candidate_qtl.gene_metadata.fetch_gene_metadata",
                    side_effect=metadata if isinstance(metadata, Exception) else None, return_value=metadata):
                value = GeneMetadataChecks(None, True).check({"reported_gene_id": "UniProt:P12345"})
                self.assertEqual((state, reason), (value["state"], value["reason"]))
                self.assertEqual("none", value["membership_effect"])

    def test_shared_offline_cache_and_major_database_namespaces(self):
        ids = ["UniProt:p12345", "TAIR:at1g12345", "NCBI Gene:12345", "GenBank:ab123456",
               "RefSeq:np_123456.1", "EnsemblPlants:AT1G12345", "Gramene:AT1G12345", "Phytozome:Gene.123"]
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / "cache.json"
            records = {":".join(reference_handle(identifier)): {"ids": [identifier], "names": ["NODE1"], "species": ["Plant alpha"]} for identifier in ids}
            cache.write_text(json.dumps({"records": [], "gene_records": records}))
            checks = GeneMetadataChecks(str(cache), False)
            with patch("candidate_qtl.gene_metadata.fetch_gene_metadata", side_effect=AssertionError("offline")):
                for identifier in ids:
                    self.assertEqual("matched", checks.check({"reported_gene": "node1", "reported_gene_id": identifier, "reported_species": "Plant alpha"})["state"])
            self.assertEqual("unconfirmed", checks.check({"reported_gene_id": "UniProt:Q00001"})["state"])
            self.assertEqual(0, checks.lookup_count)

    def test_uniprot_and_ncbi_gene_adapters_read_only_identity_fields(self):
        with patch("candidate_qtl.gene_metadata._json", return_value={"primaryAccession": "P12345",
                "organism": {"scientificName": "Plant alpha"}, "genes": [{"geneName": {"value": "NODE1"}}]}):
            self.assertEqual(["NODE1"], fetch_gene_metadata("uniprot", "P12345")["names"])
        with patch("candidate_qtl.gene_metadata._json", return_value={"result": {"12345": {
                "uid": "12345", "name": "NODE1", "otheraliases": "ALT1, ALT2", "organism": {"scientificname": "Plant alpha"}}}}):
            value = fetch_gene_metadata("geneid", "12345")
            self.assertEqual(["geneid:12345"], value["ids"])
            self.assertIn("ALT2", value["names"])
