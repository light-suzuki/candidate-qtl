"""Neutral explicit shared-gene identities, not generic accession fanout."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from candidate_qtl.core import sha256
from candidate_qtl.observation_tiers import classify_observation
from candidate_qtl.paper_evidence_table import write_paper_evidence_table
from candidate_qtl.paper_reports import bound_paper_rows
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from tests.test_resolved_paper_reference_binding import paper, record


class SharedGeneReferenceTransportTests(unittest.TestCase):
    def inputs(self, root, *, insdc=False, other=True, mutate=None):
        cache = root / "records"
        cache.mkdir()
        records = [record("Q00001"), record("Q00002")]
        if other:
            records.append(record("Q00003"))
            records[-1]["uniProtKBCrossReferences"][0]["id"] = "99999"
        for item in records:
            item["organism"]["taxonId"] = 999
        official = None
        if insdc:
            official = SeqRecord(Seq("ATG" * 100), id="AA123456.1", name="AA123456", description="Neutral gene alternative CDS")
            official.annotations.update(molecule_type="DNA", organism="Species alpha")
            official.features = [SeqFeature(SimpleLocation(0, 300), type="source",
                qualifiers={"organism": ["Species alpha"], "db_xref": ["taxon:999"], "cultivar": ["neutral-original"]}),
                SeqFeature(SimpleLocation(0, 250, strand=1), type="gene", qualifiers={"gene": ["LN1"], "locus_tag": ["neutral_locus"]})]
            for index, item in enumerate(records):
                nucleotide = "AA123456" if index < 2 else "BB123456"
                protein = f"CAA0000{index+1}.1"
                item["uniProtKBCrossReferences"] = [{"database": "EMBL", "id": nucleotide,
                    "properties": [{"key": "ProteinId", "value": protein}]}]
                if index < 2:
                    official.features.append(SeqFeature(SimpleLocation(10, 100 + index * 3, strand=1), type="CDS",
                        qualifiers={"gene": ["LN1"], "locus_tag": ["neutral_locus"], "protein_id": [protein],
                                    "note": ["alternatively spliced"], "translation": ["M" * 30]}))
        if mutate:
            mutate(records, official)
        for item in records:
            (cache / f"{item['primaryAccession']}.json").write_text(json.dumps(item), encoding="utf-8")
        if official:
            SeqIO.write(official, cache / "AA123456.embl", "embl")
        anchor = "GenBank:AA123456" if insdc else "GeneID:12345"
        original = paper(anchor)
        source = root / "mechanisms.tsv"
        pd.DataFrame([{"hypothesis_id": "neutral-reference", "case_id": "neutral", "reference_gene_name": "LN1",
            "reference_species": "Species alpha", "source_ids": "PMID:123", "proposed_stable_id": anchor,
            "paper_reports_json": json.dumps([original])}]).to_csv(source, sep="\t", index=False)
        return cache, source, original

    def resolve(self, root, *, cli=False, **kwargs):
        cache, source, original = self.inputs(root, **kwargs)
        original_hash = sha256(source)
        if cli:
            execution = subprocess.run([sys.executable, "-m", "candidate_qtl", "reference-gene-resolve",
                "--mechanisms", str(source), "--records", str(cache), "--out", str(root / "resolved"), "--json"],
                text=True, capture_output=True)
            self.assertEqual(0, execution.returncode, execution.stderr)
            result = json.loads(execution.stdout)
        else:
            result = reference_gene_resolve(str(source), str(root / "resolved"), record_cache_arg=str(cache))
        self.assertEqual(original_hash, sha256(source))
        rows = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str, keep_default_na=False).to_dict("records")
        return rows, cache, original, result

    def test_geneid_and_official_alternative_cds_transport_all_owners_only(self):
        for insdc in (False, True):
            with self.subTest(insdc=insdc), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                rows, _, original, result = self.resolve(root, insdc=insdc, cli=True)
                self.assertEqual(3, len(rows))
                for row in rows:
                    reports = json.loads(row["paper_reports_json"])
                    owner = row["proposed_stable_id"] in {"Q00001", "Q00002"}
                    expected = row["proposed_stable_id"] if owner else original["reference_gene_id"]
                    self.assertEqual([{**original, "reference_gene_id": expected}], reports)
                    classified, = bound_paper_rows(row)
                    self.assertEqual("matched" if owner else "unresolved", classified["tier_gene_binding"])
                    if owner:
                        self.assertEqual("T2", classify_observation(classified)[0])
                        proof, = json.loads(row["retrieval_provenance"])["paper_shared_reference_proofs"]
                        self.assertEqual(["Q00001", "Q00002"], proof["owners"])
                        self.assertEqual("unknown", reports[0]["paper_checklist"]["same_gene_as_focal"])
                exported = write_paper_evidence_table(rows, root / "paper-view", route="phenotype_first",
                    source_artifact=result["mechanisms_resolved"])
                table = pd.read_csv(exported["outputs"]["paper_evidence_table"], sep="\t", dtype=str)
                self.assertEqual(2, table.reference_anchor_state.eq("matched").sum())

    def test_insdc_contig_multiple_genes_loci_proteins_taxa_and_sequences_fail_closed(self):
        mutations = {
            "multiple_genes": lambda rs, gb: gb.features.append(SeqFeature(SimpleLocation(251, 290), type="gene", qualifiers={"gene": ["OTHER"]})),
            "same_name_other_locus": lambda rs, gb: gb.features[-1].qualifiers.update(locus_tag=["other_locus"]),
            "multiple_gene_ids_in_feature": lambda rs, gb: gb.features[1].qualifiers.update(gene=["LN1", "OTHER"]),
            "wrong_protein": lambda rs, gb: rs[1]["uniProtKBCrossReferences"][0]["properties"][0].update(value="CAA99999.1"),
            "wrong_source_taxid": lambda rs, gb: gb.features[0].qualifiers.update(db_xref=["taxon:998"]),
            "different_record": lambda rs, gb: setattr(gb, "id", "BB123456.1"),
            "wrong_sequence": lambda rs, gb: rs[1]["sequence"].update(value="A" * 30),
            "no_explicit_splice": lambda rs, gb: gb.features[-1].qualifiers.update(note=["unresolved product"]),
            "missing_cds_owner": lambda rs, gb: rs.pop(1),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                rows, _, original, _ = self.resolve(Path(temporary), insdc=True, mutate=mutate)
                for row in rows:
                    self.assertEqual([original], json.loads(row["paper_reports_json"]))
                    self.assertNotIn("paper_shared_reference_proofs", json.loads(row["retrieval_provenance"]))

    def test_missing_parsed_feature_location_is_unknown_for_producer_and_rejected_by_consumer(self):
        for index in (1, -1):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                rows, cache, original, _ = self.resolve(root, insdc=True)
                malformed = SeqIO.read(cache / "AA123456.embl", "embl")
                malformed.features[index].location = None
                row = next(row for row in rows if row["proposed_stable_id"] == "Q00001")
                with mock.patch("candidate_qtl.reference_shared_gene.SeqIO.read", return_value=malformed):
                    with self.assertRaisesRegex(ValueError, "location"):
                        bound_paper_rows(row)
                    result = reference_gene_resolve(str(root / "mechanisms.tsv"), str(root / "malformed-result"), record_cache_arg=str(cache))
                malformed_rows = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str, keep_default_na=False).to_dict("records")
                for candidate in malformed_rows:
                    self.assertEqual([original], json.loads(candidate["paper_reports_json"]))

    def test_shared_geneid_without_consistent_explicit_taxon_or_with_other_geneid_is_not_authorized(self):
        mutations = [lambda rs, gb: rs[1]["organism"].update(taxonId=998),
            lambda rs, gb: rs[1]["organism"].pop("taxonId"),
            lambda rs, gb: rs[1]["uniProtKBCrossReferences"].append({"database": "GeneID", "id": "77777"})]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as temporary:
                rows, _, original, _ = self.resolve(Path(temporary), mutate=mutate)
                for row in rows:
                    self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_record_enumeration_order_does_not_choose_an_isoform(self):
        observed = []
        for reverse in (False, True):
            with tempfile.TemporaryDirectory() as temporary:
                rows, _, _, _ = self.resolve(Path(temporary), insdc=True,
                    mutate=(lambda rs, gb: rs.reverse()) if reverse else None)
                observed.append([(r["proposed_stable_id"], json.loads(r["paper_reports_json"]),
                    [p["owners"] for p in json.loads(r["retrieval_provenance"]).get("paper_shared_reference_proofs", [])]) for r in rows])
        self.assertEqual(*observed)

    def test_partial_query_cannot_authorize_a_shared_gene_handle(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cache, source, original = self.inputs(root, other=False)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_reference_candidates",
                    return_value=(["Q00001", "Q00002"], "provider_unavailable", "neutral-partial-query", 2,
                                  ["incomplete_page"], [])):
                result = reference_gene_resolve(str(source), str(root / "resolved"), record_cache_arg=str(cache), allow_external_api=True)
            rows = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str, keep_default_na=False).to_dict("records")
            for row in rows:
                self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_consumers_reject_changed_raw_insdc_provider_or_resealed_inconsistent_proof(self):
        for artifact in ("AA123456.embl", "Q00002.json", "Q00003.json"):
            with self.subTest(artifact=artifact), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                rows, cache, _, result = self.resolve(root, insdc=True)
                path = cache / artifact
                path.write_bytes(path.read_bytes() + b"\n")
                row = next(row for row in rows if row["proposed_stable_id"] == "Q00001")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    bound_paper_rows(row)
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    write_paper_evidence_table(rows, root / "paper-view", route="phenotype_first", source_artifact=result["mechanisms_resolved"])
        with tempfile.TemporaryDirectory() as temporary:
            rows, _, _, _ = self.resolve(Path(temporary), insdc=True)
            row = next(row for row in rows if row["proposed_stable_id"] == "Q00001")
            provenance = json.loads(row["retrieval_provenance"])
            provenance["paper_shared_reference_proofs"][0]["owners"] = ["Q00001", "Q00003"]
            row["retrieval_provenance"] = json.dumps(provenance)
            with self.assertRaisesRegex(ValueError, "owner mismatch"):
                bound_paper_rows(row)

    def test_consumers_reject_orphan_proofs_and_deleted_shared_receipts(self):
        for change in ("delete_receipts", "delete_proofs", "orphan", "duplicate"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temporary:
                rows, _, _, _ = self.resolve(Path(temporary), insdc=True)
                row = next(row for row in rows if row["proposed_stable_id"] == "Q00001")
                provenance = json.loads(row["retrieval_provenance"])
                if change == "delete_receipts":
                    provenance.pop("paper_reference_rebindings")
                elif change == "delete_proofs":
                    provenance.pop("paper_shared_reference_proofs")
                elif change == "orphan":
                    proof = {**provenance["paper_shared_reference_proofs"][0], "handle": ["geneid", "77777"]}
                    provenance["paper_shared_reference_proofs"].append(proof)
                else:
                    provenance["paper_shared_reference_proofs"] *= 2
                row["retrieval_provenance"] = json.dumps(provenance)
                with self.assertRaisesRegex(ValueError, "proofs|proof"):
                    bound_paper_rows(row)
