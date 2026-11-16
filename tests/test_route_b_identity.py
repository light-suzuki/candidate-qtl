"""Identity/coverage regressions use synthetic sequences, no online providers."""
from __future__ import annotations

import csv
import itertools
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.route_b_identity import header_metadata, merge_header_metadata, protein_query_ids, write_identity_layer
from candidate_qtl.route_b_blast_evidence import run_local_blastp


def hsp(gene="g", ref="r", title="", taxon="", species=""):
    return "\t".join(map(str, [gene, ref, 100, 4, 4, 4, 1, 4, 1, 4,
                                1e-20, 100, "AAAA", "AAAA", title, taxon, species]))


class IdentityTests(unittest.TestCase):
    def test_uniprot_explicit_fields(self):
        r=header_metadata("sp|P1|NAME_ARATH Protein kinase OS=Arabidopsis thaliana OX=3702 GN=KIN1 PE=1 SV=2")
        self.assertEqual((r["hit_species"],r["hit_taxon_id"],r["hit_gene_symbol"],r["hit_protein_name"]),
                         ("Arabidopsis thaliana","3702","KIN1","Protein kinase"))

    def test_ncbi_header_does_not_invent_gene_symbol_or_taxon(self):
        r=header_metadata("XP_012345.1 receptor kinase [Medicago truncatula]")
        self.assertEqual(r["hit_species"],"Medicago truncatula")
        self.assertEqual(r["hit_protein_name"],"receptor kinase")
        self.assertEqual(r["hit_gene_symbol"],"unknown")
        self.assertEqual(r["hit_taxon_id"],"unknown")

    def test_suffix_alone_is_not_species_evidence(self):
        r=header_metadata("NAME_ARATH")
        self.assertEqual(r["hit_species"],"unknown")
        self.assertEqual(r["hit_protein_name"],"unknown")

    def test_database_taxonomy_without_header(self):
        r=merge_header_metadata([("", "3880", "Medicago truncatula")])
        self.assertEqual(r["hit_species"],"Medicago truncatula")
        self.assertEqual(r["hit_taxon_id"],"3880")
        self.assertEqual(r["hit_gene_symbol"],"unknown")

    def test_zero_database_taxonomy_does_not_conflict_with_header_taxonomy(self):
        title="sp|P1|NAME_ARATH Protein kinase OS=Arabidopsis thaliana OX=3702 GN=KIN1"
        r=merge_header_metadata([(title, "0", "N/A")])
        self.assertEqual(r["hit_taxon_id"],"3702")
        self.assertEqual(r["hit_metadata_state"],"complete")
        self.assertEqual(json.loads(r["hit_metadata_records_json"])[0]["database_taxon_ids"],"0")

    def test_zero_only_taxonomy_remains_unknown_without_conflict(self):
        r=merge_header_metadata([("", "0", "")])
        self.assertEqual(r["hit_taxon_id"],"unknown")
        self.assertEqual(r["hit_metadata_state"],"missing")
        self.assertEqual(json.loads(r["hit_metadata_records_json"])[0]["database_taxon_ids"],"0")

    def test_real_database_header_taxonomy_conflict_remains_conflicting(self):
        r=merge_header_metadata([("sp|P1|NAME_ARATH Protein kinase OS=Arabidopsis thaliana OX=3702 GN=KIN1", "3880", "")])
        self.assertEqual(r["hit_taxon_id"],"unknown")
        self.assertEqual(r["hit_metadata_state"],"conflicting")

    def test_zero_gene_symbol_and_protein_name_are_not_taxonomy_placeholders(self):
        r=merge_header_metadata([("sp|P1|NAME_ARATH 0 OS=Arabidopsis thaliana OX=3702 GN=0", "0", "")])
        self.assertEqual(r["hit_gene_symbol"],"0")
        self.assertEqual(r["hit_protein_name"],"0")

    def test_conflicts_preserve_all_records_without_last_writer(self):
        a=("sp|P1|X kinase OS=Species one OX=1 GN=A", "1", "Species one")
        b=("sp|P1|X kinase OS=Species two OX=2 GN=B", "2", "Species two")
        r=merge_header_metadata([a,b,a])
        self.assertEqual(r,merge_header_metadata([b,a]))
        self.assertEqual(r["hit_metadata_state"],"conflicting")
        self.assertEqual(r["hit_species"],"unknown")
        self.assertEqual(len(json.loads(r["hit_metadata_records_json"])),2)

    def test_unknown_metadata_is_not_no_hit(self):
        r=merge_header_metadata([("N/A", "N/A", "N/A")])
        self.assertEqual(r["hit_metadata_state"],"missing")

    def run_blast(self, root, lines, returncode=0):
        def run(command, **kwargs):
            Path(command[command.index("-out")+1]).write_text("\n".join(lines)+"\n")
            return subprocess.CompletedProcess(command, returncode, "", "failed" if returncode else "")
        return run_local_blastp(root/"query.fa", root/"mixed.db", root, 1, 30, .5, 1e-5,
                               run_tool=run, find_binary=lambda _: "blastp", normalise_subject=str)

    def test_all_retrieved_tied_hits_and_byte_stable_pair_metadata(self):
        lines=[hsp(ref="b",title="kinase [Species one]"), hsp(ref="a",title="kinase [Species two]"),hsp(ref="c")]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); expected=None
            for permutation in itertools.permutations(lines):
                rows,meta=self.run_blast(root,permutation)
                self.assertEqual([r["reference_gene_id"] for r in rows],["a","b","c"])
                self.assertTrue(all(r["relationship"]=="best_hit_homology" for r in rows))
                payload=Path(meta["pair_evidence_table"]).read_bytes()
                if expected is None: expected=payload
                self.assertEqual(payload,expected)
                self.assertIsNone(meta["species_quota"])

    def test_same_species_multiple_hits_are_legal(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows,_=self.run_blast(Path(tmp),[hsp(ref="a",species="Species one"),hsp(ref="b",species="Species one")])
            self.assertEqual(len(rows),2)

    def test_failure_distinct_from_successful_no_hit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            rows,failed=self.run_blast(root,[],1)
            self.assertEqual(rows,[]); self.assertEqual(failed["evaluation_state"],"provider_failure")
            _,success=self.run_blast(root,[])
            self.assertTrue(success["available"]); self.assertEqual(success["evaluation_state"],"completed")

    def test_all_gene_coverage_axes_and_missing_protein(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            hits,provider=self.run_blast(root,[hsp(ref="r",title="kinase [Species one]")])
            atlas=[{"gene_id":"g"},{"gene_id":"nohit"},{"gene_id":"missing"},{"gene_id":"family","family_id":"F"}]
            rows,meta=write_identity_layer(atlas,hits,provider,[],root,query_ids={"g","nohit"})
            byid={r["gene_id"]:r for r in rows}
            self.assertEqual(len(rows),4)
            self.assertEqual(byid["nohit"]["identity_state"],"evaluated_no_threshold_hit")
            self.assertEqual(byid["missing"]["identity_state"],"not_evaluated_no_protein_query")
            self.assertEqual(byid["family"]["unresolved_state"],"resolved_handle")
            self.assertEqual(meta["blast_identity_coverage"]["gene_count"],1)
            self.assertEqual(meta["functional_annotation_coverage"]["gene_count"],1)
            self.assertEqual(meta["source_backed_literature_coverage"]["gene_count"],0)
            self.assertEqual(meta["dark_unresolved_coverage"]["gene_count"],2)

    def test_zero_blast_taxonomy_preserves_header_identity_through_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            title="sp|P1|NAME_ARATH Protein kinase OS=Arabidopsis thaliana OX=3702 GN=KIN1"
            hits,provider=self.run_blast(root,[hsp(ref="P1",title=title,taxon="0",species="N/A")])
            self.assertEqual(hits[0]["hit_taxon_id"],"3702")
            rows,meta=write_identity_layer([{"gene_id":"g"}],hits,provider,[],root,query_ids={"g"})
            self.assertEqual(rows[0]["identity_metadata_state"],"complete")
            self.assertEqual(rows[0]["identity_reference_species"],"Arabidopsis thaliana")
            self.assertEqual(rows[0]["identity_gene_symbols"],"KIN1")
            with Path(meta["identity_hits"]["path"]).open() as handle:
                recorded=list(csv.DictReader(handle,delimiter="\t"))[0]
            self.assertEqual(recorded["hit_taxon_id"],"3702")
            self.assertEqual(json.loads(recorded["hit_metadata_records_json"])[0]["database_taxon_ids"],"0")
            self.assertEqual(recorded["rank_contribution"],"none")

    def test_reference_literature_join_does_not_require_best_hit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            hits,provider=self.run_blast(root,[hsp(ref="r1"),hsp(ref="r2")])
            _,meta=write_identity_layer([{"gene_id":"g"}],hits,provider,
                [{"reference_gene_id":"r2","literature_source_ids":"PMID:123","function_description":"kinase"}],root,query_ids={"g"})
            self.assertEqual(meta["source_backed_literature_coverage"]["gene_count"],1)
            with Path(meta["identity_hits"]["path"]).open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle,delimiter="\t"))),2)

    def test_inventory_and_hits_independent_of_input_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            hits,provider=self.run_blast(root,[hsp(gene="g1",ref="a"),hsp(gene="g2",ref="b")])
            atlas=[{"gene_id":"g1"},{"gene_id":"g2"},{"gene_id":"g3"}]
            _,one=write_identity_layer(atlas,hits,provider,[],root,query_ids={"g1","g2","g3"})
            _,two=write_identity_layer(reversed(atlas),reversed(hits),provider,[],root,query_ids={"g3","g2","g1"})
            for key in ("identity_hits","identity_inventory"):
                self.assertEqual(one[key]["sha256"],two[key]["sha256"])

    def test_duplicate_ids_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,"unique interval"):
                write_identity_layer([{"gene_id":"g"},{"gene_id":"g"}],[],{},[],root,query_ids=set())
            p=root/"q.fa"; p.write_text(">g\nAAA\n>g\nAAA\n")
            with self.assertRaisesRegex(ValueError,"unique query"):
                protein_query_ids(p)

    def test_real_prepare_consumer_with_1000_gene_inventory(self):
        # Execute the real orchestrator with synthetic atlas/search backends.
        # This does not claim installed CLI, live BLAST or biological validation.
        from candidate_qtl.route_b_prepare import route_b_prepare
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); candidates=root/"genes.tsv"; proteins=root/"q.fa"; family=root/"family.tsv"
            candidates.write_text("gene_id\tchrom\tstart\tend\n"+"".join(f"g{i:04}\tchr1\t{i*10+1}\t{i*10+9}\n" for i in range(1000)))
            proteins.write_text("".join(f">g{i:04}\nAAAA\n" for i in range(1000)))
            family.write_text("gene_id\tfamily_id\n")
            with patch("candidate_qtl.route_b_prepare._run_blastp",
                       side_effect=lambda proteins,database,out,*args:self.run_blast(out,[hsp(gene="g0000",ref="r")])):
                result=route_b_prepare(str(candidates),str(root/"out"),proteins_arg=str(proteins),
                            protein_db_arg="synthetic.db",family_index_arg=str(family),auto_local_db=False)
            with Path(result["identity_inventory"]).open() as handle:
                rows=list(csv.DictReader(handle,delimiter="\t"))
            self.assertEqual(len(rows),1000)
            self.assertEqual(rows[0]["identity_hit_count"],"1")
            manifest=json.loads(Path(result["manifest"]).read_text())
            self.assertEqual(manifest["identity_coverage"]["blast_identity_coverage"]["gene_count"],1)
            self.assertEqual(manifest["outputs"]["identity_hits"],result["identity_hits"])
            self.assertEqual(Path(result["interval_gene_atlas"]).read_bytes(),
                             (root/"out/atlas/interval_gene_atlas.tsv").read_bytes())


if __name__ == "__main__":
    unittest.main()
