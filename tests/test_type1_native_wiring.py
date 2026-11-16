"""Native Type1 wiring diagnostics; synthetic known answers, not accuracy."""
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

from candidate_qtl.relevance_tiers import T2, annotate_relevance
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for
from tests.test_relevance_tiers import bridge


class EvidenceRetentionTests(unittest.TestCase):
    def test_unknown_pipe_subjects_do_not_collapse(self):
        from candidate_qtl.route_b_prepare import _normalise_subject
        self.assertEqual(_normalise_subject("custom|one"), "custom|one")
        self.assertEqual(_normalise_subject("custom|two"), "custom|two")
        self.assertEqual(_normalise_subject("gi|123|ref|XP_123.1|"), "XP_123.1")

    def test_reannotation_preserves_all_claims(self):
        first = annotate_relevance([{"candidate_gene_id": "g"}], [bridge()])
        second = annotate_relevance(first)
        self.assertEqual(second[0]["relevance_tier"], T2)
        self.assertEqual(second[0]["relevance_evidence_json"], first[0]["relevance_evidence_json"])
        extra = bridge(); extra["source_id"] = "PMID:3"
        third = annotate_relevance(second, [extra])
        self.assertEqual(len(json.loads(third[0]["relevance_evidence_json"])), 2)

    def test_saved_claim_cannot_move_to_another_gene(self):
        row = annotate_relevance([{"candidate_gene_id": "g"}], [bridge()])[0]
        row["candidate_gene_id"] = "other"
        with self.assertRaisesRegex(ValueError, "bind"):
            annotate_relevance([row])

    def test_real_wet_keeps_multiple_claims_and_supported_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            supported = root / "supported.tsv"
            pd.DataFrame([dict(candidate_gene_id="G1", route="b", support_status="verified_semantic_support")]).to_csv(supported, sep="\t", index=False)
            atlas = root / "atlas.tsv"
            pd.DataFrame([dict(gene_id=f"G{i}", family_id="F") for i in range(1, 5)]).to_csv(atlas, sep="\t", index=False)
            records = []
            for gene in ("G1", "G2"):
                for n in (1, 2):
                    claim = bridge(); claim.update(candidate_gene_id=gene, claim_id=f"{gene}-{n}", source_id=f"PMID:{n}", support_status="unverified_route_b_review", rank_contribution="none")
                    records.append(claim)
            review = root / "review.tsv"
            claims = root / "claims.tsv"
            original = pd.DataFrame(records)
            original["member_interval_genes"] = original["candidate_gene_id"]
            original.to_csv(claims, sep="\t", index=False)
            pd.DataFrame(records).drop(columns=["mechanism", "phenotype_relation", "biological_function"]).to_csv(review, sep="\t", index=False)
            manifest = root / "review_manifest.json"
            manifest.write_text(json.dumps(dict(contract="route-review-candidates-1.0.0", stage="route_b_review_candidate_freeze", case_id="case1", output=str(review), artifact_sha256=hashlib.sha256(review.read_bytes()).hexdigest(), inputs={"claims": {"path": str(claims), "sha256": hashlib.sha256(claims.read_bytes()).hexdigest()}})))
            result = wet_candidate_handoff(str(root / "out"), route_b_candidates_arg=str(supported), route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case1")), review_candidates_arg=str(review), review_candidates_manifest_arg=str(manifest), interval_atlas_arg=str(atlas), case_id="case1", wet_budget=2, budget_policy="strict")
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(list(wet.candidate_gene_id), ["G1", "G2"])
            self.assertEqual(list(wet.relevance_tier), [T2, T2])
            for value in wet.relevance_evidence_json:
                self.assertEqual(len(json.loads(value)), 2)
            self.assertEqual(list(wet.support_tier), ["tier1_verified_support", "unverified_route_review"])
            inventory = pd.read_csv(root / "out" / "wet_relevance_inventory.tsv", sep="\t")
            self.assertEqual(len(inventory), 4)
            meta = json.loads((root / "out" / "wet_candidates_manifest.json").read_text())
            self.assertFalse(meta["relevance_selection_applied"])
            self.assertEqual(meta["wet_candidates_sha256"], hashlib.sha256(Path(result["wet_candidates"]).read_bytes()).hexdigest())


@unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "local BLAST required")
class NativeIdentityTests(unittest.TestCase):
    def test_installed_cli_mixed_db_missing_taxonomy_and_retrieval_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(751)
            sequence = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150))
            query = root / "query.fa"
            query.write_text(f">g0\n{sequence}\n")
            candidates = root / "genes.tsv"
            pd.DataFrame([dict(gene_id=f"g{i}", chrom="chr1", start=i*200+1, end=i*200+150) for i in range(1000)]).to_csv(candidates, sep="\t", index=False)
            family = root / "family.tsv"
            family.write_text("gene_id\tfamily_id\n")
            corpus = root / "corpus.tsv"
            corpus.write_text("reference_gene_id\treference_species\tdescription\tmutant_phenotype\tsource_id\n")
            env = dict(os.environ); env.pop("PYTHONPATH", None)
            sets = []
            for run, order in enumerate((range(7), reversed(range(7)))):
                fasta = root / f"db{run}.fa"
                fasta.write_text("".join(f">ref{i} protein OS=Species {i % 2} OX={100+i%2} GN=R{i}\n{sequence}\n" for i in order))
                subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "prot"], check=True, capture_output=True)
                out = root / f"out{run}"
                command = [sys.executable, "-m", "candidate_qtl", "route-b-prepare", "--candidates", str(candidates), "--proteins", str(query), "--protein-db", str(fasta), "--family-index", str(family), "--reference-function", str(corpus), "--out-dir", str(out)]
                completed = subprocess.run(command, env=env, capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                inventory = pd.read_csv(out / "route_b_identity_inventory.tsv", sep="\t")
                hits = pd.read_csv(out / "route_b_identity_hits.tsv", sep="\t")
                self.assertEqual(len(inventory), 1000)
                self.assertEqual(inventory.gene_id.nunique(), 1000)
                self.assertEqual(len(hits), 5)
                self.assertEqual(set(hits.hit_species), {"Species 0", "Species 1"})
                self.assertTrue(inventory.loc[inventory.gene_id != "g0", "identity_state"].eq("not_evaluated_no_protein_query").all())
                meta = json.loads((out / "route_b_preparation_manifest.json").read_text())
                self.assertEqual(meta["providers"]["homology_search"]["retrieval_max_target_seqs"], 5)
                self.assertEqual((out / "interval_gene_atlas.tsv").read_bytes(), (out / "atlas/interval_gene_atlas.tsv").read_bytes())
                sets.append(set(hits.hit_accession))
            # Retrieval boundary ties depend on DB order: record this limitation,
            # do not advertise post-retrieval invariance as whole-DB completeness.
            self.assertNotEqual(sets[0], sets[1])

