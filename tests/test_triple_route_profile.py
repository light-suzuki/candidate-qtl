"""Three independent routes, run as separate lanes, meeting only at the comparison.

The point of running A/B/C is that they fail differently, so the comparison must
show which characteristic error was actually cancelled per family -- not just
that the routes agreed.
"""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from candidate_qtl.cli import main
from candidate_qtl.guided_profiles import _assert_lane_isolation, _resolve_chained
from candidate_qtl.retrieval_route import triple_route_compare
from tests.bound_claim_fixtures import write_route_a_bound_claims

CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_CLV1\tchr5\t1000\t3000\tCLAVATA1 like receptor kinase acting in the meristem",
    "Psat_UNK\tchr5\t9000\t11000\thypothetical protein",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\torthogroup\trelationship\tprovider",
    "Psat_CLV1\tAtCLV1\tOG0002\tbest_hit_homology\tlocal_blastp:test",
]) + "\n"

REFERENCE_CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase restricting the stem cell pool\t"
    "fasciated stem and enlarged meristem\tmeristem_maintenance\tPMID:8402901",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed\tstarch_biosynthesis\tPMID:1978859",
]) + "\n"

# Three competing mechanisms; only CLV is actually present in the interval.
PAYLOAD = "\n".join([
    "causal_layer\tpathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
    "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
    "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
    "signaling\tmeristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
    "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
    "P30924\tArabidopsis thaliana\tCLV1\treviewed receptor kinase with a fasciation phenotype",
    "execution\thormone transport\tPIN\tauxin redistribution changes stem geometry\tPMID:1978859\t"
    "meristem lesions give the same shape\tassay auxin flux\tArabidopsis pin mutant\tfused organs\t"
    "Q9SXK1\tArabidopsis thaliana\tPIN1\tcanonical auxin efflux carrier",
    "execution\tcell wall biosynthesis\tCESA\tweakened wall changes stem shape\tPMID:12172025\t"
    "wall phenotypes rarely fasciate\tsection the stem\tArabidopsis irx mutant\tcollapsed xylem\t"
    "Q9SWW6\tArabidopsis thaliana\tCESA9\tsecondary wall cellulose synthase",
    "regulation\tmeristem transcription\tWOX\tWOX regulation changes stem-cell maintenance\t\tCLV signaling may be primary\tassay WOX expression\t\t\t\tArabidopsis thaliana\tWOX1\tindependent regulatory layer",
    "perception\tpeptide perception\tCLE receptor\tpeptide perception changes meristem size\t\ttransport may explain geometry\tmeasure peptide response\t\t\t\tArabidopsis thaliana\tCLV2\tupstream perception layer",
    "signaling\tMAP kinase signaling\tMAPK\tMAPK signaling changes division dynamics\t\twall weakness may be primary\tmeasure phosphoproteins\t\t\t\tArabidopsis thaliana\tMPK6\tsignaling complement",
    "regulation\torgan polarity\tHD-ZIP\tpolarity regulation changes stem organization\t\thormone flux may be primary\tmap polarity markers\t\t\t\tArabidopsis thaliana\tREV\tregulatory complement",
    "alternative\tlateral identity\tKANADI\tlateral identity defects alter stem geometry\t\tmeristem size may be normal\tsection organ boundaries\t\t\t\tArabidopsis thaliana\tKAN1\tregulatory alternative",
    "execution\tcell expansion\tEXPANSIN\tcell expansion defects alter stem width\t\tdivision defects may be primary\tmeasure cell dimensions\t\t\t\tArabidopsis thaliana\tEXPA1\tdownstream execution",
    "signaling\tredox signaling\tRBOH\tredox signaling changes meristem division\t\tpeptide perception may be primary\tmeasure ROS dynamics\t\t\t\tArabidopsis thaliana\tRBOHD\tindependent signal",
    "perception\treceptor signaling\tLRR-RLK\treceptor input changes stem-cell balance\t\tdownstream execution may dominate\tmeasure ligand response\t\t\t\tArabidopsis thaliana\tERECTA\tsecond perception family",
    "alternative\tmechanical feedback\tMECHANOCHANNEL\tmechanical feedback could mimic fasciation\t\tgenetic meristem signaling may be primary\tapply mechanical perturbation\t\t\t\tArabidopsis thaliana\tMSL10\tnoncanonical alternative",
]) + "\n"

# The frozen T11 legacy format verbatim (coordination/payload_freeze.tsv records
# its sha256); no stable accession exists anywhere in these seven columns.
LEGACY_PAYLOAD = "\n".join([
    "case_id\tmechanism_id\tmechanism_name\tevidence_class\tsource_log_entry\tstatus\tselection_reason",
    "pea_a\tpea_a_M1\tChalcone synthase (CHS) first committed step of anthocyanin biosynthesis\t"
    "corpus_ortholog_model_plant\tsource_log.tsv#pea_a_M1\tverified\t"
    "Flower colour and seed-coat anthocyanin phenotype text maps to anthocyanin biosynthesis; "
    "CHS is the committed first step (corpus: At CHS, Mt CHS-class).",
    "pea_a\tpea_a_M2\tLate anthocyanin pathway enzymes (DFR / LDOX/ANS)\t"
    "corpus_ortholog_model_plant\tsource_log.tsv#pea_a_M2\tverified\t"
    "Anthocyanin pigmentation loss in flower/seed coat follows late-pathway enzyme defects; "
    "DFR and LDOX corpus-verified (At).",
    "pea_a\tpea_a_M3\tMYB-bHLH-WD40 (MBW) transcriptional regulation of anthocyanin genes\t"
    "corpus_ortholog_model_plant\tsource_log.tsv#pea_a_M3\tverified\t"
    "Flower colour and seed coat anthocyanin are classically controlled by MBW complexes; "
    "TT2/TT8/TTG1 entries corpus-verified in Arabidopsis and Lotus (LjTT2a-c/LjTT8/LjTTG1).",
]) + "\n"

# Same legacy pea_a content in canonical columns (still no accession) plus one
# canonical row that carries proposed_stable_id + reference_species, so the a04
# verification separation can be asserted on both branches in one offline run.
CANONICAL_MIXED_PAYLOAD = "\n".join([
    "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
    "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
    "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
    # Issue #612: canonical-schema rows may not carry unresolvable pointer
    # tokens in source_ids; the two rows below stay unsourced on purpose so the
    # a04 no-stable-accession rejection path is still exercised.
    "anthocyanin biosynthesis\tCHS\tChalcone synthase (CHS) first committed step of anthocyanin biosynthesis\t"
    "\t\t\t\t\t\t\t\t"
    "Flower colour and seed-coat anthocyanin phenotype text maps to anthocyanin biosynthesis",
    "anthocyanin biosynthesis\tDFR\tLate anthocyanin pathway enzymes (DFR / LDOX/ANS)\t"
    "\t\t\t\t\t\t\t\t"
    "Anthocyanin pigmentation loss in flower/seed coat follows late-pathway enzyme defects",
    "meristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
    "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
    "P30924\tArabidopsis thaliana\tCLV1\treviewed receptor kinase with a fasciation phenotype",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _localization(root: Path, phrase: str, name: str = "localization.json") -> Path:
    return _write(root, name, json.dumps({
        "causal_localization": "signaling",
        "localization_confidence": "high",
        "supporting_phrase": phrase,
    }))


def _record_cache(root: Path) -> Path:
    cache = root / "records"
    cache.mkdir()
    for accession, gene_name in (("P30924", "CLV1"), ("Q9SXK1", "PIN1"), ("Q9SWW6", "CESA9")):
        (cache / f"{accession}.json").write_text(json.dumps({
            "primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Arabidopsis thaliana"},
            "sequence": {"value": "MASTPRQTVSEK" * 8},
            "genes": [{"geneName": {"value": gene_name}}],
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "8402901"}]}}],
        }), encoding="utf-8")
    return cache


def _source_cache(root: Path) -> Path:
    cache = root / "pmc"
    cache.mkdir()
    for pmid, title in (("8402901", "CLAVATA1 receptor kinase controls meristem size in Arabidopsis"),
                        ("1978859", "Starch branching enzyme SBE and the wrinkled seed character"),
                        ("12172025", "CESA cellulose synthase and secondary wall formation")):
        (cache / f"med_{pmid}.json").write_text(json.dumps({
            "cache_key": f"pmid_{pmid}", "id": pmid, "pmid": pmid, "title": title,
            "abstractText": f"Genetic loss-of-function mutants caused the phenotype. {title}",
            "journalTitle": "Plant Cell", "pubYear": "1993",
            "doi": f"10.1/{pmid}", "authorString": "A B, C D",
        }), encoding="utf-8")
    return cache


def _hit13(query: str, subject: str, identity: float, bitscore: float, title: str, *, qlen: int = 300) -> str:
    return (f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}\t{title}")


class TripleRouteProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        # These fixtures supply their own indexes or exercise missing-index
        # behavior. A workstation's reusable genome cache is unrelated input.
        for isolation in (
            patch.dict(os.environ, {"CANDIDATE_QTL_FAMILY_INDEX": ""}),
            patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]),
        ):
            isolation.start()
            self.addCleanup(isolation.stop)

    def test_hypothesis_verification_consumes_source_annotated_mechanisms(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "plan"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out),
            ]))
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            resolver = next(stage for stage in profile["stages"] if stage["command"] == "reference-gene-resolve")
            self.assertEqual(
                {"stage": "mechanism-coverage-audit", "result_key": "mechanisms"},
                resolver["inputs_from"]["mechanisms"],
            )
            source_verify = next(stage for stage in profile["stages"] if stage["command"] == "source-verify")
            self.assertEqual(
                {"stage": "reference-gene-resolve", "result_key": "mechanisms_resolved"},
                source_verify["inputs_from"]["mechanisms"],
            )
            verify = next(stage for stage in profile["stages"] if stage["command"] == "hypothesis-verify")
            self.assertEqual(
                {"stage": "source-verify", "result_key": "mechanisms_with_verified_sources"},
                verify["inputs_from"]["hypotheses"],
            )
            self.assertEqual(
                {"stage": "source-verify", "result_key": "verified_sources"},
                verify["inputs_from"]["source_records"],
            )
            freeze = next(stage for stage in profile["stages"] if stage["command"] == "bait-panel-freeze")
            self.assertEqual(
                {"stage": "source-verify", "result_key": "mechanisms_with_verified_sources"},
                freeze["inputs_from"]["hypotheses"],
            )

    def test_case_and_resources_add_the_explicit_triple_route_bridge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "plan"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                "--case", str(root / "case.yaml"), "--resources", str(root / "resources.yaml"),
            ]))
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            case_stage = profile["stages"][-1]
            self.assertEqual("case-run", case_stage["command"])
            self.assertEqual(
                ["family-to-copy", "phenotype-first-status", "route-b-prepare", "triple-route-compare"],
                sorted({value["stage"] for value in case_stage["inputs_from"].values()}),
            )
            self.assertIn("narrowing_family_copy", case_stage["inputs_from"])
            self.assertIn("narrowing_convergence", case_stage["inputs_from"])
            self.assertIn("narrowing_convergence_manifest", case_stage["inputs_from"])
            self.assertIn("narrowing_family_copy_manifest", case_stage["inputs_from"])
            self.assertEqual(
                ("triple-route-compare", "triple_route_hypothesis_comparison"),
                (case_stage["inputs_from"]["discordance_comparison"]["stage"],
                 case_stage["inputs_from"]["discordance_comparison"]["result_key"]),
            )
            self.assertEqual(
                ("phenotype-first-status", "status"),
                (case_stage["inputs_from"]["route_a_manifest"]["stage"],
                 case_stage["inputs_from"]["route_a_manifest"]["result_key"]),
            )

    def test_shared_route_is_wired_into_the_b1_overlay_plan(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "plan"
            route = root / "shared_route.tsv"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                "--shared-route", str(route),
            ]))
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            b1 = next(stage for stage in profile["stages"] if stage["command"] == "interval-reverse-lookup")
            self.assertEqual(str(route), b1["arguments"]["shared_route"])

    def _run_profile(self, root: Path, *, bridge_case_run: bool = False) -> Path:
        index_out = root / "index"
        self.assertEqual(0, main(["retrieval-index-build",
                                  "--corpus", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                                  "--database-version", "plant_ref_test", "--out-dir", str(index_out)]))
        out = root / "triple"
        candidates = _write(root, "candidates.tsv", CANDIDATES)
        payload_rows = pd.read_csv(
            _write(root, "payload.tsv", PAYLOAD), sep="\t", dtype=str
        ).fillna("").to_dict("records")
        bound_claims = write_route_a_bound_claims(root, payload_rows)
        args = [
            "guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
            "--case-id", "case_fas", "--focal-species", "Pisum sativum",
            "--phenotype", "fasciated stem with enlarged meristem",
            "--mechanism-payload", str(root / "payload.tsv"), "--payload-source", "agent",
            "--phenotype-localization-payload", str(_localization(root, "fasciated stem")),
            "--record-cache", str(_record_cache(root)), "--source-cache", str(_source_cache(root)),
            "--bound-claims", str(bound_claims),
            "--retrieval-index", str(index_out / "plant_ref_test" / "index_manifest.json"),
            "--candidates", str(candidates),
            "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
            "--reference-function", str(root / "corpus.tsv"),
        ]
        if bridge_case_run:
            case = root / "case.yaml"
            case.write_text(yaml.safe_dump({
                "case_id": "case_fas",
                "trait": {"name": "fasciated stem"},
                "source": {"source_id": "fixture", "target_interval": "chr5:1-20000"},
                "target": {"assembly": "fixture"},
                "inputs": {"candidates_resource": "candidates"},
                "ranking_policy": "deterministic_v1",
            }, sort_keys=False), encoding="utf-8")
            resources = root / "resources.yaml"
            resources.write_text(yaml.safe_dump({
                "resources": {"candidates": {"type": "candidates_tsv", "path": "candidates.tsv"}}
            }, sort_keys=False), encoding="utf-8")
            args.extend(["--case", str(case), "--resources", str(resources)])
        self.assertEqual(0, main(args))
        return out

    def test_triple_route_bridge_executes_case_run_and_preserves_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root, bridge_case_run=True)
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            case_record = profile["execution_results"][-1]
            self.assertEqual("case-run", case_record["command"])
            self.assertEqual("completed", case_record["status"])
            case_out = out / "z03_case_run"
            provider = json.loads((case_out / "narrowing_provider_state.json").read_text(encoding="utf-8"))
            self.assertEqual("evaluated", provider["status"])
            self.assertTrue((case_out / "stage_a" / "deterministic_shortlist.tsv").exists())
            self.assertTrue((case_out / "stage_a_narrowed" / "deterministic_shortlist.tsv").exists())
            options = pd.read_csv(case_out / "next_experiment" / "next_experiment_options.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(options["input_facts"].str.contains("integrated_narrowing_status=evaluated").all())
            baseline = root / "baseline_case"
            self.assertEqual(0, main([
                "case-run", "--case", str(root / "case.yaml"), "--resources", str(root / "resources.yaml"),
                "--out-dir", str(baseline),
            ]))
            self.assertEqual(
                (baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                (case_out / "rank" / "ranked_candidates.tsv").read_bytes(),
            )
            # The canonical triple-route comparison reaches the same causal
            # scenario ledger as dual-route input. Route statuses travel with
            # it, so this is an evaluable classification rather than a legacy
            # inference or a silent B-only/neither fallback.
            run_manifest = json.loads((case_out / "run_manifest.json").read_text(encoding="utf-8"))
            ledger_status = run_manifest["scenario_ledger_status"]
            self.assertEqual("present", ledger_status["discordance_input"])
            self.assertGreater(ledger_status["discordance_rows"], 0)
            self.assertEqual("evaluated", ledger_status["route_a_status"])
            self.assertEqual("evaluated", ledger_status["route_b_status"])
            ledger = pd.read_csv(case_out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(ledger.source_ids.str.contains("comparison=cmp_").any())
            self.assertTrue(ledger.rank_contribution.eq("none").all())

    def test_every_lane_runs_and_the_routes_meet_only_at_the_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root)
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            executed = {record["command"]: record["status"] for record in profile["execution_results"]}
            for command in ("phenotype-investigate", "mechanism-import", "reference-gene-resolve",
                            "source-verify", "hypothesis-verify",
                            "bait-panel-freeze", "phenotype-first-status", "retrieval-run", "retrieval-import",
                            "route-b-prepare", "interval-reverse-lookup", "reverse-literature",
                            "retrieval-interval-match", "triple-route-compare"):
                self.assertEqual("completed", executed[command], command)
            self.assertIn("route_independence", profile)
            self.assertEqual("triple-route-compare", profile["route_independence"]["meeting_point"].split()[0])

            # Route C Phase 1 never saw the interval, and the freeze precedes the match
            phase1 = json.loads((out / "c01_retrieval_run" / "route_c_local_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertFalse(phase1["interval_information_in_model_input"])
            self.assertEqual([], phase1["interval_arguments"])
            phase2 = json.loads((out / "c03_interval_match" / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertTrue(phase2["retrieval_frozen_before_interval"])
            # Route B ran with no phenotype at all
            route_b = json.loads((out / "b01_route_b" / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(route_b["phenotype_used"])
            self.assertEqual("source_backed_reverse_lookup", route_b["route_b_readiness"])
            # B1 overlay ran against the frozen B0 atlas and did not rebuild it
            b1 = json.loads((out / "b02_route_b1_overlay" / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual("b1_phenotype_overlay", b1["lane_variant"])
            self.assertTrue(b1["phenotype_used"])
            self.assertEqual(route_b["frozen_output_fingerprints"]["interval_gene_atlas"]["sha256"],
                             b1["inputs"]["atlas"]["sha256"])
            # reverse-literature chained from B0's lookup output; offline without a
            # cache it records not_evaluated_offline and still completes
            revlit = json.loads((out / "b03_reverse_literature" / "reverse_literature_manifest.json")
                                .read_text(encoding="utf-8"))
            self.assertEqual("none", revlit["provider"])
            self.assertFalse(revlit["external_api_allowed"])
            self.assertEqual(1, revlit["not_evaluated_offline"])
            # the comparison carries the four-state A/B label and the literature ledger.
            # In this fixture B0's interval hypothesis is the OG0002 orthogroup, so CLV
            # is A+C only and OG0002 is the B-only row.
            comparison = pd.read_csv(out / "z01_triple_route_compare" / "triple_route_hypothesis_comparison.tsv",
                                     sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertIn("four_state", comparison.columns)
            self.assertIn("literature_link_candidate_count", comparison.columns)
            self.assertEqual("a_only_supported", comparison.loc["CLV", "four_state"])
            self.assertEqual("b_only_supported", comparison.loc["OG0002", "four_state"])
            for family in ("PIN", "CESA"):
                self.assertEqual("a_only_supported", comparison.loc[family, "four_state"])
            # offline without a cache, reverse-literature writes no candidate rows --
            # only the status ledger records not_evaluated_offline. An absent search
            # is not an absence of literature, so the count is 0, never a miss.
            self.assertEqual("0", comparison.loc["CLV", "literature_link_candidate_count"])
            self.assertEqual("0", comparison.loc["PIN", "literature_link_candidate_count"])
            convergence = json.loads((out / "z01_triple_route_compare" / "triple_route_convergence.json")
                                     .read_text(encoding="utf-8"))
            self.assertEqual(1, convergence["route_b_literature"]["not_evaluated_offline"])
            self.assertEqual(0, convergence["route_b_literature"]["queried"])
            self.assertIn("never added to a score", convergence["route_b_literature"]["note"])

    def test_failure_modes_are_accounted_per_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root)
            comparison = pd.read_csv(out / "z01_triple_route_compare" / "triple_route_hypothesis_comparison.tsv",
                                     sep="\t", dtype=str).fillna("").set_index("family_id")
            convergence = json.loads((out / "z01_triple_route_compare" / "triple_route_convergence.json")
                                     .read_text(encoding="utf-8"))

            # CLV: reached by the phenotype-first route and by retrieval, and the frozen
            # retrieval matched a real interval gene -> A's hallucination mode is checked
            # off by interval anchoring and C's similarity mode by a source-backed mechanism.
            clv = comparison.loc["CLV"]
            self.assertEqual("a_c_convergent", clv["convergence_status"])
            self.assertEqual("True", clv["interval_anchored"])
            self.assertIn("route_c_interval_match", clv["interval_anchor_basis"])
            self.assertEqual("all_applicable_modes_neutralised", clv["failure_mode_coverage"])
            self.assertIn("route_a_hallucinated_gene_absent_from_interval", clv["neutralised_failure_modes"])
            self.assertEqual("", clv["residual_failure_modes"])

            # OG0002: the interval-first route reached the same orthogroup, and an
            # interval-blind route (C) reached it independently -> B's post-hoc
            # rationalisation mode is checked off.
            og = comparison.loc["OG0002"]
            self.assertEqual("b_c_convergent", og["convergence_status"])
            self.assertEqual("True", og["interval_blind_support"])
            self.assertIn("route_b_post_hoc_rationalisation", og["neutralised_failure_modes"])

            # PIN / CESA were proposed from the phenotype alone and are not in the
            # interval: nothing cancels Route A's characteristic error, and the row says so.
            for family in ("PIN", "CESA"):
                row = comparison.loc[family]
                self.assertEqual("phenotype_only_supported", row["convergence_status"])
                self.assertEqual("False", row["interval_anchored"])
                self.assertEqual("none_neutralised", row["failure_mode_coverage"])
                self.assertEqual("route_a_hallucinated_gene_absent_from_interval", row["residual_failure_modes"])

            cancellation = convergence["failure_mode_cancellation"]
            self.assertTrue(cancellation["route_c_interval_match_supplied"])
            self.assertIn("CLV", cancellation["fully_neutralised_families"])
            self.assertEqual({"none"}, set(comparison["rank_contribution"]))
            self.assertIn("adds no score", cancellation["residual"])

    def test_agent_authored_literature_is_verified_before_it_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root)
            manifest = json.loads((out / "a03_source_verify" / "source_verification_manifest.json")
                                  .read_text(encoding="utf-8"))
            self.assertEqual("local_record_cache", manifest["provider"])
            self.assertEqual(0, manifest["rejected_citation_count"])
            # Three citations resolve and support the mechanism claims: one for
            # CLV1, plus two for CESA9. The PIN mechanism's two papers resolve
            # but do not mention the claimed family, so they remain review-only.
            self.assertEqual(3, manifest["verified_citation_count"])
            self.assertEqual(2, manifest["review_required_citation_count"])
            validation = pd.read_csv(out / "a03_source_verify" / "source_validation.tsv", sep="\t", dtype=str).fillna("")
            unsupported = validation[validation["claim_support"] == "unsupported"]
            self.assertEqual(2, len(unsupported))
            self.assertEqual({"mech_002"}, set(unsupported["hypothesis_id"]))
            self.assertTrue(unsupported["review_reason"].eq(
                "citation_resolves_but_does_not_mention_the_claimed_subject").all())
            verified = validation[validation["verification_status"] == "verified"]
            self.assertEqual({("mech_001", "PMID:8402901"),
                              ("mech_003", "PMID:12172025"),
                              ("mech_003", "PMID:8402901")},
                             set(zip(verified["hypothesis_id"], verified["citation"])))
            self.assertTrue(verified["claim_support_status"].eq("supported_direct_genetic_test").all())
            status = pd.read_csv(out / "a03_source_verify" / "hypothesis_source_status.tsv", sep="\t", dtype=str) \
                .fillna("").set_index("hypothesis_id")
            self.assertEqual("PMID:8402901", status.loc["mech_001", "verified_source_ids"])
            self.assertEqual({"PMID:12172025", "PMID:8402901"},
                             set(status.loc["mech_003", "verified_source_ids"].split(";")))
            self.assertEqual(
                {"all_verified", "none_verified", "no_literature_citation"},
                set(status["source_verification_status"]),
            )
            mechanisms = pd.read_csv(out / "a02_mechanism_import" / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"agent"}, set(mechanisms["payload_source"]))
            outer = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("partial", outer["overall_status"])
            self.assertFalse(outer["final"])
            self.assertEqual("completed", outer["legacy_profile"]["overall_status"])
            canonical = outer["canonical_execution"]
            self.assertEqual("semantic_review_required", canonical["state"])
            self.assertTrue(canonical["required_work_remaining"])
            self.assertIsNone(canonical["wet_candidates"])
            self.assertTrue(canonical["manifest"]["sha256"])
            self.assertTrue(canonical["expected_semantic_run_dir"])


    def test_dry_run_plans_the_route_a_interval_family_match_stage(self) -> None:
        """Issue #355: a05b sits between the bait freeze and the Route A status,
        and a06 chains its falsification output as an optional input."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--phenotype", "fasciated stem"]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            stages = {stage["command"]: stage for stage in manifest["stages"]}
            self.assertIn("route-a-interval-match", stages)
            match_stage = stages["route-a-interval-match"]
            self.assertEqual("a05b_interval_family_match", Path(match_stage["arguments"]["out_dir"]).name)
            self.assertEqual(("bait-panel-freeze", "panel"),
                             (match_stage["inputs_from"]["panel"]["stage"], match_stage["inputs_from"]["panel"]["result_key"]))
            self.assertEqual(("mechanism-coverage-audit", "hypotheses"),
                             (match_stage["inputs_from"]["hypotheses"]["stage"], match_stage["inputs_from"]["hypotheses"]["result_key"]))
            self.assertIn("verified targets only", match_stage["stop_condition"])
            # a06 receives the falsification table, but never blocks on it
            a06 = stages["phenotype-first-status"]
            self.assertEqual(("route-a-interval-match", "hypothesis_interval_falsification"),
                             (a06["inputs_from"]["falsification"]["stage"], a06["inputs_from"]["falsification"]["result_key"]))
            self.assertIn("falsification", a06["optional_chain_inputs"])
            # the stage sits between a05 and a06 in the plan order
            order = [stage["command"] for stage in manifest["stages"]]
            self.assertLess(order.index("bait-panel-freeze"), order.index("route-a-interval-match"))
            self.assertLess(order.index("route-a-interval-match"), order.index("phenotype-first-status"))

    def _family_index_fixture(self, root: Path) -> Path:
        hits = _write(root, "fi_hits.tsv", "\n".join([
            _hit13("Psat_CLV1", "sp|P30924|CLV1_ARATH", 85.0, 400.0,
                   "Receptor kinase CLV1 OS=Arabidopsis thaliana OX=3702 GN=CLV1 PE=1 SV=2"),
            _hit13("Psat_UNK", "sp|P12345|UNK_ARATH", 50.0, 100.0,
                   "Hypothetical protein OS=Arabidopsis thaliana OX=3702 PE=4 SV=1"),
        ]) + "\n")
        genes = _write(root, "fi_genes.tsv", "gene_id\tchrom\tstart\tend\n"
                                              "Psat_CLV1\tchr5\t1000\t3000\n"
                                              "Psat_UNK\tchr5\t9000\t11000\n")
        index_out = root / "family_index"
        self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                  "--out-dir", str(index_out), "--reference-label", "TestRef"]))
        return index_out / "family_index.tsv"

    def test_route_a_interval_match_executes_without_an_index_as_not_evaluated(self) -> None:
        """Missing provider degrades the new stage to not_evaluated and the
        whole profile still completes safely (no family index supplied)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch.dict(os.environ, {
                "CANDIDATE_QTL_FAMILY_INDEX": "",
                "XDG_CACHE_HOME": str(root / "empty-cache"),
            }):
                out = self._run_profile(root)
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            records = {record["command"]: record for record in profile["execution_results"]}
            self.assertEqual("completed", records["route-a-interval-match"]["status"])
            self.assertEqual("completed", profile["legacy_profile"]["overall_status"])
            manifest = json.loads((out / "a05b_interval_family_match" / "route_a_interval_family_match_manifest.json")
                                  .read_text(encoding="utf-8"))
            self.assertFalse(manifest["evaluated"])
            self.assertEqual("not_evaluated_missing_reference_genome_family_index", manifest["not_evaluated_reason"])
            status = pd.read_csv(out / "a05b_interval_family_match" / "route_a_interval_family_match.tsv",
                                 sep="\t", dtype=str).fillna("")
            self.assertTrue(len(status))
            self.assertTrue(status["family_status"].eq("not_evaluated_missing_reference_genome_family_index").all())
            self.assertFalse((out / "a05b_interval_family_match" / "hypothesis_interval_falsification.tsv").exists())
            # a06 ran exactly as before, with the falsification input declared absent
            a06_status = json.loads((out / "a06_route_a_status" / "phenotype_first_status.json").read_text(encoding="utf-8"))
            self.assertIsNone(a06_status["inputs"]["falsification"])
            self.assertIsNotNone(a06_status["inputs"]["panel"])

    def test_route_a_interval_match_reaches_concrete_interval_genes_end_to_end(self) -> None:
        """Issue #355 acceptance: verified bait -> family -> interval genes ->
        phenotype-first status -> compare handoff -> family-to-copy input."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._run_profile(root)
            index = self._family_index_fixture(root)
            out2 = root / "triple2"
            args = [
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out2), "--execute",
                "--case-id", "case_fas", "--focal-species", "Pisum sativum",
                "--phenotype", "fasciated stem with enlarged meristem",
                "--mechanism-payload", str(root / "payload.tsv"), "--payload-source", "agent",
                "--phenotype-localization-payload", str(root / "localization.json"),
                "--record-cache", str(root / "records"), "--source-cache", str(root / "pmc"),
                "--bound-claims", str(root / "bound_claims.tsv"),
                "--retrieval-index", str(root / "index" / "plant_ref_test" / "index_manifest.json"),
                "--candidates", str(root / "candidates.tsv"),
                "--orthology", str(root / "orthology.tsv"),
                "--reference-function", str(root / "corpus.tsv"),
                "--family-index", str(index),
                "--interval", "chr5:1-20000",
            ]
            self.assertEqual(0, main(args))
            profile = json.loads((out2 / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            records = {record["command"]: record for record in profile["execution_results"]}
            for command in ("route-a-interval-match", "phenotype-first-status", "triple-route-compare",
                            "family-to-copy"):
                self.assertEqual("completed", records[command]["status"], command)
            self.assertEqual("completed", profile["legacy_profile"]["overall_status"])

            # a05b: the verified CLV bait forward-matches a concrete interval gene
            match_status = pd.read_csv(out2 / "a05b_interval_family_match" / "route_a_interval_family_match.tsv",
                                       sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertIn("CLV", match_status.index)
            clv = match_status.loc["CLV"]
            self.assertEqual("concrete_interval_gene_matched", clv["family_status"])
            self.assertEqual("Psat_CLV1", clv["interval_gene_ids"])
            self.assertEqual("full", clv["match_level"])
            self.assertEqual("True", clv["evaluable"])
            self.assertEqual("none", clv["rank_contribution"])
            self.assertTrue(clv["source_ids"])
            # the falsification verdicts flowed into the a06 Route A table
            a06 = pd.read_csv(out2 / "a06_route_a_status" / "phenotype_first_hypotheses.tsv",
                              sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertEqual("Psat_CLV1", a06.loc["CLV", "member_interval_genes"])
            self.assertEqual("interval_anchored", a06.loc["CLV", "anchor_status"])
            self.assertEqual("verified_bait", a06.loc["CLV", "source_backing"])
            # the comparison handoff carries family + concrete gene + provenance
            handoff = pd.read_csv(out2 / "z01_triple_route_compare" / "family_to_copy_input.tsv",
                                  sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertEqual("Psat_CLV1", handoff.loc["CLV", "member_interval_genes"])
            # family-to-copy received the same member through the atlas
            members = pd.read_csv(out2 / "z02_family_to_copy" / "family_copy_members.tsv",
                                  sep="\t", dtype=str).fillna("")
            self.assertIn("Psat_CLV1", set(members["candidate_gene_id"]))
            # the new stage wrote no rank artifact and nothing changed ranking
            self.assertFalse((out2 / "a05b_interval_family_match" / "ranked_candidates.tsv").exists())
            self.assertFalse((out2 / "a05b_interval_family_match" / "rank").exists())
            for name in ("route_a_interval_family_match.tsv",):
                frame = pd.read_csv(out2 / "a05b_interval_family_match" / name, sep="\t", dtype=str).fillna("")
                self.assertTrue(frame["rank_contribution"].eq("none").all(), name)
            comparison = pd.read_csv(out2 / "z01_triple_route_compare" / "triple_route_hypothesis_comparison.tsv",
                                     sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertTrue(comparison["rank_contribution"].eq("none").all())

    def test_rerun_in_the_same_out_dir_does_not_reuse_stale_a05b_artifacts(self) -> None:
        """A success rerun in the same out-dir must not keep a stale
        falsification artifact from an earlier run (issue #355, #302)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple_rerun"
            base = [
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
                "--case-id", "case_fas", "--focal-species", "Pisum sativum",
                "--phenotype", "fasciated stem with enlarged meristem",
                "--mechanism-payload", str(root / "payload.tsv"), "--payload-source", "agent",
                "--phenotype-localization-payload", str(root / "localization.json"),
                "--record-cache", str(root / "records"), "--source-cache", str(root / "pmc"),
                "--bound-claims", str(root / "bound_claims.tsv"),
                "--retrieval-index", str(root / "index" / "plant_ref_test" / "index_manifest.json"),
                "--candidates", str(root / "candidates.tsv"),
                "--orthology", str(root / "orthology.tsv"),
                "--reference-function", str(root / "corpus.tsv"),
            ]
            # first run: no family index -> no falsification artifact
            self._run_profile(root)
            with patch.dict(os.environ, {
                "CANDIDATE_QTL_FAMILY_INDEX": "",
                "XDG_CACHE_HOME": str(root / "empty-cache"),
            }):
                self.assertEqual(0, main(base + ["--interval", "chr5:1-20000"]))
            self.assertFalse((out / "a05b_interval_family_match" / "hypothesis_interval_falsification.tsv").exists())
            # second run in the SAME out-dir with a family index: the stale
            # 'not evaluated' state must be replaced by a real falsification
            index = self._family_index_fixture(root)
            self.assertEqual(0, main(base + ["--family-index", str(index), "--interval", "chr5:1-20000"]))
            self.assertTrue((out / "a05b_interval_family_match" / "hypothesis_interval_falsification.tsv").exists())
            status = pd.read_csv(out / "a05b_interval_family_match" / "route_a_interval_family_match.tsv",
                                 sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertEqual("concrete_interval_gene_matched", status.loc["CLV", "family_status"])
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["legacy_profile"]["final"])
            self.assertTrue(manifest["legacy_profile"]["overall_status"] == "completed")

    def test_dry_run_plans_b1_overlay_and_reverse_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--phenotype", "fasciated stem"]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            stages = {stage["command"]: stage for stage in manifest["stages"]}
            b1 = stages["interval-reverse-lookup"]
            self.assertEqual(["phenotype"], b1["required_inputs"])
            self.assertEqual([], b1["missing_inputs"])
            self.assertEqual("ready", b1["status"])
            self.assertEqual(("route-b-prepare", "interval_gene_atlas"),
                             (b1["inputs_from"]["atlas"]["stage"], b1["inputs_from"]["atlas"]["result_key"]))
            self.assertEqual(("route-b-prepare", "reference_function"),
                             (b1["inputs_from"]["reference_function"]["stage"], b1["inputs_from"]["reference_function"]["result_key"]))
            self.assertIn("atlas is not modified", b1["stop_condition"])
            self.assertIn("<from route-b-prepare.interval_gene_atlas>", b1["arguments"]["atlas"])
            revlit = stages["reverse-literature"]
            self.assertEqual(("route-b-prepare", "interval_reverse_lookup"),
                             (revlit["inputs_from"]["lookup"]["stage"], revlit["inputs_from"]["lookup"]["result_key"]))
            self.assertIn("not_evaluated_offline", revlit["stop_condition"])
            self.assertIn("none (offline cache only)", revlit["external_transmission"])
            compare = stages["triple-route-compare"]
            self.assertEqual(("reverse-literature", "literature"),
                             (compare["inputs_from"]["literature"]["stage"], compare["inputs_from"]["literature"]["result_key"]))
            self.assertEqual(("reverse-literature", "status"),
                             (compare["inputs_from"]["literature_status"]["stage"], compare["inputs_from"]["literature_status"]["result_key"]))
            self.assertIn("literature", compare["optional_chain_inputs"])
            # a phenotype-less plan marks B1 not runnable (missing input), never failed
            out2 = root / "triple_no_phenotype"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out2)]))
            manifest2 = json.loads((out2 / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            b1_no_phenotype = next(stage for stage in manifest2["stages"] if stage["command"] == "interval-reverse-lookup")
            self.assertEqual("not_runnable_missing_input", b1_no_phenotype["status"])
            self.assertEqual(["phenotype"], b1_no_phenotype["missing_inputs"])

    def test_route_a_external_verification_does_not_enable_legacy_b1_api(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_a_external_only"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                "--phenotype", "fasciated stem", "--allow-external-api",
            ]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            stages = {stage["command"]: stage for stage in manifest["stages"]}
            self.assertTrue(stages["source-verify"]["arguments"]["allow_external_api"])
            self.assertTrue(stages["hypothesis-verify"]["arguments"]["allow_external_api"])
            self.assertFalse(stages["reverse-literature"]["arguments"]["allow_external_api"])
            self.assertIn("none (offline cache only)", stages["reverse-literature"]["external_transmission"])

            opted_in = root / "route_a_and_b1_external"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(opted_in),
                "--phenotype", "fasciated stem", "--allow-external-api",
                "--allow-reverse-literature-api",
            ]))
            opted_manifest = json.loads(
                (opted_in / "guided_profile_manifest.json").read_text(encoding="utf-8")
            )
            opted_stages = {stage["command"]: stage for stage in opted_manifest["stages"]}
            self.assertTrue(opted_stages["reverse-literature"]["arguments"]["allow_external_api"])
            self.assertIn("Europe PMC", opted_stages["reverse-literature"]["external_transmission"])

    def test_compare_reports_four_state_and_literature_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_001\tCLV\tG_CLV1\tverified_bait\n"
                             "hyp_002\tWUS\t\tverified_bait\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "rev_001\tCLV\tG_CLV1\tsource_backed\n"
                             "rev_002\tPXY\tG_PXY\tsource_backed\n")
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\t"
                      "matched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank")
            route_c = _write(root, "route_c.tsv", "\n".join([
                header,
                "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1",
                "AtPXY\tArabidopsis thaliana\tPXY\tOG0002\t0.8\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t2",
                "AtFT\tArabidopsis thaliana\tFT\tOG0003\t0.7\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t3",
            ]) + "\n")
            literature = _write(root, "literature.tsv", "\n".join([
                "interval_gene_ids\treference_gene\tquery_status",
                "G_CLV1\tAtCLV1\tnot_evaluated_offline",
                "G_PXY\tAtPXY\tno_hit",
                "G_UNRELATED\tAtFT\tretrieved",
            ]) + "\n")
            lit_status = _write(root, "lit_status.tsv", "\n".join([
                "reference_gene\tquery_status",
                "AtCLV1\tnot_evaluated_offline",
                "AtPXY\tno_hit",
                "AtFT\tretrieved",
            ]) + "\n")
            out = root / "compare"
            triple_route_compare(str(route_a), str(route_b), str(route_c), str(out), case_id="demo",
                                 literature_arg=str(literature), literature_status_arg=str(lit_status))
            comparison = pd.read_csv(out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertEqual("a_and_b_supported", comparison.loc["CLV", "four_state"])
            self.assertEqual("a_only_supported", comparison.loc["WUS", "four_state"])
            self.assertEqual("b_only_supported", comparison.loc["PXY", "four_state"])
            self.assertEqual("neither_supported", comparison.loc["FT", "four_state"])
            self.assertEqual("1", comparison.loc["CLV", "literature_link_candidate_count"])
            self.assertEqual("1", comparison.loc["PXY", "literature_link_candidate_count"])
            self.assertEqual("0", comparison.loc["WUS", "literature_link_candidate_count"])
            self.assertEqual("0", comparison.loc["FT", "literature_link_candidate_count"])
            self.assertTrue(comparison["rank_contribution"].eq("none").all())
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            lit = convergence["route_b_literature"]
            self.assertEqual(1, lit["queried"])
            self.assertEqual(1, lit["no_hit"])
            self.assertEqual(1, lit["not_evaluated_offline"])
            self.assertIn("never added to a score", lit["note"])
            self.assertIn("agreement_state", comparison.columns)

    def test_plan_ends_with_family_to_copy_chained_to_handoff_and_atlas(self) -> None:
        """Issue #284: the triple-route plan gains one post-meeting stage that
        reads the canonical handoff and the frozen B0 atlas only."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "triple"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--phenotype", "fasciated stem"]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            stages = manifest["stages"]
            self.assertEqual("family-to-copy", stages[-1]["command"])
            copy_stage = stages[-1]
            self.assertEqual("z02_family_to_copy", Path(copy_stage["arguments"]["out_dir"]).name)
            self.assertEqual(("triple-route-compare", "family_to_copy_input"),
                             (copy_stage["inputs_from"]["families"]["stage"], copy_stage["inputs_from"]["families"]["result_key"]))
            self.assertEqual(("route-b-prepare", "interval_gene_atlas"),
                             (copy_stage["inputs_from"]["atlas"]["stage"], copy_stage["inputs_from"]["atlas"]["result_key"]))
            self.assertIn("<from triple-route-compare.family_to_copy_input>", copy_stage["arguments"]["families"])
            self.assertIn("score unchanged", copy_stage["stop_condition"])
            self.assertEqual(1, [stage["command"] for stage in stages].count("family-to-copy"))

    def test_execute_chains_handoff_and_atlas_into_family_to_copy(self) -> None:
        """Issue #284: --execute reaches family-to-copy with the handoff + atlas,
        unresolved families stay unresolved, and no score columns appear."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root)
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            records = {record["command"]: record for record in manifest["execution_results"]}
            self.assertEqual("completed", records["family-to-copy"]["status"])
            self.assertEqual("completed", manifest["legacy_profile"]["overall_status"])
            copy_manifest = json.loads((out / "z02_family_to_copy" / "family_copy_manifest.json").read_text(encoding="utf-8"))
            # Chained from the comparison handoff and the frozen B0 atlas; the
            # other evidence inputs were never fabricated.
            self.assertEqual(str(out / "z01_triple_route_compare" / "family_to_copy_input.tsv"),
                             copy_manifest["inputs"]["families"]["path"])
            self.assertEqual(str(out / "b01_route_b" / "interval_gene_atlas.tsv"),
                             copy_manifest["inputs"]["atlas"]["path"])
            self.assertIsNone(copy_manifest["inputs"]["orthology"])
            self.assertIsNone(copy_manifest["inputs"]["synteny"])
            self.assertIsNone(copy_manifest["inputs"]["expression"])
            self.assertIsNone(copy_manifest["inputs"]["target"])
            for name in ("family_copy_discrimination.tsv", "family_copy_shortlist.tsv", "family_copy_members.tsv"):
                columns = set(pd.read_csv(out / "z02_family_to_copy" / name, sep="\t", dtype=str).columns)
                self.assertFalse({"dry_score", "score", "dry_group"} & columns, name)
            # The handoff itself is untouched by the narrowing stage.
            handoff = pd.read_csv(out / "z01_triple_route_compare" / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(handoff["gene_family"].isin(copy_manifest["selection_basis"]).all())

    def test_family_to_copy_without_upstream_handoff_is_unavailable_not_failed(self) -> None:
        """Issue #284: when the comparison produced no handoff, the post stage is
        recorded as upstream unavailable and the profile is not a fake failure."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage = {
                "command": "family-to-copy",
                "arguments": {"families": "<from triple-route-compare.family_to_copy_input>",
                              "atlas": "<from route-b-prepare.interval_gene_atlas>", "out_dir": str(root / "z02")},
                "inputs_from": {"families": {"stage": "triple-route-compare", "result_key": "family_to_copy_input"},
                                "atlas": {"stage": "route-b-prepare", "result_key": "interval_gene_atlas"}},
            }
            # Comparison completed but produced no handoff: output_missing, not a failure.
            results = [{"command": "triple-route-compare", "status": "completed",
                        "result": {"triple_route_hypothesis_comparison": "/tmp/x.tsv"}}]
            self.assertEqual("not_runnable_upstream_output_missing:triple-route-compare.family_to_copy_input",
                             _resolve_chained(stage, results))
            # Comparison never ran: upstream unavailable, still not a failure.
            results = [{"command": "triple-route-compare", "status": "not_runnable_missing_input"}]
            self.assertEqual("not_runnable_upstream_unavailable:triple-route-compare",
                             _resolve_chained(stage, results))

    def test_lane_isolation_refuses_b_stages_reading_route_a_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a_out = root / "a03_source_verify"
            for command, input_name in (("interval-reverse-lookup", "atlas"), ("reverse-literature", "lookup")):
                stage = {"command": command, "arguments": {
                    input_name: str(a_out / "source_validation.tsv"),
                    "out_dir": str(root / "b02_route_b1_overlay"),
                }}
                with self.assertRaises(ValueError):
                    _assert_lane_isolation(stage, root)

    def test_route_a_interval_match_stays_in_lane_a_and_refuses_route_b_inputs(self) -> None:
        """The new stage is lane a: it may chain the a05 panel, but it must
        never read a lane-b artifact before the comparison (issue #355)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage = {"command": "route-a-interval-match", "arguments": {
                "panel": str(root / "a05_bait_panel_freeze" / "verified_bait_panel.tsv"),
                "out_dir": str(root / "a05b_interval_family_match"),
            }}
            _assert_lane_isolation(stage, root)  # same lane: allowed
            stage["arguments"]["family_index"] = str(root / "b01_route_b" / "family_index.tsv")
            with self.assertRaises(ValueError):
                _assert_lane_isolation(stage, root)


    def test_triple_route_legacy_7col_payload_e2e(self) -> None:
        """Issue #380: a frozen legacy 7-col payload converts through the adapter
        (rows_imported > 0, adapter map written), but the new exact-12 coverage
        gate blocks both old three-row portfolios before verification."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index_out = root / "index"
            self.assertEqual(0, main(["retrieval-index-build",
                                      "--corpus", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                                      "--database-version", "plant_ref_test", "--out-dir", str(index_out)]))
            out = root / "triple_legacy"
            self._run_payload_profile(root, _write(root, "legacy_payload.tsv", LEGACY_PAYLOAD),
                                      out, index_out, expected_exit=1)
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("partial", profile["overall_status"])
            preflight = profile["payload_schema_preflight"]
            self.assertEqual("accepted", preflight["status"])
            self.assertEqual("legacy_v1_7col", preflight["accepted_schema_version"])
            self.assertEqual(3, preflight["row_count"])
            self.assertNotIn("freeze_verification", preflight)
            a02 = json.loads((out / "a02_mechanism_import" / "mechanism_import_manifest.json")
                             .read_text(encoding="utf-8"))
            self.assertGreater(a02["rows_imported"], 0)
            self.assertEqual(3, a02["hypothesis_count"])
            adapter_map = pd.read_csv(out / "a02_mechanism_import" / "mechanism_payload_adapter_map.tsv",
                                      sep="\t", dtype=str).fillna("")
            self.assertEqual({"pea_a_M1", "pea_a_M2", "pea_a_M3"}, set(adapter_map["mechanism_id"]))
            self.assertTrue(adapter_map["mechanism_name"].str.len().gt(0).all())
            coverage = next(row for row in profile["execution_results"]
                            if row["command"] == "mechanism-coverage-audit")
            self.assertEqual("failed", coverage["status"])
            self.assertIn("one distinct family per row", coverage["error"])
            self.assertFalse((out / "a04_hypothesis_verify").exists())

            out2 = root / "triple_canonical"
            self._run_payload_profile(root, _write(root, "canonical_payload.tsv", CANONICAL_MIXED_PAYLOAD),
                                      out2, index_out, expected_exit=1)
            profile2 = json.loads((out2 / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("partial", profile2["overall_status"])
            self.assertEqual("canonical_subset", profile2["payload_schema_preflight"]["accepted_schema_version"])
            coverage2 = next(row for row in profile2["execution_results"]
                             if row["command"] == "mechanism-coverage-audit")
            self.assertEqual("failed", coverage2["status"])
            self.assertIn("invalid or empty causal_layer", coverage2["error"])
            a05b2 = json.loads((out2 / "a05b_interval_family_match" / "route_a_interval_family_match_manifest.json")
                               .read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_no_verified_bait", a05b2["not_evaluated_reason"])

    def test_triple_route_claim_extractions_reach_source_verify(self) -> None:
        """Issue #611: agent-extracted claim spans in the payload ride through
        mechanism-import into the deterministic extraction gate at a03. Offline
        the citation cannot resolve, so the gate must still reject it — wiring
        an extraction never bypasses verification."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index_out = root / "index"
            self.assertEqual(0, main(["retrieval-index-build",
                                      "--corpus", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                                      "--database-version", "plant_ref_test", "--out-dir", str(index_out)]))
            # Keep this wiring test on the valid 12-hypothesis contract. A
            # three-row legacy portfolio is now deliberately stopped by the
            # coverage gate before source verification.
            rows = [dict(zip(PAYLOAD.splitlines()[0].split("\t"), line.split("\t")))
                    for line in PAYLOAD.splitlines()[1:]]
            payload = {"hypotheses": rows,
                       "claim_extractions": [{"citation": "PMID:8402901", "claim_type": "causal",
                                              "evidence_relation": "direct_genetic_test",
                                              "quoted_support": "loss of CLV signalling enlarges the stem cell pool",
                                              "locator": "title"}]}
            out = root / "triple_extractions"
            self._run_payload_profile(root, _write(root, "payload_with_extractions.json", json.dumps(payload)),
                                      out, index_out)
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            a02 = next(item for item in profile["execution_results"]
                       if item["command"] == "mechanism-import" and item.get("status") == "completed")
            self.assertTrue(a02, "mechanism-import did not complete")
            import_manifest = json.loads((out / "a02_mechanism_import" / "mechanism_import_manifest.json")
                                         .read_text(encoding="utf-8"))
            self.assertTrue(import_manifest["claim_extractions"]["supplied"])
            a03_manifest = json.loads((out / "a03_source_verify" / "source_verification_manifest.json")
                                      .read_text(encoding="utf-8"))
            self.assertTrue(a03_manifest["claim_extraction"]["supplied"])
            self.assertEqual(1, a03_manifest["claim_extraction"]["row_count"])
            gate = pd.read_csv(out / "a03_source_verify" / "claim_extraction_gate.tsv",
                               sep="\t", dtype=str).fillna("")
            self.assertEqual("rejected_unknown_citation", gate.iloc[0]["extraction_status"])

    def test_triple_route_unknown_schema_aborts_before_stages(self) -> None:
        """Issue #380: an unknown payload schema is rejected at plan time; the
        profile aborts before any stage runs and records the rejection with the
        offending column names."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload_unknown.tsv",
                             "case_id\tmechanism_id\tunknown_schema_column\npea_a\tpea_a_M1\tx\n")
            out = root / "triple"
            self.assertEqual(1, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--execute", "--case-id", "pea_a", "--focal-species", "Pisum sativum",
                                      "--phenotype", "anthocyanin pigmentation loss",
                                      "--mechanism-payload", str(payload)]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", manifest["overall_status"])
            self.assertNotIn("execution_results", manifest)
            preflight = manifest["payload_schema_preflight"]
            self.assertEqual("rejected", preflight["status"])
            error = preflight["errors"][0]
            self.assertIn("unknown column(s)", error)
            self.assertIn("unknown_schema_column", error)
            self.assertFalse((out / "a01_phenotype_investigation").exists())
            self.assertFalse((out / "a02_mechanism_import").exists())

    def _run_payload_profile(
        self,
        root: Path,
        payload: Path,
        out: Path,
        index_out: Path,
        *,
        expected_exit: int = 0,
    ) -> None:
        args = [
            "guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
            "--case-id", "pea_a", "--focal-species", "Pisum sativum",
            "--phenotype", "anthocyanin pigmentation loss in flower and seed coat",
            "--mechanism-payload", str(payload), "--payload-source", "agent",
            "--phenotype-localization-payload", str(_localization(root, "anthocyanin pigmentation loss", "anthocyanin_localization.json")),
            "--retrieval-index", str(index_out / "plant_ref_test" / "index_manifest.json"),
            "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
            "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
            "--reference-function", str(root / "corpus.tsv"),
        ]
        self.assertEqual(expected_exit, main(args))


if __name__ == "__main__":
    unittest.main()
