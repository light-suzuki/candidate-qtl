"""Triple-route profile with the canonical family-candidate-bridge (issue #374).

The bridge normalizes Route A / B / C family evidence into one reporting table
between the comparison and family-to-copy, and family-to-copy consumes the
bridge mapping tier as a shortlist promotion gate (tier5/tier6 members are
never promoted). These E2E tests pin the wiring (plan order, chained inputs),
the executed chain (bridge -> family-to-copy -> integrated narrowing -> stage A
narrowed -> next experiment), byte-identical determinism, rank preservation
and the not_evaluated degradation vocabulary.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd
import yaml

from candidate_qtl.cli import main
from candidate_qtl.family_candidate_bridge import BRIDGE_COLUMNS, dependency_cluster_id
from tests.test_blast_evidence import MockBlastTools
from tests.bound_claim_fixtures import write_route_a_bound_claims
from tests.test_triple_route_profile import PAYLOAD, _localization

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
    "AtCLV1\tArabidopsis thaliana\tCLV\t\tCLAVATA1 receptor kinase restricting the stem cell pool\t"
    "fasciated stem and enlarged meristem\tmeristem_maintenance\tPMID:8402901",
    "AtSBE1\tArabidopsis thaliana\tSBE\t\tstarch branching enzyme I\twrinkled seed\tstarch_biosynthesis\tPMID:1978859",
]) + "\n"

# The paralog fixture: two candidate paralog families inside the interval, one
# source-backed (G_A) and one naming-only (G_W), plus a mixed family (PXY:
# G_R1 source-backed, G_R2 name-only) that cannot be narrowed to one copy.
PARALOG_CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "G_A\tchr1\t1000\t3000\tCLAVATA1 like receptor kinase",
    "G_W\tchr1\t20000\t22000\tWUSCHEL like homeobox",
    "G_R1\tchr1\t31000\t33000\tPXY like receptor",
    "G_R2\tchr1\t35000\t37000\tPXY like receptor (second copy)",
]) + "\n"

PARALOG_ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\trelationship\tpercent_identity\tquery_coverage\tprovider",
    "G_A\tAtCLV1\tone_to_one\t95\t0.9\tlocal_blastp:test",
    "G_R1\tAtPXY\tone_to_one\t90\t0.9\tlocal_blastp:test",
]) + "\n"

# The only homology layer in the paralog fixture is G_A's domain architecture;
# the PXY copies and the WOX copy stay naming-only, so the bridge marks them
# tier5 and they never reach the shortlist.
PARALOG_DOMAINS = "\n".join([
    "gene_id\tdomain_id\tdomain_name\tsource",
    "G_A\tPF00001\tLRR_receptor\tsupplied_table",
]) + "\n"

PARALOG_PAYLOAD = "\n".join([
    "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
    "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
    "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
    "meristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
    "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
    "P30924\tArabidopsis thaliana\tCLV1\treviewed receptor kinase with a fasciation phenotype",
    "organ boundary specification\tWOX\tWOX activity defines the boundary between meristem and organ\tPMID:1978859\t"
    "CLV signalling could do the same\tsection the boundary\tArabidopsis wus mutant\tstem cell loss\t"
    "\tArabidopsis thaliana\tWUS\tunverified proposal with no stable accession",
    "vascular patterning\tPXY\tPXY receptor signalling sets the vascular stem cell population\tPMID:12172025\t"
    "auxin maxima could do the same\tcount procambial cells\tArabidopsis pxy mutant\tvascular overproliferation\t"
    "\tArabidopsis thaliana\tPXY\tunverified proposal with no stable accession",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _coverage_payload(payload: str, layers: list[str]) -> str:
    """Complete a custom fixture without adding any verified target or source."""
    frame = pd.read_csv(io.StringIO(payload), sep="\t", dtype=str).fillna("")
    frame.insert(0, "causal_layer", layers)
    shared = pd.read_csv(io.StringIO(PAYLOAD), sep="\t", dtype=str).fillna("")
    existing = set(frame["gene_family"])
    additional = []
    for row in pd.concat([shared.iloc[3:], shared.iloc[:3]]).to_dict("records"):
        if row["gene_family"] in existing:
            continue
        for column in ("source_ids", "model_species_evidence", "model_species_mutant_phenotype",
                       "proposed_stable_id"):
            row[column] = ""
        row["reference_species"] = "Arabidopsis thaliana"
        row["reference_gene_name"] = f"{row['gene_family']}_reference"
        additional.append(row)
        existing.add(row["gene_family"])
        if len(frame) + len(additional) == 12:
            break
    return pd.concat([frame, pd.DataFrame(additional)]).to_csv(sep="\t", index=False, lineterminator="\n")


PARALOG_PAYLOAD = _coverage_payload(PARALOG_PAYLOAD, ["signaling", "regulation", "signaling"])


def _record_cache(root: Path) -> Path:
    cache = root / "records"
    cache.mkdir(exist_ok=True)
    for accession, gene_name in (("P30924", "CLV1"), ("Q9SXK1", "PIN1"), ("Q9SWW6", "CESA9")):
        (cache / f"{accession}.json").write_text(json.dumps({
            "primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Arabidopsis thaliana"},
            "sequence": {"value": "MASTPRQTVSEK" * 8},
            "genes": [{"geneName": {"value": gene_name}}],
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "8402901"}]}}],
        }), encoding="utf-8")
    return cache


# Sequences that pass the blast-panel row validation (no generic-motif
# trigger), so the guided Route A bait freeze -> blast re-freeze yields
# verified queries rather than review_required rows.
_BLAST_SEQUENCES = {
    "P30924": ("CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVN"
               "KQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET"),
    "Q9SXK1": ("MYTCNTWCWQTEEFKCHFAIGTHRPYTYCTNKWFDNQDNCCEEVEEAAPRMLMKADWYYRIDDDAPRQVSDLGP"
               "MAKINYNRGSIDRINHLCVWDLHQHYYWWCPFYKAFREYFAKEKQR"),
    "Q9SWW6": ("MAEKLYLLLDEEIKFCPVVGPYTCVNAKDWYNLGTNCCQEVEAAPRTLIKADWYYRVDDPAARQVSDLGPMSK"
               "INYNRGSIDRINQLCVWDLHQHYYWWCPLYKAFREYFAKEKQP"),
}


def _blast_record_cache(root: Path) -> Path:
    """Record cache with blast-panel-valid sequences (issue #374 Criterion 3)."""
    cache = root / "blast_records"
    cache.mkdir(exist_ok=True)
    gene_names = {"P30924": "CLV1", "Q9SXK1": "PIN1", "Q9SWW6": "CESA9"}
    for accession, sequence in _BLAST_SEQUENCES.items():
        (cache / f"{accession}.json").write_text(json.dumps({
            "primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Arabidopsis thaliana"},
            "sequence": {"value": sequence, "length": len(sequence)},
            "genes": [{"geneName": {"value": gene_names[accession]}}],
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "8402901"}]}}],
        }), encoding="utf-8")
    return cache


def _blast_proteins(root: Path) -> Path:
    return _write(root, "proteins.fa", "\n".join([
        ">Psat_CLV1", _BLAST_SEQUENCES["P30924"],
        ">Psat_UNK", _BLAST_SEQUENCES["Q9SXK1"],
        ">Psat_UNK2", _BLAST_SEQUENCES["Q9SWW6"],
    ]) + "\n")


def _source_cache(root: Path) -> Path:
    cache = root / "pmc"
    cache.mkdir(exist_ok=True)
    records = (
        ("8402901", "CLAVATA1 receptor kinase controls meristem size in Arabidopsis", "CLV"),
        ("1978859", "Starch branching enzyme SBE and the wrinkled seed character", "PIN"),
        ("12172025", "CESA cellulose synthase and secondary wall formation", "CESA"),
    )
    for pmid, title, claimed_family in records:
        (cache / f"med_{pmid}.json").write_text(json.dumps({
            "cache_key": f"pmid_{pmid}", "id": pmid, "pmid": pmid, "title": title,
            "abstractText": (
                f"Genetic loss-of-function mutants in {claimed_family} caused the phenotype. {title}"
            ),
            "journalTitle": "Plant Cell", "pubYear": "1993",
            "doi": f"10.1/{pmid}", "authorString": "A B, C D",
        }), encoding="utf-8")
    return cache


def _hit13(query: str, subject: str, identity: float, bitscore: float, title: str, *, qlen: int = 300) -> str:
    return (f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}\t{title}")


def _build_retrieval_index(root: Path, corpus: str, database_version: str = "plant_ref_test") -> Path:
    index_out = root / "index"
    manifest = index_out / database_version / "index_manifest.json"
    if not manifest.exists():
        assert main(["retrieval-index-build",
                     "--corpus", str(_write(root, "corpus.tsv", corpus)),
                     "--database-version", database_version, "--out-dir", str(index_out)]) == 0
    return manifest


def _family_index_fixture(root: Path, hits: list[tuple[str, str, float, float, str]], genes: str) -> Path:
    blast = "\n".join([_hit13(*hit) for hit in hits]) + "\n"
    index_out = root / "family_index"
    assert main(["family-index", "--blast-hits", str(_write(root, "fi_hits.tsv", blast)),
                 "--genes", str(_write(root, "fi_genes.tsv", genes)),
                 "--out-dir", str(index_out), "--reference-label", "TestRef"]) == 0
    return index_out / "family_index.tsv"


def _case_resources(root: Path, case_id: str, interval: str) -> tuple[Path, Path]:
    case = root / "case.yaml"
    case.write_text(yaml.safe_dump({
        "case_id": case_id,
        "trait": {"name": "fasciated stem"},
        "source": {"source_id": "fixture", "target_interval": interval},
        "target": {"assembly": "fixture"},
        "inputs": {"candidates_resource": "candidates"},
        "ranking_policy": "deterministic_v1",
    }, sort_keys=False), encoding="utf-8")
    resources = root / "resources.yaml"
    resources.write_text(yaml.safe_dump({
        "resources": {"candidates": {"type": "candidates_tsv", "path": "candidates.tsv"}}
    }, sort_keys=False), encoding="utf-8")
    return case, resources


class TripleRouteProfileBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        # Explicit fixture indexes must not be supplemented by a workstation's
        # unrelated reusable genome cache when this fixture supplies none.
        for isolation in (
            mock.patch.dict(os.environ, {"CANDIDATE_QTL_FAMILY_INDEX": ""}),
            mock.patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]),
        ):
            isolation.start()
            self.addCleanup(isolation.stop)

    def _profile_args(
        self,
        root: Path,
        out: Path,
        *,
        candidates: str = CANDIDATES,
        orthology: str = ORTHOLOGY,
        payload: str = PAYLOAD,
        domains: str | None = None,
        family_index: Path | None = None,
        interval: str = "chr5:1-20000",
        retrieval_index: Path | None = None,
        skip_retrieval: bool = False,
        case_run: bool = False,
        case_id: str = "case_fas",
        phenotype: str = "fasciated stem with enlarged meristem",
        proteins: str | None = None,
        record_cache: Path | None = None,
    ) -> list[str]:
        _write(root, "corpus.tsv", REFERENCE_CORPUS)
        if retrieval_index is None and not skip_retrieval:
            retrieval_index = _build_retrieval_index(root, REFERENCE_CORPUS)
        payload_path = _write(root, "payload.tsv", payload)
        payload_rows = pd.read_csv(payload_path, sep="\t", dtype=str).fillna("").to_dict("records")
        bound_claims = write_route_a_bound_claims(root, payload_rows)
        args = [
            "guided-profile", "--profile", "triple-route", "--out-dir", str(out), "--execute",
            "--case-id", case_id, "--focal-species", "Pisum sativum",
            "--phenotype", phenotype,
            "--phenotype-localization-payload", str(_localization(root, phenotype)),
            "--mechanism-payload", str(payload_path), "--payload-source", "agent",
            "--record-cache", str(record_cache or _record_cache(root)), "--source-cache", str(_source_cache(root)),
            "--bound-claims", str(bound_claims),
            "--candidates", str(_write(root, "candidates.tsv", candidates)),
            "--orthology", str(_write(root, "orthology.tsv", orthology)),
            "--reference-function", str(root / "corpus.tsv"),
        ]
        if retrieval_index is not None:
            args += ["--retrieval-index", str(retrieval_index)]
        if domains is not None:
            args += ["--domains", str(_write(root, "domains.tsv", domains))]
        if family_index is not None:
            args += ["--family-index", str(family_index), "--interval", interval]
        if proteins is not None:
            args += ["--proteins", proteins]
        if case_run:
            case, resources = _case_resources(root, case_id, interval)
            args += ["--case", str(case), "--resources", str(resources)]
        return args

    def _run_profile(self, root: Path, *, out_name: str = "triple", **kwargs: object) -> Path:
        out = root / out_name
        self.assertEqual(0, main(self._profile_args(root, out, **kwargs)))
        records = {record["command"]: record for record in self._manifest(out)["execution_results"]}
        for command in ("phenotype-localization-freeze", "mechanism-import", "mechanism-coverage-audit",
                        "reference-gene-resolve", "source-verify", "hypothesis-verify", "bait-panel-freeze"):
            self.assertEqual("completed", records[command]["status"], command)
        self.assertEqual("pass", records["mechanism-coverage-audit"]["result"]["status"])
        return out

    def _manifest(self, out: Path) -> dict[str, object]:
        return json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))

    # -- 1. dry-run plan -----------------------------------------------------

    def test_dry_run_plans_bridge_stage_between_compare_and_family_copy(self) -> None:
        """Issue #374: z01b sits between the comparison and family-to-copy, the
        bridge chains every route input, and family-to-copy chains the bridge
        TSV as an optional promotion gate."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "plan"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--phenotype", "fasciated stem"]))
            manifest = self._manifest(out)
            stages = manifest["stages"]
            commands = [stage["command"] for stage in stages]
            order = {command: index for index, command in enumerate(commands)}
            self.assertLess(order["triple-route-compare"], order["family-candidate-bridge"])
            self.assertLess(order["family-candidate-bridge"], order["family-to-copy"])
            self.assertEqual("family-to-copy", commands[-1])
            bridge = stages[order["family-candidate-bridge"]]
            self.assertEqual("z01b_family_candidate_bridge", Path(bridge["arguments"]["out_dir"]).name)
            self.assertEqual([], bridge["required_inputs"])
            self.assertTrue(bridge["optional_chain"])
            self.assertEqual(
                ("phenotype-first-status", "hypotheses"),
                (bridge["inputs_from"]["route_a_hypotheses"]["stage"],
                 bridge["inputs_from"]["route_a_hypotheses"]["result_key"]),
            )
            self.assertEqual(
                ("route-a-interval-match", "status"),
                (bridge["inputs_from"]["route_a_family_status"]["stage"],
                 bridge["inputs_from"]["route_a_family_status"]["result_key"]),
            )
            self.assertEqual(
                ("bait-panel-freeze", "panel"),
                (bridge["inputs_from"]["route_a_bait_panel"]["stage"],
                 bridge["inputs_from"]["route_a_bait_panel"]["result_key"]),
            )
            self.assertEqual(
                ("route-b-prepare", "orthology"),
                (bridge["inputs_from"]["route_b_orthology"]["stage"],
                 bridge["inputs_from"]["route_b_orthology"]["result_key"]),
            )
            self.assertEqual(
                ("route-b-prepare", "interval_gene_atlas"),
                (bridge["inputs_from"]["route_b_atlas"]["stage"],
                 bridge["inputs_from"]["route_b_atlas"]["result_key"]),
            )
            self.assertEqual(
                ("retrieval-import", "frozen"),
                (bridge["inputs_from"]["route_c_retrieval"]["stage"],
                 bridge["inputs_from"]["route_c_retrieval"]["result_key"]),
            )
            self.assertEqual(
                ("retrieval-interval-match", "interval_match"),
                (bridge["inputs_from"]["route_c_match"]["stage"],
                 bridge["inputs_from"]["route_c_match"]["result_key"]),
            )
            self.assertEqual(
                ("route-b-prepare", "interval_gene_atlas"),
                (bridge["inputs_from"]["genes"]["stage"],
                 bridge["inputs_from"]["genes"]["result_key"]),
            )
            self.assertEqual(
                ("route-b-prepare", "dependency_clusters"),
                (bridge["inputs_from"]["dependency_clusters"]["stage"],
                 bridge["inputs_from"]["dependency_clusters"]["result_key"]),
            )
            self.assertIn("dependency_clusters", bridge["optional_chain_inputs"])
            self.assertIn("rank_contribution stays none", bridge["stop_condition"])
            copy_stage = stages[order["family-to-copy"]]
            self.assertEqual(
                ("family-candidate-bridge", "bridge_tsv"),
                (copy_stage["inputs_from"]["bridge"]["stage"],
                 copy_stage["inputs_from"]["bridge"]["result_key"]),
            )
            self.assertIn("bridge", copy_stage["optional_chain_inputs"])
            self.assertIn("<from family-candidate-bridge.bridge_tsv>", copy_stage["arguments"]["bridge"])
            # Issue #374 route A layer: the blast-panel-search evidence is an
            # optional chained input of family-to-copy too, so a missing
            # blast-panel-search stage (proteins not supplied) reads as None.
            self.assertEqual(
                ("blast-panel-search", "evidence"),
                (copy_stage["inputs_from"]["route_a_blast"]["stage"],
                 copy_stage["inputs_from"]["route_a_blast"]["result_key"]),
            )
            self.assertIn("route_a_blast", copy_stage["optional_chain_inputs"])
            self.assertIn("<from blast-panel-search.evidence>", copy_stage["arguments"]["route_a_blast"])

    # -- 2. execute reaches the bridge and preserves rank --------------------

    def test_execute_reaches_bridge_and_preserves_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root, case_run=True)
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            for command in ("family-candidate-bridge", "family-to-copy", "case-run"):
                self.assertEqual("completed", records[command]["status"], command)
            self.assertEqual("completed", manifest["legacy_profile"]["overall_status"])
            bridge_out = out / "z01b_family_candidate_bridge"
            self.assertTrue((bridge_out / "family_candidate_bridge.tsv").exists())
            self.assertTrue((bridge_out / "family_candidate_bridge_manifest.json").exists())
            frame = pd.read_csv(bridge_out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(list(BRIDGE_COLUMNS), list(frame.columns))
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            self.assertFalse({"dry_score", "score", "dry_group", "rank"} & set(frame.columns))
            bridge_manifest = json.loads((bridge_out / "family_candidate_bridge_manifest.json")
                                         .read_text(encoding="utf-8"))
            self.assertEqual("none", bridge_manifest["rank_contribution"])
            self.assertIn("Psat_CLV1", set(frame["gene_id"]))
            # Route B dependency clusters: route-b-prepare derived one cluster
            # from the orthology fixture's OG0002, the bridge chained it in via
            # --dependency-clusters, and the Psat_CLV1 orthology row resolved
            # through it (rank_contribution stays none).
            self.assertTrue(bridge_manifest["inputs"]["dependency_clusters"]["present"])
            rbh = frame[(frame["gene_id"].eq("Psat_CLV1")) & (frame["route"].eq("b"))
                        & (frame["match_type"].eq("blastp"))]
            self.assertEqual(1, len(rbh))
            self.assertEqual(dependency_cluster_id("orthogroup:OG0002"), rbh.iloc[0]["dependency_cluster_id"])
            # The bridge TSV reached family-to-copy and the executed chain
            # reached the narrowed stage A and next-experiment outputs.
            case_out = out / "z03_case_run"
            copy_manifest = json.loads((out / "z02_family_to_copy" / "family_copy_manifest.json")
                                       .read_text(encoding="utf-8"))
            self.assertTrue(copy_manifest["bridge"]["supplied"])
            self.assertEqual(str(bridge_out / "family_candidate_bridge.tsv"),
                             copy_manifest["inputs"]["bridge"]["path"])
            self.assertTrue((case_out / "stage_a_narrowed" / "deterministic_shortlist.tsv").exists())
            options = pd.read_csv(case_out / "next_experiment" / "next_experiment_options.tsv",
                                  sep="\t", dtype=str).fillna("")
            self.assertTrue(options["input_facts"].str.contains("integrated_narrowing_status=evaluated").all())
            # rank preservation: the guided case-run's ranked table is
            # byte-identical to a plain case-run over the same inputs.
            baseline = root / "baseline_case"
            case, resources = _case_resources(root, "case_fas", "chr5:1-20000")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources),
                                      "--out-dir", str(baseline)]))
            self.assertEqual(
                (baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                (case_out / "rank" / "ranked_candidates.tsv").read_bytes(),
            )

    # -- 3. paralog E2E: only the source-backed copy is narrowed -------------

    def _paralog_family_index(self, root: Path) -> Path:
        return _family_index_fixture(
            root,
            hits=[
                ("G_A", "sp|P30924|CLV1_ARATH", 85.0, 400.0,
                 "Receptor kinase CLV1 OS=Arabidopsis thaliana OX=3702 GN=CLV1 PE=1 SV=2"),
                ("G_W", "sp|Q9FUW6|WUS_ARATH", 80.0, 350.0,
                 "Homeobox protein WUS OS=Arabidopsis thaliana OX=3702 GN=WUS PE=1 SV=1"),
                ("G_R1", "sp|Q9M2S1|PXY_ARATH", 82.0, 360.0,
                 "PXY receptor OS=Arabidopsis thaliana OX=3702 GN=PXY PE=1 SV=1"),
                ("G_R2", "sp|Q9M2S1|PXY_ARATH", 75.0, 320.0,
                 "PXY receptor OS=Arabidopsis thaliana OX=3702 GN=PXY PE=1 SV=1"),
            ],
            genes="gene_id\tchrom\tstart\tend\n"
                  "G_A\tchr1\t1000\t3000\n"
                  "G_W\tchr1\t20000\t22000\n"
                  "G_R1\tchr1\t31000\t33000\n"
                  "G_R2\tchr1\t35000\t37000\n",
        )

    def test_synthetic_paralog_e2e_narrowed_top20_keeps_unique_source_backed_only(self) -> None:
        """Issue #374 acceptance: bridge -> family-to-copy -> integrated
        narrowing -> stage A narrowed. Of the two paralog families only the
        source-backed copy (G_A, bridge tier1) enters the family copy shortlist
        and the validated N1 tier; the name-only member (G_W, bridge tier5)
        never rises to a validated-family tier."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._paralog_family_index(root)
            out = self._run_profile(root, out_name="paralog",
                                    candidates=PARALOG_CANDIDATES,
                                    orthology=PARALOG_ORTHOLOGY,
                                    payload=PARALOG_PAYLOAD,
                                    domains=PARALOG_DOMAINS,
                                    family_index=index,
                                    interval="chr1:1-50000",
                                    case_run=True,
                                    case_id="case_para",
                                    phenotype="fasciated stem with enlarged meristem")
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            for command in ("family-candidate-bridge", "family-to-copy", "case-run"):
                self.assertEqual("completed", records[command]["status"], command)

            bridge = pd.read_csv(out / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv",
                                 sep="\t", dtype=str).fillna("")
            g_a = bridge[bridge["gene_id"].eq("G_A")]
            self.assertTrue(g_a["mapping_tier"].eq("tier1_source_backed_unique").any())
            g_w = bridge[bridge["gene_id"].eq("G_W")]
            self.assertFalse(g_w.empty)
            self.assertTrue(g_w["mapping_tier"].eq("tier5_name_only").all())
            self.assertTrue(bridge["rank_contribution"].eq("none").all())

            shortlist = pd.read_csv(out / "z02_family_to_copy" / "family_copy_shortlist.tsv",
                                    sep="\t", dtype=str).fillna("")
            self.assertEqual(["G_A"], list(shortlist["candidate_gene_id"]))
            self.assertEqual("tier1_source_backed_unique",
                             shortlist.iloc[0]["bridge_mapping_tier"])
            self.assertEqual("first", shortlist.iloc[0]["first_wet_priority"])

            discrimination = pd.read_csv(out / "z02_family_to_copy" / "family_copy_discrimination.tsv",
                                         sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # The name-only member is recorded with its bridge tier5 state and
            # is never promoted; the mixed PXY family records bridge-backed
            # copy-selection loss.
            self.assertEqual("tier5_name_only", discrimination.loc["G_W", "bridge_mapping_tier"])
            self.assertEqual("covered", discrimination.loc["G_W", "bridge_missingness"])
            self.assertEqual("tier5_name_only", discrimination.loc["G_R2", "bridge_mapping_tier"])

            narrowing = pd.read_csv(out / "z03_case_run" / "integrated_narrowing" / "integrated_narrowing.tsv",
                                    sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("N1", narrowing.loc["G_A", "narrowing_tier"])
            self.assertEqual("validated_family_preferred_copy", narrowing.loc["G_A", "narrowing_reason"])
            self.assertEqual("preferred_copy", narrowing.loc["G_A", "copy_state"])
            self.assertIn(narrowing.loc["G_W", "narrowing_tier"], {"N3", "N4"})
            self.assertIn(narrowing.loc["G_R2", "narrowing_tier"], {"N3", "N4"})
            self.assertLess(int(narrowing.loc["G_A", "narrowing_priority"]),
                            int(narrowing.loc["G_W", "narrowing_priority"]))

            narrowed = pd.read_csv(out / "z03_case_run" / "stage_a_narrowed" / "deterministic_shortlist.tsv",
                                   sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertIn("G_A", narrowed.index)
            self.assertEqual("N1", narrowed.loc["G_A", "narrowing_tier"])
            self.assertEqual("first", narrowed.loc["G_A", "first_wet_priority"])
            self.assertIn(narrowed.loc["G_W", "narrowing_tier"], {"N3", "N4"})
            self.assertLess(int(narrowed.loc["G_A", "deterministic_review_order"]),
                            int(narrowed.loc["G_W", "deterministic_review_order"]))

    # -- 4. no bait / no provider / empty C degrades, never fails ------------

    def test_all_route_inputs_empty_records_bridge_as_not_runnable(self) -> None:
        """Issue #374 review: with no candidates / payload / retrieval at all,
        the bridge stage must never be recorded as failed and the profile
        overall must not be failed. Route A still emits its empty hypothesis
        table, so the bridge completes over that empty table with zero rows
        (an equivalent not-failed state); when every chained input is truly
        absent the stage is recorded as not_runnable (covered by
        test_execute_stage_guard_skips_all_empty_bridge_inputs)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "empty"
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route",
                                      "--out-dir", str(out), "--execute",
                                      "--case-id", "case_empty",
                                      "--focal-species", "Pisum sativum",
                                      "--phenotype", "fasciated stem"]))
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            self.assertNotIn("failed", manifest["overall_status"])
            self.assertNotIn("family-candidate-bridge", manifest.get("failed_stages", []))
            bridge = records["family-candidate-bridge"]
            self.assertNotEqual("failed", bridge["status"])
            if str(bridge["status"]) == "completed":
                bridge_manifest = json.loads((out / "z01b_family_candidate_bridge"
                                              / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(0, bridge_manifest["row_counts"]["total"])
                frame = pd.read_csv(out / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv",
                                    sep="\t", dtype=str)
                self.assertEqual(0, len(frame))
            else:
                self.assertTrue(str(bridge["status"]).startswith("not_runnable"), bridge)

    def test_execute_stage_guard_skips_all_empty_bridge_inputs(self) -> None:
        """Issue #374 review: when every chained route input of the bridge
        stays empty, _execute_stage must not call family_candidate_bridge()
        (which would raise "requires at least one route input" and record the
        stage as failed); it returns the not_runnable marker instead."""
        import candidate_qtl.guided_profiles as guided_profiles_module

        stage = {
            "command": "family-candidate-bridge",
            "arguments": {
                "case_id": "case_empty", "out_dir": "/tmp/nonexistent_bridge_out",
                "route_a_hypotheses": "", "route_a_family_status": "", "route_a_bait_panel": "",
                "route_b_orthology": "", "route_b_atlas": "", "route_c_retrieval": "",
                "route_c_match": "", "genes": "",
            },
        }
        result = guided_profiles_module._execute_stage(stage)
        self.assertEqual("not_runnable_no_route_input", result["guided_stage_status"])
        self.assertIn("at least one route input", result["reason"])

    def test_execute_stage_guard_keeps_a_single_route_input_runnable(self) -> None:
        """The guard must not change behavior when at least one input exists:
        the stage still dispatches to family_candidate_bridge."""
        import candidate_qtl.guided_profiles as guided_profiles_module

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = _write(root, "hypotheses.tsv", "hypothesis_id\tmember_interval_genes\nh1\t\n")
            stage = {
                "command": "family-candidate-bridge",
                "arguments": {
                    "case_id": "case_a", "out_dir": str(root / "bridge_out"),
                    "route_a_hypotheses": str(hypotheses), "route_a_family_status": "",
                    "route_a_bait_panel": "", "route_b_orthology": "", "route_b_atlas": "",
                    "route_c_retrieval": "", "route_c_match": "", "genes": "",
                },
            }
            result = guided_profiles_module._execute_stage(stage)
            self.assertNotIn("guided_stage_status", result)
            self.assertTrue(str(result.get("bridge_tsv", "")).endswith("family_candidate_bridge.tsv"))

    def test_blocked_preferred_copy_never_rises_in_narrowing(self) -> None:
        """Issue #374 review: a member whose bridge tier is tier5 (name-only)
        but whose own homology layer makes it preferred_copy keeps its
        disambiguation status in the discrimination artifact yet never enters
        the shortlist -- and the integrated narrowing chained over those
        artifacts must not turn it into an N1 (gate-blocked members are not
        promoted)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root, "families.tsv", "\n".join([
                "hypothesis_id\tgene_family\tmember_interval_genes",
                "hyp_solo\tSOLO\tG_S",
            ]) + "\n")
            _write(root, "atlas.tsv", "\n".join([
                "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\torthology_reference_gene",
                "G_S\tchr1\t1000\t3000\tSOLO\tLRR-Pkinase\tREF_S",
            ]) + "\n")
            _write(root, "orthology.tsv", "\n".join([
                "gene_id\trelationship\tprovider\tevidence_chain\treference_gene_id\torthogroup",
                "G_S\tone_to_one\torthofinder\tOF_chain\tREF_S\tOG_SOLO",
            ]) + "\n")
            bridge = root / "bridge.tsv"
            bridge.write_text("\t".join(BRIDGE_COLUMNS) + "\n"
                              + "\t".join(["demo", "G_S", "SOLO", "b", "", "", "name_only", "", "", "", "", "1",
                                           "paralog_ambiguous", "complete", "", "ff", "ff", "review_required",
                                           "dep_SOLO", "tier5_name_only", "none", ""]) + "\n", encoding="utf-8")
            copy_out = root / "copy"
            self.assertEqual(0, main(["family-to-copy", "--families", str(root / "families.tsv"),
                                      "--out-dir", str(copy_out), "--case-id", "demo",
                                      "--atlas", str(root / "atlas.tsv"),
                                      "--orthology", str(root / "orthology.tsv"),
                                      "--bridge", str(bridge)]))
            discrimination = pd.read_csv(copy_out / "family_copy_discrimination.tsv",
                                         sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", discrimination.loc["G_S", "disambiguation_status"])
            self.assertEqual("tier5_name_only", discrimination.loc["G_S", "bridge_mapping_tier"])
            self.assertEqual("blocked_tier5", discrimination.loc["G_S", "bridge_promotion_gate"])
            shortlist = pd.read_csv(copy_out / "family_copy_shortlist.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(shortlist))
            # narrowing over the same artifacts: G_S must stay off N1.
            _write(root, "ranked.tsv", "\n".join([
                "gene_id\tchrom\tstart\tend\tdry_score\tdry_rank\tdry_group\ttop_eligibility\tscreening_rank",
                "G_S\tchr1\t1000\t3000\t5.0\t1\ttop\tstrong_genetic_evidence\t1",
            ]) + "\n")
            narrow_out = root / "narrow"
            self.assertEqual(0, main(["integrated-narrowing", "--ranked", str(root / "ranked.tsv"),
                                      "--out-dir", str(narrow_out), "--case-id", "demo",
                                      "--family-copy", str(copy_out / "family_copy_discrimination.tsv"),
                                      "--family-copy-shortlist", str(copy_out / "family_copy_shortlist.tsv")]))
            narrowed = pd.read_csv(narrow_out / "integrated_narrowing.tsv",
                                   sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertNotEqual("N1", narrowed.loc["G_S", "narrowing_tier"])
            self.assertEqual("N3", narrowed.loc["G_S", "narrowing_tier"])
            self.assertEqual("base_top_reserve_rescue", narrowed.loc["G_S", "narrowing_reason"])
            self.assertEqual("preferred_copy", narrowed.loc["G_S", "copy_state"])
            self.assertEqual("", narrowed.loc["G_S", "first_wet_priority"])

    def test_no_bait_no_provider_empty_c_complete_as_not_evaluated(self) -> None:
        """No verified bait, no family index, and no retrieval index: the
        bridge runs with the available route B evidence only, its rows carry
        the not_evaluated vocabulary, and the profile overall never fails."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = "\n".join([
                "pathway_or_process\tgene_family\tmechanism_statement\tsource_ids\tcounter_hypothesis\t"
                "discriminating_observation\tmodel_species_evidence\tmodel_species_mutant_phenotype\t"
                "proposed_stable_id\treference_species\treference_gene_name\tselection_rationale",
                "meristem maintenance\tCLV\tloss of CLV signalling enlarges the stem cell pool\tPMID:8402901\t"
                "a WOX lesion could do the same\tcount meristem cell layers\tArabidopsis clv1 mutant\tfasciated stem\t"
                "\tArabidopsis thaliana\tCLV1\tno accession supplied",
                "hormone transport\tPIN\tauxin redistribution changes stem geometry\tPMID:1978859\t"
                "meristem lesions give the same shape\tassay auxin flux\tArabidopsis pin mutant\tfused organs\t"
                "\tArabidopsis thaliana\tPIN1\tno accession supplied",
                "cell wall biosynthesis\tCESA\tweakened wall changes stem shape\tPMID:12172025\t"
                "wall phenotypes rarely fasciate\tsection the stem\tArabidopsis irx mutant\tcollapsed xylem\t"
                "\tArabidopsis thaliana\tCESA9\tno accession supplied",
            ]) + "\n"
            payload = _coverage_payload(payload, ["signaling", "execution", "execution"])
            orthology = "\n".join([
                "gene_id\treference_gene_id\torthogroup\trelationship\tpercent_identity\tquery_coverage\tprovider",
                "Psat_CLV1\tAtCLV1\tOG0002\tbest_hit_homology\t0.2\t0.1\tlocal_blastp:test",
            ]) + "\n"
            out = self._run_profile(root, out_name="degraded",
                                    candidates=CANDIDATES, orthology=orthology,
                                    payload=payload, skip_retrieval=True,
                                    family_index=None, case_run=False)
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            self.assertNotIn("failed", manifest["overall_status"])
            self.assertEqual("completed", records["family-candidate-bridge"]["status"])
            # No retrieval index: route C is unavailable, never a failure.
            self.assertTrue(records["retrieval-run"]["status"].startswith("not_runnable"))
            bridge_manifest = json.loads((out / "z01b_family_candidate_bridge"
                                          / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual({"b": 2}, bridge_manifest["row_counts"]["by_route"])
            tiers = bridge_manifest["row_counts"]["by_mapping_tier"]
            self.assertGreaterEqual(tiers.get("tier6_not_evaluated", 0), 1)
            frame = pd.read_csv(out / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv",
                                sep="\t", dtype=str).fillna("")
            tier6 = frame[frame["mapping_tier"].eq("tier6_not_evaluated")]
            self.assertFalse(tier6.empty)
            self.assertTrue(tier6["reason"].eq("not_evaluated_coverage_below_threshold").all())
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            self.assertFalse(bridge_manifest["inputs"]["route_c_retrieval"]["present"])

    # -- 5. two runs are byte-identical --------------------------------------

    def test_two_runs_produce_identical_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self._run_profile(root, out_name="run1", case_run=True)
            second = self._run_profile(root, out_name="run2", case_run=True)
            for relative in (
                "z01b_family_candidate_bridge/family_candidate_bridge.tsv",
                "z02_family_to_copy/family_copy_discrimination.tsv",
                "z03_case_run/integrated_narrowing/integrated_narrowing.tsv",
                "z03_case_run/stage_a_narrowed/deterministic_shortlist.tsv",
                "z03_case_run/rank/ranked_candidates.tsv",
            ):
                with self.subTest(artifact=relative):
                    self.assertEqual(
                        hashlib.sha256((first / relative).read_bytes()).hexdigest(),
                        hashlib.sha256((second / relative).read_bytes()).hexdigest(),
                    )

    # -- 6. legacy stage A and rank unchanged --------------------------------

    def test_legacy_stage_a_and_rank_unchanged(self) -> None:
        """The bridge run's stage A and ranked tables are byte-identical to a
        plain case-run over the same inputs: the bridge adds a parallel
        reporting layer and never moves a candidate."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_profile(root, case_run=True)
            baseline = root / "baseline_case"
            case, resources = _case_resources(root, "case_fas", "chr5:1-20000")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources),
                                      "--out-dir", str(baseline)]))
            for relative in ("stage_a/deterministic_shortlist.tsv", "rank/ranked_candidates.tsv"):
                with self.subTest(artifact=relative):
                    self.assertEqual(
                        (baseline / relative).read_bytes(),
                        (out / "z03_case_run" / relative).read_bytes(),
                    )

    # -- 7. manifest bridge totals -------------------------------------------

    def test_manifest_reports_bridge_rescue_and_loss(self) -> None:
        """The family-copy manifest exposes the bridge rescue counts and
        copy-selection loss from the executed E2E run."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._paralog_family_index(root)
            out = self._run_profile(root, out_name="paralog",
                                    candidates=PARALOG_CANDIDATES,
                                    orthology=PARALOG_ORTHOLOGY,
                                    payload=PARALOG_PAYLOAD,
                                    domains=PARALOG_DOMAINS,
                                    family_index=index,
                                    interval="chr1:1-50000",
                                    case_run=True,
                                    case_id="case_para",
                                    phenotype="fasciated stem with enlarged meristem")
            copy_manifest = json.loads((out / "z02_family_to_copy" / "family_copy_manifest.json")
                                       .read_text(encoding="utf-8"))
            bridge = copy_manifest["bridge"]
            self.assertTrue(bridge["supplied"])
            self.assertIsNotNone(bridge["bridge_tsv_sha256"])
            # No member of the paralog fixture is a preferred copy on a
            # tier5/tier6 bridge tier, so nothing is blocked; the mixed PXY
            # family (source-backed G_R1 + name-only G_R2) cannot be narrowed
            # to one copy and reports bridge-backed copy-selection loss.
            self.assertEqual(0, bridge["blocked_count"])
            self.assertEqual(1, bridge["copy_selection_loss_count"])
            self.assertEqual({
                "rescue_a_only": 0, "rescue_b_only": 2, "rescue_c_only": 0, "multi_route_supported": 1,
                "not_evaluated": 0,
            }, bridge["bridge_rescue_counts"])

    # -- 8. Route A blast-panel wiring (#374 Criterion 3) -------------------

    def _blast_family_index(self, root: Path) -> Path:
        return _family_index_fixture(
            root,
            hits=[
                ("Psat_CLV1", "sp|P30924|CLV1_ARATH", 85.0, 400.0,
                 "Receptor kinase CLV1 OS=Arabidopsis thaliana OX=3702 GN=CLV1 PE=1 SV=2"),
                ("Psat_UNK", "sp|P12345|UNK_ARATH", 50.0, 100.0,
                 "Hypothetical protein OS=Arabidopsis thaliana OX=3702 PE=4 SV=1"),
            ],
            genes="gene_id\tchrom\tstart\tend\n"
                  "Psat_CLV1\tchr5\t1000\t3000\n"
                  "Psat_UNK\tchr5\t9000\t11000\n",
        )

    def _blast_raw_hits(self) -> str:
        # qseqids are the verified bait accessions (the blast panel_ids).
        return (
            "P30924\tPsat_CLV1\t99.0\t120\t120\t120\t1\t120\t1\t120\t0.0\t250.0\n"
            "Q9SXK1\tPsat_UNK\t98.0\t120\t120\t120\t1\t120\t1\t120\t1e-50\t240.0\n"
            "Q9SXK1\tPsat_UNK2\t97.0\t120\t120\t120\t1\t120\t1\t120\t1e-50\t238.0\n"
        )

    def _bridge_frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv",
                           sep="\t", dtype=str).fillna("")

    def test_dry_run_plans_blast_panel_stages_after_bait_freeze_when_proteins_supplied(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "plan"
            proteins = _blast_proteins(root)
            self.assertEqual(0, main(["guided-profile", "--profile", "triple-route", "--out-dir", str(out),
                                      "--phenotype", "fasciated stem", "--proteins", str(proteins)]))
            manifest = self._manifest(out)
            commands = [stage["command"] for stage in manifest["stages"]]
            self.assertIn("blast-panel-freeze", commands)
            self.assertIn("blast-panel-search", commands)
            self.assertLess(commands.index("bait-panel-freeze"), commands.index("blast-panel-freeze"))
            self.assertLess(commands.index("blast-panel-freeze"), commands.index("blast-panel-search"))
            self.assertLess(commands.index("blast-panel-search"), commands.index("route-a-interval-match"))
            search = next(stage for stage in manifest["stages"] if stage["command"] == "blast-panel-search")
            self.assertEqual(("blast-panel-freeze", "manifest"),
                             (search["inputs_from"]["freeze"]["stage"],
                              search["inputs_from"]["freeze"]["result_key"]))
            self.assertIn("proteins", search["required_inputs"])
            bridge = next(stage for stage in manifest["stages"] if stage["command"] == "family-candidate-bridge")
            self.assertEqual(("blast-panel-search", "evidence"),
                             (bridge["inputs_from"]["route_a_blast"]["stage"],
                              bridge["inputs_from"]["route_a_blast"]["result_key"]))
            self.assertIn("route_a_blast", bridge["optional_chain_inputs"])

    def test_route_a_blast_evidence_reaches_the_bridge_e2e(self) -> None:
        """Issue #374 Criterion 3: guided triple-route runs bait-panel-freeze ->
        blast-panel-freeze -> blast-panel-search between the bait freeze and the
        interval family match, and the blast identity / coverage / evalue are
        transcribed onto the matched Route A bridge row by the verified bait
        accession. A blast hit never promotes and never moves rank_contribution."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._blast_family_index(root)
            proteins = _blast_proteins(root)
            cache = _blast_record_cache(root)
            out = root / "blast_run"
            with MockBlastTools(raw_hits=self._blast_raw_hits()):
                self.assertEqual(0, main(self._profile_args(
                    root, out, family_index=index, interval="chr5:1-20000",
                    proteins=str(proteins), record_cache=cache)))
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            for command in ("blast-panel-freeze", "blast-panel-search", "family-candidate-bridge"):
                self.assertEqual("completed", records[command]["status"], command)
            self.assertEqual("completed", manifest["legacy_profile"]["overall_status"])
            search = json.loads((out / "a05d_blast_panel_search" / "blast_panel_search_manifest.json")
                                .read_text(encoding="utf-8"))
            self.assertEqual("evaluable", search["provider_status"])
            self.assertEqual("post_freeze_similarity_only", search["interval_usage"])
            self.assertTrue((out / "a05c_blast_panel_freeze" / "blast_panel_freeze_manifest.json").exists())
            evidence = pd.read_csv(out / "a05d_blast_panel_search" / "blast_evidence.tsv",
                                   sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(evidence))
            self.assertTrue(all(str(row["rank_contribution"]) == "none" for row in evidence.to_dict("records")))
            bridge = self._bridge_frame(out)
            clv = bridge[(bridge["gene_id"].eq("Psat_CLV1")) & bridge["route"].eq("a")].iloc[0]
            self.assertEqual("P30924", clv["reference_accession"])
            self.assertEqual("99.0", clv["identity"])
            self.assertEqual("1.0", clv["query_coverage"])
            self.assertEqual("1.0", clv["target_coverage"])
            self.assertEqual("0.0", clv["evalue"])
            self.assertEqual("resolved_single_copy", clv["ambiguity_status"])
            self.assertEqual("tier1_source_backed_unique", clv["mapping_tier"])
            self.assertTrue(bridge["rank_contribution"].eq("none").all())
            bridge_manifest = json.loads((out / "z01b_family_candidate_bridge"
                                          / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, bridge_manifest["route_a_blast"]["evidence_rows"])
            self.assertGreaterEqual(bridge_manifest["route_a_blast"]["rows_joined"], 1)

    def test_route_a_blast_unavailable_degrades_to_not_evaluated_in_bridge(self) -> None:
        """A missing BLAST+ runtime is not a miss: blast-panel-search completes
        with provider_status not_evaluated_blast_unavailable (exit 0), the
        bridge sees no transcriptable blast values, and the matched Route A row
        keeps its family-status evidence."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = self._blast_family_index(root)
            proteins = _blast_proteins(root)
            cache = _blast_record_cache(root)
            out = root / "no_blast"
            with mock.patch("candidate_qtl.blast_evidence.shutil.which", return_value=None):
                self.assertEqual(0, main(self._profile_args(
                    root, out, family_index=index, interval="chr5:1-20000",
                    proteins=str(proteins), record_cache=cache)))
            manifest = self._manifest(out)
            records = {record["command"]: record for record in manifest["execution_results"]}
            self.assertEqual("completed", records["blast-panel-search"]["status"])
            self.assertNotIn("failed", manifest["overall_status"])
            search = json.loads((out / "a05d_blast_panel_search" / "blast_panel_search_manifest.json")
                                .read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_blast_unavailable", search["provider_status"])
            bridge = self._bridge_frame(out)
            clv = bridge[(bridge["gene_id"].eq("Psat_CLV1")) & bridge["route"].eq("a")].iloc[0]
            self.assertEqual("", clv["identity"])
            self.assertEqual("", clv["query_coverage"])
            self.assertEqual("", clv["evalue"])
            self.assertEqual("tier1_source_backed_unique", clv["mapping_tier"])
            self.assertEqual("resolved_single_copy", clv["ambiguity_status"])


if __name__ == "__main__":
    unittest.main()
