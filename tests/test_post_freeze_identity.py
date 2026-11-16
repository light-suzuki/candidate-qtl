"""Post-freeze identity changes only a bound paper/member's relevance evidence."""
import copy
import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from Bio import SeqIO

from candidate_qtl.cli import main
from candidate_qtl.core import sha256
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.post_freeze_identity import (
    IDENTITY_CONTRACT, REVIEW_AXES, REVIEW_FIELDS, identity_material_bindings, load_identity_assessment,
)
from candidate_qtl.reference_sequence_bridge import build_reference_sequence_bridge
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_paper_checklist import report
from tests.test_paper_tiers import claim
from tests._issue664_helpers import manifest_for
from tests.test_reference_sequence_bridge import S, _hit, _inputs
from tests import test_route_a_review_relevance as route_a_fixtures


def descriptor(path, **extra):
    return {"path": str(path), "sha256": sha256(path), **extra}


def quote(path):
    text = Path(path).read_bytes().decode("utf-8")
    return {"start": 0, "end": len(text), "text": text}


class PostFreezeIdentityTests(unittest.TestCase):
    def fixture(self, root, *, hits=None, targets=None, state="evaluated", equivalent="yes", sequence=S,
                complete_targets=True):
        request = {"request_id": "r1", "source_id": "PMID:1234567", "reported_gene": "GeneX",
                   "reference_accession": "AA123456.1", "reference_taxid": "3888",
                   "focal_taxid": "3888", "cds_protein_id": "ABC12345.1"}
        targets = dict(targets or {"p1": sequence})
        if complete_targets:
            targets.setdefault("p2", "A" * len(sequence))
            targets.setdefault("p3", "C" * len(sequence))
        req, proteins, mapping, cache = _inputs(root, requests=[request], targets=targets,
                                               genes=[("GeneX", "ABC12345.1", sequence)])
        length = len(sequence)
        default_hits = [_hit(qlen=length, slen=length, qend=length, send=length)]
        with patch("candidate_qtl.reference_sequence_bridge._run_search",
                   return_value={"state": state, "hits": default_hits if hits is None else hits}):
            bridge = build_reference_sequence_bridge(req, proteins, mapping, root / "sequence", records_arg=cache)
        paper = {**report(phenotype_reported="yes", equivalent_target_phenotype=equivalent,
                          focal_species="yes", same_gene_as_focal="unknown"),
                 "source_id": "PMID:1234567", "reported_gene": "GeneX", "reported_species": "Plant alpha"}
        paper["source_title"] = "GeneX experimental paper"
        paper_input = root / "papers.tsv"
        pd.DataFrame([{"hypothesis_id": "h1", PAPER_FIELD: normalize_paper_reports([paper])}]).to_csv(paper_input, sep="\t", index=False)
        atlas = root / "atlas.tsv"
        pd.DataFrame({"gene_id": ["g1", "g2", "g3"]}).to_csv(atlas, sep="\t", index=False)
        canonical = root / "canonical-proteins.fa"
        canonical.write_text("".join(f">g{index}\n{value}\n" for index, value in
                                   enumerate(targets.values(), 1)), encoding="utf-8")
        source = root / "paper.json"
        source.write_text(json.dumps({"id": "1234567", "source": "MED", "pmid": "1234567",
                          "title": paper["source_title"],
                          "abstractText": "Paper 1234567 reports GeneX in Plant alpha and associates it with AA123456.1."}), encoding="utf-8")
        raw_bridge = json.loads(Path(bridge["manifest"]).read_text())
        cds = Path(raw_bridge["retrieval"][0]["raw_path"])
        item = {"proposal_id": "h1", **request, "reported_species": "Plant alpha", "sequence_reported_gene": "GeneX",
                "query_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "candidate_gene_id": "g1", "decision": "same_gene",
                "reason": "The original paper's accession, experimental taxon, complete model and competitors were reviewed.",
                "review_axes": dict.fromkeys(REVIEW_AXES, "yes"),
                "paper_source": descriptor(source, source_id="PMID:1234567"),
                "paper_quote": quote(source), "cds_quote": quote(cds)}
        review = root / "review.json"
        data = {"contract": IDENTITY_CONTRACT, "case_id": "identity-case", "focal_taxid": "3888",
                "post_freeze": True, "support_effect": "none", "atlas": descriptor(atlas),
                "candidate_proteins": descriptor(canonical),
                "paper_input": descriptor(paper_input, id_column="hypothesis_id"),
                "reference_bridge": descriptor(Path(bridge["manifest"])),
                "review": descriptor(source, reviewer_role="independent_agent", provider="fixture",
                                     model="fixture-reviewer", run_id="fixture-run"), "assessments": [item]}
        path = root / "identity.json"
        fixture = {"data": data, "item": item, "path": path, "review": review, "atlas": atlas,
                   "claims": [claim([paper], candidate_gene_id="g1", hypothesis_id="h1")],
                   "paper": paper, "bridge": bridge, "canonical": canonical}
        self.save(fixture)
        return fixture

    def save(self, fixture, *, sync_review=True):
        if sync_review:
            bridge = json.loads(Path(fixture["data"]["reference_bridge"]["path"]).read_text())
            raw_hash = bridge["retrieval"][0]["raw_sha256"]
            fixture["item"]["material_bindings"] = identity_material_bindings(fixture["data"], fixture["item"], raw_hash)
            fixture["review"].write_text(json.dumps({key: fixture["item"][key] for key in REVIEW_FIELDS}), encoding="utf-8")
            fixture["item"]["review_quote"] = quote(fixture["review"])
            fixture["data"]["review"].update(descriptor(fixture["review"]))
        fixture["path"].write_text(json.dumps(fixture["data"]), encoding="utf-8")

    def load(self, fixture, **changes):
        return load_identity_assessment(str(fixture["path"]), case_id="identity-case",
                                        atlas_arg=str(fixture["atlas"]), claims=fixture["claims"],
                                        canonical_proteins_arg=str(fixture["canonical"]), **changes)

    def test_unique_exact_review_derives_gene_bound_t1_without_mutating_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp))
            before = copy.deepcopy(f["claims"])
            hashes = {path: sha256(path) for path in Path(tmp).rglob("*") if path.is_file()}
            evidence, audit = self.load(f)
            self.assertEqual(1, audit["resolved_evidence_count"])
            self.assertEqual("none", audit["support_effect"])
            rows = annotate_relevance([{"candidate_gene_id": gene} for gene in ("g1", "g2")], f["claims"] + evidence)
            self.assertEqual(NEW_TIERS[0], rows[0]["relevance_tier"])
            self.assertNotEqual(NEW_TIERS[0], rows[1]["relevance_tier"])
            self.assertEqual(before, f["claims"])
            self.assertEqual(hashes, {path: sha256(path) for path in hashes})
            self.assertEqual("unknown", json.loads(json.loads(normalize_paper_reports([f["paper"]]))[0]["paper_checklist"])["same_gene_as_focal"])
            self.assertEqual(f["item"], evidence[0]["identity_assessment"])
            self.assertFalse(any("supported" in key for key in evidence[0]))
            original_checks = json.loads(evidence[0]["original_paper_report"]["paper_checklist"])
            new_checks = json.loads(evidence[0]["paper_checklist"])
            self.assertEqual({key: value for key, value in original_checks.items() if key != "same_gene_as_focal"},
                             {key: value for key, value in new_checks.items() if key != "same_gene_as_focal"})
            projection = json.loads(evidence[0]["homolog_projection_json"])
            self.assertEqual(("post_freeze_exact_identity", "g1", ["g1"]),
                             (projection["basis"], projection["candidate_gene_id"], projection["member_genes"]))

    def test_identity_projection_basis_requires_the_identity_contract(self):
        # The exact-identity basis is not a general bypass: a row that names it
        # without the sealed identity contract keeps the unresolved-projection T2.
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp))
            evidence, _ = self.load(f)
            forged = {key: value for key, value in evidence[0].items() if key != "post_freeze_identity_contract"}
            row = annotate_relevance([{"candidate_gene_id": "g1"}], f["claims"] + [forged])[0]
            self.assertNotEqual(NEW_TIERS[0], row["relevance_tier"])

    def test_same_gene_does_not_supply_missing_phenotype_equivalence(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), equivalent="no")
            evidence, _ = self.load(f)
            row = annotate_relevance([{"candidate_gene_id": "g1"}], f["claims"] + evidence)[0]
            self.assertNotEqual(NEW_TIERS[0], row["relevance_tier"])

    def test_repeated_identical_official_records_are_not_competing_identities(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self.fixture(Path(tmp))
            path = Path(fixture["bridge"]["manifest"])
            bridge = json.loads(path.read_text())
            first = bridge["retrieval"][0]
            duplicate = path.parent / "raw" / "repeated.txt"
            duplicate.write_bytes(Path(first["raw_path"]).read_bytes())
            bridge["retrieval"].append({**first, "raw_path": str(duplicate)})
            bridge["artifacts"][str(duplicate.relative_to(path.parent))] = sha256(duplicate)
            path.write_text(json.dumps(bridge))
            fixture["data"]["reference_bridge"] = descriptor(path)
            self.save(fixture)
            evidence, _ = self.load(fixture)
            self.assertEqual(1, len(evidence))
            duplicate.write_bytes(duplicate.read_bytes() + b"\n")
            bridge["retrieval"][-1]["raw_sha256"] = sha256(duplicate)
            bridge["artifacts"][str(duplicate.relative_to(path.parent))] = sha256(duplicate)
            path.write_text(json.dumps(bridge))
            fixture["data"]["reference_bridge"] = descriptor(path)
            self.save(fixture)
            with self.assertRaises(ValueError):
                self.load(fixture)

    def test_endpoint_mismatches_are_rejected_even_when_raw_reviewer_agrees(self):
        mutations = {"proposal_id": "other", "source_id": "PMID:7654321", "reported_gene": "GeneY",
                     "reported_species": "Plant beta", "request_id": "unknown", "candidate_gene_id": "outside",
                     "reference_taxid": "1", "focal_taxid": "1", "query_sequence_sha256": "0" * 64,
                     "reference_accession": "AA123456.2", "cds_protein_id": "OTHER.1", "sequence_reported_gene": "GeneY"}
        for field, value in mutations.items():
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                f["item"][field] = value
                self.save(f)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_hash_quote_case_and_raw_review_are_checked(self):
        for mutation in ("hash", "quote", "case", "raw_review", "source_binding", "cds_quote", "ancestor"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                if mutation == "hash":
                    f["data"]["reference_bridge"]["sha256"] = "0" * 64
                elif mutation == "quote":
                    f["item"]["paper_quote"]["text"] = "fabricated source"
                elif mutation == "case":
                    f["data"]["case_id"] = "other-case"
                elif mutation == "raw_review":
                    f["item"]["reason"] = "Changed after the reviewer ran"
                elif mutation == "source_binding":
                    f["item"]["paper_source"]["source_id"] = "PMID:9"
                elif mutation == "cds_quote":
                    f["item"]["cds_quote"]["end"] -= 1
                else:
                    Path(f["bridge"]["manifest"]).parent.joinpath("request_freeze.json").write_text("{}")
                self.save(f, sync_review=mutation != "raw_review")
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_other_copy_new_source_and_missing_review_axis_cannot_gain_t1(self):
        for mutation in ("other_copy", "new_source", "other_proposal", "no_axis", "unknown_axis"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                if mutation == "other_copy":
                    f["claims"][0]["candidate_gene_id"] = "g2"
                elif mutation == "new_source":
                    f["claims"] = [claim([{**f["paper"], "source_id": "PMID:9"}], candidate_gene_id="g1")]
                elif mutation == "other_proposal":
                    f["claims"][0]["hypothesis_id"] = "other-proposal"
                elif mutation == "no_axis":
                    f["item"]["review_axes"].pop("paper_to_accession")
                else:
                    f["item"]["review_axes"]["competitors_reviewed"] = "unknown"
                self.save(f)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_saved_paper_requires_matching_provider_record_title_and_scientific_span(self):
        for mutation in ("wrong_id", "wrong_title", "title_only"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                path = Path(f["item"]["paper_source"]["path"])
                record = json.loads(path.read_text())
                if mutation == "wrong_id":
                    record.update(id="7654321", pmid="7654321")
                elif mutation == "wrong_title":
                    record["title"] = "Another experiment"
                path.write_text(json.dumps(record), encoding="utf-8")
                f["item"]["paper_source"].update(descriptor(path))
                f["item"]["paper_quote"] = quote(path)
                if mutation == "title_only":
                    raw = path.read_text()
                    start = raw.index(record["title"])
                    f["item"]["paper_quote"] = {"start": start, "end": start + len(record["title"]), "text": record["title"]}
                self.save(f)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_materials_cannot_be_resealed_under_an_unchanged_raw_review(self):
        for material in ("paper_source", "paper_input", "atlas", "reference_bridge", "candidate_proteins"):
            with self.subTest(material=material), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                original_review = f["review"].read_bytes()
                entry = f["item"][material] if material == "paper_source" else f["data"][material]
                path = Path(entry["path"])
                if material == "paper_source":
                    value = json.loads(path.read_text())
                    value["abstractText"] = "The experiment concerns a different locus and does not connect GeneX to this accession."
                    path.write_text(json.dumps(value), encoding="utf-8")
                    f["item"]["paper_quote"] = quote(path)
                elif material == "candidate_proteins":
                    path.write_text(path.read_text() + "\n", encoding="utf-8")
                elif material in {"paper_input", "atlas"}:
                    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
                    table["replacement_material"] = "changed after review"
                    table.to_csv(path, sep="\t", index=False)
                else:
                    value = json.loads(path.read_text())
                    value["created_utc"] = "changed after review"
                    path.write_text(json.dumps(value), encoding="utf-8")
                entry.update(descriptor(path))
                self.save(f, sync_review=False)
                self.assertEqual(original_review, f["review"].read_bytes())
                with self.assertRaises(ValueError):
                    self.load(f)
                bridge = json.loads(Path(f["data"]["reference_bridge"]["path"]).read_text())
                f["item"]["material_bindings"] = identity_material_bindings(f["data"], f["item"], bridge["retrieval"][0]["raw_sha256"])
                self.save(f, sync_review=False)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_same_gene_ids_cannot_hide_different_canonical_protein_sequences(self):
        for sequence in (S + "A", "A" * len(S)):
            with self.subTest(length=len(sequence)), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                f["canonical"].write_text(f">g1\n{sequence}\n>g2\n{'A' * len(S)}\n>g3\n{'C' * len(S)}\n", encoding="utf-8")
                f["data"]["candidate_proteins"] = descriptor(f["canonical"])
                self.save(f)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_numeric_exact_hit_cannot_override_different_actual_source_and_target_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), targets={"p1": "A" * len(S)})
            with self.assertRaises(ValueError):
                self.load(f)

    def test_canonical_complete_copy_cannot_be_hidden_by_a_subset_bridge(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), complete_targets=False)
            f["canonical"].write_text(f">g1\n{S}\n>g2\n{S}\n>g3\n{'C' * len(S)}\n", encoding="utf-8")
            f["data"]["candidate_proteins"] = descriptor(f["canonical"])
            self.save(f)
            with self.assertRaises(ValueError):
                self.load(f)

    def test_actual_whole_proteome_copy_cannot_be_hidden_by_missing_alignment_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), targets={"p1": S, "p2": S}, hits=[_hit()])
            with self.assertRaises(ValueError):
                self.load(f)

    def test_partial_domain_duplicate_copy_and_unavailable_search_cannot_resolve_identity(self):
        cases = [{"hits": [_hit(slen=200)], "targets": {"p1": S * 2}},
                 {"hits": [_hit(qend=30, send=30)]},
                 {"hits": [_hit(), _hit(target="p2")], "targets": {"p1": S, "p2": S}},
                 {"hits": []}, {"hits": [], "state": "not_evaluated_blast_unavailable"}]
        for options in cases:
            with self.subTest(options=options), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp), **options)
                with self.assertRaises(ValueError):
                    self.load(f)

    def test_unknown_and_not_evaluated_are_retained_as_audit_without_positive_evidence(self):
        for decision in ("unknown", "not_evaluated"):
            with self.subTest(decision=decision), tempfile.TemporaryDirectory() as tmp:
                f = self.fixture(Path(tmp))
                f["item"]["decision"] = decision
                f["item"]["review_axes"]["locus_identity"] = "unknown"
                self.save(f)
                evidence, audit = self.load(f)
                self.assertEqual([], evidence)
                self.assertEqual(1, audit["assessment_count"])
                self.assertEqual(0, audit["resolved_evidence_count"])

    def test_distant_eligible_homolog_does_not_veto_unique_complete_exact_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), hits=[_hit(), _hit(target="p2", identity=60)],
                             targets={"p1": S, "p2": "A" * len(S)})
            evidence, _ = self.load(f)
            self.assertEqual(1, len(evidence))

    def test_identity_for_paper_a_cannot_borrow_paper_b_equivalent_phenotype(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = self.fixture(Path(tmp), equivalent="no")
            other = {**f["paper"], "source_id": "PMID:7654321"}
            checks = dict(f["paper"]["paper_checklist"], equivalent_target_phenotype="yes")
            other["paper_checklist"] = checks
            f["claims"].append(claim([other], candidate_gene_id="g1", hypothesis_id="h2"))
            evidence, _ = self.load(f)
            row = annotate_relevance([{"candidate_gene_id": "g1"}], f["claims"] + evidence)[0]
            self.assertNotEqual(NEW_TIERS[0], row["relevance_tier"])

    def test_multiple_reported_genes_in_one_annotation_model_remain_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = self.fixture(root)
            other = root / "fusion"
            other.mkdir()
            base = {"source_id": "PMID:1234567", "reference_accession": "AA123456.1",
                    "reference_taxid": "3888", "focal_taxid": "3888"}
            req, proteins, mapping, cache = _inputs(other,
                requests=[{**base, "request_id": "r1", "reported_gene": "GeneX", "cds_protein_id": "ABC12345.1"},
                          {**base, "request_id": "r2", "reported_gene": "GeneY", "cds_protein_id": "ABC12346.1"}],
                targets={"p1": S * 2}, genes=[("GeneX", "ABC12345.1", S), ("GeneY", "ABC12346.1", S)])
            with patch("candidate_qtl.reference_sequence_bridge._run_search", return_value={"state": "evaluated", "hits":
                       [_hit(slen=200), _hit(query="q000002", slen=200, sstart=101, send=200)]}):
                bridge = build_reference_sequence_bridge(req, proteins, mapping, other / "sequence", records_arg=cache)
            f["data"]["reference_bridge"] = descriptor(Path(bridge["manifest"]))
            manifest = json.loads(Path(bridge["manifest"]).read_text())
            f["item"]["cds_quote"] = quote(Path(manifest["retrieval"][0]["raw_path"]))
            self.save(f)
            correspondence = pd.read_csv(bridge["correspondence"], sep="\t")
            self.assertEqual({"multiple_reported_genes_in_model"}, set(correspondence.whole_annotation_model_identity))
            with self.assertRaises(ValueError):
                self.load(f)

    def test_real_handoff_and_cli_keep_original_generation_support_and_all_members(self):
        self._handoff_roundtrip(supported=False)

    def test_identity_update_keeps_already_verified_support_and_supported_union(self):
        self._handoff_roundtrip(supported=True)

    def _handoff_roundtrip(self, *, supported):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            identity_root, route_root = root / "identity", root / "route"
            identity_root.mkdir()
            route_root.mkdir()
            actual = {record.id: str(record.seq) for record in SeqIO.parse(
                str(route_a_fixtures.FIXTURE / "interval_proteins.fa"), "fasta")}
            f = self.fixture(identity_root, sequence=route_a_fixtures.S,
                             targets={"p1": actual["g1"], "p2": actual["g2"]})
            f["canonical"].write_text(f">g1\n{actual['g1']}\n>g2\n{actual['g2']}\n", encoding="utf-8")
            old, bridge, hypotheses, _ = route_a_fixtures.RouteAReviewRelevanceTests()._run(
                route_root, {"source_ids": "PMID:1234567", PAPER_FIELD: normalize_paper_reports([f["paper"]])})
            f["data"].update(case_id="route-a-review", atlas=descriptor(route_root / "genes.tsv"),
                              paper_input=descriptor(hypotheses, id_column="hypothesis_id"),
                              candidate_proteins=descriptor(f["canonical"]))
            self.save(f)
            frozen = {path: sha256(path) for path in root.rglob("*") if path.is_file()}
            arguments = dict(route_a_bridge_arg=bridge["bridge_tsv"], route_a_manifest_arg=bridge["manifest"],
                             review_candidates_arg=bridge["review_candidates"],
                             review_candidates_manifest_arg=bridge["review_manifest"],
                             interval_atlas_arg=str(route_root / "genes.tsv"), wet_budget=2, case_id="route-a-review")
            if supported:
                candidates = root / "supported-b.tsv"
                pd.DataFrame([{"candidate_gene_id": "g1", "route": "b", "support_status": "verified_semantic_support"}]).to_csv(candidates, sep="\t", index=False)
                arguments.update(route_b_candidates_arg=str(candidates), route_b_candidates_manifest_arg=str(
                    manifest_for(candidates, route="b", case_id="route-a-review")))
                old = wet_candidate_handoff(str(root / "wet-supported-old"), **arguments)
            new = wet_candidate_handoff(str(root / "wet-new"), identity_assessment_arg=str(f["path"]),
                identity_target_proteins_arg=f["data"]["candidate_proteins"]["path"], **arguments)
            old_rows = pd.read_csv(old["wet_candidates"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            new_rows = pd.read_csv(new["wet_candidates"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual(set(old_rows.index), set(new_rows.index))
            self.assertNotEqual(NEW_TIERS[0], old_rows.loc["g1", "relevance_tier"])
            self.assertEqual(NEW_TIERS[0], new_rows.loc["g1", "relevance_tier"])
            self.assertNotEqual(NEW_TIERS[0], new_rows.loc["g2", "relevance_tier"])
            summaries = json.loads(new_rows.loc["g1", "paper_evidence_summary_json"])
            self.assertEqual(1, len(summaries))
            self.assertEqual("T1", summaries[0]["tier"])
            self.assertEqual({"T1": 1}, json.loads(new_rows.loc["g1", "paper_tier_counts_json"]))
            for column in ("support_status", "route_a_supported", "route_b_supported", "route_c_supported"):
                if column in old_rows:
                    self.assertEqual(old_rows[column].to_dict(), new_rows[column].to_dict())
            self.assertEqual(Path(old["frozen_candidate_set"]).read_bytes(), Path(new["frozen_candidate_set"]).read_bytes())
            if supported:
                self.assertEqual("True", new_rows.loc["g1", "route_b_supported"])
            evidence = json.loads(new_rows.loc["g1", "relevance_evidence_json"])
            self.assertTrue(any(value.get("post_freeze_identity_contract") == IDENTITY_CONTRACT for value in evidence))
            self.assertEqual(frozen, {path: sha256(path) for path in frozen})
            cli_output = root / "wet-cli"
            argv = ["wet-candidate-handoff", "--out-dir", str(cli_output), "--case-id", "route-a-review",
                    "--wet-budget", "2", "--identity-assessment", str(f["path"])]
            argv.extend(["--identity-target-proteins", f["data"]["candidate_proteins"]["path"]])
            for option, key in (("--route-a-bridge", "route_a_bridge_arg"), ("--route-a-manifest", "route_a_manifest_arg"),
                                ("--review-candidates", "review_candidates_arg"),
                                ("--review-candidates-manifest", "review_candidates_manifest_arg"),
                                ("--interval-atlas", "interval_atlas_arg")):
                argv.extend([option, arguments[key]])
            if supported:
                argv.extend(["--route-b-candidates", arguments["route_b_candidates_arg"],
                             "--route-b-candidates-manifest", arguments["route_b_candidates_manifest_arg"]])
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(0, main(argv))
            cli_rows = pd.read_csv(cli_output / "wet_candidates.tsv", sep="\t", dtype=str).set_index("candidate_gene_id")
            self.assertEqual(NEW_TIERS[0], cli_rows.loc["g1", "relevance_tier"])
            self.assertEqual({"T1": 1}, json.loads(cli_rows.loc["g1", "paper_tier_counts_json"]))
            self.assertEqual(frozen, {path: sha256(path) for path in frozen})


if __name__ == "__main__":
    unittest.main()
