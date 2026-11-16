"""Phenotype-investigation layer in front of blind Route A (issue #215)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from candidate_qtl.llm_pathway import HYPOTHESIS_COLUMNS
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS
from candidate_qtl.mechanism_payload_adapter import TEMPLATE_24_COLUMNS
from tests.bound_claim_fixtures import bound_claim

MECHANISMS = [
    {"pathway_or_process": "starch biosynthesis", "gene_family": "SBE",
     "mechanism_statement": "loss of branching enzyme lowers amylopectin and gives wrinkled seed",
     "expected_role": "loss of function", "direction_consistency": "consistent with the reduced direction",
     "model_species_evidence": "Arabidopsis sbe2.1 sbe2.2 double mutant",
     "model_species_mutant_phenotype": "reduced starch branching", "source_ids": "PMID:1978859",
     "proposed_stable_id": "P30924", "reference_species": "Pisum sativum relative",
     "reference_gene_name": "SBEI", "selection_rationale": "closest characterised branching enzyme outside the focal species",
     "counter_hypothesis": "a sucrose transport lesion could give the same seed shape",
     "discriminating_observation": "measure amylopectin chain-length distribution", "confidence": "medium",
     "uncertainty": "cannot separate SBE paralogues from phenotype alone"},
    {"pathway_or_process": "sucrose transport", "gene_family": "SUT",
     "mechanism_statement": "impaired sucrose unloading reduces seed filling",
     "expected_role": "loss of function", "direction_consistency": "consistent",
     "model_species_evidence": "Medicago sut mutant", "model_species_mutant_phenotype": "reduced seed weight",
     "source_ids": "PMID:1978860", "proposed_stable_id": "Q9SXK1", "reference_species": "Arabidopsis thaliana",
     "reference_gene_name": "SUC2", "selection_rationale": "reviewed transporter with a seed-filling phenotype",
     "counter_hypothesis": "starch synthesis lesions produce a similar seed shape",
     "discriminating_observation": "assay apoplastic sucrose in developing seed", "confidence": "low",
     "uncertainty": "transport phenotypes are pleiotropic"},
    {"pathway_or_process": "cell wall biosynthesis", "gene_family": "CESA",
     "mechanism_statement": "weakened testa wall changes seed shape without changing starch",
     "expected_role": "loss of function", "direction_consistency": "consistent",
     "model_species_evidence": "Arabidopsis irx mutants", "model_species_mutant_phenotype": "collapsed xylem",
     "source_ids": "PMID:12172025", "proposed_stable_id": "Q9SWW6", "reference_species": "Arabidopsis thaliana",
     "reference_gene_name": "CESA9", "selection_rationale": "seed-coat expressed cellulose synthase",
     "counter_hypothesis": "wall phenotypes rarely give this seed shape alone",
     "discriminating_observation": "section the testa and measure wall thickness", "confidence": "low",
     "uncertainty": "little legume-specific evidence"},
    {"pathway_or_process": "storage reserve mobilization", "gene_family": "AMY",
     "mechanism_statement": "delayed reserve mobilization limits carbon supply during seed filling",
     "expected_role": "reserve mobilization", "direction_consistency": "consistent",
     "model_species_evidence": "Arabidopsis amylase mutant", "model_species_mutant_phenotype": "reduced reserve use during germination",
     "source_ids": "PMID:12172026", "proposed_stable_id": "P12345", "reference_species": "Arabidopsis thaliana",
     "reference_gene_name": "AMY3", "selection_rationale": "reserve-use enzyme with a seed-filling alternative",
     "counter_hypothesis": "the phenotype may arise from import rather than reserve mobilisation",
     "discriminating_observation": "measure reserve breakdown and carbon flux during seed filling", "confidence": "medium",
     "uncertainty": "reserve mobilisation is partly tissue and stage dependent"},
    {"pathway_or_process": "embryo growth regulation", "gene_family": "GRF",
     "mechanism_statement": "reduced embryo growth changes final seed size independently of starch production",
     "expected_role": "embryo growth regulation", "direction_consistency": "consistent",
     "model_species_evidence": "Arabidopsis growth regulator mutant", "model_species_mutant_phenotype": "smaller embryos and seeds",
     "source_ids": "PMID:12172027", "proposed_stable_id": "Q54321", "reference_species": "Arabidopsis thaliana",
     "reference_gene_name": "GRF5", "selection_rationale": "developmental growth regulator providing a tissue-level alternative",
     "counter_hypothesis": "smaller seeds may reflect reduced assimilate delivery rather than embryo growth",
     "discriminating_observation": "measure embryo cell proliferation and seed growth trajectory", "confidence": "medium",
     "uncertainty": "embryo growth and maternal supply are correlated in whole-seed assays"},
]


def _family_portfolio(count: int = 8, *, one_family: bool = False) -> list[dict]:
    rows: list[dict] = []
    for index in range(count):
        template = MECHANISMS[index % len(MECHANISMS)]
        rows.append({
            **template,
            "pathway_or_process": f"mechanism_{index + 1}",
            "gene_family": "SAME_FAMILY" if one_family else f"FAMILY_{index + 1}",
            "proposed_stable_id": f"TEST{index + 1:03d}",
        })
    return rows


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _payload(root: Path, mechanisms: list[dict] | None = None, name: str = "payload.json") -> Path:
    return _write(root, name, json.dumps({"hypotheses": mechanisms if mechanisms is not None else MECHANISMS}))


def _payload_23col(root: Path, name: str = "payload_23col.tsv") -> Path:
    columns = [column for column in MECHANISM_COLUMNS
               if column in TEMPLATE_24_COLUMNS]
    rows = []
    for record in MECHANISMS:
        row = {column: "" for column in columns}
        row.update({key: value for key, value in record.items() if key in columns})
        rows.append(row)
    path = root / name
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


LEGACY_7COL_HEADER = "case_id\tmechanism_id\tmechanism_name\tevidence_class\tsource_log_entry\tstatus\tselection_reason"
LEGACY_ROWS = [
    ["pea_a", "pea_a_M1", "Chalcone synthase (CHS) first committed step of anthocyanin biosynthesis",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M1", "verified",
     "Flower colour and seed-coat anthocyanin phenotype text maps to anthocyanin biosynthesis; "
     "CHS is the committed first step (corpus: At CHS, Mt CHS-class)."],
    ["pea_a", "pea_a_M2", "Dihydroflavonol 4-reductase (DFR) commits flux to anthocyanidin pigments",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M2", "verified",
     "Anthocyanin pigment loss maps to the late anthocyanin pathway; DFR is the committed step into anthocyanidins."],
    ["pea_a", "pea_a_M3", "R2R3-MYB transcription factor loss reduces anthocyanin pathway transcription",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M3", "verified",
     "Regulatory lesions reproduce the pigment phenotype in model species; MYB loss is the least-mutagenised step."],
]


def _legacy_payload(root: Path, rows: list[list[str]] | None = None,
                    header: str = LEGACY_7COL_HEADER, name: str = "payload_legacy.tsv") -> Path:
    lines = [header] + ["\t".join(row) for row in (rows if rows is not None else LEGACY_ROWS)]
    return _write(root, name, "\n".join(lines) + "\n")


def _freeze_manifest(root: Path, case_id: str, payload: Path, sha: str | None = None) -> Path:
    digest = sha if sha is not None else hashlib.sha256(payload.read_bytes()).hexdigest()
    return _write(root, "payload_freeze.tsv", "\n".join([
        "case_id\tmechanism_payload_sha256\tfrozen_at_utc",
        f"{case_id}\t{digest}\t2026-08-12T00:00:00Z",
    ]) + "\n")


def _investigate(root: Path, phenotype: str = "wrinkled seed with reduced starch content",
                 candidates: Path | None = None) -> Path:
    out = root / "investigation"
    argv = ["phenotype-investigate", "--out-dir", str(out), "--case-id", "case_r", "--focal-species", "Pisum sativum",
            "--phenotype", phenotype, "--direction", "reduced", "--inheritance", "recessive",
            "--tissue", "seed", "--stage", "seed filling", "--onset", "late", "--penetrance", "complete",
            "--synonym", "wrinkled seed", "--ontology-term", "TO:0000123",
            "--uncertainty", "starch content was scored on a single harvest"]
    if candidates:
        argv += ["--candidates", str(candidates)]
    assert main(argv) == 0
    return out


class InvestigationPacketTests(unittest.TestCase):
    def test_packet_records_the_full_phenotype_description_and_stays_interval_blind(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = _investigate(root)
            record = pd.read_csv(out / "phenotype_investigation.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            for column in ("direction", "inheritance", "affected_tissue", "developmental_stage",
                           "onset", "penetrance", "synonyms", "ontology_terms", "uncertainty"):
                self.assertTrue(str(record[column]).strip(), column)
            manifest = json.loads((out / "phenotype_investigation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual([], manifest["interval_blind"]["interval_arguments"])
            self.assertFalse(manifest["interval_blind"]["qtl_candidate_information_in_packet"])
            self.assertEqual(5, manifest["payload_contract"]["min_competing_mechanisms"])
            prompt = (out / "investigation_prompt.txt").read_text(encoding="utf-8")
            self.assertIn("counter-hypothesis", prompt)
            self.assertIn("discriminating observation", prompt)
            self.assertIn("target of 12 initial family-level hypotheses", prompt)
            self.assertIn("do not pad or reject", prompt)
            self.assertIn("causal-layer label", prompt)
            self.assertIn("at least five competing mechanisms", prompt)
            self.assertIn("2-3 functionally non-redundant", prompt)
            for level in ("sensing/timing", "synthesis/catabolism", "signal integration/regulation",
                          "downstream effector/transport", "developmental/tissue alternative"):
                self.assertIn(level, prompt)
            template = pd.read_csv(out / "mechanism_payload_template.tsv", sep="\t", dtype=str)
            self.assertIn("discriminating_observation", template.columns)
            self.assertIn("counter_hypothesis", template.columns)

    def test_a_packet_naming_an_interval_candidate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "candidates.tsv", "gene_id\nPsat_G0001\n")
            with self.assertRaises(ValueError) as context:
                _investigate(root, phenotype="wrinkled seed caused by Psat_G0001", candidates=candidates)
            self.assertIn("Psat_G0001", str(context.exception))

    def test_generation_stage_has_no_interval_argument(self) -> None:
        parser = build_parser()
        action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
        options = set()
        for sub_action in action.choices["phenotype-investigate"]._actions:
            options.update(sub_action.option_strings)
        for forbidden in ("--interval", "--atlas", "--gff3", "--interval-genes"):
            self.assertNotIn(forbidden, options)


class MechanismImportTests(unittest.TestCase):
    def test_pmc_namespace_alias_does_not_quarantine_a_mechanism(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sources = "PMID:1978859; PMC:PMC123456"
            record = {**MECHANISMS[0], "source_ids": sources}
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, [record])),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text())
            self.assertEqual(1, manifest["rows_imported"])
            self.assertEqual(0, manifest["rows_rejected"])
            rows = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            for field, value in record.items():
                self.assertEqual(value, rows.iloc[0][field])

    def test_generated_ids_are_disjoint_between_initial_and_gap_rounds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            initial_out = root / "initial_import"
            gap_out = root / "gap_import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, MECHANISMS[:2])),
                                      "--out-dir", str(initial_out), "--proposal-round", "1"]))
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, MECHANISMS[2:4])),
                                      "--out-dir", str(gap_out), "--proposal-round", "2"]))
            initial_ids = set(pd.read_csv(initial_out / "mechanism_hypotheses.tsv", sep="\t", dtype=str)["hypothesis_id"])
            gap_ids = set(pd.read_csv(gap_out / "mechanism_hypotheses.tsv", sep="\t", dtype=str)["hypothesis_id"])
            self.assertTrue(initial_ids.isdisjoint(gap_ids))
            self.assertTrue(all(value.startswith("mech_r2_") for value in gap_ids))

    def test_claim_extractions_in_the_payload_are_written_for_the_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {"hypotheses": MECHANISMS,
                       "claim_extractions": [{"citation": "PMID:1000003", "claim_type": "causal",
                                              "evidence_relation": "direct_genetic_test",
                                              "quoted_support": ("We list peptide-encoding loci "
                                                                 "annotated in the reference genome."),
                                              "locator": "abstract"}]}
            payload_path = _write(root, "payload_ext.json", json.dumps(payload))
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload_path),
                                      "--out-dir", str(out)]))
            extractions = pd.read_csv(out / "claim_extractions.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("PMID:1000003", extractions.iloc[0]["citation"])
            self.assertEqual("We list peptide-encoding loci annotated in the reference genome.",
                             extractions.iloc[0]["quoted_support"])
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["claim_extractions"]["supplied"])
            self.assertEqual(1, manifest["claim_extractions"]["row_count"])
            self.assertEqual(["citation", "claim_type", "evidence_relation", "quoted_support", "locator"],
                             list(extractions.columns))
            # The optional target survives the real mechanism-import writer and
            # reaches source-verify; legacy extractions above retain their schema.
            payload["hypotheses"] = [{**MECHANISMS[0], "hypothesis_id": "explicit-target"}]
            payload["claim_extractions"][0].update(
                hypothesis_id="explicit-target", citation="PMID:1978859")
            payload_path.write_text(json.dumps(payload), encoding="utf-8")
            scoped_out = root / "scoped-import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload_path),
                                      "--out-dir", str(scoped_out)]))
            scoped = pd.read_csv(scoped_out / "claim_extractions.tsv", sep="\t", dtype=str)
            self.assertEqual("explicit-target", scoped.loc[0, "hypothesis_id"])
            cache = root / "cache.json"
            cache.write_text(json.dumps([{
                "id": "1978859", "pmid": "1978859", "title": "SBEI",
                "abstractText": payload["claim_extractions"][0]["quoted_support"],
            }]), encoding="utf-8")
            verified = root / "source-verify"
            self.assertEqual(0, main([
                "source-verify", "--verification-policy", "content_v1", "--mechanisms",
                str(scoped_out / "mechanism_hypotheses.tsv"), "--out-dir", str(verified),
                "--cache", str(cache), "--claim-extraction", str(scoped_out / "claim_extractions.tsv"),
            ]))
            validation = pd.read_csv(verified / "source_validation.tsv", sep="\t", dtype=str)
            self.assertEqual("explicit-target", validation.loc[0, "hypothesis_id"])
            self.assertEqual("accepted", validation.loc[0, "extraction_status"])

    def test_a_payload_without_claim_extractions_reports_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)),
                                      "--out-dir", str(out)]))
            self.assertFalse((out / "claim_extractions.tsv").exists())
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["claim_extractions"]["supplied"])

    def test_non_verifiable_citation_tokens_are_quarantined_as_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            poisoned = [{**row, "source_ids": "input:type1_cases.tsv#q001/known_reference_function"}
                        for row in MECHANISMS]
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, poisoned)),
                                     "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text())
            self.assertEqual("not_evaluated_all_rows_rejected", manifest["import_state"])
            self.assertEqual(len(poisoned), manifest["rows_rejected"])
            self.assertEqual(0, manifest["rows_imported"])
            rejected = json.loads(Path(manifest["rejected_rows"]).read_text())
            self.assertEqual(poisoned[0], rejected[0]["raw_record"])

    def test_bad_citation_does_not_change_valid_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good = MECHANISMS[0]
            bad = {**MECHANISMS[1], "source_ids": "unsupported-db:synthetic123"}
            for name, records in (("alone", [good]), ("mixed", [good, bad])):
                main(["mechanism-import", "--payload", str(_payload(root, records)), "--out-dir", str(root / name)])
            alone = pd.read_csv(root / "alone/mechanism_hypotheses.tsv", sep="\t")
            mixed = pd.read_csv(root / "mixed/mechanism_hypotheses.tsv", sep="\t")
            pd.testing.assert_frame_equal(alone, mixed)

    def test_mixed_real_identifiers_pass_the_citation_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mixed = [{**MECHANISMS[0], "gene_family": f"MIXED_{index}",
                      "pathway_or_process": f"mixed_mechanism_{index}",
                      "proposed_stable_id": f"P3092{index}",
                      "source_ids": "PMID:1000001; UniProt:P12345"} for index in range(5)]
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, mixed)),
                                      "--out-dir", str(out)]))
            rows = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("PMID:1000001; UniProt:P12345", rows.iloc[0]["source_ids"])

    def test_llm_and_human_payloads_go_through_one_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigate(root) / "phenotype_investigation.tsv"
            llm_out = root / "import_llm"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)), "--out-dir", str(llm_out),
                                      "--investigation", str(investigation), "--payload-source", "llm",
                                      "--model", "test-model"]))
            llm = pd.read_csv(llm_out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")

            human_payload = root / "human.tsv"
            llm.drop(columns=["hypothesis_id", "parser_version", "review_status", "rank_contribution", "notes"]) \
               .assign(payload_source="human", model="").to_csv(human_payload, sep="\t", index=False)
            human_out = root / "import_human"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(human_payload), "--out-dir", str(human_out),
                                      "--investigation", str(investigation), "--payload-source", "human"]))
            human = pd.read_csv(human_out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(list(llm["pathway_or_process"]), list(human["pathway_or_process"]))
            self.assertEqual({"human"}, set(human["payload_source"]))
            for manifest_dir in (llm_out, human_out):
                manifest = json.loads((manifest_dir / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
                self.assertTrue(manifest["competing_mechanisms_met"])
                self.assertEqual(5, manifest["distinct_mechanisms"])
                self.assertEqual(5, manifest["hypotheses_with_counter_hypothesis"])
                self.assertEqual(5, manifest["hypotheses_with_discriminating_observation"])

    def test_a_single_mechanism_is_retained_with_count_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, MECHANISMS[:1])),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["competing_mechanisms_met"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("insufficient_competing_mechanisms", set(reviews["reason"]))

    def test_three_rows_are_retained_below_breadth_advisory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, MECHANISMS[:3])),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["competing_mechanisms_met"])
            self.assertEqual(5, manifest["min_competing_mechanisms"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("insufficient_competing_mechanisms", set(reviews["reason"]))

    def test_five_distinct_mechanisms_are_retained_with_breadth_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["competing_mechanisms_met"])
            self.assertEqual(5, manifest["distinct_mechanisms"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("family_node_portfolio_size_out_of_bounds", set(reviews["reason"]))

    def test_eight_distinct_family_nodes_are_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            investigation = _investigate(root) / "phenotype_investigation.tsv"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, _family_portfolio())),
                                      "--out-dir", str(out), "--investigation", str(investigation)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["competing_mechanisms_met"])
            self.assertTrue(manifest["family_node_portfolio_met"])
            self.assertEqual(8, manifest["distinct_gene_families"])

    def test_sixteen_family_nodes_are_accepted_as_oversized_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, _family_portfolio(16))),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["family_node_portfolio_met"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("family_node_portfolio_size_out_of_bounds", set(reviews["reason"]))

    def test_same_family_with_distinct_pathway_names_cannot_pad_portfolio(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, _family_portfolio(one_family=True))),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(8, manifest["distinct_mechanisms"])
            self.assertEqual(1, manifest["distinct_gene_families"])
            self.assertFalse(manifest["family_node_portfolio_met"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("duplicate_family_node", set(reviews["reason"]))

    def test_same_mechanism_gene_variants_do_not_satisfy_diversity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            variants = [{**MECHANISMS[0], "gene_family": f"SBE_variant_{index}",
                         "proposed_stable_id": f"P3092{index}"} for index in range(5)]
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, variants)),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["competing_mechanisms_met"])
            self.assertEqual(1, manifest["distinct_mechanisms"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("insufficient_competing_mechanisms", set(reviews["reason"]))

    def test_unsourced_hypotheses_are_flagged_not_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            thin = [{**MECHANISMS[0], "source_ids": "", "counter_hypothesis": "",
                     "discriminating_observation": "", "model_species_evidence": ""}, *MECHANISMS[1:]]
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, thin)), "--out-dir", str(out)]))
            reviews = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertEqual({"mechanism_without_source_id", "mechanism_without_counter_hypothesis",
                              "mechanism_without_discriminating_observation", "mechanism_without_model_species_evidence",
                              "mechanism_non_focal_target_not_evaluated"},
                             reviews)
            mechanisms = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(5, len(mechanisms))

    def test_model_prose_and_confidence_never_reach_a_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)), "--out-dir", str(out)]))
            mechanisms = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"none"}, set(mechanisms["rank_contribution"]))
            self.assertEqual({"pending_verification"}, set(mechanisms["review_status"]))
            for column in ("dry_score", "score", "weight", "rank"):
                self.assertNotIn(column, mechanisms.columns)

    def test_import_emits_the_existing_verification_chain_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)), "--out-dir", str(out)]))
            hypotheses = pd.read_csv(out / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(HYPOTHESIS_COLUMNS, list(hypotheses.columns))
            self.assertEqual({"phenotype_investigation:llm"}, set(hypotheses["proposal_source"]))

    def test_a_payload_naming_an_interval_candidate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "candidates.tsv", "gene_id\nPsat_G0001\n")
            leaky = [{**MECHANISMS[0], "reference_gene_name": "Psat_G0001"}, *MECHANISMS[1:]]
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(_payload(root, leaky)), "--out-dir", str(root / "import"),
                      "--candidates", str(candidates)])
            self.assertIn("Psat_G0001", str(context.exception))


class MechanismPayloadAdapterImportTests(unittest.TestCase):
    def test_json_source_id_list_reaches_source_verify_as_canonical_string(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_ids = [" PMID:1978859 ", "PMID:1978860", "PMID:1978859"]
            payload = _payload(root, [{**MECHANISMS[0], "source_ids": source_ids}])
            imported = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload),
                                      "--out-dir", str(imported)]))

            mechanisms = pd.read_csv(imported / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(" PMID:1978859 ;PMID:1978860;PMID:1978859",
                             mechanisms.loc[0, "source_ids"])

            cache = root / "cache"
            cache.mkdir()
            (cache / "records.json").write_text(json.dumps([
                {"id": "1978859", "pmid": "1978859", "title": "SBE starch branching enzyme",
                 "abstractText": "Loss of function in SBEI changes starch biosynthesis.",
                 "journalTitle": "Plant Journal", "pubYear": "1990", "authorString": "Example"},
                {"id": "1978860", "pmid": "1978860", "title": "SBE seed starch study",
                 "abstractText": "SBEI controls starch biosynthesis during seed filling.",
                 "journalTitle": "Plant Journal", "pubYear": "1991", "authorString": "Example"},
            ]), encoding="utf-8")
            verified = root / "source_verify"
            self.assertEqual(0, main(["source-verify", "--verification-policy", "content_v1", "--mechanisms",
                                      str(imported / "mechanism_hypotheses.tsv"),
                                      "--out-dir", str(verified), "--cache", str(cache)]))
            status = pd.read_csv(verified / "hypothesis_source_status.tsv",
                                 sep="\t", dtype=str).fillna("")
            self.assertEqual("all_verified", status.loc[0, "source_verification_status"])
            validation = pd.read_csv(verified / "source_validation.tsv",
                                     sep="\t", dtype=str).fillna("")
            self.assertEqual({"PMID:1978859", "PMID:1978860"},
                             set(validation["citation"]))

    def test_legacy_7col_payload_imports_with_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["hypothesis_count"])
            self.assertEqual("legacy_v1_7col", manifest["accepted_payload_schema"])
            self.assertEqual("legacy_v1_7col", manifest["input_schema_version"])
            self.assertEqual("mechanism-payload-adapter-1.0.0", manifest["adapter_version"])
            self.assertEqual(3, manifest["rows_read"])
            self.assertEqual(3, manifest["rows_imported"])
            self.assertEqual(0, manifest["rows_rejected"])
            self.assertIsNone(manifest["freeze_verification"])
            self.assertIsNotNone(manifest["adapter_map"])
            mechanisms_bytes = (out / "mechanism_hypotheses.tsv").read_bytes()
            self.assertEqual(hashlib.sha256(mechanisms_bytes).hexdigest(), manifest["output_sha256"])
            mechanisms = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(mechanisms))
            self.assertEqual({"pea_a"}, set(mechanisms["case_id"]))
            for column in ("mechanism_statement", "selection_rationale", "source_ids"):
                for value in mechanisms[column]:
                    self.assertTrue(str(value).strip(), column)
            self.assertEqual({"pending_verification"}, set(mechanisms["review_status"]))
            self.assertEqual({"none"}, set(mechanisms["rank_contribution"]))
            self.assertEqual({""}, set(mechanisms["pathway_or_process"]) | set(mechanisms["gene_family"]))
            adapter_map = pd.read_csv(out / "mechanism_payload_adapter_map.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(list(adapter_map["mechanism_id"]), ["pea_a_M1", "pea_a_M2", "pea_a_M3"])
            self.assertEqual(list(adapter_map["hypothesis_id"]), ["mech_001", "mech_002", "mech_003"])
            self.assertEqual({"corpus_ortholog_model_plant"}, set(adapter_map["evidence_class"]))
            self.assertEqual({"verified"}, set(adapter_map["status"]))
            self.assertEqual({"source_log.tsv#pea_a_M1", "source_log.tsv#pea_a_M2", "source_log.tsv#pea_a_M3"},
                             set(adapter_map["source_log_entry"]))
            self.assertEqual(hashlib.sha256((out / "mechanism_payload_adapter_map.tsv").read_bytes()).hexdigest(),
                             manifest["adapter_map"]["sha256"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("no_parsable_mechanism", set(reviews["reason"]))
            counts = reviews["reason"].value_counts()
            expected = {"mechanism_missing_required_field": 3, "mechanism_without_counter_hypothesis": 3,
                        "mechanism_without_discriminating_observation": 3,
                        "mechanism_without_model_species_evidence": 3, "mechanism_without_stable_accession": 3}
            for reason, count in expected.items():
                self.assertEqual(count, counts.get(reason, 0), reason)
            self.assertEqual(15, sum(counts.get(reason, 0) for reason in expected))
            missing = reviews[reviews["reason"] == "mechanism_missing_required_field"]["required_next_action"]
            for action in missing:
                self.assertIn("pathway_or_process", action)
                self.assertIn("gene_family", action)
                self.assertNotIn("mechanism_statement", action)
            self.assertEqual({"unknown: pathway_or_process, gene_family"}, set(mechanisms["notes"]))

    def test_canonical_24col_import_matches_golden(self) -> None:
        golden = Path(__file__).parent / "fixtures" / "mechanism_import_24col_golden.tsv"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload_23col(root)),
                                      "--out-dir", str(out)]))
            self.assertEqual(golden.read_bytes(), (out / "mechanism_hypotheses.tsv").read_bytes())
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("canonical_24_col_template", manifest["accepted_payload_schema"])
            self.assertEqual("none", manifest["adapter_version"])
            self.assertEqual(5, manifest["rows_read"])
            self.assertEqual(5, manifest["rows_imported"])

    def test_canonical_subset_12col_unchanged(self) -> None:
        columns = ["pathway_or_process", "gene_family", "mechanism_statement", "source_ids",
                   "counter_hypothesis", "discriminating_observation", "model_species_evidence",
                   "model_species_mutant_phenotype", "proposed_stable_id", "reference_species",
                   "reference_gene_name", "selection_rationale"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for record in MECHANISMS:
                rows.append("\t".join(str(record[column]) for column in columns))
            payload = _write(root, "payload_12col.tsv", "\t".join(columns) + "\n" + "\n".join(rows) + "\n")
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("canonical_subset", manifest["accepted_payload_schema"])
            mechanisms = pd.read_csv(out / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(5, len(mechanisms))
            self.assertTrue(manifest["competing_mechanisms_met"])

    def test_unknown_columns_fail_closed_with_column_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload_unknown.tsv", "\n".join([
                "pathway_or_process\tgene_family\tmechanism_statement\tmystery_column",
                "starch biosynthesis\tSBE\tbranching enzyme loss gives wrinkled seed\tsomething",
            ]) + "\n")
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import")])
            message = str(context.exception)
            self.assertIn("mystery_column", message)
            self.assertIn("hypothesis_id", message)
            self.assertIn("missing canonical column", message)

    def test_legacy_row_missing_mechanism_id_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            broken = [list(row) for row in LEGACY_ROWS]
            broken[1][1] = ""
            payload = _legacy_payload(root, broken)
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import")])
            message = str(context.exception)
            self.assertIn("legacy_v1_7col", message)
            self.assertIn("row 2", message)
            self.assertIn("mechanism_id is empty", message)

    def test_duplicate_mechanism_id_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            broken = [list(row) for row in LEGACY_ROWS]
            broken[2][1] = broken[0][1]
            payload = _legacy_payload(root, broken)
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import")])
            message = str(context.exception)
            self.assertIn("duplicate mechanism_id 'pea_a_M1'", message)
            self.assertIn("first at row 1", message)
            self.assertIn("row 3", message)

    def test_forbidden_truth_column_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root, header=LEGACY_7COL_HEADER + "\ttruth_gene",
                                      name="payload_leaky.tsv")
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import")])
            message = str(context.exception)
            self.assertIn("forbidden column 'truth_gene'", message)
            self.assertIn("leakage class", message)

    def test_freeze_manifest_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            freeze = _freeze_manifest(root, "pea_a", payload, sha="0" * 64)
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import"),
                      "--case-id", "pea_a", "--mechanism-freeze-manifest", str(freeze)])
            message = str(context.exception)
            self.assertIn("mechanism payload freeze mismatch", message)
            self.assertIn("expected sha256", message)

    def test_freeze_manifest_match_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            freeze = _freeze_manifest(root, "pea_a", payload)
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out),
                                      "--case-id", "pea_a", "--mechanism-freeze-manifest", str(freeze)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            expected = hashlib.sha256(payload.read_bytes()).hexdigest()
            self.assertEqual({"source": str(freeze), "expected_sha256": expected,
                              "actual_sha256": expected, "matched": True},
                             manifest["freeze_verification"])

    def test_jsonl_content_in_txt_imports_and_freeze_derives_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload.txt", "\n".join(
                json.dumps({**row, "case_id": "pea_a"}) for row in MECHANISMS
            ) + "\n")
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            freeze = _write(root, "freeze.tsv", f"case_id\tmechanism_payload_sha256\npea_a\t{digest}\n")
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out),
                                      "--mechanism-freeze-manifest", str(freeze)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("jsonl", manifest["accepted_payload_schema"])
            self.assertEqual(5, manifest["rows_imported"])
            self.assertEqual("pea_a", manifest["case_id"])
            self.assertTrue(manifest["freeze_verification"]["matched"])

    def test_empty_payload_file_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload_empty.tsv", "")
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import")])
            self.assertIn("file is empty", str(context.exception))

    def test_missing_payload_file_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError) as context:
                main(["mechanism-import", "--payload", str(root / "no_such.tsv"),
                      "--out-dir", str(root / "import")])
            self.assertIn("file not found or unreadable", str(context.exception))

    def test_headers_only_canonical_payload_stays_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _write(root, "payload_headers.tsv", "\t".join(MECHANISM_COLUMNS) + "\n")
            out = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out)]))
            manifest = json.loads((out / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("canonical_29_col", manifest["accepted_payload_schema"])
            self.assertEqual(0, manifest["hypothesis_count"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("no_parsable_mechanism", set(reviews["reason"]))


class AttributionTests(unittest.TestCase):
    def test_attribution_names_the_stage_that_stopped_the_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            imported = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, MECHANISMS[:1])),
                                      "--out-dir", str(imported)]))
            out = root / "attribution"
            self.assertEqual(0, main(["investigation-attribution", "--out-dir", str(out),
                                      "--mechanisms", str(imported / "mechanism_hypotheses.tsv")]))
            manifest = json.loads((out / "investigation_attribution.json").read_text(encoding="utf-8"))
            self.assertEqual("phenotype_research", manifest["first_incomplete_stage"])
            self.assertEqual("stopped_at:phenotype_research", manifest["outcome"])
            self.assertFalse(manifest["truth_used"])
            table = pd.read_csv(out / "investigation_attribution.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["phenotype_research", "verification", "bait_freeze", "interval_match",
                              "candidate_recovery"], list(table["stage"]))
            self.assertEqual("not_reached", table.set_index("stage").loc["verification", "status"])

    def test_a_researched_phenotype_stage_passes_and_the_next_stage_is_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            imported = root / "import"
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root)), "--out-dir", str(imported)]))
            out = root / "attribution"
            self.assertEqual(0, main(["investigation-attribution", "--out-dir", str(out),
                                      "--mechanisms", str(imported / "mechanism_hypotheses.tsv")]))
            manifest = json.loads((out / "investigation_attribution.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", manifest["stages"]["phenotype_research"]["status"])
            self.assertEqual("verification", manifest["first_incomplete_stage"])


class GuidedChainTests(unittest.TestCase):
    def _record_cache(self, root: Path) -> Path:
        cache = root / "records"
        cache.mkdir()
        for accession, name in (("P30924", "SBEI_PEA"), ("Q9SXK1", "SUC2_ARATH"), ("Q9SWW6", "CESA9_ARATH")):
            (cache / f"{accession}.json").write_text(json.dumps({
                "primaryAccession": accession, "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "uniProtkbId": name,
                "genes": [{"geneName": {"value": name.split("_")[0]}}],
                "organism": {"scientificName": "Arabidopsis thaliana"},
                "sequence": {"value": "MASTPRQTVSEK" * 8},
                "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "1978859"}]}}],
            }), encoding="utf-8")
        return cache

    def test_profile_chains_research_into_verification_and_the_frozen_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root, "genes.tsv", "gene_id\tannotation\nG_SBE\tstarch branching enzyme\nG_X\thypothetical protein\n")
            route = _write(root, "route.tsv",
                           "phenotype_token\tpathway_family\tannotation_keywords\nstarch\tstarch_biosynthesis\tstarch|branching enzyme\n")
            predict = _write(root, "predict.tsv", "\n".join([
                "case_id\tphenotype_text\tinterval\tinterval_genes",
                "case_r\twrinkled seed with reduced starch content\tchr5:1-100000\tgenes.tsv",
            ]) + "\n")
            # every reference species is Arabidopsis in the cache, so declare it
            layers = ["signaling", "regulation", "alternative", "alternative", "execution"]
            mechanisms = [
                {**record, "causal_layer": layers[index],
                 "reference_species": "Arabidopsis thaliana", "source_ids": "",
                 "proposed_stable_id": ""}
                for index, record in enumerate(MECHANISMS)
            ]
            mechanisms.extend({
                "causal_layer": "signaling", "pathway_or_process": f"complement {index}",
                "gene_family": f"ComplementFamily{index}", "mechanism_statement": f"complement mechanism {index}",
                "counter_hypothesis": "another layer may be causal",
                "discriminating_observation": f"measure complement {index}",
                "reference_species": "Arabidopsis thaliana",
                "reference_gene_name": f"ComplementRef{index}",
                "source_ids": "", "proposed_stable_id": "",
            } for index in range(1, 8))
            localization = _write(root, "localization.json", json.dumps({
                "causal_localization": "signaling", "localization_confidence": "high",
                "supporting_phrase": "wrinkled seed",
            }))
            bindings = root / "bound_claims.tsv"
            pd.DataFrame([
                bound_claim(
                    route="a", claim_id=f"mech_{index:03d}", citation="PMID:1978859",
                    subject_gene=accession, organism="Arabidopsis thaliana",
                    phenotype="wrinkled seed with reduced starch content",
                )
                for index, accession in enumerate(("P30924", "Q9SXK1", "Q9SWW6"), start=1)
            ], columns=BOUND_CLAIM_COLUMNS).to_csv(
                bindings, sep="\t", index=False, lineterminator="\n"
            )
            out = root / "guided"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "audit-predict", "--out-dir", str(out), "--execute",
                "--case-id", "case_r", "--focal-species", "Pisum sativum",
                "--phenotype", "wrinkled seed with reduced starch content",
                "--mechanism-payload", str(_payload(root, mechanisms)),
                "--phenotype-localization-payload", str(localization),
                "--record-cache", str(self._record_cache(root)),
                "--bound-claims", str(bindings),
                "--predict-manifest", str(predict), "--blind-route", str(route),
            ]))
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            executed = {record["command"]: record["status"] for record in profile["execution_results"]}
            for command in ("phenotype-investigate", "phenotype-localization-packet",
                            "phenotype-localization-freeze", "mechanism-import", "mechanism-coverage-audit",
                            "reference-gene-resolve", "hypothesis-verify",
                            "bait-panel-freeze", "audit-predict", "investigation-attribution"):
                self.assertEqual("completed", executed[command], command)
            resolution = json.loads(
                (out / "03b_reference_gene_resolution" / "reference_gene_resolution_manifest.json")
                .read_text(encoding="utf-8")
            )
            self.assertEqual(12, resolution["hypothesis_count"])
            self.assertEqual(3, resolution["resolution_origin_counts"]["provider_resolved"])
            self.assertEqual(9, resolution["resolution_origin_counts"]["not_resolved"])
            localization_manifest = out / "02c_phenotype_localization_freeze" / "phenotype_localization_freeze_manifest.json"
            import_manifest = json.loads(
                (out / "03_mechanism_import" / "mechanism_import_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                hashlib.sha256(localization_manifest.read_bytes()).hexdigest(),
                import_manifest["localization_freeze"]["sha256"],
            )
            panel = pd.read_csv(out / "05_bait_panel_freeze" / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(panel))
            self.assertEqual({"none"}, set(panel["rank_contribution"]))
            bait_manifest = json.loads((out / "05_bait_panel_freeze" / "bait_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("frozen", bait_manifest["panel_status"])
            self.assertTrue(bait_manifest["panel_sha256"])
            attribution = json.loads((out / "07_attribution" / "investigation_attribution.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", attribution["stages"]["phenotype_research"]["status"])
            self.assertEqual("complete", attribution["stages"]["bait_freeze"]["status"])
            self.assertTrue((out / "06_prediction_freeze" / "freeze_manifest.json").exists())

    def test_a_missing_payload_leaves_the_downstream_lane_unavailable_not_failed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "guided"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "audit-predict", "--out-dir", str(out), "--execute",
                "--case-id", "case_r", "--focal-species", "Pisum sativum", "--phenotype", "wrinkled seed",
            ]))
            profile = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            executed = {record["command"]: record["status"] for record in profile["execution_results"]}
            self.assertEqual("completed", executed["phenotype-investigate"])
            self.assertEqual("not_runnable_missing_input", executed["mechanism-import"])
            self.assertEqual("not_runnable_upstream_unavailable:mechanism-coverage-audit",
                             executed["reference-gene-resolve"])
            self.assertEqual("not_runnable_upstream_unavailable:reference-gene-resolve",
                             executed["bait-panel-freeze"])
            self.assertNotIn("failed", set(executed.values()))
            # attribution still runs and reports that the chain never started
            self.assertEqual("completed", executed["investigation-attribution"])
            attribution = json.loads((out / "07_attribution" / "investigation_attribution.json").read_text(encoding="utf-8"))
            self.assertEqual("phenotype_research", attribution["first_incomplete_stage"])


if __name__ == "__main__":
    unittest.main()
