from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.workflow_foundation import WorkflowState, write_foundation_outputs


class WorkflowStateTests(unittest.TestCase):
    def test_same_input_is_reused_and_resource_content_change_invalidates(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            input_path = root / "candidates.tsv"
            resource_path = root / "reference.faa"
            output_path = root / "hits.tsv"
            input_path.write_text("gene_id\nGeneA\n", encoding="utf-8")
            resource_path.write_text(">Ref\nAAAA\n", encoding="utf-8")
            output_path.write_text("GeneA\tRef\n", encoding="utf-8")

            state_path = root / "workflow_state.json"
            state = WorkflowState(state_path)
            parameters = {"provider": "blastp", "version": "2.17", "threads": 1}
            first_signature = state.signature([input_path, resource_path], parameters)
            state.record("blastp:model", first_signature, "complete", [output_path])

            reloaded = WorkflowState(state_path)
            same_signature = reloaded.signature([input_path, resource_path], parameters)
            self.assertEqual(first_signature, same_signature)
            self.assertTrue(reloaded.can_reuse("blastp:model", same_signature, [output_path]))

            # Same path and same byte length, different content: the resource hash
            # must invalidate the step rather than trusting file existence alone.
            resource_path.write_text(">Ref\nGGGG\n", encoding="utf-8")
            changed_signature = reloaded.signature([input_path, resource_path], parameters)
            self.assertNotEqual(first_signature, changed_signature)
            self.assertFalse(reloaded.can_reuse("blastp:model", changed_signature, [output_path]))


class FoundationOutputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidate_proteins = self.root / "candidate_proteins.faa"
        self.candidate_proteins.write_text(">GeneA\nMAAA\n>GeneB\nMBBB\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _panel(self, name: str, assembly: str, labels: list[str]) -> tuple[dict[str, object], Path, dict[str, Path]]:
        panel_path = self.root / f"{name}.yaml"
        panel_path.write_text(f"panel_id: {name}\n", encoding="utf-8")
        resources: dict[str, Path] = {}
        references: list[dict[str, object]] = []
        for label in labels:
            resource = self.root / f"{name}_{label}.faa"
            resource.write_text(f">{label}_record\nMAAA\n", encoding="utf-8")
            resources[label] = resource
            references.append(
                {
                    "label": label,
                    "role": "functional_reference",
                    "species": f"Species {label}",
                    "assembly_id": f"{label}_assembly",
                    "annotation_release": f"{label}_annotation_v1",
                    "resource_version": f"{label}_resource_v1",
                    "proteome": str(resource),
                    "_resolved_proteome": str(resource),
                }
            )
        return (
            {
                "native": {
                    "species": "Pisum sativum",
                    "assembly_id": assembly,
                    "annotation_release": "native_annotation_v1",
                },
                "references": references,
            },
            panel_path,
            resources,
        )

    @staticmethod
    def _candidates(*gene_ids: str) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "gene_id": gene_id,
                    "representative_transcript": f"{gene_id}.t1",
                }
                for gene_id in gene_ids
            ]
        )

    @staticmethod
    def _empty_functions() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "candidate_gene_id",
                "reference_label",
                "reference_gene_id",
                "source_id",
                "retrieval_status",
                "review_status",
            ]
        )

    @staticmethod
    def _empty_trait_matches() -> pd.DataFrame:
        return pd.DataFrame(
            columns=["gene_id", "reference_label", "reference_gene_id", "matched_terms", "source_id", "review_status"]
        )

    @staticmethod
    def _empty_variants() -> pd.DataFrame:
        return pd.DataFrame(columns=["gene_id", "chrom", "pos", "ref", "alt", "review_status"])

    @staticmethod
    def _empty_evidence() -> pd.DataFrame:
        return pd.DataFrame(columns=["gene_id", "evidence_code", "source_type", "source_id"])

    @staticmethod
    def _ranked(*gene_ids: str) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "gene_id": gene_id,
                    "automated_rank": rank,
                    "automated_score": 0.0,
                    "dry_group": "evidence_insufficient",
                    "score_breakdown": "no_scored_evidence:0",
                }
                for rank, gene_id in enumerate(gene_ids, start=1)
            ]
        )

    def _write(
        self,
        *,
        name: str,
        panel: dict[str, object],
        panel_path: Path,
        orthology: pd.DataFrame,
        provider_runs: list[dict[str, object]],
        candidates: pd.DataFrame | None = None,
        functions: pd.DataFrame | None = None,
        trait_matches: pd.DataFrame | None = None,
        evidence: pd.DataFrame | None = None,
        ranked: pd.DataFrame | None = None,
    ) -> dict[str, str]:
        out_dir = self.root / name
        out_dir.mkdir()
        candidates = candidates if candidates is not None else self._candidates("GeneA")
        return write_foundation_outputs(
            out_dir=out_dir,
            candidates=candidates,
            protein_ids=set(candidates["gene_id"]),
            panel=panel,
            panel_path=panel_path,
            orthology=orthology,
            functions=functions if functions is not None else self._empty_functions(),
            trait_matches=trait_matches if trait_matches is not None else self._empty_trait_matches(),
            variants=self._empty_variants(),
            evidence=evidence if evidence is not None else self._empty_evidence(),
            ranked=ranked if ranked is not None else self._ranked(*candidates["gene_id"]),
            provider_runs=provider_runs,
            provider_version="blastp 2.17",
            inputs={"candidate_proteins": self.candidate_proteins, "candidates": None},
            mode="local",
        )

    def test_canonical_identity_namespaces_and_one_to_many_are_preserved(self) -> None:
        panel, panel_path, resources = self._panel("identity", "native_v1", ["model_a", "model_b"])
        orthology = pd.DataFrame(
            [
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": "model_a",
                    "reference_species": "Species model_a",
                    "reference_gene_id": "SharedID",
                    "relationship": "one_to_many",
                    "orthology_confidence": "low",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "none_detected",
                    "mapping_status": "mapped",
                    "failure_reason": "",
                },
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": "model_a",
                    "reference_species": "Species model_a",
                    "reference_gene_id": "AlternativeID",
                    "relationship": "one_to_many",
                    "orthology_confidence": "low",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "none_detected",
                    "mapping_status": "mapped",
                    "failure_reason": "",
                },
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": "model_b",
                    "reference_species": "Species model_b",
                    "reference_gene_id": "SharedID",
                    "relationship": "mapped_non_rbh",
                    "orthology_confidence": "medium",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "none_detected",
                    "mapping_status": "mapped",
                    "failure_reason": "",
                },
            ]
        )
        runs = [
            {"label": label, "status": "completed", "annotation_status": "not_run", "blast": {"output": str(resource)}}
            for label, resource in resources.items()
        ]
        outputs = self._write(name="identity_v1", panel=panel, panel_path=panel_path, orthology=orthology, provider_runs=runs)
        aliases = pd.read_csv(outputs["entity_aliases"], sep="\t", dtype=str).fillna("")
        external = aliases[aliases["entity_type"].eq("database_record")]

        self.assertEqual(3, len(external))
        self.assertEqual(
            {("model_a", "SharedID"), ("model_a", "AlternativeID"), ("model_b", "SharedID")},
            set(zip(external["namespace"], external["external_id"])),
        )
        model_a = external[external["namespace"].eq("model_a")]
        self.assertEqual({"one_to_many"}, set(model_a["relationship"]))
        self.assertEqual({"one_to_many"}, set(model_a["mapping_status"]))
        self.assertEqual({"native_v1:GeneA"}, set(aliases["canonical_candidate_id"]))

        # The same bare native gene ID in another assembly must not collide with
        # the canonical identity emitted above.
        panel_v2, panel_v2_path, _ = self._panel("identity_v2", "native_v2", [])
        outputs_v2 = self._write(
            name="identity_v2_out",
            panel=panel_v2,
            panel_path=panel_v2_path,
            orthology=pd.DataFrame(columns=orthology.columns),
            provider_runs=[],
        )
        aliases_v2 = pd.read_csv(outputs_v2["entity_aliases"], sep="\t", dtype=str).fillna("")
        self.assertEqual({"native_v2:GeneA"}, set(aliases_v2["canonical_candidate_id"]))
        self.assertTrue(set(aliases["canonical_candidate_id"]).isdisjoint(aliases_v2["canonical_candidate_id"]))

    def test_coverage_distinguishes_no_hit_execution_failed_and_ambiguous(self) -> None:
        panel, panel_path, resources = self._panel("coverage", "native_v1", ["nohit", "failed", "ambiguous"])
        orthology = pd.DataFrame(
            [
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": "nohit",
                    "reference_species": "Species nohit",
                    "reference_gene_id": "",
                    "relationship": "unmapped",
                    "orthology_confidence": "",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "",
                    "mapping_status": "unmapped",
                    "failure_reason": "no_quality_passing_hit",
                },
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": "failed",
                    "reference_species": "Species failed",
                    "reference_gene_id": "",
                    "relationship": "unmapped",
                    "orthology_confidence": "",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "",
                    "mapping_status": "provider_failure",
                    "failure_reason": "blastp exited 2",
                },
                *[
                    {
                        "candidate_gene_id": "GeneA",
                        "reference_label": "ambiguous",
                        "reference_species": "Species ambiguous",
                        "reference_gene_id": reference_gene_id,
                        "relationship": "one_to_many",
                        "orthology_confidence": "low",
                        "mapping_providers": "blastp",
                        "paralog_ambiguity": "top2_near_tie",
                        "mapping_status": "mapped",
                        "failure_reason": "",
                    }
                    for reference_gene_id in ("Ref1", "Ref2")
                ],
            ]
        )
        runs = [
            {"label": "nohit", "status": "completed", "annotation_status": "not_run", "blast": {"output": str(resources["nohit"])}},
            {"label": "failed", "status": "failed", "error": "blastp exited 2", "annotation_status": "not_run"},
            {"label": "ambiguous", "status": "completed", "annotation_status": "not_run", "blast": {"output": str(resources["ambiguous"])}},
        ]
        outputs = self._write(name="coverage_out", panel=panel, panel_path=panel_path, orthology=orthology, provider_runs=runs)
        coverage = pd.read_csv(outputs["candidate_provider_coverage"], sep="\t", dtype=str).fillna("")
        status_by_provider = dict(zip(coverage["provider_name"], coverage["coverage_status"]))
        self.assertEqual("no_hit", status_by_provider["blastp:nohit"])
        self.assertEqual("execution_failed", status_by_provider["blastp:failed"])
        self.assertEqual("ambiguous", status_by_provider["blastp:ambiguous"])

    def test_categories_and_lineage_connect_raw_mapping_to_evidence_and_rank(self) -> None:
        panel, panel_path, resources = self._panel("lineage", "native_v1", ["model_a", "model_b"])
        orthology = pd.DataFrame(
            [
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": label,
                    "reference_species": f"Species {label}",
                    "reference_gene_id": "SharedID",
                    "relationship": "mapped_non_rbh",
                    "orthology_confidence": "medium",
                    "mapping_providers": "blastp",
                    "paralog_ambiguity": "none_detected",
                    "mapping_status": "mapped",
                    "failure_reason": "",
                }
                for label in ("model_a", "model_b")
            ]
        )
        functions = pd.DataFrame(
            [
                {
                    "candidate_gene_id": "GeneA",
                    "reference_label": label,
                    "reference_gene_id": "SharedID",
                    "source_id": f"PMID:{index}",
                    "retrieval_status": "source_retrieved",
                    "review_status": "unreviewed",
                    "gene_description": f"annotation from {label}",
                }
                for index, label in enumerate(("model_a", "model_b"), start=1)
            ]
        )
        trait_matches = pd.DataFrame(
            [
                {
                    "gene_id": "GeneA",
                    "reference_label": "model_a",
                    "reference_gene_id": "SharedID",
                    "matched_terms": "dwarf stem",
                    "source_id": "PMID:1",
                    "review_status": "unreviewed",
                }
            ]
        )
        evidence = pd.DataFrame(
            [
                {
                    "gene_id": "GeneA",
                    "evidence_code": "automated_candidate_function_match",
                    "source_type": "candidate_first_reference_annotation",
                    "source_id": "PMID:1",
                }
            ]
        )
        ranked = pd.DataFrame(
            [
                {
                    "gene_id": "GeneA",
                    "automated_rank": 1,
                    "automated_score": 0.625,
                    "dry_group": "reserve",
                    "score_breakdown": "automated_candidate_function_match/support:+0.625",
                }
            ]
        )
        runs = [
            {"label": label, "status": "completed", "annotation_status": "complete", "blast": {"output": str(resource)}}
            for label, resource in resources.items()
        ]
        outputs = self._write(
            name="lineage_out",
            panel=panel,
            panel_path=panel_path,
            orthology=orthology,
            provider_runs=runs,
            functions=functions,
            trait_matches=trait_matches,
            evidence=evidence,
            ranked=ranked,
        )
        observations = pd.read_csv(outputs["normalized_observations"], sep="\t", dtype=str).fillna("")
        lineage = pd.read_csv(outputs["evidence_lineage"], sep="\t", dtype=str).fillna("")

        category_by_type = observations.groupby("observation_type")["observation_category"].apply(set).to_dict()
        self.assertEqual({"mapping_robustness"}, category_by_type["identity_mapping"])
        self.assertEqual({"biological_annotation"}, category_by_type["biological_annotation"])
        self.assertEqual({"derived_result"}, category_by_type["screening_evidence"])
        self.assertEqual({"derived_result"}, category_by_type["automated_rank"])

        by_id = observations.set_index("observation_id")
        rank_id = observations.loc[observations["observation_type"].eq("automated_rank"), "observation_id"].item()
        evidence_id = observations.loc[observations["observation_type"].eq("screening_evidence"), "observation_id"].item()
        trait_id = observations.loc[observations["observation_type"].eq("trait_function_match"), "observation_id"].item()
        annotation_a_id = observations.loc[
            observations["provider_name"].eq("annotation:model_a"), "observation_id"
        ].item()
        mapping_a_id = observations.loc[
            observations["provider_name"].eq("blastp:model_a"), "observation_id"
        ].item()
        parents = lineage.groupby("observation_id")["derived_from"].apply(set).to_dict()

        self.assertIn(evidence_id, parents[rank_id])
        self.assertEqual({trait_id}, parents[evidence_id])
        self.assertEqual({annotation_a_id}, parents[trait_id])
        self.assertEqual({mapping_a_id}, parents[annotation_a_id])
        self.assertEqual("mapping_robustness", by_id.loc[mapping_a_id, "observation_category"])


if __name__ == "__main__":
    unittest.main()
