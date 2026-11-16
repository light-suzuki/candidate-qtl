from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.supported_member_union import build_supported_member_union
from tests._issue664_helpers import manifest_for


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    if "route_c" in path.name:
        manifest_for(path, route="c")
    elif "claims" in path.name:
        manifest_for(path, route="b_claims")
    elif "bridge" in path.name:
        manifest_for(path, route="a")
    return path


def test_route_c_only_admits_verified_direct_support(tmp_path: Path) -> None:
    members = _write(
        tmp_path / "route_c_members.tsv",
        "candidate_gene_id\tsupport_status\trecord_id\tsource_id\tspecies_distance\tevidence_directness\tprojection_state\tprojection_basis\n"
        "G_C1\tverified_direct_support\tR1\tPMID:1\tS1_same_species\tD1_mapped_gene_with_functional_validation\tinside_current_qtl\texact_focal_id\n"
        "G_REVIEW\treference_phenotype_review\tR2\tPMID:2\tS5_model_or_other_plant\tD0_causal_perturbation_or_rescue\tinside_current_qtl\torthology_reference_gene\n"
        "G_OUT\tverified_direct_support\tR3\tPMID:3\tS1_same_species\tD1_mapped_gene_with_functional_validation\toutside_current_qtl\texact_focal_id\n",
    )

    result = build_supported_member_union(str(tmp_path / "out"), route_c_members_arg=str(members))
    frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
    assert frame["candidate_gene_id"].tolist() == ["G_C1"]
    assert frame.loc[0, "route_support"] == "c"
    assert frame.loc[0, "route_c_supported"] == "True"
    assert frame.loc[0, "route_c_record_ids"] == "R1"
    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    assert manifest["route_c_input_status"] == "supplied"
    assert manifest["route_c_rows"] == 3
    assert manifest["route_c_verified_support_rows"] == 1
    assert manifest["route_c_nonverified_rows"] == 2
    assert manifest["route_c_supported_candidate_count"] == 1


def test_route_c_provenance_combines_deterministically_with_a_and_b(tmp_path: Path) -> None:
    bridge = _write(
        tmp_path / "bridge.tsv",
        "gene_id\troute\tmapping_tier\tfamily_id\thypothesis_id\tambiguity_status\n"
        "G1\ta\ttier1\tFA\tHA\tunique\n",
    )
    claims = _write(
        tmp_path / "claims.tsv",
        "claim_id\tverification_status\tmember_interval_genes\tfamily_id\torthogroup_id\tsemantic_unit_id\tevidence_scope\n"
        "B1\tverified\tG1;G2\tFB\tOB\tUB\twithin_species\n",
    )
    members = _write(
        tmp_path / "route_c_members.tsv",
        "candidate_gene_id\tsupport_status\trecord_ids\tsource_ids\tspecies_distances\tprojection_states\n"
        "G1\tverified_direct_support\tC1\tPMID:3\tS1_same_species\tinside_current_qtl\n"
        "G3\tverified_direct_support\tC3\tPMID:4\tS1_same_species\tinside_current_qtl\n",
    )

    result = build_supported_member_union(
        str(tmp_path / "out"),
        route_a_bridge_arg=str(bridge),
        route_b_claims_arg=str(claims),
        route_c_members_arg=str(members),
    )
    frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
    assert frame.loc["G1", "route_support"] == "a+b+c"
    assert frame.loc["G2", "route_support"] == "b"
    assert frame.loc["G3", "route_support"] == "c"
    assert frame.loc["G1", "route_c_source_ids"] == "PMID:3"
    assert frame.loc["G2", "route_c_supported"] == "False"


def test_route_c_is_optional_for_existing_a_b_union(tmp_path: Path) -> None:
    bridge = _write(
        tmp_path / "bridge.tsv",
        "gene_id\troute\tmapping_tier\nG1\ta\ttier1\n",
    )
    result = build_supported_member_union(str(tmp_path / "out"), route_a_bridge_arg=str(bridge))
    frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    assert frame.loc[0, "route_support"] == "a"
    assert frame.loc[0, "route_c_supported"] == "False"
    assert manifest["route_c_input_status"] == "not_supplied"
    assert manifest["route_c_rows"] == 0


def test_route_c_no_match_ledger_row_with_empty_gene_does_not_block_a(tmp_path: Path) -> None:
    bridge = _write(
        tmp_path / "bridge.tsv",
        "gene_id\troute\tmapping_tier\nG1\ta\ttier1\n",
    )
    members = _write(
        tmp_path / "route_c_members.tsv",
        "candidate_gene_id\tsupport_status\tprojection_state\trecord_id\n"
        "\tledger_only\tno_interval_match\tC0\n",
    )
    result = build_supported_member_union(
        str(tmp_path / "out"), route_a_bridge_arg=str(bridge), route_c_members_arg=str(members)
    )
    frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
    assert frame["candidate_gene_id"].tolist() == ["G1"]


class SupportedMemberUnionRouteCUnittestTests(unittest.TestCase):
    """Ensure the repository's unittest-discover CI executes these contracts."""

    def _run(self, function) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            function(Path(tmp))

    def test_route_c_only(self) -> None:
        self._run(test_route_c_only_admits_verified_direct_support)

    def test_abc_union(self) -> None:
        self._run(test_route_c_provenance_combines_deterministically_with_a_and_b)

    def test_optional_c(self) -> None:
        self._run(test_route_c_is_optional_for_existing_a_b_union)

    def test_empty_no_match_gene(self) -> None:
        self._run(test_route_c_no_match_ledger_row_with_empty_gene_does_not_block_a)
