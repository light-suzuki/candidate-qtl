"""Source attribution is not organism identity; generic header/consumer cases."""
from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from candidate_qtl.route_b_blast_evidence import run_local_blastp
from candidate_qtl.route_b_identity import header_metadata, merge_header_metadata
from candidate_qtl.route_b_prepare import route_b_prepare
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle


TITLE = (
    "pep chromosome:ASSEMBLY:1:10:30:1 gene:GENE1 transcript:GENE1.1 "
    "gene_biotype:protein_coding description:Example enzyme "
    "[Source:NCBI gene (formerly Entrezgene);Acc:12345]"
)


@pytest.mark.parametrize("source", ["Source:", "source:", " Source :"])
@pytest.mark.parametrize("species", ["", "N/A", "Example species"])
def test_source_bracket_keeps_raw_attribution_and_database_taxonomy(source, species):
    title = TITLE.replace("Source:", source)
    parsed = header_metadata(title, "123", species)
    assert parsed["hit_species"] == "unknown"
    assert parsed["hit_description"] == title
    assert parsed["database_taxon_ids"] == "123"
    merged = merge_header_metadata([(title, "123", species)])
    assert merged["hit_species"] == (species if species == "Example species" else "unknown")
    assert merged["hit_taxon_id"] == "123"
    assert merged["hit_metadata_state"] == "partial"
    assert json.loads(merged["hit_metadata_records_json"]) == [parsed]


def test_source_bracket_does_not_guess_species_from_accession_or_taxon():
    result = merge_header_metadata([(TITLE, "N/A", "N/A")])
    assert result["hit_species"] == result["hit_taxon_id"] == "unknown"
    assert result["hit_metadata_state"] == "partial"


def test_real_taxonomy_conflict_remains_order_independent():
    records = [(TITLE, "123", "Example species"), (TITLE, "456", "Other species")]
    result = merge_header_metadata(records)
    assert result == merge_header_metadata([records[1], records[0], records[1]])
    assert result["hit_species"] == result["hit_taxon_id"] == "unknown"
    assert result["hit_metadata_state"] == "conflicting"
    assert len(json.loads(result["hit_metadata_records_json"])) == 2


@pytest.mark.parametrize("header,taxon,species", [
    ("XP_000001.1 Example enzyme [Example species]", "", "Other species"),
    ("sp|P1|NAME_X Example enzyme OS=Example species OX=123 GN=GENE1", "456", "Example species"),
])
def test_actual_ncbi_and_uniprot_conflicts_are_not_hidden(header, taxon, species):
    assert merge_header_metadata([(header, taxon, species)])["hit_metadata_state"] == "conflicting"


def test_source_metadata_reaches_b0_and_semantic_without_changing_hit_inventory(tmp_path):
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("gene_id\tchrom\tstart\tend\ng1\tchr1\t1\t10\ng2\tchr1\t20\t30\n")
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">g1\nAAAA\n")
    family = tmp_path / "family.tsv"
    family.write_text("gene_id\tfamily_id\n")

    def fake_run(command, **kwargs):
        rows = ["\t".join(map(str, ["g1", accession, 100, 4, 4, 4, 1, 4, 1, 4,
                1e-20, 100, "AAAA", "AAAA", TITLE, "123", species]))
                for accession, species in (("REF1", "N/A"), ("REF2", "Example species"))]
        # Duplicate HSPs must not inflate hits or erase metadata.
        Path(command[command.index("-out") + 1]).write_text("\n".join(rows + rows[::-1]) + "\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    def blast(proteins, database, out, threads, identity, coverage, evalue):
        return run_local_blastp(proteins, database, out, threads, identity, coverage, evalue,
                                run_tool=fake_run, find_binary=lambda _: "blastp", normalise_subject=str)

    with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=blast):
        prepared = route_b_prepare(str(candidates), str(tmp_path / "b0"),
            proteins_arg=str(proteins), protein_db_arg="fixture.db", family_index_arg=str(family), auto_local_db=False)
    with Path(prepared["identity_hits"]).open() as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["hit_accession"] for row in hits} == {"REF1", "REF2"}
    assert {row["hit_species"] for row in hits} == {"unknown", "Example species"}
    assert all(row["hit_description"] == TITLE and row["hit_taxon_id"] == "123" for row in hits)
    assert all(float(row["query_coverage"]) == 1 and float(row["percent_identity"]) == 100 for row in hits)
    assert all(row["rank_contribution"] == "none" for row in hits)
    with Path(prepared["identity_inventory"]).open() as handle:
        inventory = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert set(inventory) == {"g1", "g2"}
    assert inventory["g1"]["identity_hit_count"] == "2"
    assert inventory["g2"]["identity_state"] == "not_evaluated_no_protein_query"
    # Identity-only units are deliberately not scheduled as semantic calls.
    # Supply frozen functional material to exercise the actual requested payload.
    lookup = tmp_path / "lookup.tsv"
    lookup.write_text("gene_id\treference_function\ng1\tExample enzyme in organ development\n")
    result = build_route_b_semantic_bundle(str(prepared["interval_gene_atlas"]), str(tmp_path / "semantic"),
        phenotype="altered organ development", identity_hits_arg=str(prepared["identity_hits"]), lookup_arg=str(lookup))
    bundle = json.loads(Path(result["bundle"]).read_text())
    assert result["interval_gene_count"] == 2
    bundle_hits = [hit for unit in bundle["units"] for gene in unit["genes"] for hit in gene["sequence_similarity_hits"]]
    assert {hit["hit_accession"] for hit in bundle_hits} == {"REF1", "REF2"}
    assert all(hit["hit_description"] == TITLE for hit in bundle_hits)
    prompt_hits = []
    for line in Path(result["prompts"]).read_text().splitlines():
        payload = json.loads(json.loads(line)["user"].split("\n\n", 1)[1])
        prompt_hits.extend(hit for gene in payload["genes"] for hit in gene["sequence_similarity_hits"])
    assert {hit["hit_species"] for hit in prompt_hits} == {"unknown", "Example species"}
    assert all(hit["interpretation"] == "sequence_similarity_handle_only" for hit in prompt_hits)
