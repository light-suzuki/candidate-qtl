"""CLI extraction keeps no-CDS models distinct from genuinely faulty phases."""
from __future__ import annotations

import contextlib
import csv
import io
import json

import pytest
from Bio import SeqIO

from candidate_qtl.cli import main


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract(tmp_path, gff, sequence="ATGAAATAA"):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(f">chr1\n{sequence}\n")
    annotation = tmp_path / "genes.gff3"
    annotation.write_text("##gff-version 3\n" + gff)
    out = tmp_path / "out"
    with contextlib.redirect_stdout(io.StringIO()):
        assert main(["extract", "--fasta", str(fasta), "--gff3", str(annotation),
                     "--interval", f"chr1:1-{len(sequence)}", "--out-dir", str(out)]) == 0
    return out


@pytest.mark.parametrize("gene_type,biotype,transcript_type", [
    ("ncRNA_gene", "lncRNA", "lnc_RNA"),
    ("ncRNA_gene", "ncRNA", "ncRNA"),
    ("pseudogene", "pseudogene", "pseudogenic_transcript"),
    ("gene", "protein_coding", "mRNA"),
])
@pytest.mark.parametrize("with_transcript", [True, False])
def test_no_cds_preserves_gene_model_and_never_invents_protein(tmp_path, gene_type, biotype, transcript_type, with_transcript):
    gff = f"chr1\ttest\t{gene_type}\t1\t9\t.\t+\t.\tID=gene:X;biotype={biotype}\n"
    if with_transcript:
        gff += (f"chr1\ttest\t{transcript_type}\t1\t9\t.\t+\t.\tID=transcript:X.1;Parent=gene:X\n"
                "chr1\ttest\texon\t1\t9\t.\t+\t.\tParent=transcript:X.1\n")
    out = extract(tmp_path, gff)
    rows = read_tsv(out / "candidates.tsv")
    assert len(rows) == 1
    row = rows[0]
    assert row["gene_id"] == "gene:X"
    assert row["gene_biotype"] == biotype
    assert row["cds_length"] == row["protein_length"] == "0"
    assert row["cds_phase_status"] == row["translation_status"] == "no_cds"
    assert row["annotation_quality"] == ("no_cds" if biotype == "protein_coding" else "noncoding_or_pseudogene")
    assert "invalid" not in row["representative_reason"]
    assert list(SeqIO.parse(out / "candidate_proteins.fa", "fasta")) == []
    assert list(SeqIO.parse(out / "candidate_cds.fa", "fasta")) == []
    assert [r.id for r in SeqIO.parse(out / "candidate_genomic.fa", "fasta")] == ["gene:X"]
    transcripts = read_tsv(out / "transcripts.tsv")
    mapping = read_tsv(out / "id_map.tsv")
    if with_transcript:
        assert row["representative_reason"] == "best available transcript; no annotated CDS"
        assert transcripts[0]["phase_status"] == "no_cds"
        assert transcripts[0]["selection_reason"] == row["representative_reason"]
        assert transcripts[0]["phase_warnings"] == ""
        assert mapping[0]["gene_id"] == "gene:X"
        assert mapping[0]["transcript_id"] == "transcript:X.1"
        assert mapping[0]["protein_id"] == mapping[0]["sequence_record_id"] == ""
    else:
        assert transcripts == mapping == []
    manifest = json.loads((out / "extract_manifest.json").read_text())
    assert manifest["candidate_count"] == 1
    assert manifest["with_protein_count"] == 0
    assert manifest["cds_phase_status_counts"] == ({"no_cds": 1} if with_transcript else {})


@pytest.mark.parametrize("sequence,segments,phase_status,quality,protein", [
    ("ATGAAATAA", [(1, 9, "0")], "valid", "complete_orf", "MK"),
    ("ATGAAATAA", [(1, 9, ".")], "missing", "invalid_cds_phase", ""),
    ("ATGAAATAA", [(1, 4, "0"), (5, 9, "0")], "inconsistent", "invalid_cds_phase", ""),
    ("AATGAAATAA", [(1, 10, "1")], "valid_partial_5prime", "partial_or_unverified_orf", "MK"),
    ("ATGAAA", [(1, 6, "0")], "valid", "partial_or_unverified_orf", "MK"),
])
def test_existing_cds_phase_qc_and_partial_translation_remain_intact(tmp_path, sequence, segments, phase_status, quality, protein):
    n = len(sequence)
    gff = (f"chr1\ttest\tgene\t1\t{n}\t.\t+\t.\tID=gene:Y;biotype=protein_coding\n"
           f"chr1\ttest\tmRNA\t1\t{n}\t.\t+\t.\tID=transcript:Y.1;Parent=gene:Y\n")
    gff += "".join(f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t{phase}\tParent=transcript:Y.1\n" for start, end, phase in segments)
    out = extract(tmp_path, gff, sequence)
    row = read_tsv(out / "candidates.tsv")[0]
    assert row["cds_phase_status"] == phase_status
    assert row["annotation_quality"] == quality
    proteins = list(SeqIO.parse(out / "candidate_proteins.fa", "fasta"))
    assert (str(proteins[0].seq) if proteins else "") == protein
    if quality == "invalid_cds_phase":
        assert row["translation_status"] == "not_translated_invalid_phase"
        assert row["representative_reason"] == "best available transcript; CDS phase invalid or missing"
