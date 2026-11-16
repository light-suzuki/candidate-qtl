from __future__ import annotations

import gzip
import json
import re
import shlex
import shutil
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path, run_tool


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _normalize_consequence(value: str) -> str:
    aliases = {"missense":"missense_variant","frameshift":"frameshift_variant","stop_gained":"stop_gained",
               "splice_donor":"splice_donor_variant","splice_acceptor":"splice_acceptor_variant",
               "5_prime_utr_variant":"5_prime_UTR_variant","3_prime_utr_variant":"3_prime_UTR_variant"}
    terms = [aliases.get(term.strip().casefold(), term.strip()) for term in re.split(r"[&,+]", value) if term.strip()]
    return "&".join(sorted(set(terms)))


def _parse_annotated_vcf(path: Path, provider: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []; csq_fields: list[str] = []
    with _open(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\">]+)", line)
                if match: csq_fields = match.group(1).strip().split("|")
                continue
            if line.startswith("#"): continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8: continue
            chrom, pos, identifier, ref, alt, _, _, info = fields[:8]
            info_map = {item.split("=", 1)[0]: item.split("=", 1)[1] if "=" in item else "" for item in info.split(";")}
            annotations = info_map.get("ANN" if provider == "snpeff" else "CSQ", "")
            for annotation in filter(None, annotations.split(",")):
                values = annotation.split("|")
                if provider == "snpeff":
                    values += [""] * (16 - len(values))
                    allele, consequence, impact, symbol, gene, feature_type, transcript = values[:7]
                    hgvsc, hgvsp = values[9], values[10]
                else:
                    mapping = dict(zip(csq_fields, values))
                    allele = mapping.get("Allele", ""); consequence = mapping.get("Consequence", ""); impact = mapping.get("IMPACT", "")
                    symbol = mapping.get("SYMBOL", ""); gene = mapping.get("Gene", ""); feature_type = mapping.get("Feature_type", ""); transcript = mapping.get("Feature", "")
                    hgvsc, hgvsp = mapping.get("HGVSc", ""), mapping.get("HGVSp", "")
                variant_id = identifier if identifier not in {"", "."} else f"{chrom}:{pos}:{ref}:{allele or alt}"
                rows.append({"variant_id":variant_id,"chrom":chrom,"pos":pos,"ref":ref,"alt":allele or alt,"provider":provider,
                             "gene_id":gene or symbol,"gene_symbol":symbol,"transcript_id":transcript,"feature_type":feature_type,
                             "consequence":_normalize_consequence(consequence),"tool_impact":impact,"hgvsc":hgvsc,"hgvsp":hgvsp,
                             "validation_status":"predicted_unvalidated","source_file":str(path)})
    return rows


def _annotation_version(path: Path, provider: str) -> str:
    patterns = ("##SnpEffVersion=", "##VEP=", "##VEP_version=")
    with _open(path) as handle:
        for line in handle:
            if not line.startswith("##"): break
            if any(line.startswith(pattern) for pattern in patterns):
                return line.strip().split("=", 1)[1].strip('"')
    return "not_recorded_in_vcf"


def variant_evidence(vcf_arg: str, fasta_arg: str, out_dir_arg: str, interval: str | None = None,
                     snpeff_vcf_arg: str | None = None, vep_vcf_arg: str | None = None,
                     execute_norm: bool = False) -> dict[str, object]:
    vcf, fasta, out = runtime_path(vcf_arg), runtime_path(fasta_arg), runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True); normalized = out / "normalized.vcf.gz"
    region = f" -r {shlex.quote(interval)}" if interval else ""
    norm_command = f"bcftools norm -f {shlex.quote(str(fasta))} -m -any --check-ref w{region} -Oz -o {shlex.quote(str(normalized))} {shlex.quote(str(vcf))}"
    versions: dict[str, str] = {}
    if execute_norm:
        if not shutil.which("bcftools"): raise ValueError("bcftools not found; run the generated command in a bcftools environment or supply a normalized VCF")
        completed = run_tool(norm_command, shell=True, text=True, capture_output=True)
        if completed.returncode: raise RuntimeError(f"bcftools norm failed: {completed.stderr.strip()}")
        run_tool(["bcftools", "index", "-f", str(normalized)], check=True)
        version = run_tool(["bcftools", "--version"], text=True, capture_output=True)
        versions["bcftools"] = version.stdout.splitlines()[0]
    source_for_fingerprint = normalized if normalized.exists() else vcf
    rows = []
    if snpeff_vcf_arg:
        snpeff_path = runtime_path(snpeff_vcf_arg); rows.extend(_parse_annotated_vcf(snpeff_path, "snpeff")); versions["snpeff"] = _annotation_version(snpeff_path, "snpeff")
    if vep_vcf_arg:
        vep_path = runtime_path(vep_vcf_arg); rows.extend(_parse_annotated_vcf(vep_path, "vep")); versions["vep"] = _annotation_version(vep_path, "vep")
    annotations = pd.DataFrame(rows, columns=["variant_id","chrom","pos","ref","alt","provider","gene_id","gene_symbol","transcript_id","feature_type","consequence","tool_impact","hgvsc","hgvsp","validation_status","source_file"])
    annotations.to_csv(out / "variant_annotations.tsv", sep="\t", index=False)
    concordance_rows = []
    if not annotations.empty:
        provider_transcripts = annotations.groupby(["variant_id","gene_id","provider"], dropna=False)["transcript_id"].apply(lambda values: set(values) - {""})
        for key, group in annotations.groupby(["variant_id","gene_id","transcript_id"], dropna=False, sort=True):
            by_provider = {provider:set("&".join(values).split("&")) - {""} for provider, values in group.groupby("provider")["consequence"]}
            snp, vep = by_provider.get("snpeff", set()), by_provider.get("vep", set())
            if snp and vep:
                status = "exact" if snp == vep else ("partial" if snp & vep else "discordant")
            else: status = "snpeff_only" if snp else "vep_only"
            snp_tx = provider_transcripts.get((key[0], key[1], "snpeff"), set()); vep_tx = provider_transcripts.get((key[0], key[1], "vep"), set())
            concordance_rows.append({"variant_id":key[0],"gene_id":key[1],"transcript_id":key[2],"snpeff_consequence":"&".join(sorted(snp)),
                                     "vep_consequence":"&".join(sorted(vep)),"concordance":status,"transcript_selection_difference":bool(snp_tx and vep_tx and snp_tx != vep_tx),
                                     "review_required":status in {"partial","discordant"},"interpretation":"annotation robustness only; not independent biological evidence"})
    concordance = pd.DataFrame(concordance_rows, columns=["variant_id","gene_id","transcript_id","snpeff_consequence","vep_consequence","concordance","transcript_selection_difference","review_required","interpretation"])
    concordance.to_csv(out / "variant_annotation_concordance.tsv", sep="\t", index=False)
    concordance[concordance["review_required"].eq(True)].to_csv(out / "review_required.tsv", sep="\t", index=False)
    evidence_rows = []
    for gene_id, group in annotations.groupby("gene_id", sort=True) if not annotations.empty else []:
        if not gene_id: continue
        evidence_rows.append({"gene_id":gene_id,"evidence_code":"predicted_variant_context","direction":"context","confidence":"low",
                              "value":";".join(sorted(set(group["consequence"]))),"source_type":"variant_annotation_concordance",
                              "source_id":str(out / "variant_annotation_concordance.tsv"),"description":"Predicted parental variant consequence; requires parental validation and is not causal evidence.",
                              "independence_group":f"variant_annotation:{source_for_fingerprint}"})
    pd.DataFrame(evidence_rows, columns=["gene_id","evidence_code","direction","confidence","value","source_type","source_id","description","independence_group"]).to_csv(out / "variant_evidence.tsv", sep="\t", index=False)
    commands = {"bcftools_norm":norm_command,"snpeff":"snpEff -v <database> normalized.vcf.gz > snpeff.vcf","vep":"vep --vcf --offline --input_file normalized.vcf.gz --output_file vep.vcf"}
    manifest = {"schema_version":"1.0","input_vcf":file_fingerprint(vcf, False),"reference_fasta":file_fingerprint(fasta, False),
                "normalized_vcf":file_fingerprint(normalized, False) if normalized.exists() else None,"normalization_executed":execute_norm,
                "commands":commands,"tool_versions":versions,"annotation_inputs":{"snpeff":snpeff_vcf_arg,"vep":vep_vcf_arg},
                "annotation_count":len(annotations),"discordance_count":int(concordance["review_required"].sum()) if not concordance.empty else 0,
                "independence_rule":"SnpEff and VEP on the same VCF are one annotation chain, not independent evidence.",
                "wet_gate":"All emitted evidence is predicted_unvalidated/context until parental genotype confirmation; Sanger and cosegregation stay separate."}
    manifest_path = out / "variant_evidence_manifest.json"; manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (out / "run_variant_annotation.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands.values()) + "\n", encoding="utf-8")
    return {"manifest":str(manifest_path),"concordance":str(out / "variant_annotation_concordance.tsv"),"evidence":str(out / "variant_evidence.tsv"),"annotation_count":len(annotations)}
