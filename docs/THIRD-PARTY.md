# Third-party licenses

This repository's own code is MIT-licensed and its documentation/templates are
CC0-1.0 (see [LICENSE](../LICENSE) and LICENSE-CC0-1.0.txt).
The components below are **not** redistributed in this repository; they are
declared as dependencies or invoked as external tools at runtime. Each remains
under its own license. Verify the exact license and version of each upstream
project before you redistribute any derived work or bundle.

## Python runtime dependencies (declared in `pyproject.toml`)

| Package | Version pin | License (indicative) |
|---|---|---|
| Biopython | >=1.83 | Biopython License Agreement (BSD-3-Clause-style) |
| pandas | >=2.0 | BSD-3-Clause |
| pysam | >=0.22 | MIT |
| PyYAML | >=6.0 | MIT |

## External tools invoked through WSL / other CLIs (not bundled)

These are called as subprocesses (directly or via `primerblast-oss` / `QTLift`).
Only tools actually used in a given run apply. Licenses are indicative and
should be confirmed against the installed version.

| Tool | Role in this workflow | License (indicative — confirm upstream) |
|---|---|---|
| NCBI BLAST+ (`blastn`, `blastp`, `makeblastdb`) | homology / specificity | Public Domain (US Gov work) |
| Primer3 (`primer3_core`) | primer design (via primerblast-oss) | GPL-2.0 |
| minimap2 | assembly-to-assembly evidence (QTLift) | MIT |
| HMMER | Pfam domain search | BSD-3-Clause |
| DIAMOND | fast protein homology | GPL-3.0 |
| MMseqs2 | fast protein homology (optional) | GPL-3.0 |
| OrthoFinder | orthogroup inference (optional) | GPL-3.0 |
| eggNOG-mapper | GO / orthology annotation (optional) | GPL / AGPL — confirm upstream |
| OMAmer | HOG assignment (optional) | confirm upstream |
| MAFFT | multiple alignment (Tier 3) | BSD-style (MAFFT license) |
| FastTree / IQ-TREE | phylogeny (Tier 3) | GPL-2.0 |
| trimAl / ClipKIT | alignment trimming | GPL-3.0 / MIT |
| InterProScan | domain / GO annotation (optional) | Apache-2.0 |
| bcftools / samtools / htslib | VCF/BAM handling (optional) | MIT/Expat |
| SnpEff / Ensembl VEP | variant annotation (optional) | confirm upstream / Apache-2.0 |
| R/qtl2 | QTL detection/remap (optional) | GPL-3.0 |
| Snakemake | multi-QTL orchestration (optional) | MIT |

## Sibling repositories by the same author

| Repository | Relationship | License |
|---|---|---|
| [primerblast-oss](https://github.com/light-suzuki/primerblast-oss) | marker/primer design backend | see repository |
| QTLift | interval coordinate lift (Step 0) | see repository |

> Note: GPL-licensed external tools (e.g. Primer3, DIAMOND, OrthoFinder) are
> invoked as separate executables and are not linked into this MIT-licensed
> code, so their copyleft terms do not extend to this repository. If you
> redistribute a bundle that includes those binaries, comply with each tool's
> license separately.
