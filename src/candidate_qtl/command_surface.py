"""Machine-readable ownership registry for the two CLI entry points (#669).

The registry is deliberately data-only: importing it must not import any feature
implementation.  This makes it safe for packaging checks and for the production
parser to decide which command modules are allowed before loading them.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

SURFACE_CLASSES = frozenset({"production_core", "production_optional", "screening", "research", "benchmark", "legacy"})

@dataclass(frozen=True)
class CommandSurface:
    command: str
    feature_module: str
    surface: str
    installed_entrypoint: str
    canonical_consumer: str
    reachability_profile: str
    deprecation_or_move_state: str = "stable"
    replacement: str = ""

    def as_dict(self) -> dict[str, str]:
        return asdict(self)


# A/C-first production (2026-11-08): Route A and Direct C form the ordinary
# workflow; interval-wide Route B research is an explicit option.
DIRECT_C_CORE_COMMANDS = frozenset({"route-c-direct-freeze", "route-c-provider-bait", "route-c-direct-project", "route-c-direct-members"})
ROUTE_B_OPTIONAL_COMMANDS = frozenset({
    "route-b-prepare", "route-b-semantic-bundle", "route-b-semantic-candidates",
    "route-b-semantic-import", "route-b-semantic-source-verify", "route-b-semantic-transfer-verify",
})


def _records() -> tuple[CommandSurface, ...]:
    # This is kept explicit so a newly added command cannot silently inherit a
    # guessed ownership.  The command modules remain the source of argparse
    # details; this table is the source of surface ownership.
    groups = {
        "production_core": "assembly assays blast_panel budget comparative diagnostics evidence extraction family homolog_linking mapping orchestration orthology package ranking scenarios tnt1 variants wet".split(),
        "production_optional": "llm budget_llm".split(),
        "screening": "discovery".split(),
        "research": "research retrieval".split(),
        "benchmark": "benchmark paired_benchmark".split(),
        "legacy": [],
    }
    records: list[CommandSurface] = []
    # Registry names are generated from the stable command declaration in each
    # module, without importing those modules.  The declarations are mirrored
    # here as command -> module lists and validated against loaded handlers.
    commands = {
      "assembly": "assembly-bridge multi-reference-locus-union", "assays": "assay-plan benchmark-reference collect-sanger-markers exon-coverage-check exon-sequencing-plan primer-plan primerblast-import",
      "blast_panel": "blast-panel-freeze blast-panel-search blast-panel-verify", "budget": "adaptive-refinement-apply adaptive-synthesis-import adaptive-synthesis-request budget-run candidate-budget-gate llm-usage supported-member-union wet-candidate-handoff",
      "budget_llm": "budget-llm-review", "comparative": "comparative-locus-benchmark comparative-locus-reconstruct external-rank propagate-terms reference-panel route-b-semantic-bundle route-b-semantic-candidates route-b-semantic-import route-b-semantic-source-verify route-b-semantic-transfer-verify semantic-claims-freeze",
      "diagnostics": "doctor provenance-ledger", "discovery": "auto-triage discover discover-promote homolog-panel-auto homolog-panel-benchmark homolog-panel-match", "evidence": "candidate-elimination direction-evidence evidence-build evidence-term-promote peak-context qtl-source-init trait-ontology-map trait-terms",
      "extraction": "annotate-plan annotation-keywords annotation-rescue candidate-flags dark-candidates extract gff-preflight liftover-markers", "family": "family-candidate-bridge", "homolog_linking": "homolog-linking-counterfactual homolog-linking-export homolog-linking-import homolog-linking-probe-export homolog-linking-score",
      "llm": "bait-panel-freeze hypothesis-verify llm-evidence-import llm-pathway-import llm-pathway-prompt llm-prompt-export llm-review-benchmark llm-run llm-shortlist-import llm-shortlist-prompt phenotype-pathway-benchmark research-packet",
      "mapping": "cosegregation qtl2-plan scaffold-rescue", "orchestration": "case-run guided-profile orchestrate", "orthology": "dual-route-compare family-index family-to-copy focal-paralog-import homology-concordance integrated-narrowing interval-atlas interval-reverse-lookup omamer-import orthofinder-import orthology-import orthology-plan orthology-tier1-plan pathway-interval-falsify phenotype-first-status reference-candidate-table reverse-literature route-a-interval-match route-b-prepare synteny-import",
      "package": "build-executor-package", "paired_benchmark": "paired-benchmark", "ranking": "candidate-pipeline dry-prioritize evidence-rank evidence-tier-benchmark evidence-tiers mapping-gate next-experiment rank ranking-qa ranking-stability-benchmark screening-shortlist weight-policy-benchmark",
      "research": "claim-verification-benchmark investigation-attribution literature mechanism-coverage-audit mechanism-import phenotype-investigate phenotype-localization-freeze phenotype-localization-packet reference-gene-cache reference-gene-resolve research source-verify",
      "retrieval": "corpus-mask phenotype-variant-benchmark reference-corpus-build reference-corpus-compare reference-corpus-fetch retrieval-hybrid-merge retrieval-import retrieval-index-build retrieval-interval-match retrieval-merge retrieval-run route-c-benchmark route-c-direct-benchmark route-c-direct-freeze route-c-direct-members route-c-direct-project route-c-direct-search route-c-provider-bait route-c-provider-comparison triple-route-compare",
      "scenarios": "scenario-create scenario-report", "tnt1": "tnt1-evidence-normalize", "variants": "allele-compare multigenome-private promoter-review sv-plan variant-evidence", "wet": "hypothesis-ledger wet-state",
      "benchmark": "audit-predict audit-score benchmark-fixture-validate benchmark-knowledge-mask benchmark-report-compare benchmark-stratify blind-benchmark challenge-benchmark-aggregate command-log-validate cross-crop-benchmark discovery-benchmark discovery-chain-benchmark knowledge-coverage masking-benchmark mechanism-coverage-benchmark next-experiment-benchmark query-eligibility reference-panel-benchmark robustness-benchmark supervised-benchmark-audit supervised-benchmark-predict-prepare supervised-benchmark-route-family-metrics supervised-benchmark-score supervised-benchmark-source-log-freeze supervised-benchmark-truth-validate supervised-benchmark-validate",
    }
    owner = {module: cls for cls, modules in groups.items() for module in modules}
    for module, names in commands.items():
        cls = owner.get(module, "legacy")
        for command in names.split():
            # A few historical feature modules contain both user and benchmark
            # commands.  Ownership is per command, never inferred from module.
            # A/C-first production: Direct C freeze/projection/members are core;
            # the interval-wide Route B research commands are optional.
            if command in DIRECT_C_CORE_COMMANDS:
                command_cls = "production_core"
            elif command in ROUTE_B_OPTIONAL_COMMANDS:
                command_cls = "production_optional"
            elif command in {"phenotype-investigate", "phenotype-localization-packet", "phenotype-localization-freeze", "mechanism-import", "mechanism-coverage-audit", "investigation-attribution", "reference-gene-cache", "reference-gene-resolve", "source-verify", "route-c-direct-search", "triple-route-compare"}:
                command_cls = "production_optional"
            elif "benchmark" in command and cls in {"production_core", "production_optional", "screening"}:
                command_cls = "benchmark"
            elif module == "retrieval":
                if command.startswith("retrieval-"):
                    command_cls = "legacy"
                elif "benchmark" in command:
                    command_cls = "benchmark"
                else:
                    command_cls = "research"
            else:
                command_cls = cls
            entry = "candidate-qtl" if command_cls in {"production_core", "production_optional", "screening"} else "candidate-qtl-bench"
            feature_module = module
            if command in {"phenotype-investigate", "phenotype-localization-packet", "phenotype-localization-freeze", "mechanism-import", "mechanism-coverage-audit", "investigation-attribution", "reference-gene-cache", "reference-gene-resolve", "source-verify"}:
                feature_module = "production_research"
            if command in {"route-c-direct-search", "route-c-direct-freeze", "route-c-provider-bait", "route-c-direct-project", "route-c-direct-members", "triple-route-compare"}:
                feature_module = "production_direct_c"
            records.append(CommandSurface(command, feature_module, command_cls, entry, "candidate_qtl", "production" if entry == "candidate-qtl" else "benchmark"))
    return tuple(records)

COMMAND_SURFACES = _records()
COMMAND_SURFACE_BY_NAME = {record.command: record for record in COMMAND_SURFACES}

if len(COMMAND_SURFACE_BY_NAME) != len(COMMAND_SURFACES) or any(r.surface not in SURFACE_CLASSES for r in COMMAND_SURFACES):
    raise RuntimeError("invalid command surface registry")

def records_for(surface: str) -> tuple[CommandSurface, ...]:
    if surface not in {"core", "all", "benchmark"}:
        raise ValueError(f"unknown CLI surface: {surface!r}")
    if surface == "all":
        return COMMAND_SURFACES
    if surface == "benchmark":
        # The bench executable is a superset: production APIs remain available
        # for reproducing a benchmark from one clean installation.
        return COMMAND_SURFACES
    return tuple(r for r in COMMAND_SURFACES if r.installed_entrypoint == "candidate-qtl")
