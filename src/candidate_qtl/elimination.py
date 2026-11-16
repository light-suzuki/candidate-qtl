"""Candidate elimination / hypothesis survival layer (issue #171).

Ranking tells you what to read first. It does not tell you what has been ruled
out, and a score subtraction cannot distinguish "the gene is out" from "the coding
hypothesis is out" from "we could not measure it". This layer manages survival
state explicitly so a candidate set can be compressed **without quietly discarding
the real answer**.

The central structural rule is that elimination strength is bound to what it is
allowed to do:

* ``hard``     - positional / genetic evidence. The only strength permitted to move
                 a gene into an ``excluded_*`` state.
* ``strong``   - mechanism-specific negative evidence at the required quality and
                 coverage. May reject or weaken **one mechanism**, never the gene.
* ``moderate`` - discriminating counterevidence. May only weaken.
* ``none``     - assay failure, partial coverage, unmeasured expression, missing
                 annotation, provider failure. Never eliminates anything; it
                 produces ``unresolved`` / ``technical_failure`` /
                 ``coverage_incomplete``.

That mapping is enforced in code, not by convention: a registry row claiming a
``moderate`` rule transitions to an excluded state is rejected when the registry
loads. Absence of evidence never becomes evidence of absence, and every state
change is appended to a history with its rule id and evidence id.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"

GENE_STATES = ("survives", "weakened", "excluded_by_genetics", "excluded_by_recombinant",
               "excluded_by_coordinate", "unresolved", "not_tested")
EXCLUDED_STATES = frozenset({"excluded_by_genetics", "excluded_by_recombinant", "excluded_by_coordinate"})
MECHANISM_STATES = ("coding_survives", "coding_weakened", "coding_rejected", "splice_survives",
                    "regulatory_survives", "sv_pav_survives", "sv_pav_rejected",
                    "expression_hypothesis_weakened", "mechanism_unresolved")
NON_EVIDENCE_STATES = ("unresolved", "technical_failure", "coverage_incomplete")
STRENGTHS = ("hard", "strong", "moderate", "none")
MECHANISM_SCOPES = ("coding", "splice", "regulatory", "sv_pav", "expression")

REGISTRY_COLUMNS = {"elimination_rule_id", "scope", "evidence_code", "minimum_quality",
                    "minimum_coverage", "strength", "default_transition", "required_review",
                    "rationale", "registry_version"}
_QUALITY_ORDER = {"any": 0, "low": 1, "medium": 2, "high": 3}

SURVIVAL_COLUMNS = [
    "candidate_gene_id", "pre_elimination_group", "current_gene_state", "surviving_hypothesis_count",
    "rejected_hypothesis_count", "hard_exclusion_count", "weakened_reason",
    "strongest_exclusion_evidence_id", "remaining_uncertainty", "recommended_next_discriminator",
    "rank_contribution", "review_status",
]
HYPOTHESIS_SURVIVAL_COLUMNS = ["candidate_gene_id", "hypothesis_id", "scope", "mechanism_state",
                               "changed_by_rule", "evidence_id", "review_status", "notes"]
EVENT_COLUMNS = ["event_id", "candidate_gene_id", "hypothesis_id", "rule_id", "evidence_id",
                 "old_state", "new_state", "scope", "strength", "quality_gate_passed",
                 "coverage_gate_passed", "review_status", "reason"]
REMAINING_COLUMNS = ["candidate_gene_id", "current_gene_state", "surviving_mechanisms",
                     "first_wet_priority", "rank_contribution", "notes"]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]


def default_rules_path() -> Path:
    """Elimination rule registry used when none is supplied.

    Authoritative copy at the repo-root ``templates/elimination_rules.tsv``, with a
    byte-identical vendored copy in the package so a clean wheel install can run
    the layer without the source tree, mirroring the weight registry.
    """
    here = Path(__file__).resolve()
    source_tree = here.parents[2] / "templates" / "elimination_rules.tsv"
    if source_tree.exists():
        return source_tree
    return here.parent / "resources" / "elimination_rules.tsv"


def load_rules(path_arg: str | None = None) -> tuple[pd.DataFrame, dict[str, object]]:
    path = runtime_path(path_arg) if path_arg else default_rules_path()
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(REGISTRY_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"elimination rule registry is missing columns: {', '.join(missing)}")
    if frame["elimination_rule_id"].duplicated().any():
        raise ValueError("elimination rule registry requires unique elimination_rule_id")
    bad_strength = sorted(set(frame["strength"]) - set(STRENGTHS))
    if bad_strength:
        raise ValueError(f"unknown strength(s): {', '.join(bad_strength)}; allowed {', '.join(STRENGTHS)}")
    for row in frame.to_dict("records"):
        transition = str(row["default_transition"])
        strength = str(row["strength"])
        # Structural guard: only positional/genetic evidence may exclude a gene, and
        # only a mechanism-scoped rule may move a mechanism state.
        if transition in EXCLUDED_STATES and strength != "hard":
            raise ValueError(f"{row['elimination_rule_id']}: only a hard rule may transition to {transition}")
        if transition in EXCLUDED_STATES and str(row["scope"]) != "gene":
            raise ValueError(f"{row['elimination_rule_id']}: gene exclusion requires scope=gene")
        if transition in MECHANISM_STATES and str(row["scope"]) not in MECHANISM_SCOPES:
            raise ValueError(f"{row['elimination_rule_id']}: mechanism transition requires a mechanism scope")
        if strength == "none" and transition not in NON_EVIDENCE_STATES:
            raise ValueError(f"{row['elimination_rule_id']}: a non-evidence rule may only transition to "
                             f"{', '.join(NON_EVIDENCE_STATES)}")
    meta = {"path": str(path), "rule_count": int(len(frame)),
            "registry_version": sorted({str(value) for value in frame["registry_version"] if str(value)}),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return frame, meta


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows)[columns] if rows else pd.DataFrame(columns=columns)


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _gate(observed: str, minimum: str) -> bool:
    return _QUALITY_ORDER.get(observed.strip().casefold(), 0) >= _QUALITY_ORDER.get(minimum.strip().casefold(), 0)


def _coverage_gate(observed: str, minimum: str) -> bool:
    try:
        return float(observed or 0) >= float(minimum or 0)
    except ValueError:
        return False


def candidate_elimination(
    candidates_arg: str,
    out_dir_arg: str,
    evidence_arg: str | None = None,
    hypotheses_arg: str | None = None,
    wet_state_arg: str | None = None,
    rules_arg: str | None = None,
    previous_history_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Apply the elimination registry and emit the surviving candidate set."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidates = _read(candidates_arg, {"gene_id"}, "candidate")
    if candidates.empty:
        raise ValueError("candidate table is empty")
    if candidates["gene_id"].duplicated().any():
        raise ValueError("candidate table requires unique gene_id")
    evidence = _read(evidence_arg, {"gene_id", "evidence_code"}, "elimination evidence")
    hypotheses = _read(hypotheses_arg, {"hypothesis_id"}, "locus hypothesis")
    wet_state = _read(wet_state_arg, {"gene_id"}, "wet state")
    rules, rules_meta = load_rules(rules_arg)
    rule_by_code: dict[str, list[dict[str, object]]] = {}
    for rule in rules.to_dict("records"):
        rule_by_code.setdefault(str(rule["evidence_code"]), []).append(rule)
    history = _read(previous_history_arg, {"event_id"}, "previous history")

    gene_state: dict[str, str] = {str(row["gene_id"]): "not_tested" for row in candidates.to_dict("records")}
    pre_group = {str(row["gene_id"]): str(row.get("dry_group", "") or row.get("disambiguation_status", ""))
                 for row in candidates.to_dict("records")}
    mechanism_state: dict[tuple[str, str], dict[str, object]] = {}
    for row in hypotheses.to_dict("records"):
        for gene_id in (value.strip() for value in str(row.get("member_gene_ids", "") or row.get("gene_id", "")).replace(",", ";").split(";") if value.strip()):
            scope = str(row.get("hypothesis_type", "") or row.get("scope", "")) or "coding"
            mechanism_state[(gene_id, str(row["hypothesis_id"]))] = {
                "scope": scope, "state": f"{scope}_survives" if f"{scope}_survives" in MECHANISM_STATES else "mechanism_unresolved",
                "rule": "", "evidence_id": "", "review": "", "notes": "",
            }
    wet_coverage = {str(row["gene_id"]): row for row in wet_state.to_dict("records")} if not wet_state.empty else {}

    events: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []
    weakened_reason: dict[str, list[str]] = {}
    hard_evidence: dict[str, str] = {}
    unresolved_notes: dict[str, list[str]] = {}

    for record in evidence.sort_values(["gene_id", "evidence_code"], kind="stable").to_dict("records"):
        gene_id = str(record["gene_id"])
        if gene_id not in gene_state:
            continue
        code = str(record["evidence_code"])
        hypothesis_id = str(record.get("hypothesis_id", ""))
        evidence_id = str(record.get("evidence_id", ""))
        quality = str(record.get("quality", "") or wet_coverage.get(gene_id, {}).get("quality", ""))
        coverage = str(record.get("coverage", "") or wet_coverage.get(gene_id, {}).get("coverage_fraction", ""))
        matched = rule_by_code.get(code, [])
        if not matched:
            review_rows.append({"reason": "unknown_evidence_code", "target": f"{gene_id}:{code}",
                                "required_next_action": "registryに存在しないevidence_code。ruleを追加するか入力を修正する",
                                "allowed_choices": "ruleを追加 | 入力を修正 | 無視して未評価のまま"})
            continue
        # Prefer the strictest rule whose gates the evidence actually passes; a rule
        # whose gates fail degrades to its non-evidence counterpart rather than firing.
        applicable = [rule for rule in matched
                      if _gate(quality, str(rule["minimum_quality"])) and _coverage_gate(coverage, str(rule["minimum_coverage"]))]
        rule = min(applicable, key=lambda item: STRENGTHS.index(str(item["strength"]))) if applicable else \
            min(matched, key=lambda item: -STRENGTHS.index(str(item["strength"])))
        quality_ok = _gate(quality, str(rule["minimum_quality"]))
        coverage_ok = _coverage_gate(coverage, str(rule["minimum_coverage"]))
        strength = str(rule["strength"])
        transition = str(rule["default_transition"])
        if not (quality_ok and coverage_ok):
            # Gate failure can never eliminate; it records the unexamined scope.
            strength, transition = "none", "coverage_incomplete" if not coverage_ok else "unresolved"

        scope = str(rule["scope"])
        if strength == "hard" and transition in EXCLUDED_STATES:
            old, new = gene_state[gene_id], transition
            gene_state[gene_id] = new
            hard_evidence[gene_id] = evidence_id or str(rule["elimination_rule_id"])
        elif strength == "strong" and transition in MECHANISM_STATES:
            key = (gene_id, hypothesis_id or f"{scope}_default")
            entry = mechanism_state.setdefault(key, {"scope": scope, "state": f"{scope}_survives", "rule": "",
                                                     "evidence_id": "", "review": "", "notes": ""})
            old, new = str(entry["state"]), transition
            entry.update({"state": new, "rule": str(rule["elimination_rule_id"]), "evidence_id": evidence_id})
            if gene_state[gene_id] == "not_tested":
                gene_state[gene_id] = "survives"
        elif strength == "moderate":
            old = gene_state[gene_id]
            if transition in MECHANISM_STATES:
                key = (gene_id, hypothesis_id or f"{scope}_default")
                entry = mechanism_state.setdefault(key, {"scope": scope, "state": f"{scope}_survives", "rule": "",
                                                         "evidence_id": "", "review": "", "notes": ""})
                old, new = str(entry["state"]), transition
                entry.update({"state": new, "rule": str(rule["elimination_rule_id"]), "evidence_id": evidence_id})
            else:
                new = "weakened" if gene_state[gene_id] not in EXCLUDED_STATES else gene_state[gene_id]
                if gene_state[gene_id] not in EXCLUDED_STATES:
                    gene_state[gene_id] = new
            weakened_reason.setdefault(gene_id, []).append(str(rule["elimination_rule_id"]))
        else:
            old = gene_state[gene_id]
            new = transition if transition in NON_EVIDENCE_STATES else "unresolved"
            if gene_state[gene_id] == "not_tested":
                gene_state[gene_id] = "unresolved"
            unresolved_notes.setdefault(gene_id, []).append(f"{rule['elimination_rule_id']}={new}")

        needs_review = str(rule["required_review"]).strip().casefold() in {"true", "1", "yes"} or not (quality_ok and coverage_ok)
        events.append({
            "event_id": f"elim_{len(events):04d}", "candidate_gene_id": gene_id, "hypothesis_id": hypothesis_id,
            "rule_id": str(rule["elimination_rule_id"]), "evidence_id": evidence_id,
            "old_state": old, "new_state": new, "scope": scope, "strength": strength,
            "quality_gate_passed": quality_ok, "coverage_gate_passed": coverage_ok,
            "review_status": "review_required" if needs_review else "applied",
            "reason": str(rule["rationale"]),
        })
        if needs_review:
            review_rows.append({"reason": f"elimination_review:{rule['elimination_rule_id']}", "target": gene_id,
                                "required_next_action": ("quality/coverage gateを満たしていないため除外に使えない。不足scopeを埋める"
                                                         if not (quality_ok and coverage_ok)
                                                         else "この遷移は人のreviewを要求するruleである"),
                                "allowed_choices": "証拠を補う | reviewして確定 | 未評価のまま保持"})

    survival_rows: list[dict[str, object]] = []
    remaining_rows: list[dict[str, object]] = []
    for gene_id in sorted(gene_state):
        state = gene_state[gene_id]
        mechanisms = {key[1]: value for key, value in mechanism_state.items() if key[0] == gene_id}
        surviving = [name for name, value in mechanisms.items() if str(value["state"]).endswith("_survives")]
        rejected = [name for name, value in mechanisms.items() if str(value["state"]).endswith("_rejected")]
        uncertainty = ";".join(unresolved_notes.get(gene_id, [])) or ("no elimination evidence supplied" if not len(evidence) else "")
        # A rejected coding hypothesis never promotes itself into a gene-level
        # exclusion: the gene survives with its other mechanisms.
        if rejected and state in {"not_tested", "survives"}:
            state = "survives"
        if state == "not_tested" and (surviving or rejected):
            state = "survives"
        gene_state[gene_id] = state
        discriminator = ("informative recombinant / cosegregation genotyping" if state in {"survives", "weakened"}
                         else ("none: already excluded" if state in EXCLUDED_STATES else "supply the missing coverage or measurement"))
        survival_rows.append({
            "candidate_gene_id": gene_id, "pre_elimination_group": pre_group.get(gene_id, ""),
            "current_gene_state": state, "surviving_hypothesis_count": len(surviving),
            "rejected_hypothesis_count": len(rejected),
            "hard_exclusion_count": 1 if state in EXCLUDED_STATES else 0,
            "weakened_reason": ";".join(weakened_reason.get(gene_id, [])),
            "strongest_exclusion_evidence_id": hard_evidence.get(gene_id, ""),
            "remaining_uncertainty": uncertainty,
            "recommended_next_discriminator": discriminator,
            "rank_contribution": "none",
            "review_status": "review_required" if (uncertainty or weakened_reason.get(gene_id)) else "context",
        })
        if state not in EXCLUDED_STATES:
            remaining_rows.append({
                "candidate_gene_id": gene_id, "current_gene_state": state,
                "surviving_mechanisms": ";".join(sorted(
                    f"{mechanisms[name]['scope']}:{name}" for name in surviving)) or "not_evaluated",
                "first_wet_priority": "high" if state == "survives" and surviving else ("low" if state == "weakened" else "medium"),
                "rank_contribution": "none",
                "notes": "survival is a triage state; it neither adds nor removes score",
            })

    if not len(evidence):
        review_rows.append({"reason": "no_elimination_evidence_supplied", "target": case_id or "case",
                            "required_next_action": "除外証拠が1件も無い。候補は全件 not_tested のまま圧縮されない",
                            "allowed_choices": "証拠を供給 | 圧縮せずに進む"})

    survival = _frame(survival_rows, SURVIVAL_COLUMNS)
    hypothesis_rows = [{
        "candidate_gene_id": key[0], "hypothesis_id": key[1], "scope": str(value["scope"]),
        "mechanism_state": str(value["state"]), "changed_by_rule": str(value["rule"]),
        "evidence_id": str(value["evidence_id"]), "review_status": "review_required" if str(value["state"]).endswith("_weakened") else "context",
        "notes": "a rejected mechanism never becomes a gene-level exclusion",
    } for key, value in sorted(mechanism_state.items())]
    hypothesis_survival = _frame(hypothesis_rows, HYPOTHESIS_SURVIVAL_COLUMNS)
    new_events = _frame(events, EVENT_COLUMNS)
    all_events = pd.concat([history[EVENT_COLUMNS] if not history.empty and set(EVENT_COLUMNS).issubset(history.columns) else pd.DataFrame(columns=EVENT_COLUMNS),
                            new_events], ignore_index=True)
    remaining = _frame(remaining_rows, REMAINING_COLUMNS)
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")

    paths = {"candidate_survival": out / "candidate_survival.tsv",
             "hypothesis_survival": out / "hypothesis_survival.tsv",
             "elimination_events": out / "elimination_events.tsv",
             "remaining_candidate_set": out / "remaining_candidate_set.tsv",
             "elimination_review_required": out / "elimination_review_required.tsv"}
    for frame, path in ((survival, paths["candidate_survival"]), (hypothesis_survival, paths["hypothesis_survival"]),
                        (all_events, paths["elimination_events"]), (remaining, paths["remaining_candidate_set"]),
                        (reviews, paths["elimination_review_required"])):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    hard_excluded = int(sum(1 for value in gene_state.values() if value in EXCLUDED_STATES))
    mechanism_only = int(sum(1 for row in survival_rows if row["rejected_hypothesis_count"] and row["current_gene_state"] not in EXCLUDED_STATES))
    unresolved = int(sum(1 for value in gene_state.values() if value in {"unresolved", "not_tested"}))
    risk_flags: list[str] = []
    if new_events.empty is False and bool((new_events["strength"] == "none").any()):
        risk_flags.append("non-evidence events present: assay failure / partial coverage / unmeasured data were recorded, not used to eliminate")
    if unresolved:
        risk_flags.append(f"{unresolved} candidate(s) remain unresolved or untested; compression is incomplete, not negative")
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "initial_candidates": int(len(candidates)), "remaining_candidates": int(len(remaining)),
        "candidate_reduction_ratio": round(len(candidates) / len(remaining), 2) if len(remaining) else None,
        "hard_excluded_count": hard_excluded, "mechanism_only_rejected_count": mechanism_only,
        "unresolved_count": unresolved,
        "gene_state_counts": {str(k): int(v) for k, v in survival["current_gene_state"].value_counts().sort_index().items()} if not survival.empty else {},
        "mechanism_state_counts": {str(k): int(v) for k, v in hypothesis_survival["mechanism_state"].value_counts().sort_index().items()} if not hypothesis_survival.empty else {},
        "false_exclusion_risk_flags": risk_flags,
        # Historical field name kept for artifact compatibility; the value is
        # the number of candidates REMOVED from the candidate space (initial
        # minus remaining), never the remaining count itself.
        "candidates_potentially_eliminated": int(len(candidates) - len(remaining)),
        "next_experiment_hint": "prefer the experiment that can discriminate the most remaining candidates, not the top-ranked one",
        "rule_registry": rules_meta,
        "strength_contract": {
            "hard": "positional / genetic evidence; the only strength allowed to exclude a gene",
            "strong": "mechanism-specific negative evidence at the required quality and coverage; rejects one mechanism, never the gene",
            "moderate": "discriminating counterevidence; may only weaken",
            "none": "assay failure / partial coverage / unmeasured / missing annotation / provider failure; never eliminates",
        },
        "inputs": {name: (file_fingerprint(runtime_path(value)) if value else None) for name, value in (
            ("candidates", candidates_arg), ("evidence", evidence_arg), ("hypotheses", hypotheses_arg),
            ("wet_state", wet_state_arg), ("previous_history", previous_history_arg))},
        "safety": [
            "absence of evidence is never converted into evidence of absence",
            "an assay failure or partial coverage never removes a candidate",
            "unmeasured expression is never treated as non-expression",
            "a missing homolog or annotation never removes a candidate",
            "a rejected coding hypothesis weakens that mechanism only; the gene keeps its regulatory and SV branches",
            "only hard positional or genetic evidence may exclude a gene, and the registry refuses any other rule that tries",
            "every state change is appended with its rule id and evidence id; nothing is overwritten silently",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "elimination_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "initial_candidates": int(len(candidates)),
            "remaining_candidates": int(len(remaining)), "hard_excluded_count": hard_excluded,
            "mechanism_only_rejected_count": mechanism_only, "unresolved_count": unresolved,
            **{name: str(path) for name, path in paths.items()}}
