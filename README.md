# candidate-qtl

**English** | [日本語](README.ja.md)

`candidate-qtl` is a research workflow and CLI that takes a plant QTL interval and a phenotype and produces an experimentally manageable candidate set for Wet validation. It uses literature and homolog knowledge from model and other plants to keep plausible causal genes, while preserving the evidence behind every decision. About 20 genes is a guide: it is not a minimum count, not a truncation cap, and not a completion condition. Mapping remains independent evidence; the tool assists map-based cloning and does not replace it or claim automatic causal-gene identification.

## The flow

```text
phenotype + QTL interval + reference FASTA/GFF3
  ├─ Route A  (required)  phenotype-first: mechanism / family / reference hypothesis
  │                       -> verified source and sequence -> deterministic projection into the interval
  ├─ Route C  (required)  Direct C, interval-blind: reported trait-gene / locus / variant literature
  │                       -> freeze -> projection into the interval
  └─ Route B  (optional)  interval-first: all interval genes -> family/orthogroup compression
                          -> semantic interpretation -> verification
        each route is frozen independently
  -> claim review -> per-paper checklist (nine answers) -> Tier T1–T6/U
  -> canonical Wet membership -> wet_candidates.tsv
       only the first Tier that overflows the review budget goes to boundary review
```

- **Route A** starts from the phenotype, proposes mechanisms and reference genes, verifies the cited source, accession and sequence, and projects verified references into the interval deterministically.
- **Route C (Direct C)** collects directly reported trait–gene/locus/variant relationships without seeing the interval or the other routes, then projects them. A searched-but-empty result is recorded and distinguished from a search that was not run.
- **Route B** is optional. It summarizes every interval gene through family/orthogroup units (an intermediate representation, not a "one gene per family" answer) and interprets them against the phenotype.
- **Independence.** A, B and C stay procedurally independent until each is frozen. Agreement afterwards is provenance, not an extra evidence count or score bonus.
- **Missing work is reported, not invented.** Unrun or failed work is `not_evaluated`, never biological negative evidence. Without Route B, Wet holds only A/C candidates and is never padded from annotation; when it is empty or small, the handoff signals `route_a_rethink` so the parent rethinks Route A once in a fresh generation.

## Tier judgement

For each paper, the original researcher (a human or an LLM agent) records the paper ID/title, the reported gene name/ID/species, nine yes/no/unknown checklist answers and short reasons. A decision table turns those answers into a Tier; the strongest applicable paper sets the gene's Tier, with no bonus for paper count.

| Tier | Meaning |
|---|---|
| T1 | Same species, same or very close reported phenotype, plus a grounded interval homolog projection |
| T2 | Reported gene–phenotype relationship for the same or a very similar cellular/developmental phenomenon |
| T3 | Reported phenotype elsewhere, linked by a concrete shared process |
| T4 | Hypothesis grounded in reported function |
| T5 | Hypothesis grounded in a specific, sourced annotation |
| T6 / U | Identity only / unresolved; kept in the full ledger, never treated as negative |

Tiers express relevance to the phenotype, not causal probability. Relevance, verified support, Wet adoption and investigation coverage are separate axes. Machines check IDs, sequences, citations, quoted spans, interval anchors and hashes; they never replace the researcher's recorded gene name/ID with a database lookup. The full rules live in [canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md) and the [LLM/machine boundary](docs/LLM-BOUNDARY.ja.md).

## Outputs

Report unresolved/pending states and actual investigation coverage together with the set.

- `wet_candidates.tsv` — the adopted Wet set (the only "final candidate set").
- `wet_relevance_inventory.tsv` — every interval gene with its relevance and assessment state.
- `wet_candidate_trace.tsv` — handoff observables per gene.

## Repository layout

| Path | Contents |
|---|---|
| `src/candidate_qtl/` | Python package and the `candidate-qtl` CLI |
| `workflow/` | Ordinary-case entry point `practical_packet_run.py` and its config template `practical_packet_run.ordinary.template.json` |
| `scripts/` | Ordinary-case helpers for the parent agent (Route B sharding, claim review session, boundary review) and local CI |
| `docs/` | Contracts, runbooks and CLI documentation (mostly Japanese) |
| `templates/`, `config/` | Input templates and example configuration |
| `tests/` | Unit and integration tests with synthetic fixtures |

## Install

- Version: `0.24.0`
- Python 3.12+. Use Linux or WSL with the `alignment` extra for FASTA/GFF processing.

```bash
git clone https://github.com/light-suzuki/candidate-qtl.git
cd candidate-qtl
python3 -m venv .venv
source .venv/bin/activate
pip install -c constraints.txt -e ".[alignment]"
candidate-qtl doctor --json
```

On Windows, `./scripts/install-local.ps1 -AddToUserPath` is also available.

## Running an ordinary case

A parent agent (or a person) owns the case and dispatches Route A, Route C and optional Route B to separate workers that start without the parent's conversation history. The wrapper below transports and checks saved responses; it does not launch researchers or provide search tools.

```bash
h=workflow/practical_packet_run.py
run=/path/to/run
cp workflow/practical_packet_run.ordinary.template.json "$run/config.json"   # then edit
python $h prepare --config "$run/config.json" --out-dir "$run/01-prepared"
# A/C-first without Route B: mark every B prompt not_evaluated, then review
python scripts/agentic_route_b.py merge --prompts "$PROMPTS" --out "$run/b-not-run.json"
python $h review --prepared "$PREPARED_MANIFEST" --sha256 "$PREPARED_SHA" \
  --responses "$run/b-not-run.json" --route-b-not-run "<reason>" --out-dir "$run/02-reviewed"
python $h finish --reviewed "$REVIEWED_MANIFEST" --sha256 "$REVIEWED_SHA" \
  --route-c-members "$C_MEMBERS" --route-c-manifest "$C_MANIFEST" --out-dir "$run/03-first-wet"
```

Each stage prints its `manifest` and `sha256` for the next stage. With Route B, merge the real shard responses instead and pass `--b-execution-manifest`.

See the [ordinary-case runbook](docs/AGENTIC-QTL-CASE.ja.md) for Route A payloads, Direct C (`route-c-direct-freeze` → `route-c-direct-project` → `route-c-direct-members`), Route B sharding with `scripts/agentic_route_b.py`, and boundary review. Command details: [CLI guide](docs/CLI.ja.md) and [CLI reference](docs/CLI-REFERENCE.md).

## Tests

```bash
python scripts/run_test_partition.py --partition core
scripts/local_ci.sh test lint
```

## Documentation

| Topic | Document |
|---|---|
| Goal and architecture | [Core overview](CORE-README.ja.md) / [core contract](docs/CORE-CONTRACT.ja.md) |
| Tier and Wet adoption policy | [Canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md) |
| Running a case | [Ordinary-case runbook](docs/AGENTIC-QTL-CASE.ja.md) / [three-route runbook](docs/TRIPLE-ROUTE-RUNBOOK.ja.md) |
| Direct C | [Direct C contract](docs/ROUTE-C-DIRECT-RESCUE.ja.md) |
| Rules for coding/research agents | [AGENTS.md](AGENTS.md) |

## Status and limits

Retention of known causal genes in development cases is diagnostic evidence only; accuracy on unseen cases is not established. Filling 20 slots does not guarantee that an unknown causal gene is retained. Benchmark datasets and evaluation tooling are not part of this public release.

## License and citation

[MIT License](LICENSE). Citation metadata: [CITATION.cff](CITATION.cff). Contributors: [CONTRIBUTORS.md](CONTRIBUTORS.md).
