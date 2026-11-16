# pea_Tl wrong-paralog regression fixture

A small, reproducible stand-in for the one V2 case (`pea_Tl`) where a wrong
paralog is the single top-scoring candidate. It fixes the **detect-and-hold**
behaviour so it cannot silently regress, and it does so *without* tuning any
weight to push the wrong candidate down.

## The two paralogs

| gene | role | evidence | dry_score | dry_group |
|------|------|----------|----------:|-----------|
| `Psat_TlB` | wrong paralog | one `curated_homolog_match` chain (two UniProt hits from a single comparative-annotation chain → **one** independent source) | 1.25 | `rescue` (held) |
| `Psat_TlA` | true gene | `peak_proximal` + `expression_detected`, both from `PMID:1` (two non-redundant channels, but **one** source) | 1.126 | `rescue` |

`Psat_TlB` leads on raw score, so it is the **sole top by score**
(`sole_top_gene_id = Psat_TlB`, `wrong_paralog_sole_top = True`).

## Why this is safe, not a failure

Neither gene has two *independent* citable sources, so the ranking gate promotes
**neither** to the Top group:

* `Psat_TlB`: `insufficient_independent_evidence` — the two UniProt hits collapse
  to one comparative-annotation chain, so it is held, not confidently promoted.
* `Psat_TlA`: `multiple_non_redundant_channels` — two channels, one source.

The missing independent signal is the same for both: **a second, independent
citable line of evidence** (an independent study, an orthogroup/synteny call, or
cosegregation) that would separate the true copy from its paralog. The right fix
is to add that signal, not to re-weight homology.

## What the regression test pins

* the true gene `Psat_TlA` is **retained** (present, `rank_max == 2`), never lost;
* the wrong paralog is **detected** as the sole top (`wrong_paralog_sole_top`);
* the system **holds** — `confident_wrong_sole_top == False`,
  `resolution_outcome == "safe_hold"` — instead of confidently emitting the wrong
  paralog;
* `Psat_TlB` still leads on score, documenting that the safe outcome comes from
  the independence gate, not from weight tuning.
