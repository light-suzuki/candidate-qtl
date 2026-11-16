# Shared source-bound claim contract

Canonical Route A/B/C support requires more than a resolved citation, matching
terms, a genetic-test cue, or a verbatim quote. The workflow must bind the
source record, subject gene, experiment or allele, asserted relation, target
phenotype/process, organism/accession, and quoted span as one claim.

The tab-separated source-bound-claim-1.0.0 artifact uses the columns exported
by candidate_qtl.claim_binding.BOUND_CLAIM_COLUMNS. A row can authorize
canonical support only when:

- binding_state is verified_complete;
- the subject identity is verified_exact or verified_mapped;
- the experiment is verified_same_subject and names that same subject;
- the submitted predicate is independently verified;
- the phenotype/process is the target_trait;
- organism and accession states are verified or explicitly not applicable;
- the quote is bound to the verified_same_claim in an abstract, full text, or
  curated database record;
- curation is deterministic or human-reviewed; and
- rank_contribution remains none.

Missing rows or axes are not_evaluated; they are not negative biological
evidence. Contradicted axes may remain explicit. Related-trait records,
secondary reviews, background mentions, ambiguous multi-gene experiments, and
curated function without a direct phenotype relation remain review-only.

Route A accepts the artifact through hypothesis-verify --bound-claims and
records its state separately from accession/species/sequence verification.
bait-panel-freeze is the canonical promotion gate and skips verified
accessions whose bound claim is incomplete.

Route B accepts it through
route-b-semantic-transfer-verify --bound-claims. Source verification remains
triage; only a complete matching bound claim plus the existing frozen
identity/transfer checks can produce verified membership.

Route C keeps its stricter sealed route-c-source-claim-binding-1.0.0 sidecar.
Its member gate adapts that sidecar through the same shared authorization API
while retaining all stable-ID, species/accession, projection, and frozen-row
hash checks.
