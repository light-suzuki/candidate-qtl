# CLI reference

Generated from the argparse registry. Run `python scripts/generate_cli_reference.py` after an intentional CLI change.

## candidate-qtl-bench

```text
usage: candidate-qtl-bench [-h] [--version] [--traceback]
                           {assembly-bridge,multi-reference-locus-union,exon-coverage-check,benchmark-reference,primerblast-import,assay-plan,primer-plan,exon-sequencing-plan,collect-sanger-markers,blast-panel-freeze,blast-panel-verify,blast-panel-search,budget-run,supported-member-union,candidate-budget-gate,wet-candidate-handoff,adaptive-synthesis-request,adaptive-synthesis-import,adaptive-refinement-apply,llm-usage,budget-llm-review,reference-panel,external-rank,propagate-terms,route-b-semantic-bundle,route-b-semantic-import,route-b-semantic-source-verify,route-b-semantic-transfer-verify,route-b-semantic-candidates,semantic-claims-freeze,comparative-locus-reconstruct,comparative-locus-benchmark,doctor,provenance-ledger,discover,discover-promote,homolog-panel-auto,homolog-panel-match,homolog-panel-benchmark,auto-triage,candidate-elimination,direction-evidence,evidence-build,evidence-term-promote,peak-context,trait-terms,trait-ontology-map,qtl-source-init,extract,gff-preflight,annotation-rescue,liftover-markers,dark-candidates,candidate-flags,annotation-keywords,annotate-plan,family-candidate-bridge,homolog-linking-export,homolog-linking-import,homolog-linking-counterfactual,homolog-linking-probe-export,homolog-linking-score,research-packet,llm-prompt-export,llm-shortlist-prompt,llm-shortlist-import,llm-review-benchmark,llm-run,llm-evidence-import,llm-pathway-prompt,llm-pathway-import,hypothesis-verify,bait-panel-freeze,phenotype-pathway-benchmark,scaffold-rescue,cosegregation,qtl2-plan,guided-profile,case-run,orchestrate,family-index,interval-atlas,family-to-copy,integrated-narrowing,phenotype-first-status,interval-reverse-lookup,route-b-prepare,reverse-literature,dual-route-compare,pathway-interval-falsify,route-a-interval-match,homology-concordance,reference-candidate-table,orthology-plan,orthology-import,orthology-tier1-plan,focal-paralog-import,orthofinder-import,omamer-import,synteny-import,build-executor-package,rank,dry-prioritize,evidence-rank,ranking-qa,ranking-stability-benchmark,screening-shortlist,evidence-tiers,evidence-tier-benchmark,weight-policy-benchmark,mapping-gate,next-experiment,candidate-pipeline,scenario-create,scenario-report,tnt1-evidence-normalize,variant-evidence,multigenome-private,sv-plan,promoter-review,allele-compare,wet-state,hypothesis-ledger,phenotype-investigate,phenotype-localization-packet,phenotype-localization-freeze,mechanism-import,mechanism-coverage-audit,investigation-attribution,reference-gene-resolve,reference-gene-cache,source-verify,route-c-direct-search,route-c-direct-freeze,route-c-provider-bait,route-c-direct-project,route-c-direct-members,triple-route-compare,claim-verification-benchmark,literature,research,reference-corpus-fetch,reference-corpus-build,reference-corpus-compare,corpus-mask,retrieval-index-build,retrieval-run,retrieval-merge,retrieval-hybrid-merge,retrieval-import,retrieval-interval-match,route-c-benchmark,route-c-provider-comparison,route-c-direct-benchmark,phenotype-variant-benchmark}
                           ...

Evidence-preserving plant QTL candidate triage and experiment decision-support workflow.

positional arguments:
  {assembly-bridge,multi-reference-locus-union,exon-coverage-check,benchmark-reference,primerblast-import,assay-plan,primer-plan,exon-sequencing-plan,collect-sanger-markers,blast-panel-freeze,blast-panel-verify,blast-panel-search,budget-run,supported-member-union,candidate-budget-gate,wet-candidate-handoff,adaptive-synthesis-request,adaptive-synthesis-import,adaptive-refinement-apply,llm-usage,budget-llm-review,reference-panel,external-rank,propagate-terms,route-b-semantic-bundle,route-b-semantic-import,route-b-semantic-source-verify,route-b-semantic-transfer-verify,route-b-semantic-candidates,semantic-claims-freeze,comparative-locus-reconstruct,comparative-locus-benchmark,doctor,provenance-ledger,discover,discover-promote,homolog-panel-auto,homolog-panel-match,homolog-panel-benchmark,auto-triage,candidate-elimination,direction-evidence,evidence-build,evidence-term-promote,peak-context,trait-terms,trait-ontology-map,qtl-source-init,extract,gff-preflight,annotation-rescue,liftover-markers,dark-candidates,candidate-flags,annotation-keywords,annotate-plan,family-candidate-bridge,homolog-linking-export,homolog-linking-import,homolog-linking-counterfactual,homolog-linking-probe-export,homolog-linking-score,research-packet,llm-prompt-export,llm-shortlist-prompt,llm-shortlist-import,llm-review-benchmark,llm-run,llm-evidence-import,llm-pathway-prompt,llm-pathway-import,hypothesis-verify,bait-panel-freeze,phenotype-pathway-benchmark,scaffold-rescue,cosegregation,qtl2-plan,guided-profile,case-run,orchestrate,family-index,interval-atlas,family-to-copy,integrated-narrowing,phenotype-first-status,interval-reverse-lookup,route-b-prepare,reverse-literature,dual-route-compare,pathway-interval-falsify,route-a-interval-match,homology-concordance,reference-candidate-table,orthology-plan,orthology-import,orthology-tier1-plan,focal-paralog-import,orthofinder-import,omamer-import,synteny-import,build-executor-package,rank,dry-prioritize,evidence-rank,ranking-qa,ranking-stability-benchmark,screening-shortlist,evidence-tiers,evidence-tier-benchmark,weight-policy-benchmark,mapping-gate,next-experiment,candidate-pipeline,scenario-create,scenario-report,tnt1-evidence-normalize,variant-evidence,multigenome-private,sv-plan,promoter-review,allele-compare,wet-state,hypothesis-ledger,phenotype-investigate,phenotype-localization-packet,phenotype-localization-freeze,mechanism-import,mechanism-coverage-audit,investigation-attribution,reference-gene-resolve,reference-gene-cache,source-verify,route-c-direct-search,route-c-direct-freeze,route-c-provider-bait,route-c-direct-project,route-c-direct-members,triple-route-compare,claim-verification-benchmark,literature,research,reference-corpus-fetch,reference-corpus-build,reference-corpus-compare,corpus-mask,retrieval-index-build,retrieval-run,retrieval-merge,retrieval-hybrid-merge,retrieval-import,retrieval-interval-match,route-c-benchmark,route-c-provider-comparison,route-c-direct-benchmark,phenotype-variant-benchmark}
    assembly-bridge     Bridge native same-species assembly IDs using optional LIS pangene candidate correspondences,
                        fam3, and explicit pairwise reciprocal sequence evidence; emit mapping and existing Route B
                        TSV adapters without symbol-only matching.
    multi-reference-locus-union
                        Build a full-ID-preserving multi-reference candidate universe and conservative locus-entity
                        union; only verified 1:1 mappings collapse duplicates.
    exon-coverage-check
                        Detect CDS windows PrimerBlast OSS results do not span and plan relaxed retries.
    benchmark-reference
                        Run known-gene virtual-QTL recovery cases on one reference without using answer IDs for
                        scoring.
    primerblast-import  Normalize PrimerBlast OSS JSON reports into an orderable marker ledger.
    assay-plan          Legacy top-candidate primerblast-oss assay plan.
    primer-plan         Generate primerblast-oss plans for direct sequencing, parental-variant markers, or a coarse
                        QTL marker grid.
    exon-sequencing-plan
                        Tile every Top-candidate representative CDS for parental PCR/Sanger.
    collect-sanger-markers
                        Collect PrimerBlast OSS tile results into an ordering table.
    blast-panel-freeze  Freeze a panel TSV into an immutable bait panel plus a manifest. Only the panel is read;
                        interval, candidate, truth and score inputs are structurally absent.
    blast-panel-verify  Re-hash a blast-panel freeze manifest against its frozen artifacts and refuse any drift.
    blast-panel-search  BLAST the verified frozen baits against the interval proteins and emit auditable blast
                        evidence. Missing tools or a provider failure degrade to not_evaluated_*, never a miss.
    budget-run          Deterministic budgeted engine: tier0 inventory, tier1 clusters, tier2 plan, handoff, resource-
                        usage ledger (issue #375).
    supported-member-union
                        Build the non-scoring union of supported Route-A bridge members and final verified Route-B
                        semantic members (#426).
    candidate-budget-gate
                        Freeze every supported candidate; report the Wet target excess without truncation or a count-
                        only stop.
    wet-candidate-handoff
                        Canonical path from frozen Route A/B artifacts to a Wet-reviewable candidate set and
                        stage_a_narrowed (#426).
    adaptive-synthesis-request
                        Emit the bounded ambiguity-resolution request for an over-budget candidate set; refuses to run
                        within budget (#426/#434).
    adaptive-synthesis-import
                        Normalize reviewer ambiguity-resolution proposals; rank_contribution stays none and no
                        proposal narrows the set (#426).
    adaptive-refinement-apply
                        Narrow a candidate set using machine-verified verdicts only; unsourced or not_evaluated
                        verdicts never remove a candidate (#426).
    llm-usage           Aggregate provider-reported token usage/call counts from llm-run raw responses without
                        estimating missing tokens or prices (#434).
    budget-llm-review   Batched, budget-gated LLM cluster review + critic + portfolio reducer over an existing budget-
                        run output directory (issue #375 T5).
    reference-panel     Validate and availability-audit a comparative reference-resource panel (native required,
                        others optional with machine-readable skip reasons).
    external-rank       Import external comparative ranks into a coverage-aware ledger (RRF consensus) without
                        touching the native dry_score.
    propagate-terms     Deterministically propagate reference-ortholog function terms onto QTL interval candidates (no
                        LLM); native-only input degrades gracefully.
    route-b-semantic-bundle
                        After B0 freeze, compress the complete interval atlas into identity-preserving
                        family/orthogroup/gene units and emit provider-neutral semantic prompts (#403).
    route-b-semantic-import
                        Import normalized llm-run responses for frozen Route-B semantic prompts. Unknown
                        sources/reference genes are unsupported; all other rows remain proposals pending relation
                        verification (#403).
    route-b-semantic-source-verify
                        Verify proposal citations/source relations with the existing source-verify engine while
                        keeping focal-family transfer as a separate unresolved axis (#403).
    route-b-semantic-transfer-verify
                        Combine source-level mechanistic support with frozen family/orthogroup/copy identity. Group-
                        level claims may verify without choosing one copy; interval_gene claims require one-to-
                        one/RBH-like transfer (#403).
    route-b-semantic-candidates
                        Freeze Route-B semantic candidate members from verified claims before comparison (#403).
    semantic-claims-freeze
                        Validate and freeze structured Route-B semantic claims; only verified source-backed rows are
                        narrowing-eligible (#403).
    comparative-locus-reconstruct
                        Reconstruct an optional candidate universe from flanking markers and precomputed multi-species
                        locus caches; never adds causal score.
    comparative-locus-benchmark
                        Evaluate frozen comparative reconstructions on synthetic/degraded-reference cases.
    doctor              Inspect required and optional WSL tools.
    provenance-ledger   Build and audit the normalized provenance DAG: nodes, edges, and shared derivation paths
                        (issue #231).
    discover            Run preliminary candidate screening from genome, interval, and phenotype inputs; promote
                        reviewed outputs to case-run for the full workflow.
    discover-promote    Create a reviewable case-run manifest from preliminary discover outputs; writes the discover
                        run's ranking defaults (min_top_score, top_n, reserve_n, rescue_n) into case.yaml.
    homolog-panel-auto  Generate a sequence-backed, pending-review homolog panel from phenotype text without
                        inspecting the QTL interval.
    homolog-panel-match
                        Match an automated homolog panel only against interval proteins and optionally freeze the
                        deterministic Stage-A 20.
    homolog-panel-benchmark
                        Blindly compare automated phenotype-to-homolog retrieval on frozen known-gene cases.
    auto-triage         Process every interval candidate against role-labelled reference proteomes without a curated
                        reference-gene panel.
    candidate-elimination
                        Manage candidate and mechanism survival state. Only hard positional/genetic evidence may
                        exclude a gene; mechanism-specific negatives reject one mechanism; assay failure, partial
                        coverage and unmeasured data never eliminate anything.
    direction-evidence  Mark phenotype-direction matches and counterevidence deterministically.
    evidence-build      Normalize annotation, expression, position, and manual evidence into one source-backed ledger.
    evidence-term-promote
                        Promote approved annotation terms into low-confidence deterministic screening evidence.
    peak-context        Import QTL peak, LOD support interval, and recombinant core into auditable position evidence.
    trait-terms         Build a controlled, approval-gated trait-term set for keyword evidence.
    trait-ontology-map  Map free trait terms to ontology IDs using a local synonym table only.
    qtl-source-init     Validate coordinate provenance and choose a safe liftover route.
    extract             Extract genes and sequences from one QTL interval.
    gff-preflight       Cross-check candidate-qtl GFF interpretation against optional AGAT/gffread outputs.
    annotation-rescue   Project known proteins only into a QTL interval and audit missing/broken native models without
                        changing rank.
    liftover-markers    Resolve pair-level marker hits into safe target interval blocks.
    dark-candidates     Export short complete no-homology ORFs as an explicit dark bucket.
    candidate-flags     Flag tandem redundancy and dark short-ORF candidates.
    annotation-keywords
                        Link best protein hits to GFF descriptions and emit auditable trait-keyword evidence.
    annotate-plan       Create candidate-only protein, GO evidence template, and local homology/domain commands.
    family-candidate-bridge
                        Normalize Route A / B / C family evidence (phenotype-first hypotheses, orthology/atlas, frozen
                        retrieval match) into one reporting table with a strict mapping tier
                        (tier1_source_backed_unique .. tier6_not_evaluated). Reporting-only: rank_contribution is
                        always none, there are no score columns, and the family_to_copy / narrowing wiring is
                        untouched. At least one route input is required; a Route B atlas is independently sufficient.
    homolog-linking-export
                        Export provider-neutral prompt package for agent-in-the-loop homolog linking.
    homolog-linking-import
                        Mechanically validate and import agent linking proposals with strict closed-world and
                        grounding substring gates.
    homolog-linking-counterfactual
                        Evaluate what truth Stage-A rank would be if linking proposals acted as support, without
                        altering actual rank.
    homolog-linking-probe-export
                        Export control memorization probe prompt with phenotype ONLY (no interval, no candidates).
    homolog-linking-score
                        Score linking results stratified by memorization probe. Probe match is three-valued (match /
                        mismatch / not_evaluated_probe_ambiguous): recall_positive only on a shared family_labels
                        alias/symbol key, recall_negative only on explicit probe abstention, and un-catalogued non-
                        matching relations are excluded from both (memory is never reported as reasoning). The alias
                        used and its provenance are recorded in the memorization_probe_score.tsv artifact.
    research-packet     Create AI-optional candidate research packets and prompts.
    llm-prompt-export   Export provider-neutral, reproducible candidate research prompts as JSONL.
    llm-shortlist-prompt
                        Export one optional Stage-B prompt that ranks the exact frozen Stage-A set.
    llm-shortlist-import
                        Validate an LLM batch ranking without changing Stage A or scoring evidence.
    llm-review-benchmark
                        Compare deterministic and optional LLM review ranks on identical frozen shortlists.
    llm-run             Run prompt JSONL through an API/LM Studio provider or normalize an offline agent artifact.
    llm-evidence-import
                        Import only human-reviewed, source-backed model claims. Each claim needs reviewed=true plus
                        reviewed_by; unresolved source IDs are marked and cannot count as independent.
    llm-pathway-prompt  Export a versioned phenotype-only prompt package asking for competing pathway/gene-family
                        hypotheses and literature-supported reference baits. The model never sees or ranks QTL
                        interval candidates.
    llm-pathway-import  Normalize raw model output into the strict pathway-proposal schema, keeping the raw payload
                        verbatim. Every proposal stays pending_verification.
    hypothesis-verify   Verify proposed accessions against an authoritative database: existence, species match,
                        retrievable sequence, source-backed function. Unverified proposals never become baits.
    bait-panel-freeze   Freeze verified proposals into an immutable bait panel, FASTA, and a query manifest consumable
                        by query-eligibility.
    phenotype-pathway-benchmark
                        Score the phenotype-only proposal lane stage by stage (pathway/family recall, accession
                        existence, verification pass, downstream causal recovery). Never pooled with the map-based
                        lane.
    scaffold-rescue     Conservatively rescue unplaced / unanchored scaffolds that are plausibly linked to the QTL
                        locus, as a separate candidate class. Being unplaced is not itself a reason; family agreement
                        alone is never high confidence.
    cosegregation       Summarize phenotype distribution by a candidate marker genotype.
    qtl2-plan           Validate F2 tables and generate or run an R/qtl2 interval scan.
    guided-profile      Start from a screening, full-case, answerless audit, or triple-route profile without hiding
                        the stable commands or review gates.
    case-run            Run/resume a deterministic manifest-driven case and emit a provenance report.
    orchestrate         Generate a laptop-safe OSS workflow for candidate annotation and primer/marker design.
    family-index        Build the cached genome-wide family index for one reference genome from a BLAST/DIAMOND hit
                        table. Run once per reference; per-case falsification is then table lookup only.
    interval-atlas      Structure every QTL-interval gene by family / domain / orthogroup / broad function without
                        using any phenotype. An exploration index for forming hypotheses from what the interval
                        actually contains; it never scores or promotes a gene.
    family-to-copy      Compare the members of a recovered gene family against each other and record what does and
                        does not discriminate between the copies. A first-Wet priority layer, never a causal claim.
    integrated-narrowing
                        Issue #303: compress a ranked candidate set into an ordinal Wet-review order by combining
                        family recovery (convergence), copy discrimination (family-to-copy), and survival state
                        (candidate-elimination). Every ranked column is passed through verbatim; narrowing_priority is
                        the global 1..N review order (tier-major), ordinal tiers N1-N4 only, never a numeric score,
                        rank_contribution stays none, and no candidate is ever excluded here (elimination alone may
                        exclude).
    phenotype-first-status
                        Route A (forward lookup): normalize the phenotype-first artifacts and detect the failure
                        signals that should trigger the interval-first route.
    interval-reverse-lookup
                        Route B (reverse lookup): start from the interval's genes, go out to their model-species
                        homologs, and derive family/process hypotheses from what is actually there. Interval-informed,
                        so never scored as blind discovery.
    route-b-prepare     Build the local Route B provider stack (homology search -> orthology table, cached source-
                        backed reference functions, interval atlas) and run the phenotype-free B0 reverse lookup.
                        Reports readiness: a broad-bin-only run is not_evaluated_for_reverse_lookup, never an
                        equivalent Route B result. The canonical atlas handoff is <out-dir>/interval_gene_atlas.tsv.
    reverse-literature  Route B stage 2: for every interval gene's reference homolog, query the literature focused by
                        the phenotype. A returned article is a link candidate, not evidence; an offline query is
                        not_evaluated, never 'no literature exists'.
    dual-route-compare  Compare the independently generated forward and reverse routes and record the convergence
                        tier: strict mechanistic / family / gene overlap / process-related / text-only, plus
                        phenotype-only / interval-only / discordant / both-weak. Convergence adds no score.
    pathway-interval-falsify
                        Falsify pathway/gene-family hypotheses against the QTL interval, coarse bucket first and the
                        fine ortholog level only for survivors. Enrichment filters hypotheses; it never scores a gene.
    route-a-interval-match
                        Route A stage (issue #355): forward-match only the verified/source-backed baits frozen by
                        bait-panel-freeze to concrete genes inside the QTL interval via a reference_genome-scoped
                        family index. Wires pathway-interval-falsify with verified targets only; zero verified baits,
                        a missing index or an unresolved interval degrade to not_evaluated, never to a miss.
    homology-concordance
                        Compare BLAST/DIAMOND/MMseqs2, domain, and orthology providers without vote inflation.
    reference-candidate-table
                        Create a review table linking QTL candidates to source-backed reference genes through
                        orthology mappings.
    orthology-plan      Generate reciprocal-BLAST refinement only for Top/Rescue candidates.
    orthology-import    Convert reciprocal-BLAST refinement results into canonical evidence rows.
    orthology-tier1-plan
                        Profile-aware Tier-1 plan: candidate query vs focal self search and a reference panel.
    focal-paralog-import
                        Flag out-of-interval paralog redundancy from a focal self-search.
    orthofinder-import  Import an OrthoFinder Orthogroups.tsv and classify 1:1/1:many/many:many.
    omamer-import       Import an OMAmer HOG placement table into orthogroup evidence.
    synteny-import      Import MCScanX/jcvi synteny anchors into synteny evidence.
    build-executor-package
                        Build a physically separated Executor package (code + executor manifests + declared reference
                        metadata only) with an allowlist package manifest and a run_capsule.json. Fails closed when
                        any truth/scoring artifact or link would enter the package.
    rank                Rank a candidate ledger with deterministic, inspectable rules.
    dry-prioritize      Rank parental-sequencing candidates using pre-Wet evidence only.
    evidence-rank       Rank candidates from the canonical evidence ledger with capped, inspectable contributions.
    ranking-qa          Audit ranking stability across eligibility, family dropout, weight perturbation, confidence,
                        provider, source-group dropout, and provider-dropout paths.
    ranking-stability-benchmark
                        Report whether frozen-case truth genes stay within Top-K under weight perturbation and
                        source/provider dropout.
    screening-shortlist
                        Freeze the non-LLM Stage-A screening shortlist before human or LLM review.
    evidence-tiers      Issue #229: report target_size-independent biological evidence tiers
                        (strong/plausible/unresolved/contradictory/knowledge_limited) side by side with the
                        configurable review_budget/assay_capacity operational cut. The operational cut never implies a
                        biological exclusion.
    evidence-tier-benchmark
                        Issue #229: Recall@budget alongside evidence-tier recall, stratified by interval-size band via
                        knowledge_coverage.benchmark_stratify. Reads frozen evidence-rank outputs only; never re-
                        ranks.
    weight-policy-benchmark
                        Compare numeric weight baselines and structural ranking gates on a frozen known-gene
                        benchmark.
    mapping-gate        Assess whether QTL mapping is adequate before candidate triage.
    next-experiment     Rank next experimental actions separately from biological candidate rank.
    candidate-pipeline  Run consistent primary evidence/ranking and optionally plan or integrate Top-candidate
                        orthology refinement.
    scenario-create     Build the canonical locus causal scenario ledger (locus_causal_scenarios.tsv) from scenario
                        rows; evidence-free rows stay unassessed and every row keeps rank_contribution=none.
    scenario-report     Render the locus causal scenario ledger as a review summary without choosing a winner.
    tnt1-evidence-normalize
                        Normalize curated Medicago Tnt1 evidence while keeping FST/insertion observations distinct
                        from causal mutant evidence and preserving cross-accession identity-bridge requirements.
    variant-evidence    Normalize a parental VCF plan and compare SnpEff/VEP transcript consequences without treating
                        predictions as Wet evidence.
    multigenome-private
                        Classify focal events against multiple control accessions as private_confirmed /
                        private_possible / shared_confirmed / not_evaluated. Keeps control callability and projection
                        denominators explicit; never adds score and never promotes a candidate.
    sv-plan             Plan assembly SV/PAV analysis and reconcile SyRI/read-based calls as unscored context.
    promoter-review     Limited non-promoting promoter/5'UTR/3'UTR/SV review route after the coding hypothesis is
                        examined. Emits strand-aware regions and intersects parental variants / SV-PAV-TE as unscored,
                        review-required context; never adds to a gene score.
    allele-compare      Call candidate variants from a pair of homologous FASTA sequences with minimap2.
    wet-state           Normalize assay coverage/results and derive append-only candidate Wet states.
    hypothesis-ledger   Track locus/mechanism hypotheses separately from gene rank.
    claim-verification-benchmark
                        Score the claim-support classifier against labelled examples (mechanistic support, associative
                        support, irrelevant mention, fabricated citation): per-class accuracy and a confusion table,
                        with cases whose record could not be retrieved reported not_evaluated instead of wrong.
    literature          Run one read-only Europe PMC literature query and save normalized JSON.
    research            Research selected candidates through read-only APIs without using an AI.
    reference-corpus-fetch
                        Download source-backed reference records (UniProtKB: function, disruption phenotype, pathway,
                        tissue, GO, evidence codes and the PubMed ids each block cites) for one or more taxa. Requires
                        --allow-external-api; there is no free-text query argument, so a case phenotype cannot shape
                        the corpus.
    reference-corpus-build
                        Build the versioned, phenotype-blind reference corpus that Routes B and C both consume.
                        Retrieved rows come from local authoritative sources; a model-proposed row is promoted only
                        when its accession appears in a retrieved source and a citation resolves. Unverified proposals
                        stay in an exclusion ledger.
    reference-corpus-compare
                        Compare two or more frozen corpora: per-field and per-species coverage (mutant phenotype,
                        pathway, GO, evidence) against a baseline label, plus optional Route C Recall@K recomputed
                        from frozen retrievals per corpus. Ranking weights are untouched and scores are never compared
                        across corpora.
    corpus-mask         Issue #237: remove the truth gene / its family / its close orthogroup from a frozen reference
                        corpus BEFORE prediction, so a benchmark can measure recovery of genuinely unfamiliar biology.
                        Mask rules come only from a benchmark-construction mask spec and are fingerprinted in the
                        manifest; no prediction command ever sees them, only the masked corpus TSV. Refuses to write a
                        masked corpus that still names the masked gene/family (self-check).
    retrieval-index-build
                        Route C: build an immutable, versioned local reference index from source-backed plant records.
                        Rows without a family/orthogroup or without a source id are excluded before Top-K and kept in
                        a ledger; every artifact is SHA-256 fingerprinted.
    retrieval-run       Route C: run a phenotype-only local retrieval against a frozen index and emit the retrieval
                        contract. The command has no interval/atlas/candidate argument, so Phase 1 cannot see the QTL
                        interval.
    retrieval-merge     Route C: merge two or more per-database retrieval TSVs by within-database rank. Raw scores are
                        never compared across databases; unannotated rows are excluded before Top-K selection.
    retrieval-hybrid-merge
                        Route C (issue #225): fuse two or more per-provider retrieval TSVs (e.g. lexical + lsa) by
                        reciprocal-rank fusion over each provider's within-provider rank. Raw retrieval_score values
                        are never compared or summed across providers -- only rank is read.
    retrieval-import    Route C: import an external phenotype/protein retrieval TSV, validate it against the
                        documented contract, and freeze it. No interval input exists on this command: the retrieval
                        stays phenotype-only and interval-blind until the separate interval-match command.
    retrieval-interval-match
                        Route C: match the frozen retrieval artifact against the QTL interval (interval atlas). Runs
                        only after retrieval-import; the frozen retrieval is never modified here.
    route-c-benchmark   Score frozen Route C retrieval artifacts on known-gene cases: Family Recall@5/@10/@20, Gene
                        Recall@20, and A-only / C-only / A+C union. Interval-blind by construction.
    route-c-provider-comparison
                        Route C (issue #225): run route-c-benchmark once per provider (e.g. lexical/lsa/hybrid) over
                        the same frozen known-gene cases and report Family/Gene Recall@K side by side. Never
                        restructures route-c-benchmark's own scoring and never pools a provider's raw retrieval_score.
    route-c-direct-benchmark
                        Score frozen Route C direct source recovery and post-freeze projection with separate
                        denominators. Contaminated/not-evaluated cases are excluded, never misses.
    phenotype-variant-benchmark
                        Score frozen per-formulation lane artifacts for phenotype wording / description-granularity
                        robustness: Route A family overlap, Route C Recall@K and retrieval overlap, Route B overlay
                        stability with an atlas-fingerprint check, A/B four-states, Stage-A stability, and stable /
                        rank-movement / catastrophic-family-dropout / stage-A-dropout classification. Never re-runs a
                        lane and never changes a score.

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --traceback           Show the Python traceback for expected input errors too (or set CANDIDATE_QTL_TRACEBACK=1).
```

## candidate-qtl-bench adaptive-refinement-apply

```text
usage: candidate-qtl-bench adaptive-refinement-apply [-h] --candidates CANDIDATES [--verification VERIFICATION]
                                                     --out-dir OUT_DIR [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --verification VERIFICATION
                        Machine verification TSV with candidate_gene_id and verdict. Omit to leave the set untouched.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench adaptive-synthesis-import

```text
usage: candidate-qtl-bench adaptive-synthesis-import [-h] [--proposals PROPOSALS] --request REQUEST
                                                     [--provider PROVIDER] [--model MODEL] --out-dir OUT_DIR
                                                     [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --proposals PROPOSALS
                        Reviewer proposals TSV. Omit to record not_evaluated_no_provider.
  --request REQUEST     adaptive_synthesis_request.json this response answers.
  --provider PROVIDER
  --model MODEL
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench adaptive-synthesis-request

```text
usage: candidate-qtl-bench adaptive-synthesis-request [-h] --candidates CANDIDATES --budget-manifest BUDGET_MANIFEST
                                                      [--ambiguity-context AMBIGUITY_CONTEXT] --out-dir OUT_DIR
                                                      [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --budget-manifest BUDGET_MANIFEST
                        candidate_budget_manifest.json; must report adaptive_refinement_required.
  --ambiguity-context AMBIGUITY_CONTEXT
                        Optional per-candidate ambiguity context.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench allele-compare

```text
usage: candidate-qtl-bench allele-compare [-h] --reference-fasta REFERENCE_FASTA --alternate-fasta ALTERNATE_FASTA
                                          [--reference-id REFERENCE_ID] [--alternate-id ALTERNATE_ID] --out-dir
                                          OUT_DIR [--minimap2 MINIMAP2] [--reference-candidates REFERENCE_CANDIDATES]
                                          [--min-coverage MIN_COVERAGE] [--sequence-type {genomic,cds,transcript}]
                                          [--json]

options:
  -h, --help            show this help message and exit
  --reference-fasta REFERENCE_FASTA
  --alternate-fasta ALTERNATE_FASTA
  --reference-id REFERENCE_ID
                        Record ID when reference FASTA contains multiple records.
  --alternate-id ALTERNATE_ID
                        Record ID when alternate FASTA contains multiple records.
  --out-dir OUT_DIR
  --minimap2 MINIMAP2   minimap2 executable name or path.
  --reference-candidates REFERENCE_CANDIDATES
                        Optional candidates.tsv used to emit genomic marker_candidates.tsv.
  --min-coverage MIN_COVERAGE
                        Minimum reference/alternate coverage for a complete comparison (0-1].
  --sequence-type {genomic,cds,transcript}
                        Declare the compared sequence type; controls frame-aware vs genomic variant classes.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench annotate-plan

```text
usage: candidate-qtl-bench annotate-plan [-h] --candidates CANDIDATES --out-dir OUT_DIR
                                         [--protein-fasta PROTEIN_FASTA] [--blastp-db BLASTP_DB]
                                         [--diamond-db DIAMOND_DB] [--hmmer-db HMMER_DB]
                                         [--eggnog-data-dir EGGNOG_DATA_DIR] [--eggnog-bin EGGNOG_BIN]
                                         [--num-threads NUM_THREADS] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate or ranked TSV containing gene_id.
  --out-dir OUT_DIR     Output directory.
  --protein-fasta PROTEIN_FASTA
                        Protein FASTA whose record IDs match gene_id.
  --blastp-db BLASTP_DB
                        BLASTP database; repeatable.
  --diamond-db DIAMOND_DB
                        DIAMOND database; repeatable.
  --hmmer-db HMMER_DB   HMMER profile database, for example a local Pfam database.
  --eggnog-data-dir EGGNOG_DATA_DIR
                        Optional eggNOG-mapper data directory for candidate-only GO/orthology annotation.
  --eggnog-bin EGGNOG_BIN
  --num-threads NUM_THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench annotation-keywords

```text
usage: candidate-qtl-bench annotation-keywords [-h] --candidates CANDIDATES --blast BLAST
                                               (--annotation-gff3 ANNOTATION_GFF3 | --subject-annotations SUBJECT_ANNOTATIONS | --use-stitle)
                                               --keyword KEYWORD --out OUT [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --blast BLAST         Candidate protein BLASTP outfmt 6 with candidate-qtl columns.
  --annotation-gff3 ANNOTATION_GFF3
                        GFF3/GFF3.GZ used to annotate BLAST subjects.
  --subject-annotations SUBJECT_ANNOTATIONS
                        Two-column subject_id/description TSV.
  --use-stitle          Read stitle from column 14 after the 13-column quality-gated BLAST fields.
  --keyword KEYWORD     Trait/function keyword or comma-separated keywords; repeatable.
  --out OUT             Canonical evidence TSV output.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench annotation-rescue

```text
usage: candidate-qtl-bench annotation-rescue [-h] --fasta FASTA --candidates CANDIDATES --query-proteins
                                             QUERY_PROTEINS --query-metadata QUERY_METADATA --interval INTERVAL
                                             --out-dir OUT_DIR [--projection-gff PROJECTION_GFF | --execute]
                                             [--projection-complete] [--miniprot MINIPROT] [--threads THREADS]
                                             [--review-limit REVIEW_LIMIT] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Target reference genome FASTA.
  --candidates CANDIDATES
                        Native candidates.tsv for the same interval.
  --query-proteins QUERY_PROTEINS
                        Source-backed homolog protein FASTA to project.
  --query-metadata QUERY_METADATA
                        TSV provenance for every query protein; source_id is mandatory.
  --interval INTERVAL   CHROM:START-END; projection is limited to this sequence.
  --out-dir OUT_DIR
  --projection-gff PROJECTION_GFF
                        Precomputed miniprot GFF3 against the generated interval-local FASTA.
  --execute             Run local miniprot now when available; otherwise emit a runnable script and review rows.
  --projection-complete
                        Assert that an imported GFF is a complete successful run; otherwise absent queries remain
                        not_evaluated.
  --miniprot MINIPROT   miniprot executable name or path.
  --threads THREADS
  --review-limit REVIEW_LIMIT
                        Maximum review appendix rows; the full ledger is always retained.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench assay-plan

```text
usage: candidate-qtl-bench assay-plan [-h] --candidates CANDIDATES --out OUT [--top-n TOP_N] [--genome GENOME]
                                      [--gff3 GFF3] [--db DB] [--primerblast-root PRIMERBLAST_ROOT] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Ranked candidate TSV.
  --out OUT             Output assay-plan TSV.
  --top-n TOP_N
  --genome GENOME       Reference FASTA for primerblast-oss assay commands.
  --gff3 GFF3           Reference GFF3 for primerblast-oss gene assays.
  --db DB               BLAST nucleotide DB; repeat for multi-reference screening.
  --primerblast-root PRIMERBLAST_ROOT
                        primerblast-oss checkout root.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench assembly-bridge

```text
usage: candidate-qtl-bench assembly-bridge [-h] --functional-evidence FUNCTIONAL_EVIDENCE --source-fam3 SOURCE_FAM3
                                           --target-fam3 TARGET_FAM3 --sequence-evidence SEQUENCE_EVIDENCE
                                           [--pangene-hash PANGENE_HASH] [--pangene-dataset PANGENE_DATASET] --out-dir
                                           OUT_DIR [--source-namespace SOURCE_NAMESPACE]
                                           [--target-namespace TARGET_NAMESPACE]
                                           [--mapping-scope {same_species_cross_assembly,same_accession_annotation_generation,cross_accession}]
                                           [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --functional-evidence FUNCTIONAL_EVIDENCE, --source-functional-evidence FUNCTIONAL_EVIDENCE, --source-functions FUNCTIONAL_EVIDENCE
                        normalized_functional_evidence.tsv from the native resource snapshot.
  --source-fam3 SOURCE_FAM3, --source-family-assignments SOURCE_FAM3
                        Source native-ID fam3 assignment TSV (native_gene_id/gene + family_id/fam3).
  --target-fam3 TARGET_FAM3, --target-family-assignments TARGET_FAM3
                        Target native-ID fam3 assignment TSV (native_gene_id/gene + family_id/fam3).
  --sequence-evidence SEQUENCE_EVIDENCE, --pairwise-sequence SEQUENCE_EVIDENCE, --sequences SEQUENCE_EVIDENCE
                        Explicit source_full_id/target_full_id pairwise sequence TSV; reciprocal state is required for
                        verification.
  --pangene-hash PANGENE_HASH
                        Optional frozen LIS pangene *.hsh.tsv[.gz] table (set ID first, gene second; named
                        gene/representative headers also accepted). Membership only proposes full-ID candidate
                        relations; it never verifies one-copy identity.
  --pangene-dataset PANGENE_DATASET
                        Exact pangene dataset generation recorded with --pangene-hash, e.g. Medicago.pan2.451K.
                        Required when --pangene-hash is supplied.
  --out-dir OUT_DIR
  --source-namespace SOURCE_NAMESPACE
                        Optional declared source namespace; otherwise inferred from each full ID.
  --target-namespace TARGET_NAMESPACE
                        Optional declared target namespace; otherwise inferred from each full ID.
  --mapping-scope {same_species_cross_assembly,same_accession_annotation_generation,cross_accession}
                        Mapping scope recorded in every mapping row and the manifest. cross_accession is for explicit
                        same-species accession-to-accession bridges such as current Medicago A17 -> R108; it does not
                        weaken identity verification.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench auto-triage

```text
usage: candidate-qtl-bench auto-triage [-h] --phenotype PHENOTYPE --resource-panel RESOURCE_PANEL --out-dir OUT_DIR
                                       [--candidates CANDIDATES] [--candidate-proteins CANDIDATE_PROTEINS]
                                       [--id-map ID_MAP] [--fasta FASTA] [--gff3 GFF3] [--interval INTERVAL]
                                       [--base-evidence BASE_EVIDENCE] [--parent-vcf PARENT_VCF]
                                       [--trait-term TRAIT_TERM]
                                       [--llm-provider {openai-responses,openai-compatible,lmstudio-cli}]
                                       [--model MODEL] [--endpoint ENDPOINT] [--allow-external-llm]
                                       [--runner {auto,native,wsl}] [--threads THREADS] [--target-size TARGET_SIZE]
                                       [--resume] [--json]

options:
  -h, --help            show this help message and exit
  --phenotype PHENOTYPE
  --resource-panel RESOURCE_PANEL
                        YAML with native and role-labelled local reference proteomes.
  --out-dir OUT_DIR
  --candidates CANDIDATES
                        Existing interval candidates TSV; requires --candidate-proteins.
  --candidate-proteins CANDIDATE_PROTEINS
                        Existing interval candidate protein FASTA; requires --candidates.
  --id-map ID_MAP       Optional native gene/transcript/protein/FASTA-record identity table from extract.
  --fasta FASTA         Native genome FASTA; use with --gff3 and --interval when candidates are not pre-extracted.
  --gff3 GFF3
  --interval INTERVAL
  --base-evidence BASE_EVIDENCE
                        Optional canonical evidence TSV to retain alongside automatic screening evidence.
  --parent-vcf PARENT_VCF
                        Optional parental VCF/VCF.GZ; absence does not block the route.
  --trait-term TRAIT_TERM
                        Additional input phenotype phrase; repeatable.
  --llm-provider {openai-responses,openai-compatible,lmstudio-cli}
                        Optional phenotype-query expansion only.
  --model MODEL
  --endpoint ENDPOINT
  --allow-external-llm
  --runner {auto,native,wsl}
  --threads THREADS
  --target-size TARGET_SIZE
  --resume
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench bait-panel-freeze

```text
usage: candidate-qtl-bench bait-panel-freeze [-h] --hypotheses HYPOTHESES --validation VALIDATION --sequences
                                             SEQUENCES [--reviewable-sequences REVIEWABLE_SEQUENCES] --out-dir OUT_DIR
                                             [--case-id CASE_ID] [--min-panel-size MIN_PANEL_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --hypotheses HYPOTHESES
  --validation VALIDATION
                        proposal_validation.tsv from hypothesis-verify.
  --sequences SEQUENCES
                        verified_sequences.json from hypothesis-verify.
  --reviewable-sequences REVIEWABLE_SEQUENCES
                        Inclusive identity/species/sequence-confirmed sequence store from hypothesis-verify, for
                        review-only projection.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --min-panel-size MIN_PANEL_SIZE
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench benchmark-reference

```text
usage: candidate-qtl-bench benchmark-reference [-h] --config CONFIG --fasta FASTA --gff3 GFF3 --protein-fasta
                                               PROTEIN_FASTA --eggnog EGGNOG --out-dir OUT_DIR [--seqid-map SEQID_MAP]
                                               [--annotation-source-id ANNOTATION_SOURCE_ID]
                                               [--interval-size INTERVAL_SIZE] [--homolog-proteins HOMOLOG_PROTEINS]
                                               [--homolog-map HOMOLOG_MAP]
                                               [--query-eligibility-manifest QUERY_ELIGIBILITY_MANIFEST]
                                               [--blastp-runner {auto,native,wsl}] [--num-threads NUM_THREADS]
                                               [--resume] [--json]

options:
  -h, --help            show this help message and exit
  --config CONFIG       TSV with case_id, trait, interval, answer_gene_id, trait_terms, and answer_source_id.
  --fasta FASTA
  --gff3 GFF3
  --protein-fasta PROTEIN_FASTA
                        Official reference protein FASTA.
  --eggnog EGGNOG       Official eggNOG-mapper annotation table for the protein FASTA.
  --out-dir OUT_DIR
  --seqid-map SEQID_MAP
  --annotation-source-id ANNOTATION_SOURCE_ID
  --interval-size INTERVAL_SIZE
                        Build an answer-centred virtual interval of this size when the config interval is blank;
                        position is not scored.
  --homolog-proteins HOMOLOG_PROTEINS
                        FASTA containing curated homolog queries; repeatable.
  --homolog-map HOMOLOG_MAP
                        TSV keyed by case_id with homolog_query_id and source provenance.
  --query-eligibility-manifest QUERY_ELIGIBILITY_MANIFEST
                        Optional query_eligibility_manifest.json to enforce fail-closed query validation.
  --blastp-runner {auto,native,wsl}
  --num-threads NUM_THREADS
  --resume              Reuse complete shared interval files and BLAST tables in out-dir.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench blast-panel-freeze

```text
usage: candidate-qtl-bench blast-panel-freeze [-h] --panel PANEL --out-dir OUT_DIR [--case-id CASE_ID]
                                              [--min-panel-size MIN_PANEL_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --panel PANEL         Panel TSV in frozen_bait_panel.tsv format.
  --out-dir OUT_DIR     Output directory for frozen_bait_panel.tsv and blast_panel_freeze_manifest.json.
  --case-id CASE_ID
  --min-panel-size MIN_PANEL_SIZE
                        Panels with fewer verified rows than this are review_required_small_panel, not rejected.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench blast-panel-search

```text
usage: candidate-qtl-bench blast-panel-search [-h] --freeze FREEZE --proteins PROTEINS --out-dir OUT_DIR
                                              [--runner {auto,native,wsl}] [--threads THREADS] [--json]

options:
  -h, --help            show this help message and exit
  --freeze FREEZE       blast_panel_freeze_manifest.json; verified before any BLAST tool runs.
  --proteins PROTEINS   Interval protein FASTA; only opened after the freeze verifies.
  --out-dir OUT_DIR     Output directory for blast_evidence.tsv and blast_panel_search_manifest.json.
  --runner {auto,native,wsl}
                        auto: native when BLAST+ is on PATH, otherwise WSL.
  --threads THREADS     blastp -num_threads.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench blast-panel-verify

```text
usage: candidate-qtl-bench blast-panel-verify [-h] --freeze FREEZE [--json]

options:
  -h, --help       show this help message and exit
  --freeze FREEZE  blast_panel_freeze_manifest.json.
  --json           Emit a stable JSON result to stdout.
```

## candidate-qtl-bench budget-llm-review

```text
usage: candidate-qtl-bench budget-llm-review [-h] --out-dir OUT_DIR [--case-id CASE_ID] --responses RESPONSES
                                             [--mechanisms MECHANISMS] [--allowed-sources ALLOWED_SOURCES]
                                             [--bridge BRIDGE]
                                             [--budget-profile {diagnostic,exhaustive,heavy,light,standard,student-budget,student_cpu,zero}]
                                             [--batch-size BATCH_SIZE] [--cache-dir CACHE_DIR] [--no-cache]
                                             [--portfolio-size PORTFOLIO_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR     Existing budget-run output directory
                        (tier0_inventory.tsv/tier1_clusters.tsv/tier2_plan.tsv/budget_handoff.tsv already present).
  --case-id CASE_ID     Case identifier recorded in every manifest.
  --responses RESPONSES
                        Canned review responses JSON for the offline FileReviewBackend; never a network call.
  --mechanisms MECHANISMS
                        mechanism_hypotheses.tsv from `mechanism-import`; enables the proposer diversity audit and the
                        mechanism-id allowlist.
  --allowed-sources ALLOWED_SOURCES
                        TSV with a source_id column of verified source ids; enables the rationale source-id allowlist.
  --bridge BRIDGE       family_candidate_bridge.tsv (from `family-candidate-bridge`); supplies the route each
                        candidate rode in on so the portfolio reducer can allocate its route-unique rescue slots.
                        Without it the reducer records no_route_information_available instead of assuming a route.
  --budget-profile {diagnostic,exhaustive,heavy,light,standard,student-budget,student_cpu,zero}
                        Budget profile preset (default: student-budget).
  --batch-size BATCH_SIZE
                        Maximum Tier1 clusters packed per review batch.
  --cache-dir CACHE_DIR
                        Content-addressed cache directory; default <out-dir>/budget_cache.
  --no-cache            Disable the content cache (recompute every batch).
  --portfolio-size PORTFOLIO_SIZE
                        Top-N portfolio size for the reducer.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench budget-run

```text
usage: candidate-qtl-bench budget-run [-h] [--case-id CASE_ID] --extract-dir EXTRACT_DIR [--out-dir OUT_DIR]
                                      [--expression EXPRESSION] [--expression_provider EXPRESSION_PROVIDER]
                                      [--domains DOMAINS] [--domains_provider DOMAINS_PROVIDER]
                                      [--orthogroup ORTHOGROUP] [--orthogroup_provider ORTHOGROUP_PROVIDER]
                                      [--family_table FAMILY_TABLE] [--family_provider FAMILY_PROVIDER]
                                      [--variant_evidence VARIANT_EVIDENCE] [--variant_provider VARIANT_PROVIDER]
                                      [--sv_evidence SV_EVIDENCE] [--sv_provider SV_PROVIDER]
                                      [--promoter_evidence PROMOTER_EVIDENCE] [--promoter_provider PROMOTER_PROVIDER]
                                      [--budget-profile {diagnostic,exhaustive,heavy,light,standard,student-budget,student_cpu,zero}]
                                      [--cache-dir CACHE_DIR] [--no-cache] [--protein-db PROTEIN_DB]
                                      [--blastp-bin BLASTP_BIN] [--flag-gene FLAG_GENE] [--exhaustive]
                                      [--bridge BRIDGE] [--portfolio-size PORTFOLIO_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --case-id CASE_ID     Case identifier recorded in every manifest.
  --extract-dir EXTRACT_DIR
                        Extract directory produced by `extract` (candidates.tsv etc.).
  --out-dir OUT_DIR     Output directory; default <extract-dir>/budget_engine.
  --expression EXPRESSION
                        expression provider table TSV.
  --expression_provider EXPRESSION_PROVIDER
                        expression provider label NAME@VERSION.
  --domains DOMAINS     domains provider table TSV.
  --domains_provider DOMAINS_PROVIDER
                        domains provider label NAME@VERSION.
  --orthogroup ORTHOGROUP
                        orthogroup provider table TSV.
  --orthogroup_provider ORTHOGROUP_PROVIDER
                        orthogroup provider label NAME@VERSION.
  --family_table FAMILY_TABLE
                        family provider table TSV.
  --family_provider FAMILY_PROVIDER
                        family provider label NAME@VERSION.
  --variant_evidence VARIANT_EVIDENCE
                        variant provider table TSV.
  --variant_provider VARIANT_PROVIDER
                        variant provider label NAME@VERSION.
  --sv_evidence SV_EVIDENCE
                        sv provider table TSV.
  --sv_provider SV_PROVIDER
                        sv provider label NAME@VERSION.
  --promoter_evidence PROMOTER_EVIDENCE
                        promoter provider table TSV.
  --promoter_provider PROMOTER_PROVIDER
                        promoter provider label NAME@VERSION.
  --budget-profile {diagnostic,exhaustive,heavy,light,standard,student-budget,student_cpu,zero}
                        Budget profile preset.
  --cache-dir CACHE_DIR
                        Content-addressed cache directory; default <out-dir>/budget_cache.
  --no-cache            Disable the content cache (recompute every query).
  --protein-db PROTEIN_DB
                        Proteome database file or directory searched by blastp; absent means provider missing.
  --blastp-bin BLASTP_BIN
                        blastp executable name or path.
  --flag-gene FLAG_GENE
                        Force the Tier 2 expansion trigger for a gene; repeatable.
  --exhaustive          Also run the unlimited exhaustive arm and the gating comparison.
  --bridge BRIDGE       family_candidate_bridge.tsv (from `family-candidate-bridge`); supplies the route each
                        candidate rode in on so the portfolio reducer can allocate its route-unique rescue slots.
                        Without it the reducer records no_route_information_available instead of assuming a route.
  --portfolio-size PORTFOLIO_SIZE
                        Top-N portfolio size for the deterministic/LLM-assisted reducer (issue #375 T5).
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench build-executor-package

```text
usage: candidate-qtl-bench build-executor-package [-h] --out-dir OUT_DIR [--repo-root REPO_ROOT] [--include INCLUDE]
                                                  [--reference-assets REFERENCE_ASSETS] [--run-command RUN_COMMAND]
                                                  [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR     New output directory outside the repository root.
  --repo-root REPO_ROOT
                        Repository root; default: detected from the candidate-qtl pyproject.toml.
  --include INCLUDE     Repo-relative glob to include; repeatable. Default: the built-in allowlist patterns.
  --reference-assets REFERENCE_ASSETS
                        Optional TSV (path_wsl required; resource_key/provider/version/retrieved_at_utc optional) of
                        external reference assets fingerprinted into the capsule.
  --run-command RUN_COMMAND
                        Quoted Executor command line (e.g. "candidate-qtl audit-predict --manifest m.tsv") recorded in
                        the capsule as the declared execution command.
  --full-hash           Force full SHA-256 for large reference assets.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench candidate-budget-gate

```text
usage: candidate-qtl-bench candidate-budget-gate [-h] --candidates CANDIDATES --out-dir OUT_DIR
                                                 [--wet-budget WET_BUDGET] [--budget-policy {target,strict}]
                                                 [--adaptive-round ADAPTIVE_ROUND]
                                                 [--max-adaptive-rounds MAX_ADAPTIVE_ROUNDS] [--case-id CASE_ID]
                                                 [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Supported member TSV with candidate_gene_id or gene_id. Membership is never truncated.
  --out-dir OUT_DIR
  --wet-budget WET_BUDGET
                        Operational review target, not a scientific cutoff.
  --budget-policy {target,strict}
                        Use strict only for an explicitly fixed resource capacity.
  --adaptive-round ADAPTIVE_ROUND
  --max-adaptive-rounds MAX_ADAPTIVE_ROUNDS
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench candidate-elimination

```text
usage: candidate-qtl-bench candidate-elimination [-h] --candidates CANDIDATES --out-dir OUT_DIR [--evidence EVIDENCE]
                                                 [--hypotheses HYPOTHESES] [--wet-state WET_STATE] [--rules RULES]
                                                 [--previous-history PREVIOUS_HISTORY] [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate/ranked/shortlist TSV with gene_id.
  --out-dir OUT_DIR
  --evidence EVIDENCE   Elimination evidence TSV (gene_id/evidence_code[/hypothesis_id/evidence_id/quality/coverage]).
  --hypotheses HYPOTHESES
                        locus_hypotheses.tsv, for mechanism-level state.
  --wet-state WET_STATE
                        candidate_wet_state.tsv, for coverage and quality context.
  --rules RULES         Override the versioned elimination rule registry (templates/elimination_rules.tsv).
  --previous-history PREVIOUS_HISTORY
                        elimination_events.tsv from an earlier run; history is append-only.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench candidate-flags

```text
usage: candidate-qtl-bench candidate-flags [-h] --candidates CANDIDATES --out OUT [--tandem-distance TANDEM_DISTANCE]
                                           [--max-aa MAX_AA] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --out OUT
  --tandem-distance TANDEM_DISTANCE
  --max-aa MAX_AA
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench candidate-pipeline

```text
usage: candidate-qtl-bench candidate-pipeline [-h] --candidates CANDIDATES --out-dir OUT_DIR [--blast BLAST]
                                              [--diamond DIAMOND] [--hmmer HMMER] [--eggnog EGGNOG]
                                              [--expression EXPRESSION] [--target-tissue TARGET_TISSUE]
                                              [--target-stage TARGET_STAGE] [--min-expression MIN_EXPRESSION]
                                              [--target-fold TARGET_FOLD] [--manual-evidence MANUAL_EVIDENCE]
                                              [--refined-evidence REFINED_EVIDENCE] [--peak PEAK]
                                              [--interval INTERVAL] [--peak-window PEAK_WINDOW]
                                              [--boundary-window BOUNDARY_WINDOW] [--weights WEIGHTS] [--top-n TOP_N]
                                              [--reserve-n RESERVE_N] [--rescue-n RESCUE_N]
                                              [--min-top-score MIN_TOP_SCORE] [--protein-fasta PROTEIN_FASTA]
                                              [--target-proteome TARGET_PROTEOME]
                                              [--relative-proteome RELATIVE_PROTEOME] [--num-threads NUM_THREADS]
                                              [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --out-dir OUT_DIR
  --blast BLAST
  --diamond DIAMOND
  --hmmer HMMER
  --eggnog EGGNOG
  --expression EXPRESSION
  --target-tissue TARGET_TISSUE
  --target-stage TARGET_STAGE
  --min-expression MIN_EXPRESSION
  --target-fold TARGET_FOLD
  --manual-evidence MANUAL_EVIDENCE
  --refined-evidence REFINED_EVIDENCE
                        Canonical evidence produced by orthology-import.
  --peak PEAK
  --interval INTERVAL
  --peak-window PEAK_WINDOW
  --boundary-window BOUNDARY_WINDOW
  --weights WEIGHTS
  --top-n TOP_N
  --reserve-n RESERVE_N
  --rescue-n RESCUE_N
  --min-top-score MIN_TOP_SCORE
  --protein-fasta PROTEIN_FASTA
                        Candidate proteins normalized to gene_id.
  --target-proteome TARGET_PROTEOME
  --relative-proteome RELATIVE_PROTEOME
                        LABEL=FASTA; repeatable.
  --num-threads NUM_THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench case-run

```text
usage: candidate-qtl-bench case-run [-h] --case CASE --resources RESOURCES --out-dir OUT_DIR
                                    [--mode {deterministic,llm-assisted}] [--resume] [--strict-fingerprints] [--json]

options:
  -h, --help            show this help message and exit
  --case CASE           Case manifest (YAML). Ranking falls back to top_n=3, reserve_n=3, rescue_n=3,
                        min_top_score=2.0 when the manifest declares none.
  --resources RESOURCES
  --out-dir OUT_DIR
  --mode {deterministic,llm-assisted}
  --resume
  --strict-fingerprints
                        Hash every input by full content and record fingerprint strength.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench claim-verification-benchmark

```text
usage: candidate-qtl-bench claim-verification-benchmark [-h] --cases CASES --out-dir OUT_DIR [--cache CACHE]
                                                        [--cache-is-complete] [--allow-external-api] [--full-hash]
                                                        [--json]

options:
  -h, --help            show this help message and exit
  --cases CASES         TSV: case_id/citation/expected_claim_support_status + optional
                        reference_gene_name/gene_family/pathway_or_process/expected_citation_verified.
  --out-dir OUT_DIR
  --cache CACHE         Directory (or file) of cached Europe PMC records; without it every resolvable case is
                        not_evaluated rather than wrong.
  --cache-is-complete   Assert that the cache holds every record that exists for these citations, so a cache miss is a
                        genuine non-resolution. Required to evaluate fabricated-citation cases offline.
  --allow-external-api
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench collect-sanger-markers

```text
usage: candidate-qtl-bench collect-sanger-markers [-h] --plan PLAN --out OUT [--run] [--json]

options:
  -h, --help   show this help message and exit
  --plan PLAN  exon_sequencing_plan.tsv
  --out OUT    Output candidate_sanger_markers.tsv
  --run        Run the generated PrimerBlast OSS shell script first.
  --json       Emit a stable JSON result to stdout.
```

## candidate-qtl-bench comparative-locus-benchmark

```text
usage: candidate-qtl-bench comparative-locus-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR [--full-hash]
                                                       [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  TSV: case_id/degradation_type/truth_orthogroup/reconstructed_universe_tsv/reconstruction_manife
                       st_json.
  --out-dir OUT_DIR
  --full-hash
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench comparative-locus-reconstruct

```text
usage: candidate-qtl-bench comparative-locus-reconstruct [-h] --markers MARKERS --reference-panel REFERENCE_PANEL
                                                         --marker-projections MARKER_PROJECTIONS
                                                         --reference-locus-genes REFERENCE_LOCUS_GENES
                                                         [--focal-recovery FOCAL_RECOVERY]
                                                         [--native-genes NATIVE_GENES]
                                                         [--absence-evidence ABSENCE_EVIDENCE]
                                                         [--min-reference-species MIN_REFERENCE_SPECIES]
                                                         [--core-support-fraction CORE_SUPPORT_FRACTION]
                                                         [--case-id CASE_ID] --out-dir OUT_DIR [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --markers MARKERS     FASTA with at least two flanking marker sequences; native chromosome coordinates are not
                        required.
  --reference-panel REFERENCE_PANEL
                        TSV: species/role/available/cache_version (+ optional phylogenetic_distance and cache paths).
  --marker-projections MARKER_PROJECTIONS
                        Precomputed marker hits: marker_id/species/seqid/start (+ optional orientation/hit_count).
  --reference-locus-genes REFERENCE_LOCUS_GENES
                        Cached projected-block genes: species/gene_id/orthogroup/order_index.
  --focal-recovery FOCAL_RECOVERY
                        Precomputed focal assembly/scaffold/transcriptome hits: orthogroup/focal_gene_id/source_class.
  --native-genes NATIVE_GENES
                        Optional native annotation genes: gene_id/orthogroup; retained as the primary candidate class.
  --absence-evidence ABSENCE_EVIDENCE
                        Optional callable absence evidence: orthogroup/evidence_status/source_id.
  --min-reference-species MIN_REFERENCE_SPECIES
  --core-support-fraction CORE_SUPPORT_FRACTION
  --case-id CASE_ID
  --out-dir OUT_DIR
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench corpus-mask

```text
usage: candidate-qtl-bench corpus-mask [-h] --corpus CORPUS --mask-spec MASK_SPEC --tier
                                       {leave_gene_out,leave_family_out,leave_close_orthogroup_out} --out-dir OUT_DIR
                                       [--case-id CASE_ID] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --corpus CORPUS       Frozen reference corpus TSV from reference-corpus-build (or a prior corpus-mask output, for
                        stacking tiers).
  --mask-spec MASK_SPEC
                        TSV: case_id/truth_gene_id (+ optional truth_family_id, truth_orthogroup_id, aliases),
                        generated from truth in the benchmark-construction phase. May hold every case; --case-id
                        selects one row.
  --tier {leave_gene_out,leave_family_out,leave_close_orthogroup_out}
                        leave_gene_out removes the exact reference gene; leave_family_out also removes its
                        family/subfamily and clearly-equivalent aliases (requires truth_family_id);
                        leave_close_orthogroup_out additionally removes the whole orthogroup (requires
                        truth_orthogroup_id too). Each tier is a superset of the one before it.
  --out-dir OUT_DIR
  --case-id CASE_ID     Select one row of --mask-spec; optional when the spec has exactly one row.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench cosegregation

```text
usage: candidate-qtl-bench cosegregation [-h] --input INPUT --out OUT [--genotype-column GENOTYPE_COLUMN]
                                         [--phenotype-column PHENOTYPE_COLUMN] [--json]

options:
  -h, --help            show this help message and exit
  --input INPUT         TSV based on templates/f2_genotype_phenotype.tsv.
  --out OUT             Output summary TSV.
  --genotype-column GENOTYPE_COLUMN
  --phenotype-column PHENOTYPE_COLUMN
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench dark-candidates

```text
usage: candidate-qtl-bench dark-candidates [-h] --candidates CANDIDATES --out OUT [--max-aa MAX_AA] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --out OUT
  --max-aa MAX_AA
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench direction-evidence

```text
usage: candidate-qtl-bench direction-evidence [-h] --evidence EVIDENCE --expected {increase,decrease,gain,loss} --out
                                              OUT [--json]

options:
  -h, --help            show this help message and exit
  --evidence EVIDENCE
  --expected {increase,decrease,gain,loss}
  --out OUT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench discover

```text
usage: candidate-qtl-bench discover [-h] --fasta FASTA --gff3 GFF3 --interval INTERVAL --phenotype PHENOTYPE
                                    [--phenotype-term PHENOTYPE_TERM] [--approved-terms APPROVED_TERMS]
                                    [--direction {increase,decrease,gain,loss,unknown}]
                                    [--inheritance {recessive,dominant,semidominant,unknown}] [--evidence EVIDENCE]
                                    [--protein-db PROTEIN_DB] [--no-auto-local-db] [--offline | --allow-external-api]
                                    [--dry-run] [--external-provider {ncbi-blast}]
                                    [--external-max-candidates EXTERNAL_MAX_CANDIDATES] [--email EMAIL]
                                    [--top-n TOP_N] [--threads THREADS] [--flank FLANK] [--seqid-map SEQID_MAP]
                                    [--gff-qc-report GFF_QC_REPORT] [--qtl-source-manifest QTL_SOURCE_MANIFEST]
                                    [--full-hash] --out-dir OUT_DIR [--json]

options:
  -h, --help            show this help message and exit
  --fasta FASTA
  --gff3 GFF3
  --interval INTERVAL   CHROM:START-END; repeat for split blocks.
  --phenotype PHENOTYPE
  --phenotype-term PHENOTYPE_TERM
  --approved-terms APPROVED_TERMS
                        Approved trait-term TSV from `trait-terms build`; replaces auto token splitting.
  --direction {increase,decrease,gain,loss,unknown}
  --inheritance {recessive,dominant,semidominant,unknown}
                        Expected segregation model. It affects evidence only when a source-backed row supplies
                        inheritance_model and molecular_effect.
  --evidence EVIDENCE   Optional canonical source-backed evidence TSV; repeatable. Rows without source_id are
                        quarantined.
  --protein-db PROTEIN_DB
                        Preferred local protein BLAST database.
  --no-auto-local-db    Disable discovery of the configured local protein DB.
  --offline             Prohibit all external sequence submission; use GFF3/local databases only.
  --allow-external-api  Allow candidate protein sequences to be transmitted to the external NCBI BLAST API when no
                        local protein database is selected.
  --dry-run             Plan only: report whether a real run would transmit candidate sequences externally and write
                        the plan to the manifest, without running BLAST or submitting anything.
  --external-provider {ncbi-blast}
  --external-max-candidates EXTERNAL_MAX_CANDIDATES
  --email EMAIL
  --top-n TOP_N         Dry ranking batch size; also caps Reserve and Rescue slots (default: 10). The Dry score gate
                        is min_top_score=1.0.
  --threads THREADS
  --flank FLANK
  --seqid-map SEQID_MAP
                        GFF_SEQID=FASTA_SEQID pairs or a two-column TSV, for references whose GFF and FASTA seqids
                        differ.
  --gff-qc-report GFF_QC_REPORT
                        Optional GenomeTools/AGAT validator report to record in the manifest.
  --qtl-source-manifest QTL_SOURCE_MANIFEST
                        qtl_source_preflight.json coordinate-provenance gate (same gate as extract).
  --full-hash           SHA-256 all inputs, including multi-gigabyte FASTA files.
  --out-dir OUT_DIR
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench discover-promote

```text
usage: candidate-qtl-bench discover-promote [-h] --discover-dir DISCOVER_DIR --out-dir OUT_DIR [--case-id CASE_ID]
                                            [--json]

options:
  -h, --help            show this help message and exit
  --discover-dir DISCOVER_DIR
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench doctor

```text
usage: candidate-qtl-bench doctor [-h] [--json]

options:
  -h, --help  show this help message and exit
  --json      Emit a stable JSON result to stdout.
```

## candidate-qtl-bench dry-prioritize

```text
usage: candidate-qtl-bench dry-prioritize [-h] --candidates CANDIDATES --out OUT [--top-n TOP_N]
                                          [--reserve-n RESERVE_N] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        TSV based on templates/dry_candidate_prioritization.tsv.
  --out OUT             Output TSV with Dry score, group, rank, and score breakdown.
  --top-n TOP_N         Number of non-rescue candidates in the first sequencing batch.
  --reserve-n RESERVE_N
                        Number of non-rescue candidates retained as the next batch.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench dual-route-compare

```text
usage: candidate-qtl-bench dual-route-compare [-h] --route-a ROUTE_A --route-b ROUTE_B --out-dir OUT_DIR
                                              [--case-id CASE_ID] [--route-b-manifest ROUTE_B_MANIFEST]
                                              [--benchmark-mode] [--json]

options:
  -h, --help            show this help message and exit
  --route-a ROUTE_A     phenotype_first_hypotheses.tsv
  --route-b ROUTE_B     family_process_hypotheses.tsv
  --out-dir OUT_DIR
  --case-id CASE_ID
  --route-b-manifest ROUTE_B_MANIFEST
                        interval_first_status.json from interval-reverse-lookup; records whether Route B saw the
                        phenotype (B0 phenotype-blind = false, B1 phenotype overlay = true).
  --benchmark-mode      Mark the artifact as a benchmark-mode comparison (B0 phenotype-blind Route B).
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench evidence-build

```text
usage: candidate-qtl-bench evidence-build [-h] --candidates CANDIDATES --out-dir OUT_DIR [--blast BLAST]
                                          [--diamond DIAMOND] [--hmmer HMMER] [--eggnog EGGNOG]
                                          [--expression EXPRESSION] [--target-tissue TARGET_TISSUE]
                                          [--target-stage TARGET_STAGE] [--target-condition TARGET_CONDITION]
                                          [--min-expression MIN_EXPRESSION] [--target-fold TARGET_FOLD]
                                          [--manual-evidence MANUAL_EVIDENCE] [--existing-evidence EXISTING_EVIDENCE]
                                          [--peak PEAK] [--interval INTERVAL] [--peak-window PEAK_WINDOW]
                                          [--boundary-window BOUNDARY_WINDOW] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Transcript-aware candidates.tsv from extract.
  --out-dir OUT_DIR
  --blast BLAST         BLASTP outfmt 6 file; repeatable.
  --diamond DIAMOND     DIAMOND outfmt 6 file; repeatable.
  --hmmer HMMER         HMMER --tblout file; repeatable.
  --eggnog EGGNOG       eggNOG-mapper annotations file; repeatable.
  --expression EXPRESSION
                        Long-form TSV with gene_id, tissue, value, and optional is_target.
  --target-tissue TARGET_TISSUE
                        Target tissue name; repeatable unless expression TSV has is_target.
  --target-stage TARGET_STAGE
                        Target stage; combined with --target-tissue unless is_target is supplied.
  --target-condition TARGET_CONDITION
                        Restrict target expression rows to these condition values; repeatable.
  --min-expression MIN_EXPRESSION
  --target-fold TARGET_FOLD
  --manual-evidence MANUAL_EVIDENCE
                        Source-backed TSV using the canonical evidence schema.
  --existing-evidence EXISTING_EVIDENCE
                        Existing canonical evidence TSV, for example orthology refinement results.
  --peak PEAK           QTL peak as CHROM:POS or POS.
  --interval INTERVAL   QTL interval; inferred from candidate table when unique.
  --peak-window PEAK_WINDOW
                        Analyst-supplied peak uncertainty half-width for scored peak_proximal support; 0 scores exact
                        overlap only and reports the 5%-interval heuristic as unscored context.
  --boundary-window BOUNDARY_WINDOW
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench evidence-rank

```text
usage: candidate-qtl-bench evidence-rank [-h] --candidates CANDIDATES --evidence EVIDENCE --out-dir OUT_DIR
                                         [--weights WEIGHTS] [--top-n TOP_N] [--reserve-n RESERVE_N]
                                         [--rescue-n RESCUE_N] [--min-top-score MIN_TOP_SCORE] [--case-id CASE_ID]
                                         [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --evidence EVIDENCE
  --out-dir OUT_DIR
  --weights WEIGHTS     Optional TSV with evidence_code and weight.
  --top-n TOP_N
  --reserve-n RESERVE_N
  --rescue-n RESCUE_N   Maximum Rescue-group candidates retained.
  --min-top-score MIN_TOP_SCORE
  --case-id CASE_ID     Optional case identifier recorded in the decision artifact.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench evidence-term-promote

```text
usage: candidate-qtl-bench evidence-term-promote [-h] --evidence EVIDENCE --out-dir OUT_DIR --term TERM
                                                 [--demote-unapproved-existing] [--json]

options:
  -h, --help            show this help message and exit
  --evidence EVIDENCE
  --out-dir OUT_DIR
  --term TERM           Approved exact token/phrase; repeatable.
  --demote-unapproved-existing
                        Convert pre-existing trait matches outside the approved vocabulary to unscored context.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench evidence-tier-benchmark

```text
usage: candidate-qtl-bench evidence-tier-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR
                                                   [--review-budget REVIEW_BUDGET] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST   TSV: one row per case (case_id/ranked_tsv[/answer_gene_id/knowledge_coverage_tsv]).
  --out-dir OUT_DIR
  --review-budget REVIEW_BUDGET
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench evidence-tiers

```text
usage: candidate-qtl-bench evidence-tiers [-h] --ranked RANKED --out-dir OUT_DIR [--review-budget REVIEW_BUDGET]
                                          [--reserve-size RESERVE_SIZE] [--knowledge-coverage KNOWLEDGE_COVERAGE]
                                          [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED       ranked_candidates.tsv from evidence-rank or candidate-pipeline (needs top_eligibility).
  --out-dir OUT_DIR
  --review-budget REVIEW_BUDGET
                        Assay-capacity cut (same meaning/default as screening-shortlist's --target-size); changes
                        operational_state only, never evidence_tier.
  --reserve-size RESERVE_SIZE
                        Extra rank window right after review_budget marked operational_state=reserve; 0 disables it.
  --knowledge-coverage KNOWLEDGE_COVERAGE
                        Optional knowledge_coverage.tsv from knowledge-coverage; enables the knowledge_limited tier.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench exon-coverage-check

```text
usage: candidate-qtl-bench exon-coverage-check [-h] --plan PLAN --out-dir OUT_DIR [--relax-flank RELAX_FLANK]
                                               [--relax-product-extra RELAX_PRODUCT_EXTRA] [--json]

options:
  -h, --help            show this help message and exit
  --plan PLAN           exon_sequencing_plan.tsv from exon-sequencing-plan.
  --out-dir OUT_DIR
  --relax-flank RELAX_FLANK
                        Suggested wider flank for retry windows.
  --relax-product-extra RELAX_PRODUCT_EXTRA
                        Extra product-size span suggested for retries.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench exon-sequencing-plan

```text
usage: candidate-qtl-bench exon-sequencing-plan [-h] --candidates CANDIDATES --gff3 GFF3 --out-dir OUT_DIR
                                                --primerblast-root PRIMERBLAST_ROOT --genome GENOME --db DB
                                                [--num-threads NUM_THREADS] [--top-n TOP_N] [--flank FLANK]
                                                [--max-target-span MAX_TARGET_SPAN] [--overlap OVERLAP]
                                                [--seqid-map SEQID_MAP] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --gff3 GFF3
  --out-dir OUT_DIR
  --primerblast-root PRIMERBLAST_ROOT
                        primerblast-oss checkout root.
  --genome GENOME       Reference FASTA for primerblast-oss.
  --db DB               BLAST nucleotide DB; repeat for multi-reference screening.
  --num-threads NUM_THREADS
                        Laptop-safe worker count for BLAST-backed primer checks.
  --top-n TOP_N
  --flank FLANK
  --max-target-span MAX_TARGET_SPAN
  --overlap OVERLAP
  --seqid-map SEQID_MAP
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench external-rank

```text
usage: candidate-qtl-bench external-rank [-h] --native-rank NATIVE_RANK --external-ranks EXTERNAL_RANKS --out-dir
                                         OUT_DIR [--top-k TOP_K] [--rrf-k RRF_K]
                                         [--discordance-spread DISCORDANCE_SPREAD] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --native-rank NATIVE_RANK
                        Native ranked_candidates.tsv (source of the candidate universe and native rank).
  --external-ranks EXTERNAL_RANKS
                        Long-format TSV: one row per (ranker_id, gene_id) with coverage.
  --out-dir OUT_DIR
  --top-k TOP_K         Top-K threshold used to compare native and consensus prioritization.
  --rrf-k RRF_K         Reciprocal-rank-fusion constant.
  --discordance-spread DISCORDANCE_SPREAD
                        Normalized cross-family rank spread that flags disagreement.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench extract

```text
usage: candidate-qtl-bench extract [-h] --fasta FASTA --gff3 GFF3 --interval INTERVAL --out-dir OUT_DIR
                                   [--seqid-map SEQID_MAP] [--full-hash] [--gff-qc-report GFF_QC_REPORT]
                                   [--flank FLANK] [--qtl-source-manifest QTL_SOURCE_MANIFEST] [--json]

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Reference FASTA; Windows paths are accepted through the PowerShell wrapper.
  --gff3 GFF3           Reference GFF3 or GFF3.GZ.
  --interval INTERVAL   CHROM:START-END; repeat for split/non-collinear blocks.
  --out-dir OUT_DIR     Output directory.
  --seqid-map SEQID_MAP
                        Optional GFF_SEQID=FASTA_SEQID pairs (comma-separated) or a two-column TSV map.
  --full-hash           Read and SHA-256 all inputs, including multi-gigabyte FASTA files.
  --gff-qc-report GFF_QC_REPORT
                        Optional GenomeTools/AGAT validator report to record in the manifest.
  --flank FLANK         Add gene flanks to candidate_genomic.fa for primer design.
  --qtl-source-manifest QTL_SOURCE_MANIFEST
                        qtl_source_preflight.json provenance gate.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench family-candidate-bridge

```text
usage: candidate-qtl-bench family-candidate-bridge [-h] [--route-a-hypotheses ROUTE_A_HYPOTHESES]
                                                   [--route-a-family-status ROUTE_A_FAMILY_STATUS]
                                                   [--route-a-bait-panel ROUTE_A_BAIT_PANEL]
                                                   [--route-a-projection-panel ROUTE_A_PROJECTION_PANEL]
                                                   [--route-a-blast ROUTE_A_BLAST]
                                                   [--route-a-blast-hits ROUTE_A_BLAST_HITS]
                                                   [--route-b-orthology ROUTE_B_ORTHOLOGY]
                                                   [--route-b-atlas ROUTE_B_ATLAS]
                                                   [--route-b-family-index ROUTE_B_FAMILY_INDEX]
                                                   [--route-c-retrieval ROUTE_C_RETRIEVAL]
                                                   [--route-c-match ROUTE_C_MATCH]
                                                   [--dependency-clusters DEPENDENCY_CLUSTERS] [--genes GENES]
                                                   --out-dir OUT_DIR --case-id CASE_ID [--json]

options:
  -h, --help            show this help message and exit
  --route-a-hypotheses ROUTE_A_HYPOTHESES
                        phenotype_first_hypotheses.tsv (hypothesis_id +
                        member_interval_genes[/gene_family/proposed_stable_id]).
  --route-a-family-status ROUTE_A_FAMILY_STATUS
                        route_a_interval_family_match.tsv (hypothesis_id +
                        family_status[/interval_gene_ids/source_ids]); per-hypothesis family status from the Route A
                        interval match stage.
  --route-a-bait-panel ROUTE_A_BAIT_PANEL
                        verified_bait_panel.tsv from bait-panel-freeze (hypothesis_id +
                        function_source_ids[/source_provider]); hypotheses outside the panel are unverified bait,
                        never a match.
  --route-a-projection-panel ROUTE_A_PROJECTION_PANEL
                        projection_bait_panel.tsv for identity-confirmed review-only BLAST members; never contributes
                        verified support.
  --route-a-blast ROUTE_A_BLAST
                        blast_evidence.tsv from blast-panel-search (panel_id + top_hit_identity/top_hit_query_coverage
                        /top_hit_target_coverage/top_hit_evalue/paralog_count/ambiguity_status/provider_status);
                        joined to Route A matched rows by reference_accession=panel_id (the verified bait accession).
                        A blast hit is similarity evidence only for an existing matched family row; an eligible hit
                        may create one sequence-match row only when family-index status explicitly falsifies that
                        family (#530).
  --route-a-blast-hits ROUTE_A_BLAST_HITS, --route-a-blast-eligible ROUTE_A_BLAST_HITS
                        eligible_interval_hits.tsv from blast-panel-search; all long-form threshold-passing interval
                        hits consumed by Route A family/copy review (#532).
  --route-b-orthology ROUTE_B_ORTHOLOGY
                        Orthology table (gene_id/reference_gene_id +
                        relationship/orthogroup/percent_identity/query_coverage/evalue/provider[/source_id]);
                        best_hit_homology is a best hit, not a reciprocal one.
  --route-b-atlas ROUTE_B_ATLAS
                        interval_gene_atlas.tsv (gene_id + evolutionary_family_id/orthogroup/domain_architecture/domai
                        n_ids/domain_source/family_id/family_label/focal_paralog_count/tandem_expansion); also emits
                        one row per atlas evidence axis.
  --route-b-family-index ROUTE_B_FAMILY_INDEX
                        Route B0 family_index.tsv. An interval_candidates index is consumed only for positive Route A
                        family projection after both lanes freeze: every matching interval copy is retained; no-hit
                        never proves absence.
  --route-c-retrieval ROUTE_C_RETRIEVAL
                        route_c_retrieval_merged.tsv (reference_gene_id/rank + database_version/model_version); joined
                        to the interval match on (reference_gene_id, rank).
  --route-c-match ROUTE_C_MATCH
                        route_c_interval_match.tsv (reference_gene_id/rank/matched_interval_gene +
                        match_basis/family_id/orthogroup_id); rows without any family id are skipped and counted.
  --dependency-clusters DEPENDENCY_CLUSTERS
                        evidence_dependency_clusters.tsv
                        (dependency_cluster/dependency_status/root_source_id/source_ids); overrides the route-default
                        dependency roots.
  --genes GENES         One-column gene_id TSV; genes with no bridge row are reported in the manifest.
  --out-dir OUT_DIR
  --case-id CASE_ID     Non-empty case identifier recorded in every bridge row and the manifest.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench family-index

```text
usage: candidate-qtl-bench family-index [-h] --blast-hits BLAST_HITS --out-dir OUT_DIR [--genes GENES]
                                        [--id-map ID_MAP] [--family-map FAMILY_MAP]
                                        [--reference-label REFERENCE_LABEL] [--reference-species REFERENCE_SPECIES]
                                        [--family-hits FAMILY_HITS] [--min-identity MIN_IDENTITY]
                                        [--min-query-coverage MIN_QUERY_COVERAGE] [--max-evalue MAX_EVALUE]
                                        [--case-id CASE_ID] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --blast-hits BLAST_HITS
                        Reference proteome vs model/relative proteome, 12/13-column outfmt 6.
  --out-dir OUT_DIR
  --genes GENES         Gene coordinate table (gene_id/chrom/start/end); genes without a qualifying hit are kept as an
                        explicit blind spot.
  --id-map ID_MAP       Protein/transcript to gene id map from extract.
  --family-map FAMILY_MAP
                        Optional curated family labels (subject_id/family_id[/symbol]); overrides the symbol-stem
                        heuristic.
  --reference-label REFERENCE_LABEL
                        Focal reference genome label recorded in the manifest.
  --reference-species REFERENCE_SPECIES
  --family-hits FAMILY_HITS
                        Top hits kept per gene for the coarse family bucket.
  --min-identity MIN_IDENTITY
  --min-query-coverage MIN_QUERY_COVERAGE
  --max-evalue MAX_EVALUE
  --case-id CASE_ID     Optional case identifier recorded in the manifest.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench family-to-copy

```text
usage: candidate-qtl-bench family-to-copy [-h] --families FAMILIES --out-dir OUT_DIR [--atlas ATLAS]
                                          [--orthology ORTHOLOGY] [--synteny SYNTENY] [--expression EXPRESSION]
                                          [--position POSITION] [--subfamily-map SUBFAMILY_MAP]
                                          [--family-target FAMILY_TARGET] [--bridge BRIDGE]
                                          [--route-a-blast ROUTE_A_BLAST] [--route-a-blast-hits ROUTE_A_BLAST_HITS]
                                          [--tandem-distance TANDEM_DISTANCE] [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --families FAMILIES   family_process_hypotheses.tsv (or any TSV with
                        hypothesis_id/gene_family/member_interval_genes).
  --out-dir OUT_DIR
  --atlas ATLAS         interval_gene_atlas.tsv, for coordinates and domain architecture (canonical <out-
                        dir>/interval_gene_atlas.tsv from route-b-prepare or the root-level interval-atlas file).
  --orthology ORTHOLOGY
                        Orthology table
                        (gene_id[/relationship/provider/evidence_chain/reference_gene_id/orthogroup/paralog_count]).
  --synteny SYNTENY     Synteny table (gene_id/synteny_status[/provider]).
  --expression EXPRESSION
                        Expression table (gene_id/relevant_expression_status); unmeasured is never a negative.
  --position POSITION   Positional evidence table (gene_id/position_class); broad interval membership does not
                        discriminate.
  --subfamily-map SUBFAMILY_MAP
                        Curated subfamily labels (subfamily + gene_id/reference_gene_id/family_id); otherwise an
                        orthogroup is used, or the answer stays subfamily_unresolved.
  --family-target FAMILY_TARGET
                        Verified target specification (target_id + optional hypothesis_id/expected_family_id/expected_
                        orthogroup/reference_gene_id/expected_subfamily/expected_domain_architecture/expected_synteny_
                        anchor/expected_expression_context/source_ids/verification_status). Copy selection is scored
                        against this target; without it the family-internal modal subfamily/domain are descriptive
                        statistics and outlier QC only and never pick the preferred copy. verification_status other
                        than verified/source_backed/curated is treated as pending_review and can never be the sole
                        basis of a preferred copy.
  --bridge BRIDGE       family_candidate_bridge.tsv (canonical 22-column bridge, e.g. from the family-candidate-bridge
                        command). Looked up per (gene_id, family_id) pair, never by gene_id alone, so a gene's bridge
                        row under a different family never leaks into this family's gate. Supplies mapping tier /
                        dependency cluster / ambiguity / route support: members with a tier5/tier6 bridge tier are
                        never promoted into the shortlist (recorded as
                        bridge_promotion_gate=blocked_tier5/blocked_tier6, disambiguation status unchanged), and
                        family-level route rescue (route_rescue_state:
                        rescue_a_only/rescue_b_only/rescue_c_only/multi_route_supported/not_evaluated, the last when
                        every bridge row a family has is tier6_not_evaluated) and bridge-backed copy-selection loss
                        (copy_selection_loss) are recorded as separate columns.
  --route-a-blast ROUTE_A_BLAST
                        blast_evidence.tsv from blast-panel-search. Only evaluable verified rows with
                        ambiguity_status=unique and threshold-passing top hits become a route_a_blastp homology layer
                        for the hit member gene; provenance is recorded in the discrimination column
                        copy_set_evidence_source (present only when this input is supplied).
  --route-a-blast-hits ROUTE_A_BLAST_HITS, --route-a-blast-eligible ROUTE_A_BLAST_HITS
                        eligible_interval_hits.tsv from blast-panel-search; all threshold-passing interval targets
                        retained for copy review (#532).
  --tandem-distance TANDEM_DISTANCE
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench focal-paralog-import

```text
usage: candidate-qtl-bench focal-paralog-import [-h] --focal-blast FOCAL_BLAST --candidates CANDIDATES --out OUT
                                                [--close-ratio CLOSE_RATIO] [--json]

options:
  -h, --help            show this help message and exit
  --focal-blast FOCAL_BLAST
                        Candidate-vs-focal-proteome BLAST/DIAMOND outfmt 6 result.
  --candidates CANDIDATES
                        Candidate table defining the QTL-interval gene set.
  --out OUT             Canonical evidence TSV output.
  --close-ratio CLOSE_RATIO
                        Paralog bitscore fraction of the self hit to count as close.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench gff-preflight

```text
usage: candidate-qtl-bench gff-preflight [-h] --fasta FASTA --gff3 GFF3 --interval INTERVAL --out-dir OUT_DIR
                                         [--seqid-map SEQID_MAP] [--normalized-gff NORMALIZED_GFF]
                                         [--agat-cds AGAT_CDS] [--gffread-cds GFFREAD_CDS] [--execute] [--json]

options:
  -h, --help            show this help message and exit
  --fasta FASTA
  --gff3 GFF3
  --interval INTERVAL
  --out-dir OUT_DIR
  --seqid-map SEQID_MAP
  --normalized-gff NORMALIZED_GFF
  --agat-cds AGAT_CDS
  --gffread-cds GFFREAD_CDS
  --execute
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench guided-profile

```text
usage: candidate-qtl-bench guided-profile [-h] --profile
                                          {screening,full-case,audit-benchmark,audit-predict,triple-route} --out-dir
                                          OUT_DIR [--execute] [--fasta FASTA] [--gff3 GFF3] [--interval INTERVAL]
                                          [--phenotype PHENOTYPE] [--phenotype-context PHENOTYPE_CONTEXT]
                                          [--allow-external-api] [--allow-reverse-literature-api]
                                          [--protein-db PROTEIN_DB] [--no-auto-local-db]
                                          [--screening-dir SCREENING_DIR] [--case CASE] [--resources RESOURCES]
                                          [--case-id CASE_ID] [--query-manifest QUERY_MANIFEST]
                                          [--query-fasta QUERY_FASTA] [--blind-manifest BLIND_MANIFEST]
                                          [--blind-route BLIND_ROUTE] [--predict-manifest PREDICT_MANIFEST]
                                          [--focal-species FOCAL_SPECIES] [--direction DIRECTION]
                                          [--inheritance INHERITANCE] [--mechanism-payload MECHANISM_PAYLOAD]
                                          [--phenotype-localization-payload PHENOTYPE_LOCALIZATION_PAYLOAD]
                                          [--gap-fill-mechanisms GAP_FILL_MECHANISMS]
                                          [--gap-fill-mechanism-import-manifest GAP_FILL_MECHANISM_IMPORT_MANIFEST]
                                          [--gap-request-audit GAP_REQUEST_AUDIT]
                                          [--mechanism-freeze-manifest MECHANISM_FREEZE_MANIFEST]
                                          [--payload-source {llm,agent,human}] [--model MODEL]
                                          [--record-cache RECORD_CACHE] [--source-cache SOURCE_CACHE]
                                          [--bound-claims BOUND_CLAIMS] [--agent-claim-review AGENT_CLAIM_REVIEW]
                                          [--route-b-agent-claim-review ROUTE_B_AGENT_CLAIM_REVIEW]
                                          [--route-b-source-review-records ROUTE_B_SOURCE_REVIEW_RECORDS]
                                          [--reverse-lit-cache REVERSE_LIT_CACHE] [--shared-route SHARED_ROUTE]
                                          [--retrieval-index RETRIEVAL_INDEX] [--proteins PROTEINS]
                                          [--reference-function REFERENCE_FUNCTION] [--orthology ORTHOLOGY]
                                          [--domains DOMAINS] [--family-index FAMILY_INDEX]
                                          [--annotations ANNOTATIONS] [--copy-orthology COPY_ORTHOLOGY]
                                          [--copy-synteny COPY_SYNTENY] [--copy-expression COPY_EXPRESSION]
                                          [--copy-position COPY_POSITION] [--copy-subfamily-map COPY_SUBFAMILY_MAP]
                                          [--copy-target COPY_TARGET] [--annotation-rescue ANNOTATION_RESCUE]
                                          [--benchmark-mode] [--deterministic-subfamily-refinement]
                                          [--route-c-direct-frozen ROUTE_C_DIRECT_FROZEN]
                                          [--route-c-direct-freeze-manifest ROUTE_C_DIRECT_FREEZE_MANIFEST]
                                          [--route-c-direct-qtl-id ROUTE_C_DIRECT_QTL_ID]
                                          [--route-c-direct-eligible-hits ROUTE_C_DIRECT_ELIGIBLE_HITS]
                                          [--candidates CANDIDATES] [--evidence EVIDENCE] [--causal-gene CAUSAL_GENE]
                                          [--json]

options:
  -h, --help            show this help message and exit
  --profile {screening,full-case,audit-benchmark,audit-predict,triple-route}
  --out-dir OUT_DIR
  --execute             Execute runnable stages. Default is a dry-run manifest.
  --fasta FASTA
  --gff3 GFF3
  --interval INTERVAL
  --phenotype PHENOTYPE
  --phenotype-context PHENOTYPE_CONTEXT
                        Optional phenotype_context.yaml (case_id, developmental_stage, tissue,
                        environment_or_treatment, ontology term ...). Recorded for audit in the Route A/C manifests;
                        the free text is never replaced or inferred, and no structured field is inserted into any
                        query.
  --allow-external-api  Permit profile-specific external providers: NCBI BLAST for screening, and UniProt/Europe PMC
                        verification or retrieval for audit-predict/triple-route. Review the profile manifest's
                        external_transmission fields before execution.
  --allow-reverse-literature-api
                        Triple-route only: separately opt in to external Europe PMC lookups for the legacy diagnostic
                        Route B1 reverse-literature stage. Canonical Route B does not require this.
  --protein-db PROTEIN_DB
                        Screening only: preferred local protein BLAST database (same meaning as `discover --protein-
                        db`).
  --no-auto-local-db    Screening only: disable discovery of the configured local protein DB.
  --screening-dir SCREENING_DIR
  --case CASE
  --resources RESOURCES
  --case-id CASE_ID
  --query-manifest QUERY_MANIFEST
  --query-fasta QUERY_FASTA
  --blind-manifest BLIND_MANIFEST
  --blind-route BLIND_ROUTE
  --predict-manifest PREDICT_MANIFEST
                        audit-benchmark / audit-predict: answerless prediction manifest
                        (case_id/phenotype_text/interval/interval_genes). An answer/truth column is rejected.
  --focal-species FOCAL_SPECIES
                        focal species for the phenotype-investigation and verification stages (audit-predict / triple-
                        route).
  --direction DIRECTION
  --inheritance INHERITANCE
  --mechanism-payload MECHANISM_PAYLOAD
                        audit-predict / triple-route: LLM or hand-authored mechanism payload for mechanism-import
                        (canonical 29/24-column TSV, canonical-subset TSV, JSON/JSONL, or the frozen legacy 7-column
                        TSV, which is converted with an explicit adapter).
  --phenotype-localization-payload PHENOTYPE_LOCALIZATION_PAYLOAD
                        audit-predict / triple-route: phenotype-only localization JSON, created before the mechanism
                        payload and frozen by the profile.
  --gap-fill-mechanisms GAP_FILL_MECHANISMS
                        Optional already-imported blind round-2 mechanism_hypotheses.tsv responding to a prior
                        coverage request.
  --gap-fill-mechanism-import-manifest GAP_FILL_MECHANISM_IMPORT_MANIFEST
                        Required with --gap-fill-mechanisms; matching blind round-2 mechanism_import_manifest.json.
  --gap-request-audit GAP_REQUEST_AUDIT
                        Required with --gap-fill-mechanisms; prior mechanism_coverage_audit.json that requested the
                        exact labels.
  --mechanism-freeze-manifest MECHANISM_FREEZE_MANIFEST
                        TSV with case_id/mechanism_payload_sha256 (coordination/payload_freeze.tsv format); the
                        payload bytes must match the frozen hash or the import fails closed.
  --payload-source {llm,agent,human}
                        Who produced the mechanism payload: an API completion (llm), this agent driving the flow with
                        its own literature tools (agent), or a person (human). All three pass the same gate.
  --model MODEL         Model identifier recorded with an LLM payload.
  --record-cache RECORD_CACHE, --records RECORD_CACHE
                        audit-predict / triple-route: local UniProt-shaped record cache for offline verification.
  --source-cache SOURCE_CACHE
                        triple-route only: cached Europe PMC records for offline source-verify.
  --bound-claims BOUND_CLAIMS
                        audit-predict / triple-route: shared source-bound claim TSV required before Route-A bait
                        promotion.
  --agent-claim-review AGENT_CLAIM_REVIEW
                        triple-route: pre-reviewed independent-agent packet for the exact frozen source/hypothesis
                        inputs; alternative to --bound-claims.
  --route-b-agent-claim-review ROUTE_B_AGENT_CLAIM_REVIEW
                        Route-B-only independent review packet; never reused by Route A.
  --route-b-source-review-records ROUTE_B_SOURCE_REVIEW_RECORDS
                        Route-B-only source snapshot matching the semantic claims.
  --reverse-lit-cache REVERSE_LIT_CACHE
                        triple-route only: query-keyed cache directory for reverse-literature (offline use).
  --shared-route SHARED_ROUTE
                        triple-route only: pre-fixed phenotype-to-pathway route used by the non-scoring B1 overlay.
  --retrieval-index RETRIEVAL_INDEX
                        triple-route only: index_manifest.json from retrieval-index-build (Route C Phase 1).
  --proteins PROTEINS   triple-route only: candidate protein FASTA for Route B homology and, when --family-index is
                        absent, an interval-scoped fallback family index.
  --reference-function REFERENCE_FUNCTION
                        triple-route only: locally cached source-backed reference records for Route B.
  --orthology ORTHOLOGY
  --domains DOMAINS
  --family-index FAMILY_INDEX
                        Preferred reusable genome-wide family_index.tsv. For triple-route, omission triggers lookup in
                        resources.yaml, CANDIDATE_QTL_FAMILY_INDEX, then the per-user cache; Route B may still
                        generate an interval-scoped fallback from --proteins + --protein-db if none is found.
  --annotations ANNOTATIONS
  --copy-orthology COPY_ORTHOLOGY
                        Triple-route canonical copy stage only: explicit native family-to-copy orthology evidence
                        file. Omission is not negative evidence.
  --copy-synteny COPY_SYNTENY
                        Triple-route canonical copy stage only: explicit native family-to-copy synteny evidence file.
                        Omission is not negative evidence.
  --copy-expression COPY_EXPRESSION
                        Triple-route canonical copy stage only: explicit native family-to-copy expression evidence
                        file. Omission is not negative evidence.
  --copy-position COPY_POSITION
                        Triple-route canonical copy stage only: explicit native family-to-copy position evidence file.
                        Omission is not negative evidence.
  --copy-subfamily-map COPY_SUBFAMILY_MAP
                        Triple-route canonical copy stage only: explicit native family-to-copy subfamily_map evidence
                        file. Omission is not negative evidence.
  --copy-target COPY_TARGET
                        Triple-route canonical copy stage only: explicit native family-to-copy target evidence file.
                        Omission is not negative evidence.
  --annotation-rescue ANNOTATION_RESCUE
                        triple-route only: annotation_rescue.tsv from annotation-rescue, passed to the Route A
                        interval family match stage; a rescue provider that did not complete holds absence rather than
                        falsifying it.
  --benchmark-mode      Mark verification as benchmark mode (rejects focal-lineage / answer-backmapped markers).
  --deterministic-subfamily-refinement
                        triple-route --execute only: when the supported member union exceeds the Wet budget, apply one
                        machine-verified curated cross-subfamily refinement pass and re-gate instead of stopping at an
                        adaptive synthesis request.
  --route-c-direct-frozen ROUTE_C_DIRECT_FROZEN
                        Triple-route only: pre-interval frozen Direct-C rescue TSV. Search/freeze must have completed
                        before this profile opens the interval.
  --route-c-direct-freeze-manifest ROUTE_C_DIRECT_FREEZE_MANIFEST
                        Required with --route-c-direct-frozen; binds the interval-blind frozen bytes.
  --route-c-direct-qtl-id ROUTE_C_DIRECT_QTL_ID
                        Required with --route-c-direct-frozen; QTL identifier recorded during post-freeze projection.
  --route-c-direct-eligible-hits ROUTE_C_DIRECT_ELIGIBLE_HITS
                        Optional precomputed eligible BLAST hits used only after Direct-C freeze verification.
  --candidates CANDIDATES
  --evidence EVIDENCE
  --causal-gene CAUSAL_GENE
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-linking-counterfactual

```text
usage: candidate-qtl-bench homolog-linking-counterfactual [-h] --ranked RANKED --proposals PROPOSALS --truth-genes
                                                          TRUTH_GENES [--out-tsv OUT_TSV] [--bonus BONUS] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED       ranked_candidates.tsv from discovery.
  --proposals PROPOSALS
                        homolog_linking_proposals.tsv.
  --truth-genes TRUTH_GENES
                        Semicolon-separated truth gene IDs.
  --out-tsv OUT_TSV     Optional output path for counterfactual ranking TSV.
  --bonus BONUS         Counterfactual support score bonus (default: 0.5).
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-linking-export

```text
usage: candidate-qtl-bench homolog-linking-export [-h] --out-dir OUT_DIR --case-id CASE_ID --evidence EVIDENCE
                                                  [--phenotype PHENOTYPE] [--discover-manifest DISCOVER_MANIFEST]
                                                  [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR     Output directory for prompt package.
  --case-id CASE_ID     Case identifier.
  --evidence EVIDENCE   04_evidence.tsv from discovery.
  --phenotype PHENOTYPE
                        Phenotype description text.
  --discover-manifest DISCOVER_MANIFEST
                        discover_manifest.json from discovery.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-linking-import

```text
usage: candidate-qtl-bench homolog-linking-import [-h] --prompt PROMPT --response RESPONSE --out-dir OUT_DIR
                                                  --reviewed-by REVIEWED_BY [--payload-source PAYLOAD_SOURCE]
                                                  [--provider PROVIDER] [--model MODEL]
                                                  [--max-proposals MAX_PROPOSALS] [--json]

options:
  -h, --help            show this help message and exit
  --prompt PROMPT       homolog_linking_prompt.json or export manifest.
  --response RESPONSE   Agent linking response JSON file.
  --out-dir OUT_DIR     Output directory for proposal artifacts.
  --reviewed-by REVIEWED_BY
                        Reviewer identity (separate from model proposer).
  --payload-source PAYLOAD_SOURCE
                        Proposer category (default: agent).
  --provider PROVIDER   Provider / engine identifier.
  --model MODEL         Model name / version.
  --max-proposals MAX_PROPOSALS
                        Maximum number of candidate proposals allowed (budget).
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-linking-probe-export

```text
usage: candidate-qtl-bench homolog-linking-probe-export [-h] --out-dir OUT_DIR --case-id CASE_ID --phenotype PHENOTYPE
                                                        [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR     Output directory.
  --case-id CASE_ID     Case identifier.
  --phenotype PHENOTYPE
                        Phenotype description text.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-linking-score

```text
usage: candidate-qtl-bench homolog-linking-score [-h] --linking-manifest LINKING_MANIFEST
                                                 [--probe-response PROBE_RESPONSE] --truth-genes TRUTH_GENES
                                                 --truth-family TRUTH_FAMILY [--json]

options:
  -h, --help            show this help message and exit
  --linking-manifest LINKING_MANIFEST
                        homolog_linking_manifest.json.
  --probe-response PROBE_RESPONSE
                        Memorization probe response JSON (optional).
  --truth-genes TRUTH_GENES
                        Semicolon-separated truth gene IDs.
  --truth-family TRUTH_FAMILY
                        Truth family or gene symbol for probe matching.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-panel-auto

```text
usage: candidate-qtl-bench homolog-panel-auto [-h] --phenotype PHENOTYPE --species SPECIES --out-dir OUT_DIR
                                              [--direction DIRECTION] [--inheritance INHERITANCE] [--term TERM]
                                              [--llm-provider {openai-responses,openai-compatible,lmstudio-cli}]
                                              [--model MODEL] [--endpoint ENDPOINT] [--allow-external-api]
                                              [--allow-external-llm] [--panel-size PANEL_SIZE] [--per-query PER_QUERY]
                                              [--json]

options:
  -h, --help            show this help message and exit
  --phenotype PHENOTYPE
                        Free-text phenotype supplied before candidate inspection.
  --species SPECIES     Focal species name, used to exclude exact-species reference entries.
  --out-dir OUT_DIR
  --direction DIRECTION
  --inheritance INHERITANCE
  --term TERM           Additional user-supplied phenotype phrase; repeatable.
  --llm-provider {openai-responses,openai-compatible,lmstudio-cli}
                        Optional provider used only to propose search hypotheses.
  --model MODEL         Model name required with --llm-provider.
  --endpoint ENDPOINT   Optional LLM HTTP endpoint, or lms executable path for lmstudio-cli.
  --allow-external-api  Allow phenotype text to be sent to the UniProt REST API.
  --allow-external-llm  Allow phenotype text to be sent to a non-loopback LLM endpoint.
  --panel-size PANEL_SIZE
  --per-query PER_QUERY
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-panel-benchmark

```text
usage: candidate-qtl-bench homolog-panel-benchmark [-h] --reference-manifest REFERENCE_MANIFEST --out-dir OUT_DIR
                                                   [--species SPECIES] [--case CASE]
                                                   [--llm-provider {openai-responses,openai-compatible,lmstudio-cli}]
                                                   [--model MODEL] [--endpoint ENDPOINT] [--allow-external-api]
                                                   [--allow-external-llm] [--panel-size PANEL_SIZE]
                                                   [--per-query PER_QUERY] [--target-size TARGET_SIZE]
                                                   [--runner {auto,native,wsl}] [--threads THREADS] [--resume]
                                                   [--benchmark-curated-terms] [--json]

options:
  -h, --help            show this help message and exit
  --reference-manifest REFERENCE_MANIFEST
                        benchmark-reference manifest containing frozen phenotype-only cases.
  --out-dir OUT_DIR
  --species SPECIES
  --case CASE           Locus ID to benchmark; repeatable. Default: all phenotype-only cases.
  --llm-provider {openai-responses,openai-compatible,lmstudio-cli}
  --model MODEL
  --endpoint ENDPOINT
  --allow-external-api
  --allow-external-llm
  --panel-size PANEL_SIZE
  --per-query PER_QUERY
  --target-size TARGET_SIZE
  --runner {auto,native,wsl}
  --threads THREADS
  --resume
  --benchmark-curated-terms
                        Ablation only: include pre-curated phenotype terms. Primary benchmark uses raw trait text
                        only.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homolog-panel-match

```text
usage: candidate-qtl-bench homolog-panel-match [-h] --candidates CANDIDATES --candidate-proteins CANDIDATE_PROTEINS
                                               --panel PANEL --panel-fasta PANEL_FASTA --out-dir OUT_DIR
                                               [--base-evidence BASE_EVIDENCE] [--target-size TARGET_SIZE]
                                               [--runner {auto,native,wsl}] [--threads THREADS] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --candidate-proteins CANDIDATE_PROTEINS
  --panel PANEL
  --panel-fasta PANEL_FASTA
  --out-dir OUT_DIR
  --base-evidence BASE_EVIDENCE
                        Existing canonical evidence TSV; when supplied, rank and freeze Stage A.
  --target-size TARGET_SIZE
  --runner {auto,native,wsl}
  --threads THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench homology-concordance

```text
usage: candidate-qtl-bench homology-concordance [-h] --hit-table HIT_TABLE --out-dir OUT_DIR [--top-n TOP_N]
                                                [--domain-table DOMAIN_TABLE] [--orthology-table ORTHOLOGY_TABLE]
                                                [--candidates CANDIDATES] [--provider-metadata PROVIDER_METADATA]
                                                [--json]

options:
  -h, --help            show this help message and exit
  --hit-table HIT_TABLE
                        PROVIDER=12/13-column outfmt TSV; repeatable
  --out-dir OUT_DIR
  --top-n TOP_N
  --domain-table DOMAIN_TABLE
                        PROVIDER=TSV with gene_id/domain_id/start/end
  --orthology-table ORTHOLOGY_TABLE
                        PROVIDER=TSV with gene_id/orthogroup/relationship
  --candidates CANDIDATES
  --provider-metadata PROVIDER_METADATA
                        PROVIDER=JSON with tool/database version and command
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench hypothesis-ledger

```text
usage: candidate-qtl-bench hypothesis-ledger [-h] --out-dir OUT_DIR [--hypotheses HYPOTHESES] [--case-id CASE_ID]
                                             [--evidence EVIDENCE] [--wet-results WET_RESULTS] [--sv-events SV_EVENTS]
                                             [--wet-candidate-state WET_CANDIDATE_STATE]
                                             [--annotation-rescue ANNOTATION_RESCUE]
                                             [--dual-route-comparison DUAL_ROUTE_COMPARISON]
                                             [--previous-history PREVIOUS_HISTORY] [--ranked RANKED] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --hypotheses HYPOTHESES
                        Hypothesis proposal/current-state TSV.
  --case-id CASE_ID
  --evidence EVIDENCE   Canonical evidence table with observation_id or evidence_id.
  --wet-results WET_RESULTS
                        Normalized wet_results.tsv with result_id.
  --sv-events SV_EVENTS
                        Canonical SV table with sv_event_id, call_id, or cluster_id.
  --wet-candidate-state WET_CANDIDATE_STATE
                        Normalized candidate_wet_state.tsv.
  --annotation-rescue ANNOTATION_RESCUE
                        annotation-rescue audit TSV; native_discordant/reference_rescue/ambiguous_overlap rows seed
                        annotation_missing hypotheses.
  --dual-route-comparison DUAL_ROUTE_COMPARISON
                        dual-route-compare comparison TSV; A/B disagreement seeds linked_qtl or
                        unknown_locus_mechanism hypotheses.
  --previous-history PREVIOUS_HISTORY
                        Previous hypothesis_history.tsv for append-only transitions.
  --ranked RANKED       Optional deterministic rank input; fingerprinted and never modified.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench hypothesis-verify

```text
usage: candidate-qtl-bench hypothesis-verify [-h] --hypotheses HYPOTHESES --out-dir OUT_DIR --focal-species
                                             FOCAL_SPECIES [--record-cache RECORD_CACHE]
                                             [--source-records SOURCE_RECORDS] [--bound-claims BOUND_CLAIMS]
                                             [--agent-claim-review AGENT_CLAIM_REVIEW]
                                             [--source-review-records SOURCE_REVIEW_RECORDS] [--allow-external-api]
                                             [--benchmark-mode] [--min-panel-size MIN_PANEL_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --hypotheses HYPOTHESES
                        llm_pathway_hypotheses.tsv
  --out-dir OUT_DIR
  --focal-species FOCAL_SPECIES
                        Focal species; a reference gene from it is rejected.
  --record-cache RECORD_CACHE
                        Local UniProt-shaped JSON record file or directory for offline verification.
  --source-records SOURCE_RECORDS
                        verified_sources.json from source-verify; required to validate publication annotations in
                        canonical Route A.
  --bound-claims BOUND_CLAIMS
                        Shared bound-claim TSV joining citation, gene, experiment and phenotype. Without it,
                        source/cue evidence remains review-only.
  --agent-claim-review AGENT_CLAIM_REVIEW
                        Independent agent review packet with frozen raw response, prompt and separate executor
                        provenance; alternative to --bound-claims.
  --source-review-records SOURCE_REVIEW_RECORDS
                        source_review_records.json from source-verify; required with --agent-claim-review for exact
                        text/ID/species checks.
  --allow-external-api  Allow accession lookups against the UniProt REST API.
  --benchmark-mode      Additionally reject focal-lineage / answer-backmapped markers.
  --min-panel-size MIN_PANEL_SIZE
                        Panels below this size are review-required, not rejected.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench integrated-narrowing

```text
usage: candidate-qtl-bench integrated-narrowing [-h] --ranked RANKED --out-dir OUT_DIR
                                                [--family-recovery FAMILY_RECOVERY] [--convergence CONVERGENCE]
                                                [--family-copy FAMILY_COPY]
                                                [--family-copy-shortlist FAMILY_COPY_SHORTLIST]
                                                [--remaining REMAINING] [--case-id CASE_ID]
                                                [--family-copy-manifest FAMILY_COPY_MANIFEST] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED       ranked_candidates.tsv from evidence-rank or candidate-pipeline; all columns are passed through
                        verbatim.
  --out-dir OUT_DIR
  --family-recovery FAMILY_RECOVERY
                        family_to_copy_input.tsv (canonical triple-route handoff:
                        hypothesis_id/gene_family/member_interval_genes) resolving gene -> family.
  --convergence CONVERGENCE
                        triple_route_hypothesis_comparison.tsv (full per-family recovery states) or
                        triple_route_convergence.json (degraded: convergent families only).
  --family-copy FAMILY_COPY
                        family_copy_discrimination.tsv; disambiguation_status is passed through as copy_state.
  --family-copy-shortlist FAMILY_COPY_SHORTLIST
                        family_copy_shortlist.tsv; first_wet_priority is passed through when present.
  --remaining REMAINING
                        remaining_candidate_set.tsv from candidate-elimination; survival-state columns only.
  --case-id CASE_ID
  --family-copy-manifest FAMILY_COPY_MANIFEST
                        family_copy_manifest.json from family-to-copy.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench interval-atlas

```text
usage: candidate-qtl-bench interval-atlas [-h] --candidates CANDIDATES --out-dir OUT_DIR [--family-index FAMILY_INDEX]
                                          [--domains DOMAINS] [--orthology ORTHOLOGY] [--annotations ANNOTATIONS]
                                          [--flags FLAGS] [--function-bins FUNCTION_BINS]
                                          [--tandem-distance TANDEM_DISTANCE] [--dark-max-aa DARK_MAX_AA]
                                          [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Interval candidates TSV (gene_id + optional chrom/start/end/protein_length).
  --out-dir OUT_DIR
  --family-index FAMILY_INDEX
                        family_index.tsv; reused rather than recomputed.
  --domains DOMAINS     Domain table (gene_id/domain_id[/domain_name/start/source]).
  --orthology ORTHOLOGY
                        Orthology table (gene_id[/orthogroup/relationship/reference_gene_id/provider]).
  --annotations ANNOTATIONS
                        Reference annotation table
                        (gene_id[/description/function_description/go_terms/protein_family]).
  --flags FLAGS         candidate_flags.tsv, for dark_state and tandem context.
  --function-bins FUNCTION_BINS
                        Override the versioned broad-function bin map (templates/broad_function_bins.tsv).
  --tandem-distance TANDEM_DISTANCE
  --dark-max-aa DARK_MAX_AA
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench interval-reverse-lookup

```text
usage: candidate-qtl-bench interval-reverse-lookup [-h] --atlas ATLAS --out-dir OUT_DIR
                                                   [--reference-function REFERENCE_FUNCTION]
                                                   [--llm-proposals LLM_PROPOSALS] [--phenotype PHENOTYPE]
                                                   [--shared-route SHARED_ROUTE] [--min-members MIN_MEMBERS]
                                                   [--min-semantic-coverage MIN_SEMANTIC_COVERAGE] [--case-id CASE_ID]
                                                   [--json]

options:
  -h, --help            show this help message and exit
  --atlas ATLAS         interval_gene_atlas.tsv: the canonical <out-dir>/interval_gene_atlas.tsv from route-b-prepare,
                        or the root-level file from interval-atlas.
  --out-dir OUT_DIR
  --reference-function REFERENCE_FUNCTION
                        Reference gene function table
                        (reference_gene_id[/function_description/pathway/mutant_phenotype/source_id]).
  --llm-proposals LLM_PROPOSALS
                        llm_pathway_hypotheses.tsv imported from an LLM run of the emitted reverse-lookup prompt.
  --phenotype PHENOTYPE
                        Optional phenotype overlay; recorded in the status, and it makes the lane non-blind.
  --shared-route SHARED_ROUTE
                        Optional pre-fixed phenotype->pathway route TSV
                        (phenotype_token/pathway_family/annotation_keywords). Its pathway family and keywords form a
                        separately reported, non-scoring B1 overlay vocabulary.
  --min-members MIN_MEMBERS
  --min-semantic-coverage MIN_SEMANTIC_COVERAGE
                        Minimum fraction (0-1) of interval genes with family/orthogroup/domain/broad-function/source-
                        backed-function annotation. At least one annotated gene is always required for semantic
                        evaluability.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench investigation-attribution

```text
usage: candidate-qtl-bench investigation-attribution [-h] --out-dir OUT_DIR [--mechanisms MECHANISMS]
                                                     [--validation VALIDATION] [--panel PANEL]
                                                     [--interval-match INTERVAL_MATCH] [--ranked RANKED]
                                                     [--case-id CASE_ID] [--min-mechanisms MIN_MECHANISMS] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --mechanisms MECHANISMS
  --validation VALIDATION
  --panel PANEL
  --interval-match INTERVAL_MATCH
  --ranked RANKED
  --case-id CASE_ID
  --min-mechanisms MIN_MECHANISMS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench liftover-markers

```text
usage: candidate-qtl-bench liftover-markers [-h] (--hits HITS | --marker-fasta MARKER_FASTA) [--db DB]
                                            [--threads THREADS] --out-dir OUT_DIR [--max-product MAX_PRODUCT]
                                            [--max-block-gap MAX_BLOCK_GAP] [--json]

options:
  -h, --help            show this help message and exit
  --hits HITS           Pair table or 10-column BLAST TSV.
  --marker-fasta MARKER_FASTA
                        Primer FASTA with IDs ending _F/_R; runs blastn-short.
  --db DB               Nucleotide BLAST DB for --marker-fasta.
  --threads THREADS
  --out-dir OUT_DIR
  --max-product MAX_PRODUCT
  --max-block-gap MAX_BLOCK_GAP
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench literature

```text
usage: candidate-qtl-bench literature [-h] --query QUERY --out OUT [--limit LIMIT] [--json]

options:
  -h, --help     show this help message and exit
  --query QUERY  Europe PMC search query.
  --out OUT      Output JSON path.
  --limit LIMIT
  --json         Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-evidence-import

```text
usage: candidate-qtl-bench llm-evidence-import [-h] --candidates CANDIDATES --input INPUT --out OUT
                                               [--resolve-against PATH [PATH ...]] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --input INPUT
  --out OUT
  --resolve-against PATH [PATH ...]
                        Local corpus files (reference proteome FASTA, annotation TSV, citation list) a source_id must
                        be found in. Offline and format-agnostic. Without this, accepted rows are marked unverified
                        and share one independence group, so a well-formed but non-existent accession cannot reach
                        independent_experimental_support.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-pathway-import

```text
usage: candidate-qtl-bench llm-pathway-import [-h] --raw RAW --out-dir OUT_DIR [--input-manifest INPUT_MANIFEST]
                                              [--model MODEL] [--proposal-source PROPOSAL_SOURCE] [--case-id CASE_ID]
                                              [--execution-mode {blind,assisted}] [--json]

options:
  -h, --help            show this help message and exit
  --raw RAW             Raw model output (.json/.jsonl/.txt), or a manually written proposal file.
  --out-dir OUT_DIR
  --input-manifest INPUT_MANIFEST
                        llm_input_manifest.json from llm-pathway-prompt.
  --model MODEL         Model/provider version that produced the output.
  --proposal-source PROPOSAL_SOURCE
                        Proposal source, e.g. a provider name or `manual`.
  --case-id CASE_ID
  --execution-mode {blind,assisted}
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-pathway-prompt

```text
usage: candidate-qtl-bench llm-pathway-prompt [-h] --out-dir OUT_DIR --case-id CASE_ID --focal-species FOCAL_SPECIES
                                              --phenotype PHENOTYPE
                                              [--direction {increase,decrease,gain,loss,unknown}]
                                              [--inheritance {recessive,dominant,semidominant,unknown}]
                                              [--organ ORGAN] [--stage STAGE] [--mode {blind,assisted}]
                                              [--constraint CONSTRAINTS] [--exclusion EXCLUSIONS]
                                              [--assisted-context ASSISTED_CONTEXT] [--candidates CANDIDATES] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --case-id CASE_ID
  --focal-species FOCAL_SPECIES
  --phenotype PHENOTYPE
  --direction {increase,decrease,gain,loss,unknown}
  --inheritance {recessive,dominant,semidominant,unknown}
  --organ ORGAN         Affected organ or tissue.
  --stage STAGE         Developmental stage.
  --mode {blind,assisted}
                        blind: phenotype only. assisted: researcher context allowed, recorded separately in the
                        manifest.
  --constraint CONSTRAINTS
                        Researcher constraint; repeatable (assisted or blind).
  --exclusion EXCLUSIONS
                        Known exclusion; repeatable.
  --assisted-context ASSISTED_CONTEXT
                        Assisted mode only: file with researcher hypotheses / verified literature / ontology terms.
  --candidates CANDIDATES
                        Interval candidates TSV used ONLY to audit that no candidate id leaked into a blind prompt.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-prompt-export

```text
usage: candidate-qtl-bench llm-prompt-export [-h] --packet PACKET --out OUT [--json]

options:
  -h, --help       show this help message and exit
  --packet PACKET  research_packet.json
  --out OUT        Prompt JSONL output.
  --json           Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-review-benchmark

```text
usage: candidate-qtl-bench llm-review-benchmark [-h] --config CONFIG --out-dir OUT_DIR [--json]

options:
  -h, --help         show this help message and exit
  --config CONFIG    TSV with case, truth, shortlist, and two review ranking paths.
  --out-dir OUT_DIR
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-run

```text
usage: candidate-qtl-bench llm-run [-h] --prompts PROMPTS --out-dir OUT_DIR --provider
                                   {openai-responses,openai-compatible,lmstudio-cli,codex-agent} --model MODEL
                                   [--endpoint ENDPOINT] [--agent-file AGENT_FILE] [--api-key-env API_KEY_ENV]
                                   [--allow-external-api] [--max-output-tokens MAX_OUTPUT_TOKENS] [--timeout TIMEOUT]
                                   [--json]

options:
  -h, --help            show this help message and exit
  --prompts PROMPTS
  --out-dir OUT_DIR
  --provider {openai-responses,openai-compatible,lmstudio-cli,codex-agent}
  --model MODEL
  --endpoint ENDPOINT   HTTP endpoint, or lms executable path for lmstudio-cli.
  --agent-file AGENT_FILE
                        Offline Codex/OpenCode JSON/JSONL artifact; use with --provider codex-agent.
  --api-key-env API_KEY_ENV
  --allow-external-api
  --max-output-tokens MAX_OUTPUT_TOKENS
  --timeout TIMEOUT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-shortlist-import

```text
usage: candidate-qtl-bench llm-shortlist-import [-h] --shortlist SHORTLIST --responses RESPONSES --out-dir OUT_DIR
                                                [--json]

options:
  -h, --help            show this help message and exit
  --shortlist SHORTLIST
  --responses RESPONSES
                        llm_candidate_responses.json from llm-run.
  --out-dir OUT_DIR
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-shortlist-prompt

```text
usage: candidate-qtl-bench llm-shortlist-prompt [-h] --packet PACKET --out OUT [--json]

options:
  -h, --help       show this help message and exit
  --packet PACKET  research_packet.json for the deterministic shortlist.
  --out OUT        Output prompt JSONL.
  --json           Emit a stable JSON result to stdout.
```

## candidate-qtl-bench llm-usage

```text
usage: candidate-qtl-bench llm-usage [-h] --raw-responses RAW_RESPONSES [--llm-manifest LLM_MANIFEST] --out-dir
                                     OUT_DIR [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --raw-responses RAW_RESPONSES
                        llm_raw_responses.jsonl from llm-run.
  --llm-manifest LLM_MANIFEST
                        Optional llm_run_manifest.json for provider/model/prompt-count validation.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench mapping-gate

```text
usage: candidate-qtl-bench mapping-gate [-h] --metrics METRICS --out-dir OUT_DIR [--thresholds THRESHOLDS] [--json]

options:
  -h, --help            show this help message and exit
  --metrics METRICS     TSV with metric, value, and source columns.
  --out-dir OUT_DIR
  --thresholds THRESHOLDS
                        Optional TSV with threshold and value columns.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench mechanism-coverage-audit

```text
usage: candidate-qtl-bench mechanism-coverage-audit [-h] --mechanisms MECHANISMS --mechanism-import-manifest
                                                    MECHANISM_IMPORT_MANIFEST --localization-freeze-manifest
                                                    LOCALIZATION_FREEZE_MANIFEST --investigation INVESTIGATION
                                                    --out-dir OUT_DIR [--candidates CANDIDATES]
                                                    [--gap-request-audit GAP_REQUEST_AUDIT] [--json]

options:
  -h, --help            show this help message and exit
  --mechanisms MECHANISMS
  --mechanism-import-manifest MECHANISM_IMPORT_MANIFEST
  --localization-freeze-manifest LOCALIZATION_FREEZE_MANIFEST
  --investigation INVESTIGATION
  --out-dir OUT_DIR
  --candidates CANDIDATES
  --gap-request-audit GAP_REQUEST_AUDIT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench mechanism-import

```text
usage: candidate-qtl-bench mechanism-import [-h] --payload PAYLOAD [--out-dir OUT_DIR] [--preflight]
                                            [--investigation INVESTIGATION] [--payload-source {llm,agent,human}]
                                            [--model MODEL] [--case-id CASE_ID]
                                            [--retrieval-provenance RETRIEVAL_PROVENANCE] [--candidates CANDIDATES]
                                            [--min-mechanisms MIN_MECHANISMS] [--proposal-round PROPOSAL_ROUND]
                                            [--interval-informed]
                                            [--mechanism-freeze-manifest MECHANISM_FREEZE_MANIFEST]
                                            [--localization-freeze-manifest LOCALIZATION_FREEZE_MANIFEST]
                                            [--benchmark-mode] [--json]

options:
  -h, --help            show this help message and exit
  --payload PAYLOAD
  --out-dir OUT_DIR
  --preflight           Validate the original payload and paper reports without creating import artifacts.
  --investigation INVESTIGATION
  --payload-source {llm,agent,human}
  --model MODEL
  --case-id CASE_ID
  --retrieval-provenance RETRIEVAL_PROVENANCE
  --candidates CANDIDATES
  --min-mechanisms MIN_MECHANISMS
  --proposal-round PROPOSAL_ROUND
  --interval-informed
  --mechanism-freeze-manifest MECHANISM_FREEZE_MANIFEST
  --localization-freeze-manifest LOCALIZATION_FREEZE_MANIFEST
  --benchmark-mode
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench multi-reference-locus-union

```text
usage: candidate-qtl-bench multi-reference-locus-union [-h] --primary-candidates PRIMARY_CANDIDATES
                                                       [--alternate-candidates LABEL=PATH] --mapping MAPPING --out-dir
                                                       OUT_DIR [--primary-label PRIMARY_LABEL] [--case-id CASE_ID]
                                                       [--json]

options:
  -h, --help            show this help message and exit
  --primary-candidates PRIMARY_CANDIDATES, --primary PRIMARY_CANDIDATES
                        Primary/run-target candidate table.
  --alternate-candidates LABEL=PATH, --alternate LABEL=PATH
                        Alternate candidate table; repeat LABEL=PATH. Alternate-only IDs are retained.
  --mapping MAPPING, --assembly-mapping MAPPING
                        Full assembly_bridge_mapping.tsv; all states are audited, only safe links collapse.
  --out-dir OUT_DIR
  --primary-label PRIMARY_LABEL
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench multigenome-private

```text
usage: candidate-qtl-bench multigenome-private [-h] --out-dir OUT_DIR --focal-vcf FOCAL_VCF --gff3 GFF3 [--case CASE]
                                               [--control-vcf CONTROL_VCF] [--callable-bed CALLABLE_BED]
                                               [--projection-ledger PROJECTION_LEDGER] [--native-genes NATIVE_GENES]
                                               [--reference-consequence REFERENCE_CONSEQUENCE]
                                               [--native-consequence NATIVE_CONSEQUENCE]
                                               [--focal-accession FOCAL_ACCESSION]
                                               [--phenotype-control PHENOTYPE_CONTROLS]
                                               [--comparison-accession COMPARISON_ACCESSIONS] [--case-id CASE_ID]
                                               [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --focal-vcf FOCAL_VCF
                        Normalized focal accession VCF/VCF.GZ.
  --gff3 GFF3           Reference GFF3 used for gene assignment and CDS denominators.
  --case CASE           Case manifest containing private_variant_model (focal/controls/callability/loci).
  --control-vcf CONTROL_VCF
                        ACCESSION=path to a control VCF; repeatable.
  --callable-bed CALLABLE_BED
                        ACCESSION=path to a callable-regions BED; repeatable.
  --projection-ledger PROJECTION_LEDGER
                        Per reference gene projection state TSV (gene_id/projection_state[/projection_run_complete]).
  --native-genes NATIVE_GENES
                        Optional native annotation gene table, reported as a denominator.
  --reference-consequence REFERENCE_CONSEQUENCE
                        Reference-frame consequence TSV (chrom/pos/ref/alt/consequence[/gene_id]).
  --native-consequence NATIVE_CONSEQUENCE
                        Focal native-frame consequence TSV with the same columns.
  --focal-accession FOCAL_ACCESSION
                        Focal accession; overrides the case manifest.
  --phenotype-control PHENOTYPE_CONTROLS
                        Required phenotype control accession; repeatable.
  --comparison-accession COMPARISON_ACCESSIONS
                        Additional comparison accession kept in separate columns; repeatable.
  --case-id CASE_ID
  --full-hash           SHA-256 all inputs, including large VCFs.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench next-experiment

```text
usage: candidate-qtl-bench next-experiment [-h] --ranked RANKED --case CASE --out-dir OUT_DIR
                                           [--mapping-gate MAPPING_GATE]
                                           [--candidate-path-decision CANDIDATE_PATH_DECISION]
                                           [--candidate-flags CANDIDATE_FLAGS] [--assay-outcomes ASSAY_OUTCOMES]
                                           [--exon-coverage EXON_COVERAGE]
                                           [--informative-recombinants INFORMATIVE_RECOMBINANTS]
                                           [--wet-candidate-state WET_CANDIDATE_STATE] [--hypotheses HYPOTHESES]
                                           [--remaining-candidates REMAINING_CANDIDATES] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED       Deterministic ranked candidate TSV; read-only biological priority input.
  --case CASE           Case manifest with decision_support constraints.
  --out-dir OUT_DIR
  --mapping-gate MAPPING_GATE
  --candidate-path-decision CANDIDATE_PATH_DECISION
  --candidate-flags CANDIDATE_FLAGS
  --assay-outcomes ASSAY_OUTCOMES
  --exon-coverage EXON_COVERAGE
  --informative-recombinants INFORMATIVE_RECOMBINANTS
  --wet-candidate-state WET_CANDIDATE_STATE
                        Normalized candidate_wet_state.tsv from wet-state.
  --hypotheses HYPOTHESES
                        Locus hypothesis ledger; read-only action context, never biological rank input.
  --remaining-candidates REMAINING_CANDIDATES
                        Canonical remaining_candidate_set.tsv from candidate-elimination; post-elimination candidate
                        counts are used when supplied, never ranked reordering.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench omamer-import

```text
usage: candidate-qtl-bench omamer-import [-h] --omamer OMAMER --candidates CANDIDATES --out OUT [--json]

options:
  -h, --help            show this help message and exit
  --omamer OMAMER       OMAmer search output TSV.
  --candidates CANDIDATES
  --out OUT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench orchestrate

```text
usage: candidate-qtl-bench orchestrate [-h] --candidates CANDIDATES --out-dir OUT_DIR [--protein-fasta PROTEIN_FASTA]
                                       [--blastp-db BLASTP_DB] [--diamond-db DIAMOND_DB] [--hmmer-db HMMER_DB]
                                       [--eggnog-data-dir EGGNOG_DATA_DIR] [--primerblast-root PRIMERBLAST_ROOT]
                                       [--genome GENOME] [--gff3 GFF3] [--nucleotide-db NUCLEOTIDE_DB]
                                       [--candidate-genomic-fasta CANDIDATE_GENOMIC_FASTA] [--interval INTERVAL]
                                       [--variants VARIANTS] [--vcf VCF] [--top-n TOP_N] [--num-threads NUM_THREADS]
                                       [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Dry candidate table or ranked TSV.
  --out-dir OUT_DIR
  --protein-fasta PROTEIN_FASTA
  --blastp-db BLASTP_DB
  --diamond-db DIAMOND_DB
  --hmmer-db HMMER_DB
  --eggnog-data-dir EGGNOG_DATA_DIR
  --primerblast-root PRIMERBLAST_ROOT
  --genome GENOME
  --gff3 GFF3
  --nucleotide-db NUCLEOTIDE_DB
                        BLAST nucleotide DB for primerblast-oss; repeatable.
  --candidate-genomic-fasta CANDIDATE_GENOMIC_FASTA
  --interval INTERVAL   Optional QTL interval for coarse grid-marker design.
  --variants VARIANTS   Optional confirmed parental variant TSV.
  --vcf VCF
  --top-n TOP_N
  --num-threads NUM_THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench orthofinder-import

```text
usage: candidate-qtl-bench orthofinder-import [-h] --orthogroups ORTHOGROUPS --candidates CANDIDATES --focal-species
                                              FOCAL_SPECIES --out OUT [--json]

options:
  -h, --help            show this help message and exit
  --orthogroups ORTHOGROUPS
                        OrthoFinder Orthogroups.tsv.
  --candidates CANDIDATES
  --focal-species FOCAL_SPECIES
                        Column name for the focal species in Orthogroups.tsv.
  --out OUT             Canonical evidence TSV output.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench orthology-import

```text
usage: candidate-qtl-bench orthology-import [-h] --manifest MANIFEST --out OUT [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  orthology_plan.json after running its shell script.
  --out OUT
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench orthology-plan

```text
usage: candidate-qtl-bench orthology-plan [-h] --ranked RANKED --protein-fasta PROTEIN_FASTA --target-proteome
                                          TARGET_PROTEOME --relative-proteome RELATIVE_PROTEOME --out-dir OUT_DIR
                                          [--top-n TOP_N] [--num-threads NUM_THREADS] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED
  --protein-fasta PROTEIN_FASTA
                        Candidate proteins normalized to gene_id.
  --target-proteome TARGET_PROTEOME
                        Complete target-species proteome with gene-compatible IDs.
  --relative-proteome RELATIVE_PROTEOME
                        LABEL=FASTA; repeatable for near relatives.
  --out-dir OUT_DIR
  --top-n TOP_N
  --num-threads NUM_THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench orthology-tier1-plan

```text
usage: candidate-qtl-bench orthology-tier1-plan [-h] --candidates CANDIDATES --protein-fasta PROTEIN_FASTA
                                                --focal-proteome FOCAL_PROTEOME
                                                [--reference-proteome REFERENCE_PROTEOME] --out-dir OUT_DIR
                                                [--profile {accuracy64,laptop16}] [--profiles-table PROFILES_TABLE]
                                                [--aligner {diamond,blastp}] [--threads THREADS] [--max-hits MAX_HITS]
                                                [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate table; only its gene_id proteins are used as query.
  --protein-fasta PROTEIN_FASTA
                        Candidate proteins normalized to gene_id.
  --focal-proteome FOCAL_PROTEOME
                        Whole focal-species proteome (search target, reused across QTLs).
  --reference-proteome REFERENCE_PROTEOME
                        LABEL=FASTA reference panel proteome; repeatable.
  --out-dir OUT_DIR
  --profile {accuracy64,laptop16}
  --profiles-table PROFILES_TABLE
                        Optional runtime_profiles.tsv overriding the built-in defaults.
  --aligner {diamond,blastp}
  --threads THREADS     Override the profile thread count.
  --max-hits MAX_HITS   Override the profile per-query hit cap.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench pathway-interval-falsify

```text
usage: candidate-qtl-bench pathway-interval-falsify [-h] --family-index FAMILY_INDEX --hypotheses HYPOTHESES --out-dir
                                                    OUT_DIR [--annotation-rescue ANNOTATION_RESCUE]
                                                    [--interval INTERVAL] [--candidates CANDIDATES]
                                                    [--tandem-distance TANDEM_DISTANCE] [--max-q MAX_Q]
                                                    [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --family-index FAMILY_INDEX
                        family_index.tsv from family-index.
  --hypotheses HYPOTHESES
                        llm_pathway_hypotheses.tsv, or any TSV with hypothesis_id +
                        gene_family/reference_gene_name/proposed_stable_id.
  --out-dir OUT_DIR
  --annotation-rescue ANNOTATION_RESCUE
                        annotation_rescue.tsv (miniprot projection of the reference proteins onto the interval
                        sequence). Without it an interval miss stays absent_annotation_only -- a hold -- because 'not
                        in the annotated gene list' is not 'absent'.
  --interval INTERVAL   CHROM:START-END resolved against the index coordinates; repeatable.
  --candidates CANDIDATES
                        Interval candidates TSV, as an alternative to --interval.
  --tandem-distance TANDEM_DISTANCE
                        Same-family genes closer than this collapse into one block.
  --max-q MAX_Q         BH-adjusted threshold separating enriched from uninformative survivors.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench peak-context

```text
usage: candidate-qtl-bench peak-context [-h] [--candidates CANDIDATES] --out-dir OUT_DIR [--lod-profile LOD_PROFILE]
                                        [--lod-drop LOD_DROP] [--peak-metadata PEAK_METADATA] [--peak PEAK]
                                        [--interval INTERVAL] [--peak-window PEAK_WINDOW] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidates.tsv with gene_id/chrom/start/end; position evidence is emitted per gene when
                        supplied.
  --out-dir OUT_DIR
  --lod-profile LOD_PROFILE
                        TSV with chrom, position, and lod columns for a LOD-drop support interval.
  --lod-drop LOD_DROP   LOD drop for the support interval (default 1.5).
  --peak-metadata PEAK_METADATA
                        Single-QTL TSV following templates/qtl_peak_metadata.tsv.
  --peak PEAK           QTL peak as CHROM:POS or POS when no profile/metadata is available.
  --interval INTERVAL   QTL interval CHROM:START-END; used only as context, never as positive evidence.
  --peak-window PEAK_WINDOW
                        Analyst-supplied uncertainty half-width for peak_proximal support.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-first-status

```text
usage: candidate-qtl-bench phenotype-first-status [-h] --out-dir OUT_DIR [--hypotheses HYPOTHESES]
                                                  [--validation VALIDATION] [--panel PANEL]
                                                  [--falsification FALSIFICATION] [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --hypotheses HYPOTHESES
                        llm_pathway_hypotheses.tsv from the phenotype-first proposal route.
  --validation VALIDATION
                        proposal_validation.tsv from hypothesis-verify.
  --panel PANEL         verified_bait_panel.tsv from bait-panel-freeze.
  --falsification FALSIFICATION
                        hypothesis_interval_falsification.tsv from pathway-interval-falsify.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-investigate

```text
usage: candidate-qtl-bench phenotype-investigate [-h] --out-dir OUT_DIR --case-id CASE_ID --focal-species
                                                 FOCAL_SPECIES --phenotype PHENOTYPE [--direction DIRECTION]
                                                 [--inheritance INHERITANCE] [--tissue TISSUE] [--stage STAGE]
                                                 [--onset ONSET] [--penetrance PENETRANCE] [--synonym SYNONYM]
                                                 [--ontology-term ONTOLOGY_TERM] [--exclusion EXCLUSION]
                                                 [--constraint CONSTRAINT] [--uncertainty UNCERTAINTY]
                                                 [--candidates CANDIDATES] [--benchmark-mode] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --case-id CASE_ID
  --focal-species FOCAL_SPECIES
  --phenotype PHENOTYPE
  --direction DIRECTION
  --inheritance INHERITANCE
  --tissue TISSUE
  --stage STAGE
  --onset ONSET
  --penetrance PENETRANCE
  --synonym SYNONYM
  --ontology-term ONTOLOGY_TERM
  --exclusion EXCLUSION
  --constraint CONSTRAINT
  --uncertainty UNCERTAINTY
  --candidates CANDIDATES
  --benchmark-mode
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-localization-freeze

```text
usage: candidate-qtl-bench phenotype-localization-freeze [-h] --payload PAYLOAD --investigation INVESTIGATION
                                                         --out-dir OUT_DIR [--json]

options:
  -h, --help            show this help message and exit
  --payload PAYLOAD
  --investigation INVESTIGATION
  --out-dir OUT_DIR
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-localization-packet

```text
usage: candidate-qtl-bench phenotype-localization-packet [-h] --investigation INVESTIGATION --out-dir OUT_DIR
                                                         [--candidates CANDIDATES] [--json]

options:
  -h, --help            show this help message and exit
  --investigation INVESTIGATION
  --out-dir OUT_DIR
  --candidates CANDIDATES
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-pathway-benchmark

```text
usage: candidate-qtl-bench phenotype-pathway-benchmark [-h] --config CONFIG --out-dir OUT_DIR [--top-n TOP_N] [--json]

options:
  -h, --help         show this help message and exit
  --config CONFIG    TSV: case_id/expected_pathway/expected_gene_family[/causal_gene_id/hypotheses/validation/panel/ra
                     nked].
  --out-dir OUT_DIR
  --top-n TOP_N
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench phenotype-variant-benchmark

```text
usage: candidate-qtl-bench phenotype-variant-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR [--full-hash]
                                                       [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  TSV: one row per case (case_id/truth_family_id/variants_tsv + optional
                       truth_gene_id/base_phenotype and per-variant artifact patterns route_a_families_tsv/route_c_ret
                       rieval_tsv/route_b_overlay_status_json/route_b_atlas_manifest_json/route_b_families_tsv/ranked_
                       tsv; {variant_id} is substituted, else <stem>_<variant_id>.<ext>). variants.tsv rows:
                       variant_id/variant_kind/phenotype_text.
  --out-dir OUT_DIR
  --full-hash
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench primer-plan

```text
usage: candidate-qtl-bench primer-plan [-h] --stage {direct-sequencing,variant-marker,qtl-grid} --out-dir OUT_DIR
                                       --primerblast-root PRIMERBLAST_ROOT --genome GENOME --db DB
                                       [--num-threads NUM_THREADS] [--candidates CANDIDATES] [--top-n TOP_N]
                                       [--gff3 GFF3] [--candidate-genomic-fasta CANDIDATE_GENOMIC_FASTA]
                                       [--variants VARIANTS] [--vcf VCF] [--interval INTERVAL] [--n-markers N_MARKERS]
                                       [--spacing SPACING] [--json]

options:
  -h, --help            show this help message and exit
  --stage {direct-sequencing,variant-marker,qtl-grid}
  --out-dir OUT_DIR
  --primerblast-root PRIMERBLAST_ROOT
                        primerblast-oss checkout root.
  --genome GENOME       Reference FASTA for primerblast-oss.
  --db DB               BLAST nucleotide DB; repeat for multi-reference screening.
  --num-threads NUM_THREADS
                        Laptop-safe worker count for BLAST-backed primer checks.
  --candidates CANDIDATES
                        Dry-ranked candidate TSV; required for direct-sequencing.
  --top-n TOP_N
  --gff3 GFF3           Reference GFF3; required for direct-sequencing.
  --candidate-genomic-fasta CANDIDATE_GENOMIC_FASTA
                        Optional candidate genomic FASTA for tiled Sanger amplicons.
  --variants VARIANTS   Confirmed parental variant TSV for variant-marker stage.
  --vcf VCF             Optional parental VCF supplied to primerblast-oss.
  --interval INTERVAL   QTL interval for qtl-grid, CHROM:START-END.
  --n-markers N_MARKERS
  --spacing SPACING
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench primerblast-import

```text
usage: candidate-qtl-bench primerblast-import [-h] --results RESULTS --out-dir OUT_DIR [--candidates CANDIDATES]
                                              [--export {kasp,asp,hrm}] [--json]

options:
  -h, --help            show this help message and exit
  --results RESULTS     PrimerBlast OSS JSON report (design or CAPS/InDel/assay); repeatable.
  --out-dir OUT_DIR
  --candidates CANDIDATES
                        Optional candidate ledger to join gene IDs by gene_id/transcript/name.
  --export {kasp,asp,hrm}
                        Emit an allele-specific ordering scaffold; repeatable.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench promoter-review

```text
usage: candidate-qtl-bench promoter-review [-h] --genes GENES --gff3 GFF3 --out-dir OUT_DIR
                                           [--parental-variants PARENTAL_VARIANTS] [--sv-calls SV_CALLS]
                                           [--coverage COVERAGE] [--motif-hits MOTIF_HITS] [--hypotheses HYPOTHESES]
                                           [--promoter-window PROMOTER_WINDOW] [--case-id CASE_ID] [--force] [--json]

options:
  -h, --help            show this help message and exit
  --genes GENES         Target genes TSV (gene_id + optional representative_transcript).
  --gff3 GFF3           Reference GFF3 with gene/mRNA/exon/UTR features.
  --out-dir OUT_DIR
  --parental-variants PARENTAL_VARIANTS
                        Parental SNV/indel TSV: chrom/start/end (+ optional flags/feature_id/source_id).
  --sv-calls SV_CALLS   SV/PAV/TE call TSV: chrom/start/end (+ optional observation_class/flags/feature_id/source_id).
  --coverage COVERAGE   Coding exon coverage / Wet-state TSV; gates route start and records unexamined scope.
  --motif-hits MOTIF_HITS
                        Optional reference-only motif hit TSV; recorded as context, never causal.
  --hypotheses HYPOTHESES
                        Optional locus_hypotheses.tsv; matching typed observations are linked by observation ID.
  --promoter-window PROMOTER_WINDOW
                        Promoter window in bp upstream of TSS; repeatable for multiple scopes (default 2000).
  --case-id CASE_ID     Locus/case id recorded in the manifest.
  --force               Run the review scope even if coding coverage has not been supplied/examined.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench propagate-terms

```text
usage: candidate-qtl-bench propagate-terms [-h] --candidates CANDIDATES [--orthology ORTHOLOGY]
                                           [--reference-annotation REFERENCE_ANNOTATION] [--trait-terms TRAIT_TERMS]
                                           --out-dir OUT_DIR [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Interval candidates.tsv (gene_id universe).
  --orthology ORTHOLOGY
                        TSV mapping candidate_gene_id -> reference_gene_id (+ optional
                        relationship/source_genome/reference_role/reference_species/confidence).
  --reference-annotation REFERENCE_ANNOTATION
                        TSV of reference_gene_id -> term (+ optional ontology_id/ontology_name/source).
  --trait-terms TRAIT_TERMS
                        Optional analyst trait terms (one per line or a term column) for trait_relevance tagging.
  --out-dir OUT_DIR
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench provenance-ledger

```text
usage: candidate-qtl-bench provenance-ledger [-h] --entries ENTRIES --ledger-dir LEDGER_DIR --out-dir OUT_DIR
                                             [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --entries ENTRIES     TSV with node_kind, source_id, label, optional parent_node; optional
                        node_id/source_ref/edge_kind columns.
  --ledger-dir LEDGER_DIR
                        Directory receiving normalized nodes.tsv / edges.tsv.
  --out-dir OUT_DIR     Directory receiving provenance_graph.tsv and the manifest JSON.
  --case-id CASE_ID     Case identifier recorded into the manifest.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench qtl-source-init

```text
usage: candidate-qtl-bench qtl-source-init [-h] --manifest MANIFEST --out-dir OUT_DIR [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  YAML/JSON QTL source manifest.
  --out-dir OUT_DIR
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench qtl2-plan

```text
usage: candidate-qtl-bench qtl2-plan [-h] --genotypes GENOTYPES --phenotypes PHENOTYPES --markers MARKERS --trait
                                     TRAIT [--covariate COVARIATE] --out-dir OUT_DIR [--n-perm N_PERM] [--alpha ALPHA]
                                     [--error-prob ERROR_PROB] [--execute] [--rscript RSCRIPT] [--json]

options:
  -h, --help            show this help message and exit
  --genotypes GENOTYPES
                        Wide TSV: individual_id then AA/AB/BB marker columns.
  --phenotypes PHENOTYPES
                        TSV with individual_id, trait, and optional covariates.
  --markers MARKERS     TSV with marker_id, chrom, pos_bp, and optional pos_cm.
  --trait TRAIT
  --covariate COVARIATE
                        Phenotype covariate; repeatable.
  --out-dir OUT_DIR
  --n-perm N_PERM
  --alpha ALPHA
  --error-prob ERROR_PROB
  --execute             Run the generated script; requires R/qtl2.
  --rscript RSCRIPT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench rank

```text
usage: candidate-qtl-bench rank [-h] --ledger LEDGER --out OUT [--research-summary RESEARCH_SUMMARY]
                                [--go-evidence GO_EVIDENCE] [--json]

options:
  -h, --help            show this help message and exit
  --ledger LEDGER       TSV based on templates/candidate_ledger.tsv.
  --out OUT             Output ranked TSV.
  --research-summary RESEARCH_SUMMARY
                        Optional research_summary.tsv from candidate-qtl research.
  --go-evidence GO_EVIDENCE
                        Optional gene-level GO evidence TSV from annotate-plan.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench ranking-qa

```text
usage: candidate-qtl-bench ranking-qa [-h] --candidates CANDIDATES --evidence EVIDENCE --out-dir OUT_DIR
                                      [--weights WEIGHTS] [--top-k TOP_K] [--perturbation PERTURBATION]
                                      [--source-dropout SOURCE_DROPOUT] [--provider-dropout PROVIDER_DROPOUT]
                                      [--benchmark-mode] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --evidence EVIDENCE
  --out-dir OUT_DIR
  --weights WEIGHTS
  --top-k TOP_K
  --perturbation PERTURBATION
                        Fractional weight perturbation; repeatable (defaults: 0.10, 0.25, 0.50).
  --source-dropout SOURCE_DROPOUT
                        Resolved independence group to drop (e.g. comparative_annotation_chain, PMID:12345);
                        repeatable. Defaults to every non-empty group found in the evidence.
  --provider-dropout PROVIDER_DROPOUT
                        source_type to drop (e.g. blastp, eggnog, qtl_position); repeatable. Defaults to every non-
                        empty source_type found in the evidence.
  --benchmark-mode      Mark this run as calibration/benchmark rather than a normal case.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench ranking-stability-benchmark

```text
usage: candidate-qtl-bench ranking-stability-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR [--top-k TOP_K]
                                                       [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  JSON manifest with cases: case_id, candidates, evidence, ranked_tsv, answer_gene_id.
  --out-dir OUT_DIR
  --top-k TOP_K
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-candidate-table

```text
usage: candidate-qtl-bench reference-candidate-table [-h] --candidates CANDIDATES --reference-genes REFERENCE_GENES
                                                     --orthology ORTHOLOGY [--annotation ANNOTATION] --out-dir OUT_DIR
                                                     [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate or ranked TSV with gene_id.
  --reference-genes REFERENCE_GENES
                        Source-backed TSV following templates/reference_gene_panel.tsv.
  --orthology ORTHOLOGY
                        PROVIDER=TSV following templates/candidate_reference_orthology.tsv; repeatable.
  --annotation ANNOTATION
                        Optional TSV with gene_id and annotation columns.
  --out-dir OUT_DIR
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-corpus-build

```text
usage: candidate-qtl-bench reference-corpus-build [-h] [--from-uniprot UNIPROT] [--from-blastdb BLASTDB]
                                                  [--from-dump DUMP] [--from-payload PAYLOAD]
                                                  [--orthogroup-map ORTHOGROUP_MAP] --out-dir OUT_DIR --corpus-version
                                                  CORPUS_VERSION [--citation-cache CITATION_CACHE] [--benchmark-mode]
                                                  [--preset-acceptance PRESET_ACCEPTANCE] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --from-uniprot UNIPROT
                        Downloaded UniProtKB table from reference-corpus-fetch (or an equivalent export): curated
                        function, disruption phenotype, pathway, tissue/stage, GO and evidence codes. Repeatable.
  --from-blastdb BLASTDB
                        Local protein BLAST database prefix; every entry is read with blastdbcmd and its UniProt-style
                        defline parsed. Repeatable.
  --from-dump DUMP      Exported reference table (TSV with reference_gene_id ...). Repeatable.
  --from-payload PAYLOAD
                        Model/agent-proposed rows (.tsv/.json/.jsonl) adding function and mutant-phenotype text.
                        Repeatable.
  --orthogroup-map ORTHOGROUP_MAP
                        TSV mapping reference_gene_id to orthogroup_id/family_id from a source-backed orthology
                        resource. Repeatable.
  --out-dir OUT_DIR
  --corpus-version CORPUS_VERSION
                        Version label; the corpus is frozen under it and referenced by that name downstream.
  --citation-cache CITATION_CACHE
                        Cached Europe PMC records used to resolve the citations of proposed rows offline.
  --benchmark-mode      Bar focal-lineage / answer-shaped rows and record that the corpus was built for a benchmark.
  --preset-acceptance PRESET_ACCEPTANCE
                        Preset id from the packaged refcorpus_presets.json: the build verifies the preset's required
                        gene symbols are present in the corpus as usable gene symbols (casefold; fail-closed: one
                        missing required term aborts before anything is written, so a previously frozen corpus cannot
                        be silently replaced) and records the per-term result in the manifest. Optional terms are
                        recorded only.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-corpus-compare

```text
usage: candidate-qtl-bench reference-corpus-compare [-h] --corpus CORPUS --out-dir OUT_DIR
                                                    [--baseline-label BASELINE_LABEL]
                                                    [--retrieval-manifest RETRIEVAL_MANIFEST] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --corpus CORPUS       LABEL=PATH for a reference corpus TSV. Repeat at least twice.
  --out-dir OUT_DIR
  --baseline-label BASELINE_LABEL
                        Label the deltas are measured against; defaults to the first --corpus.
  --retrieval-manifest RETRIEVAL_MANIFEST
                        Optional TSV (corpus_label/case_id/retrieval_tsv/truth_family_id[/truth_gene_id]) of frozen
                        Route C retrievals produced from each corpus for the same cases.
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-corpus-fetch

```text
usage: candidate-qtl-bench reference-corpus-fetch [-h] [--taxon TAXON] [--preset PRESET] --out-dir OUT_DIR
                                                  [--allow-external-api] [--include-unreviewed] [--max-rows MAX_ROWS]
                                                  [--page-size PAGE_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --taxon TAXON         NCBI taxon id (3702 Arabidopsis, 3880 Medicago, 3847 soybean, 4530 rice, 4577 maize, 4081
                        tomato ...). Repeatable; may be combined with --preset (union; an explicit --taxon wins over
                        the same id from the preset).
  --preset PRESET       Versioned corpus preset from the packaged refcorpus_presets.json (e.g. legume_v1): expands to
                        its taxa and per-taxon reviewed flags, so a taxon that needs TrEMBL entries is fetched with
                        them. May be combined with --taxon.
  --out-dir OUT_DIR
  --allow-external-api  Required. Only the taxon id and the field list are transmitted.
  --include-unreviewed  Include unreviewed (TrEMBL) entries; reviewed-only by default because curated comment blocks
                        are what the routes need.
  --max-rows MAX_ROWS   Stop after this many rows per taxon (0 = all).
  --page-size PAGE_SIZE
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-gene-cache

```text
usage: candidate-qtl-bench reference-gene-cache [-h] --mechanisms MECHANISMS --out OUT [--allow-external-api]
                                                [--timeout TIMEOUT] [--max-queries MAX_QUERIES]
                                                [--max-records MAX_RECORDS]
                                                [--max-elapsed-seconds MAX_ELAPSED_SECONDS]
                                                [--resume-manifest RESUME_MANIFEST] [--json]

options:
  -h, --help            show this help message and exit
  --mechanisms MECHANISMS
  --out OUT
  --allow-external-api
  --timeout TIMEOUT
  --max-queries MAX_QUERIES
  --max-records MAX_RECORDS
  --max-elapsed-seconds MAX_ELAPSED_SECONDS
  --resume-manifest RESUME_MANIFEST
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-gene-resolve

```text
usage: candidate-qtl-bench reference-gene-resolve [-h] [--mechanisms MECHANISMS] --out OUT [--allow-external-api]
                                                  [--records RECORDS] [--timeout TIMEOUT]
                                                  [--focal-species FOCAL_SPECIES] [--identity-review IDENTITY_REVIEW]
                                                  [--sequence-requests SEQUENCE_REQUESTS] [--post-freeze]
                                                  [--proteome PROTEOME] [--protein-gene-map PROTEIN_GENE_MAP]
                                                  [--runner {auto,native,wsl}] [--threads THREADS] [--json]

options:
  -h, --help            show this help message and exit
  --mechanisms MECHANISMS
  --out OUT
  --allow-external-api
  --records RECORDS
  --timeout TIMEOUT
  --focal-species FOCAL_SPECIES
  --identity-review IDENTITY_REVIEW
                        Identity review ledger TSV (reference_gene_name, reference_species, accession,
                        record_gene_name, verdict, evidence_source_ids, reviewer, rationale). A same_gene row with a
                        citation confirms a provisional declared accession whose record carries record_gene_name.
  --sequence-requests SEQUENCE_REQUESTS
                        Post-freeze public INSDC accession requests; separate from the blind mechanism resolver.
  --post-freeze
  --proteome PROTEOME
  --protein-gene-map PROTEIN_GENE_MAP
  --runner {auto,native,wsl}
  --threads THREADS
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reference-panel

```text
usage: candidate-qtl-bench reference-panel [-h] --panel PANEL --out-dir OUT_DIR [--full-hash] [--json]

options:
  -h, --help         show this help message and exit
  --panel PANEL      TSV declaring resource_id/role/resource_type/taxon/source/uri_or_path per resource.
  --out-dir OUT_DIR
  --full-hash
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench research

```text
usage: candidate-qtl-bench research [-h] --candidates CANDIDATES --out-dir OUT_DIR --species SPECIES --trait TRAIT
                                    [--top-n TOP_N] [--literature-limit LITERATURE_LIMIT] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate or ranked TSV.
  --out-dir OUT_DIR     Output directory.
  --species SPECIES
  --trait TRAIT
  --top-n TOP_N
  --literature-limit LITERATURE_LIMIT
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench research-packet

```text
usage: candidate-qtl-bench research-packet [-h] --candidates CANDIDATES --out-dir OUT_DIR --species SPECIES --trait
                                           TRAIT [--top-n TOP_N] [--fetch-literature]
                                           [--literature-limit LITERATURE_LIMIT] [--evidence EVIDENCE]
                                           [--contributions CONTRIBUTIONS] [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Candidate or ranked TSV.
  --out-dir OUT_DIR     Output directory.
  --species SPECIES
  --trait TRAIT
  --top-n TOP_N
  --fetch-literature    Read Europe PMC for each selected candidate and cache JSON results.
  --literature-limit LITERATURE_LIMIT
  --evidence EVIDENCE   Canonical evidence TSV to embed.
  --contributions CONTRIBUTIONS
                        Score contributions TSV to embed.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-hybrid-merge

```text
usage: candidate-qtl-bench retrieval-hybrid-merge [-h] --retrieval RETRIEVAL --out-dir OUT_DIR [--top-k TOP_K]
                                                  [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --retrieval RETRIEVAL
                        Per-provider retrieval TSV. Repeat at least twice; each input must carry exactly one
                        model_name/model_version pair (one provider).
  --out-dir OUT_DIR
  --top-k TOP_K
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-import

```text
usage: candidate-qtl-bench retrieval-import [-h] --retrieval RETRIEVAL --out-dir OUT_DIR [--phenotype PHENOTYPE]
                                            [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --retrieval RETRIEVAL
                        External retrieval TSV following templates/route_c_retrieval.tsv (reference_gene_id/reference_
                        species/family_id|orthogroup_id/retrieval_score/matched_modality/matched_description_or_phenot
                        ype/model_name/model_version/database_version/rank). Columns outside the contract are
                        rejected.
  --out-dir OUT_DIR
  --phenotype PHENOTYPE
                        Phenotype text recorded in the status; not sent to any model by this command.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-index-build

```text
usage: candidate-qtl-bench retrieval-index-build [-h] --corpus CORPUS --out-dir OUT_DIR --database-version
                                                 DATABASE_VERSION [--corpus-label CORPUS_LABEL]
                                                 [--provider {lexical,lsa,neural}] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --corpus CORPUS       Reference corpus TSV (reference_gene_id/reference_species/source_id + family_id and/or
                        orthogroup_id + description/function_description/mutant_phenotype/pathway/go_terms/keywords).
                        Repeatable.
  --out-dir OUT_DIR
  --database-version DATABASE_VERSION
                        Version label for this index; retrieval scores are only comparable within one
                        database_version.
  --corpus-label CORPUS_LABEL
                        Human-readable corpus name recorded in the manifest.
  --provider {lexical,lsa,neural}
                        Retrieval provider (issue #225): 'lexical' (default, deterministic TF-IDF, issue #217), 'lsa'
                        (latent-semantic: TF-IDF + truncated SVD, no new dependency), or 'neural'
                        (sentence_transformers adapter; reports not_evaluated_provider_unavailable when the backend is
                        not installed, as it is not on this workstation).
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-interval-match

```text
usage: candidate-qtl-bench retrieval-interval-match [-h] --retrieval RETRIEVAL --atlas ATLAS --out-dir OUT_DIR
                                                    [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --retrieval RETRIEVAL
                        route_c_retrieval_frozen.tsv produced by retrieval-import.
  --atlas ATLAS         interval_gene_atlas.tsv (gene_id/family_id/orthogroup/orthology_reference_gene): the canonical
                        <out-dir>/interval_gene_atlas.tsv from route-b-prepare or the root-level interval-atlas file.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-merge

```text
usage: candidate-qtl-bench retrieval-merge [-h] --retrieval RETRIEVAL --out-dir OUT_DIR [--top-k TOP_K]
                                           [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --retrieval RETRIEVAL
                        Per-database retrieval TSV. Repeat at least twice; each input must carry one unique
                        database_version.
  --out-dir OUT_DIR
  --top-k TOP_K
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench retrieval-run

```text
usage: candidate-qtl-bench retrieval-run [-h] --index INDEX --phenotype PHENOTYPE --out-dir OUT_DIR [--top-k TOP_K]
                                         [--case-id CASE_ID] [--phenotype-context PHENOTYPE_CONTEXT]
                                         [--provider {lexical,lsa,neural}] [--json]

options:
  -h, --help            show this help message and exit
  --index INDEX         index_manifest.json (or the index directory) from retrieval-index-build. Its artifacts are re-
                        hashed before retrieval.
  --phenotype PHENOTYPE
                        Phenotype text. This is the only query input; it never leaves the machine.
  --out-dir OUT_DIR
  --top-k TOP_K
  --case-id CASE_ID
  --phenotype-context PHENOTYPE_CONTEXT
                        Optional phenotype_context.yaml (case_id, developmental_stage, tissue,
                        environment_or_treatment, ontology term ...). Recorded for audit in the retrieval status; the
                        retrieval query stays the phenotype free text and no structured field is appended to it.
                        Supported for the default lexical provider.
  --provider {lexical,lsa,neural}
                        Must match the provider the --index was built with ('lexical' default / 'lsa' / 'neural').
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench reverse-literature

```text
usage: candidate-qtl-bench reverse-literature [-h] --lookup LOOKUP --out-dir OUT_DIR [--phenotype PHENOTYPE]
                                              [--cache CACHE] [--allow-external-api] [--limit LIMIT]
                                              [--case-id CASE_ID] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --lookup LOOKUP       interval_reverse_lookup.tsv from route-b-prepare
                        (gene_id/reference_gene[/reference_function]).
  --out-dir OUT_DIR
  --phenotype PHENOTYPE
                        Focuses each query. Only generic phenotype words are transmitted.
  --cache CACHE         Query-keyed cached Europe PMC results, so a rerun is byte-identical offline.
  --allow-external-api  Permit read-only Europe PMC queries.
  --limit LIMIT         Articles kept per reference gene.
  --case-id CASE_ID
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-a-interval-match

```text
usage: candidate-qtl-bench route-a-interval-match [-h] --panel PANEL --hypotheses HYPOTHESES --out-dir OUT_DIR
                                                  [--family-index FAMILY_INDEX] [--interval INTERVAL]
                                                  [--candidates CANDIDATES] [--annotation-rescue ANNOTATION_RESCUE]
                                                  [--tandem-distance TANDEM_DISTANCE] [--max-q MAX_Q]
                                                  [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --panel PANEL         verified_bait_panel.tsv from bait-panel-freeze; only its hypothesis_ids are matched.
  --hypotheses HYPOTHESES
                        llm_pathway_hypotheses.tsv from mechanism-import.
  --out-dir OUT_DIR
  --family-index FAMILY_INDEX
                        reference_genome-scoped family_index.tsv; absent ->
                        not_evaluated_missing_reference_genome_family_index.
  --interval INTERVAL   CHROM:START-END resolved against the index coordinates; repeatable.
  --candidates CANDIDATES
                        Interval candidates TSV, as an alternative to --interval.
  --annotation-rescue ANNOTATION_RESCUE
                        annotation_rescue.tsv; a rescue provider that did not complete holds absence rather than
                        falsifying it.
  --tandem-distance TANDEM_DISTANCE
  --max-q MAX_Q
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-prepare

```text
usage: candidate-qtl-bench route-b-prepare [-h] --candidates CANDIDATES --out-dir OUT_DIR [--proteins PROTEINS]
                                           [--protein-db PROTEIN_DB] [--reference-function REFERENCE_FUNCTION]
                                           [--mutant-functional-evidence MUTANT_FUNCTIONAL_EVIDENCE]
                                           [--cross-accession-map CROSS_ACCESSION_MAP] [--orthology ORTHOLOGY]
                                           [--domains DOMAINS] [--family-index FAMILY_INDEX]
                                           [--annotations ANNOTATIONS] [--case-id CASE_ID] [--threads THREADS]
                                           [--min-identity MIN_IDENTITY] [--min-coverage MIN_COVERAGE]
                                           [--max-evalue MAX_EVALUE] [--no-auto-local-db] [--phenotype PHENOTYPE]
                                           [--json]

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
                        Interval candidate table (gene_id ...), e.g. candidates.tsv from extract.
  --out-dir OUT_DIR
  --proteins PROTEINS   Candidate protein FASTA (candidate_proteins.fa) for the local homology search. When --family-
                        index is absent, the same standard 13-column BLAST output generates an explicitly interval-
                        scoped fallback family index.
  --protein-db PROTEIN_DB
                        Local BLASTP database used when interval reference/annotation material is incomplete; defaults
                        to CANDIDATE_QTL_PROTEIN_DB. Complete usable material skips the search.
  --reference-function REFERENCE_FUNCTION
                        Local reference TSV or directory of UniProt JSON records; repeatable. Annotation and citation
                        links remain investigation material, not verified claims. Defaults to
                        CANDIDATE_QTL_REFERENCE_FUNCTION.
  --mutant-functional-evidence MUTANT_FUNCTIONAL_EVIDENCE
                        normalized_functional_evidence.tsv from tnt1-evidence-normalize or the LORE1 adapter
                        (accession-native FUNCTION_COLUMNS rows). Projected onto the reference-function contract via
                        --cross-accession-map; repeatable. Rows without published/genotyped mutant support are skipped
                        and counted, never turned into negatives.
  --cross-accession-map CROSS_ACCESSION_MAP
                        Explicit cross-accession mapping table
                        (source_full_id/target_full_id/mapping_basis/relationship/mapping_state ...). mapping_basis
                        must be machine-checkable
                        (explicit_id_map/pangene_family/synteny/reciprocal_sequence/composite); symbol-only identity
                        transfer is rejected. Repeatable; without it every emitted row keeps its native ID with
                        mapping_state=not_evaluated.
  --orthology ORTHOLOGY
                        Pre-computed orthology table. Usable interval mappings are retained; local homology search
                        fills missing material when a protein query and database are available. Empty/off-interval
                        rows do not suppress fallback.
  --domains DOMAINS     Domain table (gene_id/domain_id[/domain_name/source/start]).
  --family-index FAMILY_INDEX
                        Reusable genome-wide family_index.tsv from family-index. Preferred over the interval-scoped
                        fallback generated from --proteins + --protein-db.
  --annotations ANNOTATIONS
                        Per-gene annotation table (gene_id/description/go_terms...).
  --case-id CASE_ID
  --threads THREADS
  --min-identity MIN_IDENTITY
                        Minimum percent identity for a usable local homology hit.
  --min-coverage MIN_COVERAGE
                        Minimum query coverage (0-1) for a usable local homology hit.
  --max-evalue MAX_EVALUE
  --no-auto-local-db    Disable discovery of the configured local protein DB.
  --phenotype PHENOTYPE
                        Marks which entries of the (phenotype-independent) interval inventory look related. The atlas
                        is always built without it, so enumeration precedes filtering.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-semantic-bundle

```text
usage: candidate-qtl-bench route-b-semantic-bundle [-h] [--atlas ATLAS] [--phenotype PHENOTYPE]
                                                   [--extend-bundle EXTEND_BUNDLE]
                                                   [--parent-bundle-sha256 PARENT_BUNDLE_SHA256]
                                                   [--research-discoveries RESEARCH_DISCOVERIES]
                                                   [--focal-species FOCAL_SPECIES]
                                                   [--b-review-membership-policy {initial_bibliography_v1,source_function_v1}]
                                                   [--literature LITERATURE] [--literature-status LITERATURE_STATUS]
                                                   [--lookup LOOKUP] [--identity-hits IDENTITY_HITS] --out-dir OUT_DIR
                                                   [--case-id CASE_ID] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --atlas ATLAS         Frozen interval_gene_atlas.tsv; required for an initial bundle.
  --phenotype PHENOTYPE
                        Phenotype introduced only after B0; required for an initial bundle.
  --extend-bundle EXTEND_BUNDLE
                        Frozen parent bundle to extend with newly retrieved sources; never edits the parent.
  --parent-bundle-sha256 PARENT_BUNDLE_SHA256
                        Coordinator-pinned parent hash for --extend-bundle.
  --research-discoveries RESEARCH_DISCOVERIES
                        Worker JSON with case_id, parent_bundle_sha256 and discovered_sources; not claim responses.
  --focal-species FOCAL_SPECIES
  --b-review-membership-policy {initial_bibliography_v1,source_function_v1}
  --literature LITERATURE
                        Optional reverse_literature.tsv link candidates; retrieval is not evidence.
  --literature-status LITERATURE_STATUS
                        Optional reverse_literature_status.tsv preserving no-hit/not-evaluated states.
  --lookup LOOKUP       Optional interval_reverse_lookup.tsv carrying source-backed B0 function/pathway/source fields
                        into the frozen semantic unit.
  --identity-hits IDENTITY_HITS
                        Optional route_b_identity_hits.tsv preserving all B0 sequence-similarity retrieval handles for
                        each interval gene; these are not orthology or trait support.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-semantic-candidates

```text
usage: candidate-qtl-bench route-b-semantic-candidates [-h] --claims CLAIMS --out-dir OUT_DIR [--case-id CASE_ID]
                                                       [--json]

options:
  -h, --help         show this help message and exit
  --claims CLAIMS    semantic_claims.tsv (or semantic_claims_frozen.tsv).
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-semantic-import

```text
usage: candidate-qtl-bench route-b-semantic-import [-h] --bundle BUNDLE --prompts PROMPTS --responses RESPONSES
                                                   [--llm-manifest LLM_MANIFEST] --out-dir OUT_DIR [--case-id CASE_ID]
                                                   [--json]

options:
  -h, --help            show this help message and exit
  --bundle BUNDLE       route_b_semantic_bundle.json
  --prompts PROMPTS     route_b_semantic_prompts.jsonl used by llm-run
  --responses RESPONSES
                        llm_candidate_responses.json from llm-run
  --llm-manifest LLM_MANIFEST
                        Optional llm_run_manifest.json supplying provider/model metadata.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-semantic-source-verify

```text
usage: candidate-qtl-bench route-b-semantic-source-verify [-h] --proposals PROPOSALS [--bundle BUNDLE] [--cache CACHE]
                                                          [--allow-external-api] --out-dir OUT_DIR [--case-id CASE_ID]
                                                          [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --proposals PROPOSALS
                        semantic_claim_proposals.tsv from route-b-semantic-import.
  --bundle BUNDLE       Optional frozen route_b_semantic_bundle.json. Enables source-backed curated-function
                        verification for database accessions grounded in the exact B0 unit; this is weaker than
                        mechanistic literature and does not verify transfer.
  --cache CACHE         Optional local literature-record cache used by source-verify.
  --allow-external-api  Allow source-verify to resolve citation identifiers through Europe PMC with an official NCBI
                        PubMed/PMC E-utilities fallback. Only citation identifiers are transmitted.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-b-semantic-transfer-verify

```text
usage: candidate-qtl-bench route-b-semantic-transfer-verify [-h] --proposals PROPOSALS --source-verification
                                                            SOURCE_VERIFICATION --bundle BUNDLE
                                                            [--bound-claims BOUND_CLAIMS]
                                                            [--agent-claim-review AGENT_CLAIM_REVIEW]
                                                            [--source-review-records SOURCE_REVIEW_RECORDS]
                                                            [--record-cache RECORD_CACHE] --out-dir OUT_DIR
                                                            [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --proposals PROPOSALS
                        semantic_claim_proposals.tsv
  --source-verification SOURCE_VERIFICATION
                        semantic_claim_source_verification.tsv
  --bundle BUNDLE       route_b_semantic_bundle.json
  --bound-claims BOUND_CLAIMS
                        Shared bound-claim TSV. Without a complete matching row, source/cue evidence cannot authorize
                        canonical membership.
  --agent-claim-review AGENT_CLAIM_REVIEW
                        Independent Route-B review packet for the exact semantic proposal/source snapshot.
  --source-review-records SOURCE_REVIEW_RECORDS
                        source_review_records.json matching the exact Route-B proposal citations.
  --record-cache RECORD_CACHE, --records RECORD_CACHE
                        Frozen UniProt-shaped provider records for exact accession/species checks.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-benchmark

```text
usage: candidate-qtl-bench route-c-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR [--full-hash] [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  TSV: one row per case (case_id/retrieval_tsv/truth_family_id[/truth_gene_id/route_a_tsv]).
  --out-dir OUT_DIR
  --full-hash
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-direct-benchmark

```text
usage: candidate-qtl-bench route-c-direct-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR [--full-hash] [--json]

options:
  -h, --help           show this help message and exit
  --manifest MANIFEST  TSV following templates/route_c_direct_benchmark.tsv; truth is opened only after freeze.
  --out-dir OUT_DIR
  --full-hash
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-direct-freeze

```text
usage: candidate-qtl-bench route-c-direct-freeze [-h] --proposals PROPOSALS --out-dir OUT_DIR
                                                 [--citation-cache CITATION_CACHE] [--allow-external-api]
                                                 [--claim-extraction CLAIM_EXTRACTION] [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --proposals PROPOSALS
  --out-dir OUT_DIR
  --citation-cache CITATION_CACHE
  --allow-external-api
  --claim-extraction CLAIM_EXTRACTION
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-direct-members

```text
usage: candidate-qtl-bench route-c-direct-members [-h] --frozen FROZEN --freeze-manifest FREEZE_MANIFEST --projection
                                                  PROJECTION --projection-manifest PROJECTION_MANIFEST --atlas ATLAS
                                                  --out-dir OUT_DIR [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --frozen FROZEN
  --freeze-manifest FREEZE_MANIFEST
  --projection PROJECTION
  --projection-manifest PROJECTION_MANIFEST
  --atlas ATLAS
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-direct-project

```text
usage: candidate-qtl-bench route-c-direct-project [-h] --frozen FROZEN --freeze-manifest FREEZE_MANIFEST --atlas ATLAS
                                                  --qtl-id QTL_ID [--eligible-hits ELIGIBLE_HITS]
                                                  [--identity-mapping IDENTITY_MAPPING]
                                                  [--reference-bridge REFERENCE_BRIDGE] --out-dir OUT_DIR
                                                  [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --frozen FROZEN
  --freeze-manifest FREEZE_MANIFEST
  --atlas ATLAS
  --qtl-id QTL_ID
  --eligible-hits ELIGIBLE_HITS
  --identity-mapping IDENTITY_MAPPING
  --reference-bridge REFERENCE_BRIDGE
                        Hash-verified public sequence bridge; keeps all matched copies review-only.
  --out-dir OUT_DIR
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-direct-search

```text
usage: candidate-qtl-bench route-c-direct-search [-h] --phenotype PHENOTYPE --focal-species FOCAL_SPECIES
                                                 --species-scope SPECIES_SCOPE --out-dir OUT_DIR [--cache CACHE]
                                                 [--allow-external-api] [--limit LIMIT] [--case-id CASE_ID]
                                                 [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --phenotype PHENOTYPE
  --focal-species FOCAL_SPECIES
  --species-scope SPECIES_SCOPE
  --out-dir OUT_DIR
  --cache CACHE
  --allow-external-api
  --limit LIMIT
  --case-id CASE_ID
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-provider-bait

```text
usage: candidate-qtl-bench route-c-provider-bait [-h] --frozen FROZEN --freeze-manifest FREEZE_MANIFEST --out-dir
                                                 OUT_DIR [--cache CACHE] [--allow-external-api] [--json]

options:
  -h, --help            show this help message and exit
  --frozen FROZEN
  --freeze-manifest FREEZE_MANIFEST
  --out-dir OUT_DIR
  --cache CACHE         Directory of cached UniProt literature-search JSON (uniprot_lit_<PMID>.json).
  --allow-external-api  Query UniProt by public PubMed ID when the cache lacks it.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench route-c-provider-comparison

```text
usage: candidate-qtl-bench route-c-provider-comparison [-h] --provider-manifest PROVIDER_MANIFEST --out-dir OUT_DIR
                                                       [--full-hash] [--json]

options:
  -h, --help            show this help message and exit
  --provider-manifest PROVIDER_MANIFEST
                        LABEL=PATH for a route-c-benchmark manifest whose retrieval_tsv values are that provider's own
                        retrieval for the same cases. Repeat at least twice.
  --out-dir OUT_DIR
  --full-hash
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench scaffold-rescue

```text
usage: candidate-qtl-bench scaffold-rescue [-h] --scaffolds SCAFFOLDS --out-dir OUT_DIR [--marker-hits MARKER_HITS]
                                           [--synteny SYNTENY] [--orthogroup-context ORTHOGROUP_CONTEXT]
                                           [--placement-support PLACEMENT_SUPPORT] [--scaffold-genes SCAFFOLD_GENES]
                                           [--max-marker-scaffolds MAX_MARKER_SCAFFOLDS] [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --scaffolds SCAFFOLDS
                        Scaffold table (scaffold_id[/length]).
  --out-dir OUT_DIR
  --marker-hits MARKER_HITS
                        Marker/flanking alignment (marker_id/scaffold_id[/pair_id/identity]); paired markers beat lone
                        hits and multi-mapping markers are rejected.
  --synteny SYNTENY     Synteny bridge (scaffold_id[/block_id/anchor_species]).
  --orthogroup-context ORTHOGROUP_CONTEXT
                        Gene-content context (scaffold_id[/orthogroup]); weak support only.
  --placement-support PLACEMENT_SUPPORT
                        Hi-C / linkage map / optical map / assembler placement (scaffold_id/support_type[/source]).
  --scaffold-genes SCAFFOLD_GENES
                        Genes on scaffolds (gene_id/scaffold_id).
  --max-marker-scaffolds MAX_MARKER_SCAFFOLDS
                        A marker hitting more scaffolds than this is treated as a repeat.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench scenario-create

```text
usage: candidate-qtl-bench scenario-create [-h] --input INPUT [--case CASE] --out-dir OUT_DIR [--json]

options:
  -h, --help         show this help message and exit
  --input INPUT      Scenario rows as TSV or YAML (list of mappings).
  --case CASE        Optional case_id overriding/completing every row's case_id.
  --out-dir OUT_DIR  Directory receiving locus_causal_scenarios.tsv.
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench scenario-report

```text
usage: candidate-qtl-bench scenario-report [-h] --ledger LEDGER [--json]

options:
  -h, --help       show this help message and exit
  --ledger LEDGER  locus_causal_scenarios.tsv from scenario-create.
  --json           Emit a stable JSON result to stdout.
```

## candidate-qtl-bench screening-shortlist

```text
usage: candidate-qtl-bench screening-shortlist [-h] --ranked RANKED --out-dir OUT_DIR [--target-size TARGET_SIZE]
                                               [--rescue-size RESCUE_SIZE] [--json]

options:
  -h, --help            show this help message and exit
  --ranked RANKED       ranked_candidates.tsv from evidence-rank or candidate-pipeline.
  --out-dir OUT_DIR
  --target-size TARGET_SIZE
                        Primary deterministic review budget.
  --rescue-size RESCUE_SIZE
                        Separate Rescue appendix budget; not counted in the primary shortlist.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench semantic-claims-freeze

```text
usage: candidate-qtl-bench semantic-claims-freeze [-h] --claims CLAIMS --out-dir OUT_DIR [--case-id CASE_ID]
                                                  [--route ROUTE] [--json]

options:
  -h, --help         show this help message and exit
  --claims CLAIMS    Structured semantic claims TSV following candidate_qtl.semantic_claims.CLAIM_COLUMNS.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --route ROUTE      Expected route label for every claim (default interval_first).
  --json             Emit a stable JSON result to stdout.
```

## candidate-qtl-bench source-verify

```text
usage: candidate-qtl-bench source-verify [-h] --mechanisms MECHANISMS --out-dir OUT_DIR [--cache CACHE]
                                         [--allow-external-api] [--id-column ID_COLUMN]
                                         [--source-column SOURCE_COLUMN]
                                         [--verification-policy {initial_bibliography_v1,content_v1}]
                                         [--focal-species FOCAL_SPECIES] [--claim-extraction CLAIM_EXTRACTION]
                                         [--full-hash] [--reference-bridge REFERENCE_BRIDGE] [--json]

options:
  -h, --help            show this help message and exit
  --mechanisms MECHANISMS
  --out-dir OUT_DIR
  --cache CACHE
  --allow-external-api
  --id-column ID_COLUMN
  --source-column SOURCE_COLUMN
  --verification-policy {initial_bibliography_v1,content_v1}
                        Default: ID/title metadata checks for the retained evidence papers. Explicit content_v1
                        preserves strict content verification.
  --focal-species FOCAL_SPECIES
  --claim-extraction CLAIM_EXTRACTION
  --full-hash
  --reference-bridge REFERENCE_BRIDGE
                        Hash-verified sequence correspondence for post-freeze parent review; never rewrites paper
                        facts or promotes support.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench supported-member-union

```text
usage: candidate-qtl-bench supported-member-union [-h] [--route-a-bridge ROUTE_A_BRIDGE]
                                                  [--route-b-claims ROUTE_B_CLAIMS]
                                                  [--route-b-candidates ROUTE_B_CANDIDATES]
                                                  [--route-c-members ROUTE_C_MEMBERS]
                                                  [--route-a-manifest ROUTE_A_MANIFEST]
                                                  [--route-b-claims-manifest ROUTE_B_CLAIMS_MANIFEST]
                                                  [--route-b-candidates-manifest ROUTE_B_CANDIDATES_MANIFEST]
                                                  [--route-c-manifest ROUTE_C_MANIFEST] --out-dir OUT_DIR
                                                  [--case-id CASE_ID] [--json]

options:
  -h, --help            show this help message and exit
  --route-a-bridge ROUTE_A_BRIDGE
                        Frozen family_candidate_bridge.tsv; only route=a rows that are not tier6/not_evaluated
                        contribute.
  --route-b-claims ROUTE_B_CLAIMS
                        Final semantic claims TSV; only verification_status=verified rows contribute all
                        member_interval_genes.
  --route-b-candidates ROUTE_B_CANDIDATES
                        Frozen route_b_verified_candidates.tsv; only route=b and
                        support_status=verified_semantic_support rows contribute.
  --route-c-members ROUTE_C_MEMBERS
                        Frozen Direct-C member ledger; only support_status=verified_direct_support rows contribute.
  --route-a-manifest ROUTE_A_MANIFEST
                        Producer manifest for Route A bridge (required for canonical input).
  --route-b-claims-manifest ROUTE_B_CLAIMS_MANIFEST
                        Producer manifest for Route B claims.
  --route-b-candidates-manifest ROUTE_B_CANDIDATES_MANIFEST
                        Producer manifest for Route B frozen candidates.
  --route-c-manifest ROUTE_C_MANIFEST
                        route_c_direct_members_manifest.json.
  --out-dir OUT_DIR
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench sv-plan

```text
usage: candidate-qtl-bench sv-plan [-h] --reference REFERENCE --query QUERY --qtl-interval QTL_INTERVAL --out-dir
                                   OUT_DIR [--candidates CANDIDATES] [--syri-out SYRI_OUT] [--call-table CALL_TABLE]
                                   [--breakpoint-tolerance BREAKPOINT_TOLERANCE] [--decision DECISION] [--json]

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
  --query QUERY
  --qtl-interval QTL_INTERVAL
  --out-dir OUT_DIR
  --candidates CANDIDATES
  --syri-out SYRI_OUT
  --call-table CALL_TABLE
                        PROVIDER=canonical.tsv; repeatable
  --breakpoint-tolerance BREAKPOINT_TOLERANCE
  --decision DECISION   candidate_path_decision.json; must request broaden_to_regulatory_sv_te_cnv
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench synteny-import

```text
usage: candidate-qtl-bench synteny-import [-h] --anchors ANCHORS --candidates CANDIDATES --out OUT
                                          [--format {anchors,collinearity}] [--json]

options:
  -h, --help            show this help message and exit
  --anchors ANCHORS     jcvi .anchors or MCScanX .collinearity file.
  --candidates CANDIDATES
  --out OUT
  --format {anchors,collinearity}
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench tnt1-evidence-normalize

```text
usage: candidate-qtl-bench tnt1-evidence-normalize [-h] --evidence EVIDENCE --out-dir OUT_DIR [--registry REGISTRY]
                                                   [--case-id CASE_ID] [--json]

options:
  -h, --help           show this help message and exit
  --evidence EVIDENCE  Curated Tnt1 TSV following templates/tnt1_evidence.tsv.
  --out-dir OUT_DIR
  --registry REGISTRY  Optional reference-resource registry used only to label exact native namespaces as
                       registered/unregistered; it does not infer identity.
  --case-id CASE_ID
  --json               Emit a stable JSON result to stdout.
```

## candidate-qtl-bench trait-ontology-map

```text
usage: candidate-qtl-bench trait-ontology-map [-h] --trait-terms TRAIT_TERMS --ontology-map ONTOLOGY_MAP --out OUT
                                              [--json]

options:
  -h, --help            show this help message and exit
  --trait-terms TRAIT_TERMS
                        TSV following templates/trait_terms.tsv.
  --ontology-map ONTOLOGY_MAP
                        Local TSV with term, ontology_id, and optional ontology_name/match_type/confidence.
  --out OUT             Normalized trait-ontology TSV.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench trait-terms

```text
usage: candidate-qtl-bench trait-terms [-h] {build} ...

positional arguments:
  {build}
    build     Classify, approve, and export trait terms; reject answer-leakage in benchmark mode.

options:
  -h, --help  show this help message and exit
```

## candidate-qtl-bench triple-route-compare

```text
usage: candidate-qtl-bench triple-route-compare [-h] --route-a ROUTE_A --route-b ROUTE_B --route-c ROUTE_C --out-dir
                                                OUT_DIR [--case-id CASE_ID]
                                                [--route-c-interval-match ROUTE_C_INTERVAL_MATCH] [--json]

options:
  -h, --help            show this help message and exit
  --route-a ROUTE_A
  --route-b ROUTE_B
  --route-c ROUTE_C
  --out-dir OUT_DIR
  --case-id CASE_ID
  --route-c-interval-match ROUTE_C_INTERVAL_MATCH
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench variant-evidence

```text
usage: candidate-qtl-bench variant-evidence [-h] --vcf VCF --fasta FASTA --out-dir OUT_DIR [--interval INTERVAL]
                                            [--snpeff-vcf SNPEFF_VCF] [--vep-vcf VEP_VCF] [--execute-norm] [--json]

options:
  -h, --help            show this help message and exit
  --vcf VCF
  --fasta FASTA
  --out-dir OUT_DIR
  --interval INTERVAL
  --snpeff-vcf SNPEFF_VCF
  --vep-vcf VEP_VCF
  --execute-norm
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench weight-policy-benchmark

```text
usage: candidate-qtl-bench weight-policy-benchmark [-h] --manifest MANIFEST --out-dir OUT_DIR
                                                   [--blocked-term BLOCKED_TERM]
                                                   [--bootstrap-iterations BOOTSTRAP_ITERATIONS] [--seed SEED]
                                                   [--crop CROP] [--json]

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST   benchmark_manifest.json from benchmark-reference.
  --out-dir OUT_DIR
  --blocked-term BLOCKED_TERM
                        Known answer/locus term excluded from upstream benchmark vocabulary; repeatable and recorded.
  --bootstrap-iterations BOOTSTRAP_ITERATIONS
  --seed SEED
  --crop CROP           Crop label when cases do not carry one.
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench wet-candidate-handoff

```text
usage: candidate-qtl-bench wet-candidate-handoff [-h] [--route-a-bridge ROUTE_A_BRIDGE]
                                                 [--route-b-claims ROUTE_B_CLAIMS]
                                                 [--route-b-candidates ROUTE_B_CANDIDATES]
                                                 [--route-c-members ROUTE_C_MEMBERS]
                                                 [--route-a-manifest ROUTE_A_MANIFEST]
                                                 [--route-b-claims-manifest ROUTE_B_CLAIMS_MANIFEST]
                                                 [--route-b-candidates-manifest ROUTE_B_CANDIDATES_MANIFEST]
                                                 [--route-c-manifest ROUTE_C_MANIFEST]
                                                 [--refined-candidates REFINED_CANDIDATES]
                                                 [--refinement-manifest REFINEMENT_MANIFEST]
                                                 [--refinement-original-candidates REFINEMENT_ORIGINAL_CANDIDATES]
                                                 [--refinement-original-union-manifest REFINEMENT_ORIGINAL_UNION_MANIFEST]
                                                 [--ranked RANKED] [--family-recovery FAMILY_RECOVERY]
                                                 [--convergence CONVERGENCE] [--family-copy FAMILY_COPY]
                                                 [--family-copy-shortlist FAMILY_COPY_SHORTLIST]
                                                 [--ambiguity-context AMBIGUITY_CONTEXT]
                                                 [--interval-atlas INTERVAL_ATLAS]
                                                 [--review-candidates REVIEW_CANDIDATES]
                                                 [--review-candidates-manifest REVIEW_CANDIDATES_MANIFEST]
                                                 [--identity-assessment IDENTITY_ASSESSMENT]
                                                 [--identity-target-proteins IDENTITY_TARGET_PROTEINS]
                                                 [--reference-function REFERENCE_FUNCTION]
                                                 [--wet-fill-mode {tiered_interval_fill,supported_only}]
                                                 [--deterministic-subfamily-refinement] --out-dir OUT_DIR
                                                 [--wet-budget WET_BUDGET] [--budget-policy {target,strict}]
                                                 [--adaptive-round ADAPTIVE_ROUND]
                                                 [--max-adaptive-rounds MAX_ADAPTIVE_ROUNDS]
                                                 [--min-route-candidates MIN_ROUTE_CANDIDATES] [--case-id CASE_ID]
                                                 [--json]

options:
  -h, --help            show this help message and exit
  --route-a-bridge ROUTE_A_BRIDGE
                        Frozen family_candidate_bridge.tsv.
  --route-b-claims ROUTE_B_CLAIMS
                        Final semantic claims TSV.
  --route-b-candidates ROUTE_B_CANDIDATES
                        Frozen route_b_verified_candidates.tsv.
  --route-c-members ROUTE_C_MEMBERS
                        Frozen Direct-C member ledger. Verified direct support enters the union; review-only rows may
                        fill tier2 Wet slots.
  --route-a-manifest ROUTE_A_MANIFEST
                        Producer manifest for Route A bridge.
  --route-b-claims-manifest ROUTE_B_CLAIMS_MANIFEST
                        Producer manifest for Route B claims.
  --route-b-candidates-manifest ROUTE_B_CANDIDATES_MANIFEST
                        Producer manifest for Route B frozen candidates.
  --route-c-manifest ROUTE_C_MANIFEST
                        route_c_direct_members_manifest.json.
  --refined-candidates REFINED_CANDIDATES
                        Machine-verified refined candidate set; re-entry requires adaptive-round >=1 and an original
                        canonical union above the Wet budget.
  --refinement-manifest REFINEMENT_MANIFEST
                        adaptive_refinement_manifest.json binding the verification and refined output.
  --refinement-original-candidates REFINEMENT_ORIGINAL_CANDIDATES
                        Original supported member union consumed by adaptive-refinement-apply.
  --refinement-original-union-manifest REFINEMENT_ORIGINAL_UNION_MANIFEST
                        Canonical supported_member_union_manifest.json binding the original candidate set.
  --ranked RANKED       Deterministic ranked_candidates.tsv. Read-only ordering inside the complete frozen Wet
                        membership.
  --family-recovery FAMILY_RECOVERY
                        Optional family recovery handoff for narrowing.
  --convergence CONVERGENCE
                        Optional dual/triple-route convergence table for narrowing.
  --family-copy FAMILY_COPY
                        Optional family-to-copy table for narrowing.
  --family-copy-shortlist FAMILY_COPY_SHORTLIST
                        Optional family-copy shortlist for narrowing.
  --ambiguity-context AMBIGUITY_CONTEXT
                        Optional per-candidate ambiguity context handed to the adaptive synthesis request.
  --interval-atlas INTERVAL_ATLAS
                        Frozen interval_gene_atlas.tsv. Enables tiered review-slot filling when the supported set is
                        under the Wet budget.
  --review-candidates REVIEW_CANDIDATES
                        Sealed A/B review candidates TSV; repeat for multiple routes. Retained ahead of generic fill
                        without verified support.
  --review-candidates-manifest REVIEW_CANDIDATES_MANIFEST
                        Matching review producer manifest; repeat in the same order as --review-candidates.
  --identity-assessment IDENTITY_ASSESSMENT
                        Explicit hash-bound post-freeze parent/independent identity review; derives member-bound Tier
                        evidence without changing support or frozen papers.
  --identity-target-proteins IDENTITY_TARGET_PROTEINS
                        Canonical interval candidate protein FASTA; required with --identity-assessment to verify the
                        exact target sequence.
  --reference-function REFERENCE_FUNCTION
                        Optional reference function/mutant-phenotype corpus; lets a filled candidate be tiered as
                        tier2_reference_phenotype.
  --wet-fill-mode {tiered_interval_fill,supported_only}
                        How unused Wet review slots are filled. Filled candidates carry no support and no score.
  --deterministic-subfamily-refinement
                        When the gate reports over-budget, apply one machine-verified curated cross-subfamily
                        refinement pass (requires --family-copy) and re-gate once instead of only emitting an adaptive
                        synthesis request.
  --out-dir OUT_DIR
  --wet-budget WET_BUDGET
                        Operational review target, not a scientific cutoff.
  --budget-policy {target,strict}
                        Use strict only for an explicitly fixed resource capacity.
  --adaptive-round ADAPTIVE_ROUND
  --max-adaptive-rounds MAX_ADAPTIVE_ROUNDS
  --min-route-candidates MIN_ROUTE_CANDIDATES
                        Operational A/C-first signal: a final Wet set below this asks the parent to rethink Route A in
                        a fresh generation. Never pads or trims membership.
  --case-id CASE_ID
  --json                Emit a stable JSON result to stdout.
```

## candidate-qtl-bench wet-state

```text
usage: candidate-qtl-bench wet-state [-h] --out-dir OUT_DIR [--assays ASSAYS] [--coverage COVERAGE]
                                     [--results RESULTS] [--previous-history PREVIOUS_HISTORY]
                                     [--primer-results PRIMER_RESULTS] [--exon-coverage EXON_COVERAGE]
                                     [--min-coverage MIN_COVERAGE]
                                     [--required-read-direction {bidirectional,forward,reverse,any}] [--json]

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR
  --assays ASSAYS       Canonical assay registry TSV.
  --coverage COVERAGE   Actual Wet target coverage TSV; in-silico design coverage does not pass the Wet gate.
  --results RESULTS     Parent/F1/individual-level Wet result TSV.
  --previous-history PREVIOUS_HISTORY
                        Previous wet_state_history.tsv to append without overwriting.
  --primer-results PRIMER_RESULTS
                        PrimerBlast import TSV; imported as design state only.
  --exon-coverage EXON_COVERAGE
                        exon coverage_windows.tsv; imported as design coverage only.
  --min-coverage MIN_COVERAGE
  --required-read-direction {bidirectional,forward,reverse,any}
  --json                Emit a stable JSON result to stdout.
```
