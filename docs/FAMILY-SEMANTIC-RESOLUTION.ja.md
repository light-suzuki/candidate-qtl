# Family / member / copy semantic resolution

`family-to-copy`のdeterministic evidenceでcopyを決められない場合だけ、LLM/agentを**semantic proposal layer**として使います。LLMはfinal judgeではありません。

```text
supported_member_union
    ↓
canonical family recovery
    ↓
deterministic family-to-copy
    ├─ preferred / co-preferred → STOP（LLM不要）
    │
    └─ unresolved family/copy
          ↓
      bounded semantic prompt
          ↓
      LLM / agent proposal
          ↓
      proposed pair verification
      (interval member + supplied reference + strong orthology)
          ↓
      machine competitor scan
      (same reference × all eligible copies)
          ↓
      preferred / co-preferred ordering
          ↓
      Wet membershipは不変
```

## Canonical provenance gate

semantic fallbackは任意の`family_copy`を受けません。production wrapperは先に#514のcanonical handoffを実行し、次を確認します。

- `family_copy_producer_class = canonical_supported_family_copy`
- `family_recovery_source = supported_member_union`
- `legacy_triple_route_family_context_consumed = false`

この条件を満たさないlegacy/diagnostic family-copy artifactしか無い場合、semantic layerは`not_evaluated_noncanonical_family_copy_provenance`として停止します。Wet membership自体はcanonical supported unionからそのまま維持されます。

## Cheap-first rule

semantic fallbackへ送る前にdeterministic resolutionを優先します。

1. exact stable ID
2. deterministic family / orthology
3. expanded symbol / synonym match
4. deterministic `family-to-copy`
5. **それでも未解決familyだけ** semantic resolution
6. proposed pairのmachine verification
7. 同じreferenceを持つ全copyのmachine competitor scan
8. 未解決なら複数copyをそのまま保持

`preferred_copy` / `co_preferred_tie`が既に得られたfamilyはsemantic promptを作りません。

## Prompt対象

現在のfallback対象は:

- `family_supported_copy_unresolved`
- `paralog_ambiguous`
- `tandem_ambiguous`
- `insufficient_discriminating_evidence`
- `family_resolution_required`

`not_in_validated_copy_set`は「別copy set / family identity問題」の可能性があるため、このcopy promotion fallbackでは自動昇格対象にしません。

1 familyあたり最大20 memberに制限し、それを超えるfamilyはreviewへ残します。

## LLM / agentに渡すもの

frozen artifactから次をbundleします。

- family ID
- interval member gene IDs
- current deterministic copy state
- orthology class / relationship
- supplied reference gene IDs
- subfamily / target-relative subfamily state
- domain architecture / target-relative domain state
- orthogroup
- annotation / description
- canonical family bridgeのroute / mapping tier / ambiguity/provider state

候補集合はLLMが追加できません。reference geneもpromptに供給されたもの以外はimport時にunsupportedになります。

## LLMがしてよいこと

semantic relationを**提案**できます。

- `same_family`
- `same_subfamily`
- `likely_reference_relative`
- `distinct_subfamily`
- `unresolved`

また、どのinterval copyがどの supplied reference gene に対応すると解釈するかを提案できます。

ここで重要なのは、LLMが「A1がAtXに対応しそう」と提案しても、**A1が唯一のcopyであるとはまだ決まらない**ことです。LLM proposalは「このreference関係をmachineで検査する」という仮説です。

## LLMがしてはいけないこと

- numeric gene score
- causal probability / confidence score
- candidate追加・削除
- supplied evidenceにないaccession/reference/sourceの発明
- proposalだけで`preferred_copy`を成立させること
- 同等のparalogを黙って省略することでsingle-copy扱いにすること

`score` / `causal_score` / `probability` / `causal_probability` / `confidence_score`をmodel JSONに含めるとfail-closedで拒否します。

## Machine verification: 2段階

### 1. proposed pair verification

proposalのcandidate-reference pairについて:

1. candidate geneが元のunresolved family memberである
2. proposed reference geneがそのcandidateのfrozen atlasに既に存在する
3. machine artifactのorthology relationが1:1 / one-to-one / RBH相当である
4. semantic relationが`same_subfamily`または`likely_reference_relative`

を確認します。

ここで通ったものは**preliminary `verified_copy_mapping`**です。まだsingle-copy確定ではありません。

### 2. machine competitor scan

LLMが1 copyだけを提案しても、そのreferenceに対して別copyが同等のmachine relationを持つ可能性があります。そこでfinal gateは、提案されたreferenceごとに**family内の全eligible unresolved copyを再走査**します。

- supplied referenceが同じ
- frozen atlas上でそのreferenceにlink済み
- orthologyが1:1 / one-to-one / RBH相当

を満たすcopyをmachine側で全て拾います。

結果:

- machine-equivalent copyが1つだけ → `preferred_copy`
- 2つ以上 → **LLMが省略していても全て `co_preferred_tie`**
- competitor scanを通らないもの → 元のdeterministic stateを保持

つまり、LLM proposalはreference hypothesisを提示できますが、**uniquenessを決めるのはmachine competitor scan**です。

## Wet membershipとの境界

このfallbackは**membershipを変更しません**。

`wet_candidate_handoff`へ渡るsupported A/B member unionとcanonical family recoveryはsemantic copy fallbackより先に確定します。semantic verificationは、すでにsupportされたWet candidate set内のorderingに使う`family_copy` / shortlistを補助するだけです。

したがって:

- semantic responseが無い → original deterministic copy artifactsと#514のWet resultをそのまま使う
- provider failure → negativeにしない
- verification失敗 → original unresolved stateを維持
- A1がverified preferredになってもA2をWet membershipから削除しない
- A1だけ提案されてもA2がmachine-equivalentならA1/A2をco-tieに戻す
- semantic ordering後にWet membership setまたはsupported-union SHAが変わればfail-closed

## Artifact flow

canonical profileでは:

```text
<profile>/b04_canonical_semantic/
  07_supported_union_for_family/
  08_canonical_family_recovery/
  09_canonical_family_to_copy/
  10_family_semantic/
    01_request/
      family_semantic_request.json
      family_semantic_prompts.jsonl
    02_llm_run/
      llm_candidate_responses.json
      llm_run_manifest.json        # optional
    03_import/
      family_semantic_proposals.tsv
    04_preliminary_verify/
      family_semantic_copy_verification.tsv
      family_semantic_copy_shortlist.tsv
    05_competitor_gate/
      family_semantic_copy_final.tsv
      family_semantic_copy_shortlist_final.tsv
    canonical_family_semantic_wet_manifest.json
```

verified semantic orderingが成立した場合だけ、final ordering handoffを:

```text
<profile>/z05_canonical_family_semantic_wet_handoff/
```

へ出します。semantic responseが無い、またはsemantic promotionが無い場合は#514の既存Wet handoffを再生成しません。

production orderingへ渡すのは`05_competitor_gate`のfinal artifactsです。`04_preliminary_verify`はLLM-proposed pairの検証中間物であり、single-copy uniquenessの正本ではありません。

## Design principle

> LLMは人間研究者なら名称・文脈・reference biologyを読んで考える部分を担当する。candidate membershipやcopy promotion、特に「他に同等copyが無い」という主張は、機械検証できる事実でゲートする。

この分離により、semantic flexibilityを使いながら`not_evaluated != negative`、paralog保持、no-score、source/identity verificationの不変条件を維持します。
