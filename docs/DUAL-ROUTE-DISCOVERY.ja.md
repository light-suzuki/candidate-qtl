# Route A / Route B discovery contract

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

この文書は、phenotype-first Route Aとinterval-first Route Bの役割分担、freeze independence、production合流点を説明します。

Direct Route CはA/B discoveryそのものではありません。ただし、事前freeze済みDirect-C recordをpost-B0 projection/member gateへ通した場合、eligible memberは同じdownstream `supported_member_union`へoptionalに合流できます。

historical `dual-route-compare` / `triple-route-compare`の旧family handoffはproduction正本ではありません。

文献・claim確認が未完のA/B提案も、機械で区間memberへの対応を確認できれば、
verified supportとは別のreview候補としてWetへ保持します。汎用のinterval fillより先に扱い、
20件という運用目標による早期cutoffを避けます。`supported_member_union`の意味は変えません。

---

## 1. A/Bの対称構造

```text
Route A（phenotype-first）             Route B（interval-first）
   phenotype                            QTL区間・全gene inventory
       ↓                                          ↓
 LLM：機構・参照geneを仮説化          B0：family/orthogroup圧縮
       ↓                                          ↓
 機械：参照ID/配列・区間投影           LLM：機能から表現型との橋渡し
       ↓                                          ↓
 書誌・member照合                        書誌・member照合
       └──────────────────┬───────────────────────┘
                          ↓
                   A/Bの独立結果を合流
                 ┌────────┴─────────┐
                 ↓                  ↓
          verified support     未支持でも適格なreview
                 ↓                  ↓
        supported_member_union   review候補（別入力）
                 └────────┬─────────┘
                          ↓
                family/copyとWet handoff
                          ↓
                   wet_candidates.tsv
```

**A-only/B-onlyも適格ならWet候補に残せます。** 新規`initial_bibliography_v1`では、初期candidate admissionを独立全文レビューの完了に依存させません。review候補が`supported_member_union`に昇格したことにはならず、verified supportには既存の厳格な確認が必要です。旧bundle/凍結世代は旧policyで再生します。

- Route A = `LLM/agent -> machine`
- Route B = `machine -> LLM/agent -> machine verification`
- familyはfinal answerではなく中間表現
- raw interval geneをfamily/orthogroup/domain単位へ圧縮する
- A/Bはfreeze前に互いのresult artifactを見ない
- agreementはprovenanceでありscore bonusではない

---

## 2. Route A — phenotype-first

```text
phenotype
  -> mechanism / pathway / reference gene / family hypothesis
  -> source / accession / species / sequence verification
  -> verified baitまたはidentity確認済みのreview hypothesis
  -> deterministic interval projection
  -> verified support ＋ 適格なA-review interval gene（別属性）
```

LLM/agentがしてよいこと:

- biological mechanism / pathway hypothesis
- gene family / reference gene proposal
- source-backed search terms
- competing explanation

machineが確定するもの:

- stable accession / species
- source relation
- sequence
- interval projection
- 区間内gene/copyと凍結参照のidentity・投影対応（review候補を含む。verified supportとは別）

LLM proseだけでinterval geneをcause扱いしません。

---

## 3. Route B — interval-first

```text
all interval genes
  -> B0 phenotype-blind inventory freeze
  -> family / orthogroup / subfamily / domain compression
  -> LLM/agentによる既知機能調査・表現型との照合
  -> 新規通常：paper ID/title・reference/memberの機械照合
       -> 初期Tier（T1〜T6/U）＋適格なB-review候補
  -> 厳格support・旧policyのみ：独立原文確認＋transfer verification
       -> route_b_verified_candidates.tsv
  -> eligible route_b_review_candidates.tsv（supported unionとは別にWetへ）
```

区間の数百geneをrawでLLMへ投げません。

historical phenotype-aware B1 reverse lookupはdiagnostic/compatibility用途であり、production Route-B membershipには入りません。

互換性用旧bundleのB-only reviewは`source_function_v1`に従い、出典付き参照機能と区間対応・具体的な条件付き推論を確認します。対象種の因果論文・同一表現型は不要ですが、出典IDや汎用機能だけでは足りません。正本は [canonical Route B](CANONICAL-ROUTE-B.ja.md) と [Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md) です。

---

## 4. Freeze independence

A/Bが同じbiologyや同じpublic sourceへ到達すること自体は禁止しません。ただしfreeze前に互いのresult artifactを見せません。

目的はfailure modeを分離することです。

- Route A miss: phenotype -> mechanism/reference hypothesisを外した
- Route B miss: annotation/family/reference coverageまたはsemantic transferを外した

freeze後にagreement/disagreementを比較できますが、agreementでnumeric scoreを加算しません。

---

## 5. Optional Direct Cとの境界

Direct Cは、既報のdirect trait–gene/locus/variant factを回収するknown-answer rescueです。

```text
pre-interval Direct-C search/freeze
  -> post-B0 projection
  -> route-c-direct-members
  -> verified_direct_supportだけoptional unionへ
```

Direct CはA/B correctnessの必須dependencyではありません。

current mainではA/B/C provenance mergeまで配線済みです。#661のbinding gateは実装済み・close済みですが、#663/#662の残課題があるためDirect-C automatic support全体はprovisionalです。

- #661: source gene/species/trait/experimentの封印済みbindingを要求する実装済みgate。既存 [producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md) を利用
- #663: projection order invariance / full provenance
- #662: 検索局所失敗の隔離・取得済み結果の再利用。strict lineageは任意

S0/S1 D2、cross-species、orthology/BLAST-onlyはreview-only、outside/no-hit/unresolvedはledger-onlyです。

legacy lexical/LSA/neural Route C retrievalはdiagnostic/historicalであり、Direct Cとは別です。

---

## 6. Production member merge

production合流はhistorical comparison tableではなく:

```text
Route A supported members
+
Route B verified members
+
optional Direct-C verified_direct_support members
        ↓
supported_member_union.tsv
```

union rule:

- A-only / B-only / C-only -> retain
- A+B / A+C / B+C / A+B+C -> retain + provenance
- ambiguous supported paralogs -> retain multiple
- provider/not-evaluated/unverified -> negativeへ変換しない
- Direct-C review-only / ledger-only -> support unionへ入れない

`supported_member_union.tsv`はrankではありません。

---

## 7. Canonical family recovery

member unionの後にfamily-level production handoffを作ります。

```text
supported_member_union
  -> canonical family recovery
```

cheap-first identity hierarchy:

1. exact raw family identity
2. deterministic canonical family key
3. mutually unique expanded symbol/synonym match
4. ambiguityが残れば別family hypothesisを保持
5. genuine semantic ambiguityだけ後段agentへ

expanded matchはstrict convergenceでもscore bonusでもありません。

---

## 8. Family -> copy

canonical family memberだけをdeterministic `family_to_copy()`へ渡します。

使えるevidence:

- orthogroup / orthology class
- verified reference relation
- subfamily
- domain architecture
- synteny
- expression
- positional evidence

family-level supportだけでsingle copyを強制しません。

代表state:

```text
preferred_copy
co_preferred_tie
family_supported_copy_unresolved
paralog_ambiguous
tandem_ambiguous
insufficient_discriminating_evidence
family_resolution_required
not_in_validated_copy_set
```

family/copyのsemantic resolutionは、canonical `family_semantic_resolution.py` と
`family_semantic_competitor_gate.py` を経由します。

---

## 9. Optional semantic copy fallback

決定論的evidenceでcopyを決めきれないfamilyだけagent/LLMへproposalを依頼できます。

LLMはcandidate-reference mappingを提案します。machineが:

- original family memberか
- frozen referenceに存在するか
- 1:1/RBH相当か
- equivalent competitor copyが無いか

を再確認します。

1 copyなら`preferred_copy`、複数なら`co_preferred_tie`として全copyを戻します。

---

## 10. Historical compare commands

`dual-route-compare` / `triple-route-compare`は次に使えます。

- A/B/C agreement/disagreement report
- historical benchmark reproduction
- source-dependence audit
- diagnostic convergence view

ただし旧:

```text
comparison table
  -> family_to_copy_input.tsv
  -> production family-to-copy
```

はnon-canonicalです。新しいproduction consumerを追加しません。

---

## 11. Final Wet boundary

production membership正本は`wet_candidates.tsv`です。

通常targetは#817の`upper_evidence_scoped_v3`です。全gene台帳と元の適格性を保存し、累積25を初めて超えるTierだけを比較します。上位は保持、後続は台帳に残し、元の調査担当者の理由付き判断を`original_assessment_v1`として機械照合します。20未満だけで追加探索せず、25で強制切断しません。旧`evidence_scoped_v2`・版なしscopeは旧policyで再生します。詳細は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)に従います。以下はlegacy strict/`supported_only`/明示refinementの互換規則です。

```text
supported <= wet_budget（通常の目標は20件）
  -> 全supportedを tier1_verified_support として保持
  -> 空きreview枠を tier2_reference_phenotype で埋める
  -> 次に tier3_family_or_domain で埋める
  -> tier4_dark はfillしない

supported > wet_budget
  -> 全supportedを保持して通常のWet handoff
  -> 必要な場合だけverified refinement
```

fill rowは`membership_basis=interval_review_coverage`であり、support、score、route evidenceを持ちません。

Direct-C `reference_phenotype_review`はreview coverageとして検討できますが、verified supportへ自動転移しません。

---

## 12. Summary

```text
A: phenotype -> semantic hypothesis -> machine-verified member
B: interval -> machine compression -> semantic interpretation -> machine-verified member
C: optional pre-frozen direct fact -> post-B0 projection -> narrow member gate
                         ↓
                supported member union
                         ↓
                canonical family recovery
                         ↓
                deterministic family-to-copy
                         ↓
                optional semantic fallback
                         ↓
                   Wet membership/order
```

> LLMへ任せるのは人間なら意味を読む部分です。member identity、source、orthology、copy uniqueness、membershipはmachine-verifiable contractで固定します。
