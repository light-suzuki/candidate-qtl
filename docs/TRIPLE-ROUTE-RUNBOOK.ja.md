# triple-route 実行runbook — A/B + optional Direct C -> family/copy -> Wet

このrunbookは、current mainで実案件を同じ順序で実行するための手順です。

Wetはverified supportに加え、区間memberへの対応を確認したA/B提案を明示的なreview候補として保持します。
文献・claim確認が未完の提案を汎用fillより先に残し、20件は早期cutoffではなく後段精査の運用目標とします。
review候補は`supported_member_union`へ昇格しません。

`triple-route`という名前には歴史的なlexical retrieval laneも残っていますが、production membershipは比較表やlegacy `family_to_copy_input.tsv`では決まりません。

```text
               phenotype + QTL interval
                          │
            ┌─────────────┴─────────────┐
            ↓                           ↓
 Route A：phenotype-first     Route B：interval-first
 LLM→機械で参照/区間照合       B0→LLM→機械でmember照合
            └─────────────┬─────────────┘
                          ↓
                   A/B候補を統合
             ┌────────────┴───────────┐
             ↓                        ↓
       verified support          適格review候補
             ↓                        ↓
    supported_member_union     review入力（別経路）
             └───────────┬────────────┘
                         ↓
              family/copy・Wet handoff
                         ↑
           optional Route C（最後のsafety net）
           QTL-blind search/freeze
             → post-B0 projection / member gate
             → verified support / review / ledger
                         ↓
             wet_candidates.tsv＋trace/report
```

A/B候補の支持確認とWet採用は別です。Cは任意の最終セーフティネットとして従来どおり扱い、Cの起動条件・promotion gate・実装をこの図の修正で変えません。Cのsupportだけが`supported_member_union`へ入り、reviewは別入力、no-hit等は台帳に残ります。

family/orthogroupは最終回答ではなく中間表現です。

---

## 1. Core rules

1. Route A、Route B、Direct Cはfreeze前に互いのcandidate resultを見ない。
2. Route Aはphenotype-first: LLM/agentがmechanism/reference/familyを提案し、machineがsource/accession/sequence/interval memberを検証する。
3. Route Bはinterval-first: machineがB0 inventoryをfreezeし、family/orthogroup/domain単位へ圧縮してからsemantic interpretationする。
4. Direct Cはoptional: 別contextで表現型・対象種から自由に既報gene/locus/variantを検索し、A/Bと区間を見ずにfreezeする。B0 atlas後にだけprojectする。
5. LLM outputだけでcandidate membershipを変更しない。
6. provider failure / not measured / no responseは`not_evaluated`でありnegativeではない。
7. A/B/C agreementはprovenanceでありnumeric score bonusではない。
8. final membershipの正本は`wet_candidates.tsv`。
9. 通常targetの採用・補充・境界精査は[canonical membership](CANONICAL-WET-MEMBERSHIP.ja.md)に従う。旧support tierとdark除外を通常targetの選定規則へ流用しない。
10. fill rowは`membership_basis=interval_review_coverage`でありsupportではない。
11. supported > `wet_budget`を任意件数へ切らない。件数だけで生物学的失敗とは扱わない。
12. Direct CはA/B correctnessの必須dependencyではない。
13. 受理済みDirect-C救出候補はWetへ直接保持し、合計25件超・低いTierを理由に外さない。境界精査とfinishでも保護し、support/review属性は変えない。

### Direct Cのcurrent safety boundary

Direct Cのcanonical handoffと#661のbinding gateは実装済みで、#661はclose済みです。#663/#662の残課題があるためautomatic support全体はprovisionalです。

- #661: 封印済みsource gene/species/trait/experiment bindingを要求する実装済みgate。既存 [producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md) を使う
- #663: projection order invariance / provenance completeness
- #662: 検索局所失敗の隔離・取得済み結果の再利用。strict search lineageは任意

current runでDirect-C tier1 rowが出た場合はsource bindingとprojection basisを人間が確認します。安全側ではreview-onlyとして扱います。

---

## 2. 必要な入力

### Route A/B共通

- phenotype text
- focal species
- QTL interval
- reference FASTA / GFF3
- interval proteins
- 既存の参照対応・family/domain・機能注釈があれば再利用。不足時は設定済みprotein DBで補う（[参照fallback](INVESTIGATION-MATERIAL.ja.md)）
- optional domains / orthology / reference function cache

### Direct Cを使う場合

QTLを開く前に用意する:

- `route_c_direct_rescue_frozen.tsv`
- `route_c_direct_rescue_freeze_manifest.json`
- QTL ID（projection時に指定）
- optional `eligible_interval_hits.tsv`

`triple-route` profileはDirect-C proposal extractionもBLAST executionも行いません。

---

## 3. Route A — interval-blind hypothesis generation

現行Route Aの概念:

```text
phenotype-investigate
  -> phenotype-localization-packet
  -> external fresh localization generation
  -> phenotype-localization-freeze
  -> external fresh 12-hypothesis portfolio
  -> mechanism-import
  -> mechanism-coverage-audit
       pass -> reference resolution / source verification
       gap  -> at most 3 fresh complementary hypotheses
  -> verified bait freeze
  -> interval projection
```

Route A initial generationへ渡してはいけないもの:

- interval gene IDs
- Route B output
- Direct-C result
- truth gene/family/answer paper
- prior benchmark score

`guided-profile`自体はLLMを内部実行しません。必要なpayloadが無ければreview/request artifactを出して停止し、外部agent/LLMのresponseを同じcontractで戻します。

---

## 4. Route B — B0を先にfreezeする

production Route B:

```text
interval_gene_atlas.tsv
  -> family / orthogroup / domain / reference compression
  -> semantic bundle + prompts
  -> external LLM/agent response（新資料は再freeze）
  -> 新規通常：書誌・参照/member照合、初期Tier計算
       -> eligible route_b_review_candidates.tsv → Wet候補へ
  -> 厳格support・旧source_function_v1：原文/機能review、transfer verification
       -> route_b_verified_candidates.tsv → supported_member_union
```

数百geneをrawでLLMへ投げません。B0 atlasはphenotype-freeで、manifest/hashにより固定します。

semantic responseが無い場合はbiology missにせず、例えば:

```text
state = semantic_review_required
provider_state = not_evaluated_no_semantic_response
```

として止めます。

標準`llm-run`を使う例:

```bash
candidate-qtl llm-run \
  --prompts <OUT>/b04_canonical_semantic/01_bundle/route_b_semantic_prompts.jsonl \
  --out-dir <OUT>/b04_canonical_semantic/02_llm_run \
  --provider <PROVIDER> \
  --model <MODEL>
```

responseを置いた後、同じ`guided-profile --profile triple-route --execute`を再実行します。

legacy phenotype-aware reverse lookupはdiagnostic/compatibility用途であり、canonical Route-B membership入力ではありません。

---

## 5. Direct C — searchとfreezeをQTL前に完了する

### 5.1 Search — 通常はハーネス標準のWeb検索

C担当は、A/Bとは別の新しいcontextに**対象種と表現型だけ**を渡して起動します。問いは「この植物でこの形質に関わると報告された遺伝子、古典座位、変異は何か？」です。A/B回答、QTL/区間遺伝子、事前gene/family候補、truth、親の調査履歴は渡しません。複数の報告が同じ遺伝子でも、QTLに入らない報告でも、pre-freezeでは落としません。CはAの結果を参照して差分探索する経路ではありません。

Codex/OpenCode等の担当ハーネスで、実際に利用可能な標準Web検索・取得ツールを使います。生成済みベンチ入力を変更せず通常用の短い指示を得る場合:

```bash
python -c "from candidate_qtl.route_c_direct_pipeline import render_harness_native_direct_prompt; print(render_harness_native_direct_prompt())"
```

必要な場合だけ既存の専用検索も使えます（必須ではない）:

```bash
candidate-qtl route-c-direct-search \
  --phenotype "<PHENOTYPE>" \
  --focal-species "<SPECIES>" \
  --species-scope templates/route_c_species_scope.tsv \
  --out-dir <C_SEARCH_OUT> \
  --case-id <CASE> \
  --allow-external-api
```

実際の検索語、読んだ出典URL・本文範囲、PMID/DOI、gene/locus/variant名、実験種、状態を保存し、**元の17列`route_c_direct_proposals.tsv`契約**へ出力します。自由探索でも最終的な提案の書式とID/種/形質照合は維持します。検索hitはcandidate linkであり、evidenceではありません。検索ツールが使えなかった場合は未評価として扱い、検索完了していないのに`completed_empty`を報告しません。

### 5.2 Freeze

```bash
candidate-qtl route-c-direct-freeze \
  --proposals <route_c_direct_proposals.tsv> \
  --out-dir <C_FREEZE_OUT> \
  ...source cache/provider options...
```

freeze artifactにQTL/candidate/truth fieldを入れません。

generic verificationはcitation/cue triageであり、それだけでauto-supportにはしません。#661のRoute-C固有binding gateは実装済みです。[binding producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md) で原典に基づく人間の確認を封印し、bound manifestをprojection/memberへ渡します。gene、source species/accession、experiment、phenotypeが同じclaimへ結び付いている必要があります。

---

## 6. `guided-profile triple-route` 入口

通常入口:

```bash
candidate-qtl guided-profile \
  --profile triple-route \
  --execute \
  --out-dir <OUT> \
  --case-id <CASE> \
  --focal-species "<SPECIES>" \
  --phenotype "<PHENOTYPE>" \
  ...Route A/B inputs...
```

Direct Cを接続する場合:

```bash
candidate-qtl guided-profile \
  --profile triple-route \
  --execute \
  --out-dir <OUT> \
  --case-id <CASE> \
  --focal-species "<SPECIES>" \
  --phenotype "<PHENOTYPE>" \
  --route-c-direct-frozen <C_FREEZE_OUT>/route_c_direct_rescue_frozen.tsv \
  --route-c-direct-freeze-manifest <C_FREEZE_OUT>/route_c_direct_rescue_freeze_manifest.json \
  --route-c-direct-qtl-id <QTL_ID> \
  --route-c-direct-eligible-hits <eligible_interval_hits.tsv> \
  ...Route A/B inputs...
```

3点セット:

- `--route-c-direct-frozen`
- `--route-c-direct-freeze-manifest`
- `--route-c-direct-qtl-id`

のどれかだけを渡すとfail closedします。`--route-c-direct-eligible-hits`はoptionalです。

profile内部では:

```text
B0 interval atlas generation
  -> verify Direct-C freeze
  -> route-c-direct-project
  -> route-c-direct-members
  -> canonical post-B0 family/semantic/Wet extension
```

の順です。

---

## 7. Direct-C projectionの読み方

projection states:

```text
inside_current_qtl
outside_current_qtl
no_interval_match
mapping_unresolved
not_projected
```

重要:

- `no_interval_match`はcurrent intervalでeligible matchが無いこと
- `outside_current_qtl`はfocal genome上の位置を実際に解決しQTL外と確認したこと

interval BLAST no-hitをoutside-QTL presenceへ変換しません。現行standard projectorはwhole-genome outside resolverを持たないため、`outside_current_qtl`自動確定は未実装です。

#663完了前は、exact/orthology/BLAST evidenceが入力順で上書きされる可能性があるため、selected `projection_basis`を人間が確認します。

---

## 8. Direct-C member state

### `verified_direct_support`

current codeが要求するもの:

```text
S0/S1
+ D0/D1
+ source verified
+ eligible direct record
+ verified gene/locus state
+ direct target/related trait state
+ inside_current_qtl
+ exact focal ID / assembly alias mapping
```

このrowだけがC supportとしてunionへ入ります。#661の封印済みbindingも必要です。#663/#662の残課題がある間はprovisionalとして扱います。

### `reference_phenotype_review`

- S0/S1 D2
- cross-species D0-D2
- orthology/BLAST-only projection
- direct factはあるがfocal copy supportへ安全に転移できない

### `ledger_only`

- outside
- no interval match
- unresolved
- not projected
- unverified / context-only

review-onlyとledger-onlyをverified supportへ混ぜません。

---

## 9. Production merge

production merge input:

```text
Route A supported members
+
Route B source/transfer verified members
+
optional Direct-C verified_direct_support rows
        ↓
supported_member_union.tsv
```

union rule:

- A-only / B-only / C-onlyを保持
- A+B / A+C / B+C / A+B+Cをprovenanceとして保持
- agreementでscoreを増やさない
- ambiguous copiesをcollapseしない
- not-evaluated/unverified/review-onlyをsupportへ変換しない

historical `dual-route-compare` / `triple-route-compare`の旧family handoffはproduction正本ではありません。

---

## 10. Canonical family recovery

`supported_member_union.tsv`からcanonical family recoveryを作ります。

cheap-first identity hierarchy:

1. exact raw family identity
2. deterministic canonical family key
3. mutually unique expanded symbol/synonym match
4. ambiguityが残れば別family hypothesisとして保持
5. genuinely semantic ambiguityだけ後段agentへ

family labelがなくorthogroupだけなら`orthogroup:<ID>`として保持し、family名を捏造しません。

---

## 11. Deterministic family-to-copy

canonical family memberだけを`family_to_copy()`へ渡します。

代表state:

```text
preferred_copy
co_preferred_tie
family_supported_copy_unresolved
paralog_ambiguous
tandem_ambiguous
copy_weakened
insufficient_discriminating_evidence
family_resolution_required
not_in_validated_copy_set
```

`preferred_copy`は最初にWetするcopyでありcause断定ではありません。

family supportだけで1 copyを強制しません。

---

## 12. Optional semantic copy fallback

決定論的に解けないfamilyだけagent/LLMへrequestを出せます。

対象例:

- `family_supported_copy_unresolved`
- `paralog_ambiguous`
- `tandem_ambiguous`
- `insufficient_discriminating_evidence`
- `family_resolution_required`

agentはcandidate-reference mappingを提案するだけです。machineがmember identity、frozen reference、orthology、competitor copyを再確認します。

responseが無くてもWet membershipは止めず、deterministic unresolved stateを保持します。

family/copyのsemantic resolutionはcanonical semantic pathで実行します。

---

## 13. Wet budget gate

通常targetは#817の`upper_evidence_scoped_v3`です。全gene台帳と元の適格性を保ち、累積25を初めて超えるTierだけを比較します。上位は保持、後続は台帳に残し、元の調査担当者の理由付き判断を`original_assessment_v1`として機械照合します。採用/比較待ち/理由付き保留を区別し、20との差分だけで追加調査せず、25で強制切断しません。旧`evidence_scoped_v2`・版なしscopeは旧policyで再生します。以下はlegacy strict/`supported_only`/明示refinementの互換規則です。詳細は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)に従います。

### supported <= `wet_budget`（通常の目標は20件）

`wet_candidates.tsv`へ全supportedを`tier1_verified_support`としてfreezeします。

空きreview枠:

1. `tier2_reference_phenotype`
2. `tier3_family_or_domain`
3. `tier4_dark`は使わない

fill rowは`membership_basis=interval_review_coverage`です。support、score、route flagを持ちません。

Direct-C cross-species/review-only rowは「見る理由」であり、transferred supportではありません。

### supported > `wet_budget`

全supportedを保持して通常のWet handoffへ進みます。20件は目安で、15〜30件程度も普通の結果です。
件数だけで停止・警告・理由説明を要求せず、この幅もhard limitにはしません。

```text
全supportedを保持 -> wet_candidates.tsv
必要時だけverified refinement
```

proposal単独でgeneを削除しません。固定容量を明示した`--budget-policy strict`の場合だけ、超過時のbounded refinement / reviewを要求します。

---

## 14. Ordering boundary

membership freeze後だけorderingを行います。

```text
wet_candidates.tsv
  -> ordering-only narrowing
  -> same gene set, new review order
```

ordering前後でgene setが変わればfail closedします。

---

## 15. Provenance checks

production artifactでは少なくとも次を確認します。

```text
family_recovery_source = supported_member_union
family_copy_producer_class = canonical_supported_family_copy
supported_union_identity_verified = true
rank_contribution = none
```

Direct Cを使った場合は追加で:

- Direct-C freeze SHA-256
- projection manifestとatlas fingerprint
- `route_c_direct_members` policy version
- C record/source IDs
- support_status
- species distance / directness
- selected projection basis

を確認します。

#662完了後はsearch cassette / proposal / freeze lineageも必須にします。

---

## 16. Legacy retrievalとの区別

`guided-profile triple-route`にはlegacy lexical/LSA/neural retrieval stageも残っています。これはdiagnostic/historicalです。

```text
legacy retrieval
  -> similarity-based phenotype retrieval
  -> historical compare/provenance
  -> no Direct-C support semantics
```

```text
Direct C
  -> source-backed direct report
  -> strict freeze/projection/member gate
  -> optional supported union
```

同じ「Route C」という古い名前が残っていても、この2つを混同しません。

---

## 17. Current blockers / stop conditions

次の場合はDirect-C support promotionを信頼せずreviewへ戻します。

- sourceにstable ID/species/trait bindingが無い
- symbolだけ一致し別stable IDがproposalに supplied
- real quoteはあるがexperimented geneとphenotypeの関係が未確認
- variant-only proposalでverification pathが不完全
- exactとBLASTが異なるgeneを指す
- exact basisがBLAST row orderで上書きされた疑い
- method not evaluatedとno hitを区別できない
- search/proposal/freeze lineageがstrict benchmarkでsealされていない
- whole-genome evidenceなしでoutside-QTLを主張している

---

## 18. Validation

文書・contract変更:

```bash
python -m unittest tests.test_agents_contract_consistency -v
```

Direct-C behavior変更ではさらに:

- `tests.test_route_c_direct_rescue`
- `tests.test_route_c_direct_pipeline`
- `tests.test_route_c_direct_members`
- `tests.test_supported_member_union_route_c`
- relevant reachability tests

を実行します。

現行配線と#661のbinding gateが実装済みであることと、#663/#662を含むscientific/reproducibility guaranteeが完成したことを分けて報告します。
