# Direct Route C — trait–gene/locus/variant rescue

> Parent Issue: [#640]()
> Phase 1 schema/validator: PR #641 merged
> Current main: search -> freeze -> projection -> `route-c-direct-members` -> optional supported-union/family-copy/Wet handoff implemented
> Change class: `EXPERIMENTAL` / safety fixes tracked as `CORE_FIX`

Direct Route Cは、A/Bと競争する「第3の予測器」ではありません。

**通常のC調査は、独立したLLMに「この植物のこの表現型に関わると報告された遺伝子・古典座位・変異には何があるか？」と普通に尋ねることから始めます。** 入力は対象種と表現型だけ。Aの回答、Bの台帳、QTL位置、候補gene/family一覧、truthを与えず、予め遺伝子名・座位・区間内hitを指定しません。Aと同じ報告を拾ってもよく、将来のQTL区間と一致しない報告も検索段階で捨てません。検索結果は実出典に対応させ、17列の既存proposalへ整理してからfreezeします。

ハーネスが提供するWeb検索・取得ツールを研究担当が直接使う方法を通常入口とします。**`route-c-direct-search`（Europe PMC + species-scope表）は使用可能な追加手段であり、先に実行する義務はありません。** この変更は調査入口の軽量化で、後段のsource/claim/identity・freeze・projection・Wet保護条件の緩和ではありません。

> **対象形質について、同種・近縁種・モデル植物で直接報告されたtrait–gene / trait–locus / trait–variant関係を回収し、QTLを見る前後を分離した証拠鎖として残すoptional known-answer rescue lane**

です。

---

## 1. 現在の実装状態

current mainでは次が実装されています。

```text
独立C担当のハーネス標準検索・取得（通常）
  -> 直接報告を既存17列proposalへ記録
  | optional: route-c-direct-search + species-scopeの文献検索

route-c-direct-freeze
  -> generic source_verifyを再利用
  -> direct record schemaへnormalize
  -> QTL/intervalを見ないままSHA-256 freeze

既存claim-binding producer prepare/seal
  -> 人または独立agentが原典を確認したbindingとbound freeze manifest

route-c-provider-bait（任意、記号だけのrecord向け）
  -> recordの報告種のUniProtKB entryのうち、recordの引用PMIDを参考文献に持ち、
     構造化gene名がrecordの報告記号と一致するものを1件だけbaitに選ぶ
  -> 0件・別配列の複数件はfail closed。同名でほぼ同一配列の複数件（変異アレル等）は
     reviewed > mRNA > 長さ > accession の順で代表を選ぶ。区間は見ない
  -> panel.tsv（blast-panel-freeze入力）とprovider_identity_edges.tsv

blast-panel-freeze -> blast-panel-search
  -> panel_id = record_id のeligible_interval_hits.tsv

route-c-direct-project
  -> frozen bytesを再検証
  -> interval atlas / orthology reference / precomputed eligible BLAST hitsへ投影
  -> --identity-mapping provider_identity_edges.tsv を渡すと、BLAST hitのsource
     identityを報告記号として記録する（state=verified_provider_sequence_reference）。
     読込時にprovider entryのbytes、frozen record、凍結bait panelの配列SHA、
     search manifestのhits SHAを再検証する。review-onlyでsupportは作らない

route-c-direct-members
  -> verified_direct_support
   | reference_phenotype_review
   | ledger_only
  を分離

supported-member-union
  -> optional A/B/C provenance merge
  -> canonical family/copy
  -> wet-candidate-handoff

route-c-direct-benchmark
  -> source recoveryとprojection successを別denominatorで評価
```

`guided-profile --profile triple-route --execute`には、pre-frozen Direct CをB0 atlas生成後に投影し、member adapterからexisting supported union/Wetへ渡す配線があります。

したがって、Direct Cを「search/projectionだけ実装、orchestrator/Wet未接続」と説明するのは古い状態です。

### ただし、完成したtrusted support laneではない

配線とRoute-C固有binding gateは実装済みです。実装済み範囲と残課題を分けます。

1. [#661]() — **実装済み・close済み**。citation、identity、species/accession、experiment、phenotypeを封印済みsidecarへ結び、member昇格前に検査する。人間または独立agentのレビューから作る手順は [binding producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md)
2. [#663]() — exact/orthology/BLAST projectionのorder invarianceとprovenance completeness
3. [#662]() — close済み。検索局所失敗の隔離と取得済み結果の再利用。strict search lineageは任意

A/C優先（2026-11-08）では、機械検査を通った完全なbindingとprojectionの一致で`verified_direct_support`を認め、人間の別確認を必須にしません。競合・未認可・不完全なrowは`reference_phenotype_review`以下に留めます。binding欠落をcitation/cueだけで補ってsupportへ上げません。#663の残受入（alias・複数copy・順序不変性）は継続して追跡します。

---

## 2. A / B / C の責務

```text
C = 既知を落とさない
    known-answer missを防ぐ

A = phenotypeからbiologyを考える
    semantic / mechanistic missを防ぐ

B = intervalを全部見る
    hypothesis missを防ぐ
```

| lane | 起点 | 問い | 強み | 弱点 |
|---|---|---|---|---|
| Direct C | phenotype + species/cultivarだけ | この形質に関連する既報gene/locus/variantは何か（自由探索） | geneリスト・QTL・A/Bに縛られず、古典座位・変異も回収 | 未知biologyは作れない。出典とidentityの検証は後段 |
| Route A | phenotype | どのmechanism/pathway/family/reference geneが妥当か | 類似形質・機構・他種biologyから外挿 | 正しい機構を思いつかないと落ちる |
| Route B | QTL interval | 区間に何があり何をしているか | 仮説に依存しない区間網羅 | annotation/reference coverageが弱いと意味を読めない |

A/B/Cは3票ではありません。agreementはprovenanceであり、独立したscore bonusではありません。

---

## 3. C と A の境界

境界はspecies distanceではなく**claim type**です。

### Direct Cに入る

```text
Arabidopsis gene A loss-of-function
  -> target phenotypeを直接示す
  -> PMID / DOI verified
```

### Route Aに入る

```text
A is a receptor kinase
B belongs to a related signaling family
therefore B may affect the phenotype
```

- Direct C = already reported factのretrieval + verification
- Route A = new biological hypothesisのinference

Arabidopsisでもdirect mutant relationならCです。Pisum論文でもgeneric functionから形質関与を推論しただけならCではありません。

---

## 4. 検索範囲とevidence directness

検索範囲は広く、主張は厳しくします。

### species distance

```text
S0_same_cultivar_or_accession
S1_same_species
S2_same_genus
S3_close_clade
S4_broader_legume
S5_model_or_other_plant
```

### evidence directness

```text
D0_causal_perturbation_or_rescue
D1_mapped_gene_with_functional_validation
D2_direct_trait_gene_association
D3_context_only
```

species distanceとevidence directnessは別軸です。

- S5/D0 Arabidopsis knockout: 遠いが強いreference seed
- S1/D3 same-species expression: 近いがcontext-only

D3をprojection-eligible direct recordやsupportへ昇格しません。S/Dを1つのopaque numerical scoreにしません。

---

## 5. Direct-C frozen record

Template / implementation:

- [`templates/route_c_direct_rescue.tsv`](../templates/route_c_direct_rescue.tsv)
- [`src/candidate_qtl/route_c_direct_rescue.py`](../src/candidate_qtl/route_c_direct_rescue.py)

主要列:

```text
record_id
query_phenotype
focal_species
focal_cultivar_or_accession
reported_phenotype
reference_species
reference_cultivar_or_accession
reported_gene_id
reported_gene_symbol
reported_locus
reported_variant
species_distance
evidence_directness
evidence_type
source_id / PMID / PMCID / DOI
source_verification_state
gene_identity_state
trait_relation_state
reference_sequence_id
direct_rescue_eligibility
rank_contribution = none
```

freeze artifactには次を入れません。

```text
qtl_id
interval_gene_id
focal_gene_id
candidate rank
truth / answer
```

`focal_species` / cultivarはspecies-distanceを監査するquery contextであり、QTL情報ではありません。

### `eligible_direct_record`

少なくとも以下を要求します。

- source verified
- gene/locus identity stateが許可された状態
- direct target/related trait relation
- D0–D2

これは**projectionへ進む資格**であり、Wet supportへ自動昇格する資格ではありません。

---

## 6. Current search / freeze boundary

```text
phenotype + focal species/cultivar
  -> independent agent + native harness search/fetch (ordinary)
      [optional: route-c-direct-search for literature links]
  -> source-grounded gene/locus/variant proposal (17 columns)
  -> route-c-direct-freeze
  -> route_c_direct_rescue_frozen.tsv
  |  interval-blind up to here
  v
post-freeze projection
```

`route-c-direct-search`はspecies-scope tableに基づく専用の任意検索です。通常Cは、実際に利用できるハーネス標準のWeb検索・取得を使い、別contextで自由に調査します。どちらも検索hitは候補linkであり、evidenceではありません。未実行・検索障害は`not_evaluated`、本当に検索して該当がなければ`completed_empty`と区別します。

**重複・区間外・gene名未同定を入口で除外しない。** 既報の`reported_locus`や`reported_variant`だけでもproposalに記録でき、同じgene/locusの別論文を別`record_id`として保持します。QTLとの関係はfreeze後に評価し、`no_interval_match`を`outside_current_qtl`と同一視しません。A/Bとの重複はfreeze後の統合で処理し、Cが独立調査だった事実を消しません。

現行のfreezeはgeneric `source_verify()`を再利用します。generic cue/term matchingはcitation triageであって、gene/species/trait/experimentが同じclaimに結合した証明ではありません。その後、実装済み#661のbinding sidecarを [既存producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md) で検証・封印し、bound manifestをprojection/memberへ渡します。元freezeを上書きしません。

### #661で実装されたmembership前のbinding

Direct-C membership前には、最低でも別軸で次を確認します。

```text
citation_resolution_state
source_organism_state
source_accession_or_cultivar_state
reported_identity_type/value
reported_identity_verification_state
identity_mapping_basis
experimental_target_gene_state
reported_variant_state
phenotype_match_state
perturbation_to_phenotype_binding_state
source_scope
```

- symbol occurrenceだけで別のsubmitted stable IDをverifyしない
- locus occurrenceだけでgene IDへ移らない
- real quote + agent-supplied relation labelだけでsupportを作らない
- S0/S1はproposal自己申告だけでなくsource metadataからverify/computeする
- variant-only recordは明示的なstate machineが完成するまでreview-only

---

## 7. FreezeしてからQTLへ投影する

```text
route_c_direct_rescue_frozen.tsv
  -> freeze hash verification
  -> exact ID / assembly mapping / orthology / eligible BLAST evidence
  -> route_c_direct_projection.tsv
  -> route_c_direct_evidence_chain.tsv
```

### projection states

```text
inside_current_qtl
outside_current_qtl
no_interval_match
mapping_unresolved
not_projected
```

#### `inside_current_qtl`

focal geneへ解決され、current interval atlas内にある。

#### `outside_current_qtl`

focal genome上のgene/coordinateへ実際に解決され、その位置がcurrent QTL外にある。

#### `no_interval_match`

current intervalに対する評価可能なsearchを行ったが、eligible matchが無い。`focal_gene_id`を主張しない。

#### `mapping_unresolved`

source identityをfocal genomeへ安全に解決できない。

#### `not_projected`

provider/input未実行、またはrecordがprojection-eligibleでない。negativeではない。

### `no_interval_match != outside_current_qtl`

interval-only BLAST no-hitから、homologがgenomeの別位置に存在すると推定してはいけません。

現行standard projectorはinterval atlas中心です。whole-genome/assembly resolverによる`outside_current_qtl`自動判定は未実装です。outside stateはschemaとして存在しますが、standard interval-only runが自動的に生成する保証ではありません。

---

## 8. Projection evidenceの現在のblocker

現行projectorは、同じfocal geneに対するexact、orthology、BLAST evidenceを1つのdictionary slotへ順に書きます。このため後から来たBLASTがexact basisを上書きし、support eligibilityを変える可能性があります。

[#663]()のtarget contract:

```text
append-only mapping evidence ledger
  -> exact / alias / orthology / BLASTを全件保持
  -> deterministic, versioned reducer
  -> selected_projection_basis
     + all_projection_bases
     + ambiguity/conflict state
```

必須条件:

- input row orderで結果を変えない
- exact/aliasをBLASTが消さない
- one-to-manyを保持する
- method not evaluatedとevaluated no-hitを分ける
- corroborationとconflictを明示する

`projection_basis`の採否は既存のmember gate（complete binding＋projection authorization、競合はreview）が機械的に決めます。#663の残受入は継続して追跡します。

---

## 9. `route-c-direct-members` とcurrent promotion rule

[`src/candidate_qtl/route_c_direct_members.py`](../src/candidate_qtl/route_c_direct_members.py) は、projection ledgerをcanonical candidate handoff用に分類します。

### `verified_direct_support`

現行コードでは全て必要です。

```text
species_distance = S0 or S1
evidence_directness = D0 or D1
source_verification_state = verified
direct_rescue_eligibility = eligible_direct_record
gene_identity_state = verified_gene or verified_locus_only
trait_relation_state = direct_target_trait or direct_related_trait
projection_state = inside_current_qtl
projection_basis = exact_focal_id or assembly_or_alias_mapping
focal_gene_id non-empty
```

このrowだけが`supported_member_union`へ入ります。

### `reference_phenotype_review`

代表例:

- S0/S1 D2
- cross-species D0–D2
- orthology/BLAST-only projection
- source direct factはあるがsingle-copy supportへ転移できない

これは`tier2_reference_phenotype`相当の「見る理由」であり、transferred causal supportではありません。

### Source gene とprotein endpointの型付き対応（#663、`CORE_FIX`）

`route-c-direct-project --identity-mapping <TSV>` は任意の追加入力です。
source bindingのgene/locus IDとBLAST ledgerのprotein accessionが異なる場合、
原典に存在する対応を検査し、`members -> read_direct_c_review_candidates -> Wet`へ渡します。
元のgene、protein、native GFF IDは置換せず、BLAST-onlyはreviewのままです。
この入力だけで`verified_direct_support`を作りません。

TSVの列は以下の完全一致です。

```text
record_id
source_identity_field
source_identity_value
target_identity_field
target_identity_value
reference_species
source_namespace
mapping_basis
source_record_path
source_record_sha256
quoted_mapping
reviewed_by
rank_contribution
```

source fieldは凍結recordの`reported_gene_id`または`reported_locus`、
target fieldは同recordの`reference_sequence_id`です。IDは完全一致で検査します。
`mapping_basis`は既存binding契約の`curated_id_map`または`source_record_namespace`、
`rank_contribution=none`です。原典の相対pathはTSVの親directory基準です。
原典SHA-256と実在するquoteを検査し、TSV全体もprojection manifestでhash-boundにします。
members生成前とWet consumerで原典まで再検査します。

現行の原典adapterはUniProt JSONです。`primaryAccession`がprotein endpoint、
`organism.scientificName`がreference species、`uniProtKBCrossReferences`中の明示的な
`database`/`id`がsource gene endpointであることを検査します。
gene namespaceは`TAIR`と`GeneID`を扱い、任意のxref、transcript、自由文中のID共起、
reviewer名だけで対応を認可しません。bare accessionまたは明示的な`UniProt:` namespaceを受理します。

別gene/種/protein、未知namespace、対応の競合、hash/quoteの改変は拒否します。
mapping未供給は従来の欠測扱いで、gene prefixの除去やID推測をしません。
同じsourceに対する複数のBLAST候補は既存ledgerとcopy ambiguityで保持します。
異なるmapping methodの競合は既存review consumerの拒否規則を維持します。
native GFF IDは既存BLAST投影から取得し、一般的なassembly alias producerや他providerの
typed mapping adapterはこの修正では未実装です。#663全受入の完了とは区別します。

### `ledger_only`

- `outside_current_qtl`
- `no_interval_match`
- `mapping_unresolved`
- `not_projected`
- unverified / context-only / otherwise ineligible

candidate universeへ入れず、known-reference ledgerには残します。

### Current safety interpretation

上記は**実装されている分類規則**です。#661のbinding gateは実装済みで、A/C優先ではその機械検査の結果で`verified_direct_support`を認めます。#663の残受入があるため、配線とbindingの実装を完成したscientific guaranteeとは扱いません。

---

## 10. A/B/C merge とWet

current mainではDirect Cのverified memberをexisting `supported_member_union`へ渡せます。

```text
Route A supported members
+
Route B verified members
+
optional Direct-C verified_direct_support rows
        ↓
supported_member_union.tsv
        ↓
canonical family recovery
        ↓
family-to-copy
        ↓
wet-candidate-handoff
        ↓
wet_candidates.tsv
```

unionはrankではありません。

- A-only / B-only / C-onlyを保持
- A+B / A+C / B+C / A+B+Cをprovenanceとして保持
- convergenceでnumeric bonusを加えない
- ambiguous copiesをsingle geneへ強制しない
- review-only / ledger-only C rowをverified supportへ混ぜない

Wet rule:

```text
supported <= wet_budget
  -> all supported = tier1_verified_support
  -> fill with tier2_reference_phenotype
  -> then tier3_family_or_domain
  -> never tier4_dark

supported > wet_budget
  -> bounded verified refinement
  -> re-gate
  -> review/mapping required if still over budget
```

fill rowsは`membership_basis=interval_review_coverage`です。Direct Cに別Top-20 producerはありません。

---

## 11. Guided-profile wiring

Direct-C search/freezeはQTLを開く前に別stageで完了させます。`triple-route`はそのfrozen artifactを受けます。

```bash
candidate-qtl guided-profile \
  --profile triple-route --execute \
  --out-dir <OUT> \
  --case-id <CASE> \
  --focal-species "<SPECIES>" \
  --phenotype "<PHENOTYPE>" \
  --route-c-direct-frozen <route_c_direct_rescue_frozen.tsv> \
  --route-c-direct-freeze-manifest <route_c_direct_rescue_freeze_manifest.json> \
  --route-c-direct-qtl-id <QTL_ID> \
  --route-c-direct-eligible-hits <eligible_interval_hits.tsv> \
  ...Route A/B inputs...
```

規則:

- `--route-c-direct-frozen`
- `--route-c-direct-freeze-manifest`
- `--route-c-direct-qtl-id`

は3点セットです。`--route-c-direct-eligible-hits`はoptionalです。

profileはB0 atlas生成後に`route-c-direct-project`と`route-c-direct-members`を呼び、canonical post-B0 extensionへC member artifactを渡します。

profileはproposal extractionもBLAST executionも行いません。

---

## 12. 既存部品の流用

| need | component |
|---|---|
| 通常の文献探索 | C担当が使えるハーネス標準search/fetch。専用検索コマンドは任意 | 
| PMID/PMCID/DOI verification | `source_verify.py::source_verify()`（A/Bと同じ既存実装） |
| Direct-C schema/freeze | `route_c_direct_rescue.py`, `route_c_direct_pipeline.py` |
| reference bait freeze | `blast_panel.py::freeze_blast_panel()` |
| BLAST execution | existing `blast_evidence.py::run_panel_blastp()` outside Direct-C profile |
| threshold-passing hits | `eligible_interval_hits.tsv` |
| interval inventory | `interval_atlas.py` |
| member gate | `route_c_direct_members.py` |
| family/copy | `canonical_family_recovery.py`, `family_to_copy.py` |
| supported union | `supported_member_union.py` |
| Wet | `wet_candidate_handoff.py` |

新C専用threshold、score、Top-20、family/copy pipelineを作りません。

---

## 13. Literature reproducibility blocker

現行search cacheは完全なimmutable evidence cassetteではありません。

[#662]()が要求するもの:

- provider/request contractを含むcontent-addressed cassette
- raw responseとSHA-256
- original provider retrieval timeとcurrent run timeの分離
- per-query cache/live provenance
- provider failureをquery単位で保持
- live→offline replayのbyte stability
- Direct-C search candidate -> proposal -> freezeのstrict cryptographic lineage
- ordinary manual proposalとstrict benchmark proposalの区別

retrieval hitはlink discoveryであり、evidenceやrank scoreにはなりません。

---

## 14. 旧Route C retrieval

既存のlexical / LSA / neural phenotype retrievalは残します。

```text
reference corpus
  -> retrieval index
  -> phenotype-only retrieval
  -> freeze
  -> interval match
```

用途:

- historical benchmark reproduction
- vocabulary mismatch diagnostic
- semantic retrieval experiment
- retrieval baseline

これは**legacy Route C retrieval**です。current Direct Cのprimary meaningではありません。similarity scoreをsource-backed direct evidenceへ変換しません。

詳細: [`ROUTE-C-RETRIEVAL.ja.md`](ROUTE-C-RETRIEVAL.ja.md)

---

## 15. 現在未実装の範囲

- 内蔵proposal-extraction LLM/agent runner
- 汎用的な自動semantic relation proof、variant-onlyの完全binding（#661で実装したRoute-C固有の封印sidecarとは別）
- Direct-C profile内でのBLAST execution
- whole-genome mappingによる`outside_current_qtl`自動確定
- order-invariant projection evidence ledger/reducer
- 検索局所失敗隔離・取得済み結果再利用の残受入（#662。strict search lineageは任意）
- cross-species evidenceのautomatic support promotion
- variant-only proposalのcomplete verification/projection state machine
- #663/#662を含むrealistic same-species / cross-species / provider-degraded E2E受入の残範囲

schemaやCLIが存在することを理由に、これらを実装済みと呼びません。

---

## 16. Tests

関連tests:

- `tests/test_route_c_direct_rescue.py`
  - schema/template同期
  - S/D vocabulary
  - pre-freeze QTL field拒否
  - D3/unverified promotion拒否
  - projection state discipline
  - `no_interval_match != outside_current_qtl`

- `tests/test_route_c_direct_pipeline.py`
  - multi-species search
  - generic source verification
  - freeze/projection hash改変検知
  - one-to-many projection
  - no-hit semantics
  - source/projection分離benchmark

- `tests/test_supported_member_union_route_c.py`
  - C-only verified support admission
  - A/B/C provenance merge
  - review/outside rowsのunion除外
  - C未供給でもA/B behavior維持

- `tests/test_agents_contract_consistency.py`
  - top-level docsがDirect-C current wiringを未接続と誤記しないこと
  - #661の実装済みbindingと#662/#663の残課題を区別すること

claim bindingの回帰は`tests/test_issue_661_direct_c_binding_gate.py`と`tests/test_route_c_claim_binding_producer.py`を参照します。projectionの順序・provenanceの残受入は#663で管理します。fixture通過を実案件の科学性能とは扱いません。

---

## 17. 一文でのcurrent definition

> **Direct Route Cは、同種・近縁種・モデル植物で直接報告されたtrait–gene/locus/variant関係をsource-backedに回収し、QTLを見ないままfreezeした後、interval atlas・ID・orthology・precomputed BLAST evidenceでcurrent QTLへ投影し、optional member adapterからexisting supported-union/Wetへ接続するlaneである。配線はcurrent mainに存在し、A/C優先ではAと並ぶ通常経路である。機械検査を通ったbindingとprojectionでsupportを認め、cross-species/BLAST-onlyはreview evidenceのまま保持する。**
