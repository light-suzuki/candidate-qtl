# Core contract — canonical workflow と変更境界

#1002の標準T1は、同種の同一・同等または非常に近い具体的表現型の既報と、適格な区間内相同投影に基づく第一段候補です。理由付きの他種配列proxyと複数copyを許容し、全長100%・一意1:1・同一座位確定は要求しません。`same_gene_as_focal=unknown/no` とreview supportを保ち、厳密な同一性assessmentは任意の後段に置きます。既存の参照・配列・source/member/hash検証を保ち、一般domain・短い局所一致・無根拠proxyは新T1にしません。判定の正本は[canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md)です。ユーザー指定によりTier方針の版は増やさず、現行mainの標準とします。

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

この文書は `candidate-qtl` のarchitecture preflightです。実装前に、今回の変更がCORE、optional integration、diagnostic、governanceのどれかを判定してください。

Issue #469 は **closed historical stabilization roadmap** です。現在の作業順序を古いIssue DAGから復元せず、current `main`、`CORE-README.ja.md`、このcontract、active Issue/PRの依存関係を優先します。一方、#469で確立した「core correctnessを新機能より優先する」「truth-aware tuningをproductionへ持ち込まない」というdisciplineはstanding guardrailとして維持します。

---

## 1. このrepoの目的

目的は、QTL区間から原因候補を透明・再現可能に絞り、次のWet実験または追加mappingを決めやすくすることです。

対象にはfine mapping前の生のQTL区間を含める。目的はQTL mappingを加速するため、候補原因遺伝子を実験可能なWet集合へ絞ること。coding variantだけでなく、SV・調節・発現経路で原因となる遺伝子も保持し、同一区間に原因遺伝子が複数あればすべて残す。親系統variantや発現情報は必須入力にしない。追加mapping・シーケンスへ進む中間工程として総費用・時間を改善し、mapping主体の絞り込みとの比較を目指すが、達成済みの性能・費用保証ではない。
候補数は入力QTL境界から導く。合成N=1,000や任意の物理窓へ合わせて区間を切り詰めない。

通常解析は、モデル植物・他植物に蓄積された機能や表現型の知見を対象区間へ橋渡しできることを前提にします。全区間geneの機能解明を完了させる必要はありません。研究の乏しい遺伝子は未解決として台帳に残してよく、調査した上での根拠不足と未読・未実行を区別します。知識不足を非原因の判定にはしません。判断と採否の正本は[canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md)です。

目的ではありません:

- benchmark scoreを最大化すること
- LLMに最終原因遺伝子を断定させること
- familyを必ず1遺伝子へcollapseすること
- 任意のRoute BをA/C correctnessの必須dependencyにすること
- 20件へ任意truncateすること
- route / metric / taxonomyを増やすこと自体

---

## 2. canonical A/C core と optional Route B

### Route A — phenotype-first

`phenotype -> literature/biology -> mechanism/family/reference bait -> verification/freeze -> deterministic interval projection`

- initial hypothesisはinterval-blind。
- source/accession/species/sequenceは検証してからtrusted artifactへ昇格。
- 通常のWeb調査では、一次文献に裏付けられた対象種と同種のreferenceも利用できる。同種という理由だけでは除外しない。明示的な`--benchmark-mode`ではnon-focal制約を適用し、そのpolicyをimport・review session・verifyまで保持する。
- familyから区間内copyへの投影はmachine側で行う。
- 複数copyが残ってよい。
- Route B outputとbenchmark truthはfreeze前に見せない。

Issue #688の承認済み方式として、Route Aは`independent_agent_review`も受理する。
提案者と別実行のagentが同じgene・experiment・phenotype・organismの関係をレビューし、
machineがfreeze済み原文/入力/回答hash、引用位置、providerのID対応と種を再照合する。
これは意味の決定的証明でも人間確認でもない。方式名だけのTSVでは昇格できず、
既存support predicate・accession/sequence・投影・copy gateは維持する。
Direct Cの確認方式やcross-species tier1条件をこの変更で拡張しない。

### Route B — interval-first

`QTL interval -> all genes -> machine annotation/family/orthogroup compression -> source-backed interpretation -> transfer verification -> interval members`

- B0 inventory/compressionはphenotype-free。
- raw hundreds of genesをそのままLLMへ渡さない。
- Route A outputをfreeze前に見せない。
- LLM proseはnumerical scoreにしない。
- group-level supportはplausible interval copiesへ機械的に展開する。

入力容量を超える場合は、全依頼対象を候補・根拠・出典の対応を維持して分割する。分割数・サイズは実行上の設定であり、科学的な候補除外条件ではない。
各分割の提案を凍結して機械的に合流し、分割ごとのTop-Kや固定枠で候補を落とさない。同一候補の複数claim・参照・出典は所属の重複整理後も保持する。
全unitの処理状態を記録し、未処理・失敗・材料不足はnot_evaluated、実際に解釈したうえでclaimなしは別状態とする。部分結果を完全なB実行として扱わない。
圧縮率が低いことだけで実行を止めず、実行・検索・人手・追加実験を含む費用を評価する。

### 実行主体と中間artifact

production coreは単一LLMを前提にしない。Codex/OpenCode等の親AgentがA/B/C、Route B unit、boundary candidateをsubagentへ分配し、結果を既存consumerへ戻してよい。

中間TSV/JSON/manifestは型付きhandoff・監査・checkpointとしてproductionで利用できる。これらの存在を「手作業」や失敗とみなさない。一方、Agentが自由形式の表を作って既存hash/ID contractを迂回することは認めない。

Agentic executionでも次を守る。

- 親Agentが全assigned rosterを所有し、missing/failedは`not_evaluated`として明示する。
- deterministic merge/validatorがunknown/duplicate ID、case mismatch、source/input hash mismatchを拒否する。
- subagent/query/shard単位の障害は局所化し、健全なartifactを保存する。
- executor/model/source/reviewer provenanceとcheckpoint/resume位置を保持する。
- final adopted Wetとreportが正本であり、個々のsubagent回答や中間TSVをcausal truthとして扱わない。

再開時の「元回答・資料を使う」は、追加文献取得禁止や1round終了を意味しない。固定入力での再生はユーザーまたは割当benchmark policyが指定した範囲に限り、benchmarkという呼称や検証器のoffline設定から研究側の検索禁止を推測しない。許可された追加資料は新世代へ凍結し、旧回答・旧freezeを保持する。停止済みrunを自動再開する許可ではない。

検索範囲・round budget・途中終了は親Agentの実行判断として、根拠となる今回の指示、割当roundと実施round、残る具体的な問いを既存plan/receiptへ記録する。予算未消化の終了をbudget exhausted、未取得を外部資料不存在、部分調査を文献調査完了と報告しない。全unknownの追加探索は要求しない。

T1が0の場合、代表の`tier_evidence_json`だけでなく全paper reportsと対象座位の対応を確認し、未探索・同種原報ありだが同一性未確定・Tier算出誤りを分ける。代理参照の種と原報の報告遺伝子の種を混同しない。Type1やtruth一致を理由にT1へ昇格せず、曖昧copy・部分配列の限界を保持する。

研究者へ1000 geneの手分割やA/B/Cの手結合を要求しないことは運用上の目標だが、one-shot CLIや中間ファイルゼロはcore correctnessの条件ではない。

### Wetの確認候補とverified supportの分離

早期の文献・claim検証は、Wetの確認候補をすべて排除するための条件ではありません。
Aは機械で参照名・種・配列を確認して区間へ投影し、Bはfreeze済みの区間memberへの対応を確認します。
そのうえで証拠確認が未完の提案を、検証状態と理由を伴うreview候補として保持します。
この候補は`supported_member_union`へ昇格せず、Wetでは汎用のinterval review fillより先に残します。
20件は運用目標です。候補を失う早期cutoffにはせず、多すぎる場合は後段の精査で絞ります。

取得済みの参照注釈はBの調査材料として保持し、検証済みclaimと区別します。Aの有意な部分一致も、凍結配列・検索原本・区間所属を検証して要レビュー候補に残します。全長coverage未達を支持へ昇格させません。共通基準と1,000遺伝子での確認範囲は [調査材料の保持](INVESTIGATION-MATERIAL.ja.md) を参照してください。
投影に使う参照のID不整合・配列欠損、区間外member、改変されたfreezeを許容する変更ではありません。

論文の意味、表現型の等価性、報告遺伝子と配列参照の生物学的な対応は元研究役のLLMが判断し、理由を記録します。模式例では、強い表現型報告のあるMtAの配列を取得できない場合も、元研究役が資料と生物学的知識からMtAと相同参照AtBの対応を判断し、AtBの配列を通じて区間内PsCを候補化できます。通常の初期候補化にMtA配列の入手、AtB自身の表現型論文、別の本文reviewerを必須にしません。原報のMtA・種・表現型・引用は保存し、AtBを配列参照、PsCを対象遺伝子として区別します。機械は選択したAtBのidentity・配列とAtBからPsCへの区間対応を検査し、この橋渡しをverified supportへ自動昇格させません。資料の扱いは[調査材料の保持](INVESTIGATION-MATERIAL.ja.md)と[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)に従います。

### Route A / B comparison

A/Bは個別freeze後に初めて比較します。

- same interval geneだけならagreement provenance。
- convergenceは追加の独立evidenceでもscore bonusでもない。
- disagreementは片方を自動削除する理由ではない。

### Direct Route C — required independent reported-fact search

Direct Cは、**A/Bの回答にもQTL区間にも依存しない、対象種と表現型だけからの自由な既報探索**を担当する必須laneです。同種を起点に、必要なら近縁種・モデル植物まで調べ、直接報告されたtrait–gene/locus/variant関係（古典座位、variant-onlyを含む）を回収します。Aとの重複・後で区間外となる報告・未同定のgene IDを事前除外しません。現行の形式化は既存17列proposalと後段のsource/identity/freeze/projectionが担当します。研究担当はハーネス標準のWeb検索・取得を優先し、専用`route-c-direct-search`は任意です。新規の一般機構推論はA、区間内遺伝子の機能解釈は任意のBに残します。

current mainの実装:

```text
native harness research (or optional route-c-direct-search)
  -> 17-column proposal -> route-c-direct-freeze
  -> route-c-direct-project
  -> route-c-direct-members
  -> optional supported-member-union
  -> family/copy
  -> wet-candidate-handoff
```

`guided-profile --profile triple-route --execute`は、事前freeze済みDirect-C artifactをB0 atlas生成後に投影してこのhandoffへ接続できます。つまりDirect Cは未接続ではありません。

ただし、Direct Cのautomatic support promotionは次のblockerが閉じるまで**experimental / provisional**です。

- [#661](): Route-C固有のstructured binding gateは実装済み・close済み。source identity・実験・形質・projectionを結ぶ封印済みbindingを要求する。手順は [binding producer](ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md)
- [#663](): exact / orthology / BLAST evidenceをorder-invariantかつprovenance-completeにreduceする
- [#662](): 検索の局所失敗を隔離し、取得済み結果を再利用する。strict search lineageは任意

現行コード上のmember policy:

```text
S0/S1 + D0/D1
+ source verified / direct record eligible
+ sealed source-claim binding authorizes projection
+ inside_current_qtl
+ exact focal ID または assembly/alias mapping
  -> verified_direct_support

S0/S1 D2
または cross-species D0-D2
または orthology/BLAST-only projection
  -> reference_phenotype_review

outside / no interval match / unresolved / not projected
  -> ledger_only
```

契約上の運用:

- #661のbinding producerは`curated_human`または記名した`independent_agent_review`を受け付け、同じquote・citation・frozen identity・retrieved scope検査を行う。
- 完全bindingと認可されたprojectionは追加human reviewなしでsupportへ進む。矛盾・未認可・未完行はreview-onlyとする。
- 通常はA/Cを実施し、検索済み空結果とprovider failure・未実行を区別する。旧A/B凍結世代は元の契約で再生する。
- Direct Cに別score、別Top-20、別family/copy pipelineを作らない。
- legacy lexical/LSA/neural retrievalはdiagnostic/historicalであり、Direct Cのevidenceへ読み替えない。

### Direct Cで現在未実装または外部責務

- 候補論文からproposalを抽出する内蔵LLM/agent runner
- BLAST executionそのもの
- whole-genome resolverによる`outside_current_qtl`自動確定
- cross-species evidenceの自動tier1昇格
- variant-only recordの完全なverification/projection contract

---

## 3. Canonical Wet boundary

base A/Cと任意Bのverified supportを、既存member gateを通して`supported_member_union`へ統合します。一方、A/Bの**初期適格review候補**とCの適格reviewはsupportへ昇格させず別入力でWetへ合流します。family/copyを扱った後のproduction membership正本は`wet_candidates.tsv`です。

通常の`tiered_interval_fill` + `target` handoffは#817の件数非依存policyです。B0 atlasを含む全interval rosterのrelevance台帳と元の適格supported/review gene/copyを保存し、今回のWet集合は以下のv3規則で定めます。20は目安であり、最低件数・切断上限・調査完了の条件ではありません。
旧 `evidence_scoped_v2` の凍結scopeは、最初に25を超えるTierと後続Tierの全適格候補を同じimmutable boundary scopeに含め、独立review consumerを使う元の規則で再生します。新規通常handoffには後述のv3を使います。
8/12件でも件数不足によるlower-tier全件調査を要求しません。調査完了は依頼済み必要作業の実状態で決めます。旧versionなしscopeは旧whole-tier/continuationとして再生し、旧freezeを書き換えません。

`--allow-provisional-wet`は旧whole-tier policyの互換例外です。hash固定済み全区間screenがあり、超過境界が比較根拠のないT4/Uだけの場合に具体的T1–T3と保護Cを暫定採用し、残りを`interval_unresolved.tsv`へ保存します。新通常policyの20未満はそれだけでprovisionalにしません。

strict capacity、`supported_only`、明示的refinementは互換性を保つlegacy diagnostic pathです。
これらは通常のrelevance selectionを置き換えず、T1を落とす機械的shortcutにもなりません。

fill rowは`membership_basis=interval_review_coverage`であり、verified supportではありません。

初期仮説12件・Wet20件は運用上の目安で、件数だけを科学的合否に使いません。
最終15〜30件程度も通常結果であり、ずれの理由説明は不要です。この幅もhard limitではありません。
実験設備などの固定容量を明示的に指定したときだけ`--budget-policy strict`を使います。

新通常handoffは`upper_evidence_scoped_v3`を使います。適格な報告根拠の新T1〜T3（旧T1/T2も対象）があれば上位集合とprotected Cを対象にし、下位は元の適格性を保って台帳へ残します。旧T3は新T3の報告根拠と同一視しません。上位0なら従来の適格poolを使います。累積25を初めて超えるTierだけを既存boundary比較へ渡し、それより上は保持、後続Tierは台帳に残して追加比較しません。active集合は保持prefix・この境界・protected Cです。25は切断上限ではなく、同Tierの順位付けや無理な縮小を要求しません。既存の根拠で比較し、別agentによるTier/本文再審査を追加しません。元eligible集合・active集合・policy版を封印し、完了/未完による下位復帰や追加roundは起動しません。旧`evidence_scoped_v2`・版なしscopeは旧挙動を再生します。request生成はLLM実行完了ではありません。
v3の封印済み境界は元担当の理由付き判断を`original_assessment_v1`としてhash/roster/identity照合で適用し、reviewer欄を空に保ちます。未評価・欠落は保持して履歴未完とし、全未評価を完了としません。旧v2/legacy/genericの独立review契約は維持します。
`claim_review_session.finish_wet_review`はsource hashとscopeを検証し、完了結果を次世代Wetへ接続します。
完了・未完・round budget exhaustedを区別し、第1段の凍結集合を保持します。

- scoreを新しく持たせない
- route flagを持たせない
- 他種referenceの観察事実と、相同性・機構からfocal geneへ橋渡しした候補仮説を区別する。後者を実証済みのtrait transferとは断定しないが、推論であることを理由に候補化を禁止しない
- reviewable tier2/tier3が不足したら`unfilled_budget_slots`を明示する
- dark geneで件数を作らない

`--wet-fill-mode supported_only`はdiagnostic/comparison用です。

`wet_candidate_trace.tsv`はunion/review在籍、bound claim数、relevance、first-stage dispositionを記録する監査台帳です。
`upstream_llm_assessment=not_evaluated_by_handoff`を含み、traceだけからLLM性能やcausal truthを判定しません。
Type3のtrue causal gene全件の最終Wet retention、actual unique gene count、mappingとのcost/time比較を別々に評価します。Type1/2の回収は先行診断、expanded-set retentionとfirst-stage retentionも診断列として分けて報告します。

互換出力には `tier1_verified_support`、`unverified_route_review`、`tier2_reference_phenotype`、
`tier3_family_or_domain`、`tier4_dark` を残します。これらはsupport/membership vocabularyであり、
新規A/BのT1〜T6/U relevanceにも、旧凍結世代のT1〜T4/U relevanceにも同一視しません。

---

## 4. Evidence invariants

以下はarchitecture invariantです。変更には明示的なIssue、理由、tests、docsが必要です。

- `not_evaluated != negative`
- provider/search/measurement failure != biological absence
- annotation absence != sequence absence
- generic homology / conservation = context。trait-specific supportへ自動昇格しない
- `paralog_ambiguous` = epistemic/copy uncertainty。automatic counterevidenceにしない
- source不明・未検証claim = trusted supportにしない
- LLM/model memory = citable sourceではない
- LLM prose = numerical scoreではない
- route convergence = score bonusではない
- benchmark truth = production rule inputではない
- interval no-hit != outside-QTL mapping
- cross-species direct fact + homology != focal-gene causal support
- symbol/cue occurrence != stable-ID/species/trait/experiment binding
- a real quoted span != automatically verified relation label

### Direct C claim-binding invariant

source verificationのcitation existenceとgeneric cue triageだけでcanonical membershipを作りません。

最低でも別軸として保持すべきもの:

```text
citation resolution
source organism/accession
reported identity type/value
identity verification/mapping basis
experimental target gene
reported variant
phenotype match
perturbation -> phenotype binding
source scope
```

1つのsymbol occurrenceが別に入力されたstable IDをauthorizeしてはいけません。locus occurrenceからgene IDへ移る場合もexplicit mappingが必要です。

### Direct C projection invariant

projectionはmethodごとのevidenceを上書きせず保持し、deterministic reducerでdecisionを作ります。

- verified exact/alias evidenceを後続BLASTが消してはいけない
- exact、orthology、BLASTの全provenanceを保持する
- one-to-manyをsingle copyへ潰さない
- method not evaluatedとevaluated no-hitを分ける
- input row orderでartifact bytesやmembershipが変わらない

これは#663完了後の必須contractです。現行実装を完成済み保証として記述しません。

---

## 5. Lane coverage — 実行したものだけを主張する

production base coreはRoute A / Direct Cとfreeze後のmerge/Wet handoffです（A/C-first）。Route Bは任意です。未実行Bを候補やcoverageの入力へ含めません。

責務laneを実行していないのに、そのlaneが必要なUse Case Typeの性能を主張してはいけません。単独の`auto-triage`やscreening shortlistをA/B/C全体の性能として報告しません。

### `auto-triage` の位置づけ

`auto-triage`はpanel-freeなscreening入口です。

- Route Aのmechanism/source/bait freezeを行わない
- Route Bのcanonical interval-atlas→semantic→transfer verificationを完遂しない
- Direct Cのdirect-report search/freeze/bindingを行わない
- A/B/C freeze後mergeを代替しない

### Use Case Type と責務lane

| Type | 回収責務lane | 必須/主要seam | `auto-triage`単独で評価可能か |
|---|---|---|---|
| Type 1: 調査漏れ、またはidentifier / assembly mismatch | Direct C、Route A、identity・liftover・homology照合 | freeze/projectionと既存verification | 不可 |
| Type 2: 表現型記述または観察スケール（micro/macro）の不一致 | Route A multi-mechanism、Route B semantic interpretation。retrievalは診断 | mechanism import、route-b-semantic-* | 不可 |
| Type 3: 対象植物の原因geneは未知だが、homolog機能や他植物の知見から推論可能 | Route A transfer、Route B family/orthogroup、Direct C direct factはreview seed | mechanism import、route-b-semantic-*、optional Direct-C ledger | 不可 |
| Type 4: 真に種特異的で新規 | mapping / dark / annotation rescue / Wet | case dependent | 保持保証なし |

Direct Cを実行した場合も、current experimental promotionとreview-only evidenceを別々に報告します。

---

## 6. Family / copy contract

family、orthogroup、subfamilyは中間表現です。

```text
family recovered
  -> validated copy set
  -> orthology / subfamily / domain / synteny / expression / position
  -> preferred copy OR multiple retained copies
```

- family-level supportだけでsingle copyを強制しない
- one-to-many / many-to-many / unresolvedは集合として保持する
- `not_in_validated_copy_set`を単なるmissingnessとして自動昇格しない
- family/copy resolutionはcanonical semantic pathを使い、deterministicなcopy判定を維持する

---

## 7. Change class

PR作成前に必ず1つ選びます。

### `CORE_FIX`

existing canonical/optional integrated behaviorがcontractどおり動かないbugを直す。

例: #661 claim binding、#663 projection overwrite。

### `CORE_INTEGRATION`

既存部品をcanonical pathへ接続する。新Route、score、evidence meaningを作らない。

### `BENCH_ONLY`

fixture、scorer、benchmark report、audit、frozen runだけを変更する。

### `EXPERIMENTAL`

optional provider、new retrieval experiment、adaptive research stage等。A/C core要件にしない。

### `GOVERNANCE`

docs、agent contract、repository organization、change-control、contract consistency testだけを変更する。production behaviorを変えない。

---

## 8. Stabilization guardrails

active architecture Issueがない限り、次を追加・変更しません。

- new Route
- new top-level CLI command family
- new benchmark family / headline metric
- new adaptive stage
- ranking weight / threshold tuning
- benchmark truthを見て追加したsynonym / alias / special case
- evidence direction / trait-specificityの意味
- Direct CをA/Bの必須dependencyへすること
- familyを常に1 geneへcollapseするpolicy
- Direct Cのblockerをdocsだけで解消済み扱いすること

---

## 9. Agent stop conditions

次のどれかなら、隣接機能を勝手に実装せずscope conflictをIssue/PRへ明記します。

1. active Issueのacceptance criteriaを越える変更が必要。
2. new route / CLI / metric / taxonomy / provider / adaptive stageが必要。
3. weight / threshold / prompt / synonym / family alias変更が必要。
4. benchmark truthを見た後でproduction behaviorを変えたくなった。
5. missing dataをnegativeとして扱わないと実装できない。
6. source-backedでないtextをsupportにしないと成果が出ない。
7. core bugをfixture側special caseで隠したくなった。
8. unit testは通るがreal CLI/consumerから到達できない。
9. Direct Cのcitation/keyword triageをclaim bindingとして扱わないとpromotionできない。
10. projection evidenceを上書きしないとschemaへ収まらない。

---

## 10. Minimal-diff rule

core PRは原則として:

- 1つのroot causeを直す
- acceptance criteriaに必要なfilesだけ触る
- regression testを追加する
- opportunistic refactorを混ぜない
- benchmark headlineを同じPRで更新しない

GOVERNANCE PRは、production behaviorを変えず、矛盾する上位文書をまとめて同期し、再発しやすい主張をcontract-consistency testで固定します。

---

## 11. Current priority

現在の優先判断:

1. **A/C core correctness / reachability**
2. **Direct C safety blockers** — #663 projection determinism（#661のRoute-C binding gateは実装済み）
3. **Direct C reproducibility** — #662 局所失敗隔離・取得済み資料再利用（strict lineageは任意）
4. **realistic E2E validation** — same-species support、cross-species review-only、provider-degraded
5. **outside-QTL/global resolver decision**
6. **core integration debt / legacy isolation**
7. **repository cleanup**
8. **new benchmark generation after contract freeze**

#469、#426等はhistorical contextであり、current work controllerではありません。

---

## 12. Validation before merge

- targeted unit/contract tests
- `tests/test_agents_contract_consistency.py` for top-level docs/AGENTS changes
- generated CLI reference check when CLI changes
- source and installed reachability gates for integration changes
- missing-provider/non-negative states
- input/freeze/output hashes
- no unrelated diff

Hosted GitHub Actionsのavailabilityはbiology/code correctnessそのものではありません。実行できたvalidationを正確に記録します。
