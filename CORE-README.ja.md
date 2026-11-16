# CORE-README — このrepositoryが目指すもの

この文書は `candidate-qtl` の目的と全体構成の説明です。Tier・Wet採否は[canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md)（`docs/CANONICAL-WET-MEMBERSHIP.ja.md`）、変更境界は[core contract](docs/CORE-CONTRACT.ja.md)を正本とします。通常の新規A/Bは`initial_bibliography_v1`を使い、旧凍結世代は元の規則で再生します。

## 1. 一言でいうと

> **QTL区間と表現型から、真の原因候補を失わずに、実験可能な約20遺伝子のWet検証候補集合へ整理するflow。**

**支持確定をWet採用の必須条件にしません。** 初回mappingの生のQTL区間も受け取り、追加mapping・シーケンス・Wetへ進む前の候補選定を支援します。mappingは独立した証拠として保ち、研究全体の時間・費用を減らすことが目標です。達成済みの性能保証ではありません。

**探索の前提は、モデル植物や他植物で研究された遺伝子の知見を、対象植物の相同遺伝子へ活用できることです。** 全geneを台帳に保持し、初期判断の範囲を記録しますが、全遺伝子の機能解明を完了条件にはしません。研究の乏しい遺伝子は根拠不足を明示して未解決のまま残せます。

**20件は目安であり、最低件数・切断上限・調査完了の条件ではありません。20件を超える既存候補は数字だけで切らない。** 必要な確認を終えて12件なら12件でよく、20との差分だけで追加探索を要求しません。既に挙がった弱いreview候補・別機構・曖昧copyは広く残し、本命と同じ直接証拠や長い機構説明を一律に要求しません。

主な評価軸は、**実際に採用したWet集合での正解保持と実遺伝子数**です。以下は考え方の例で、実測結果ではありません。

- 本命5件に正解がなくても、補充15件に正解が残った → Wet保持として成功。supportとは読み替えない。
- supported 3件だけを返して正解を落とした → 支持確認ができても候補保持には失敗。
- 正解を含む3 copyが残った → 保持として成功。3遺伝子と数え、単一copy確定は要求しない。
- 区間の全遺伝子を残した → 保持はしても、小さな実験集合にする目標を満たさない。

### 1.1 base productionの範囲

base productionは、**参照genome/GFF3、QTL区間、表現型だけで開始できるWet前の候補整理**です。map-based cloningを置き換えず、Wetで直接調べられる約20–25遺伝子へ区間を圧縮したところを一つの終了点とします。

次はbase productionの必須条件にしません。親系統の配列/VCF・親間変異、SV/PAV/TE/調節解析、発現データ、組換体の追加やfine mappingの反復、全geneの厳密な1:1 orthology、単一copy確定、annotation rescue、複数原因遺伝子モデルの自動確定。これらは候補を得た後に研究者が必要に応じて進む後段です。

Route Bの第1段では、各区間geneへ低コストな**identity handle**（複数種のprotein DBへのBLASTで「何に似ているか」）を付けることを優先します。BLAST類似度はidentity/投影の材料であり、形質の支持ではありません。参照検索は不足する材料を補うfallbackで、既存の参照対応・family/domain・機能注釈を先に使います。詳しくは[調査材料の保持](docs/INVESTIGATION-MATERIAL.ja.md)。

### 1.2 relevance Tier

Tierは support確認、mappingの確かさ、copyの曖昧さではなく、**今回の形質との生物学的関連度**だけを表します。原調査担当が論文ごとに9項目のyes/no/unknownと理由を記録し、決定表がTierを計算します。適用できる最も強い論文が遺伝子のTierを決め、論文数での加点はしません。

- T1/T2/T3：具体的な遺伝子–表現型の既報（T1は同種の同一・非常に近い表現型 + 妥当な区間内相同投影）
- T4/T5：報告機能・出典付き注釈に基づく仮説
- T6/U：identityのみ・未解決。全遺伝子台帳に保持し、非原因とはしない

次はTierとは別の軸です。

- support / verification：`verified / review / context / not_evaluated`
- mapping / identity quality：`family_candidate_bridge.mapping_tier`
- family resolution：`superfamily -> family -> subfamily/orthogroup -> specific homolog -> focal copy set`
- copy state：`preferred_copy / co_preferred_tie / paralog_ambiguous / tandem_ambiguous / ...`
- route provenance：A / C / 任意のB

unique/ambiguousをTierの強弱へ変換しません。T2 familyに区間内3 copyがあり安全に区別できなければ、3 geneともT2として残し得ます。条件の詳細は[決定表](docs/CANONICAL-WET-MEMBERSHIP.ja.md#初回abのチェック表と軽いid照合paper_checklist_v1)に集約します。

### 1.3 実行モデル — Agent orchestration

実行主体は単一のLLM呼出しに限りません。**親Agentがcaseを所有し、Route A/B/C、Route B shard、Wet境界レビューをboundedなsubagent taskへ分配してよい**ものとします。中間TSV/JSON/manifestは削除対象ではなく、Agent間の型付き通信・監査・checkpoint・resumeのartifactです。研究者がそれらを手で作る必要はありません。

必須なのは次の境界です。

- 全assigned IDをrosterとして保持し、未回答・失敗を`not_evaluated`として残す。
- unknown/duplicate ID、case不一致、source/input hash不一致を機械側で拒否する。
- executor/model/source/proposal/reviewの来歴を保存する。
- 局所失敗で健全な結果を捨てず、有限budget後もcheckpoint/resume位置を残す。
- 最終的に採用Wet集合と読める報告へ到達する。

科学担当は親の会話履歴を引き継がない新しいcontextで起動し、許可した入力artifactを記録します。運用の正本は[通常案件の手順](docs/AGENTIC-QTL-CASE.ja.md)、入口は`workflow/practical_packet_run.py`（`prepare` → `review` → `finish`）です。

## 2. canonical production path

```text
表現型 + QTL区間
  ├─ A：表現型から機構・参照遺伝子を考え、区間へ投影（必須）
  ├─ C：既報のtrait–gene/locus/variantを区間を見ずに検索・投影（必須）
  └─ B：区間の全遺伝子を整理して表現型との関係を考える（任意）
  -> 各経路を独立にfreeze
  -> claim review -> 論文ごとの9回答 -> Tier T1〜T6/U
  -> 支持済み部分集合（supported_member_union）と明示reviewを区別
  -> 全gene台帳を保存 / T1と既存の適格supported・reviewを引き継ぐ
  -> wet_candidates.tsvをfreeze
       <=25：境界レビューなしでWetへ
       >25：累積で最初に超えたTierだけを境界レビューへ、後続Tierは台帳
```

A/Cがbase coreで、Bは任意です。Cは検索を実施し、検索済みの空結果と未実行を区別します。A/Bの明示reviewはsupported unionへ昇格せず、同じWetへ別入力として合流します。source identity・種・引用と投影台帳を検査できたDirect Cのreview候補は、件数や低いTierを理由に落としません。

Bを使わない場合、WetにはA/Cの候補だけが入り、区間注釈で数を埋めません。候補が0件または既定3件未満のとき、handoffは`route_a_rethink`を記録し、親がAを新しい世代で一度考え直します。この合図は採用集合を増減させません。

**LLMに期待する中心的な仕事は、他種の文献と対象種の遺伝子を橋渡しする仮説形成です。** Aでは「他種のこの遺伝子の報告が対象の表現型に似ている → 対象種の相同遺伝子が候補」、Bでは「この区間遺伝子は他種のこの遺伝子・familyに対応する → 報告された作用から今回の表現型を説明できそうだ」と考えます。報告された事実と対象種への推論を分けて記録します。特定のモデル種・形質に専用の選抜規則は設けません。

模式例として、原論文の報告遺伝子X → 配列を使う代理参照Y → 対象区間の遺伝子Zという橋渡しも通常の仮説形成に含みます。XとYの対応と表現型の意味は元の調査担当が理由付きで判断し、Yのidentity・配列とZへの投影は既存の機械検査に通します。詳しくは[報告遺伝子・代理参照・対象座位](docs/CANONICAL-WET-MEMBERSHIP.ja.md#報告遺伝子代理参照対象座位を分ける)。

計算量は主に区間内の遺伝子数に依存し、全ゲノム資産は前計算して再利用します。実QTLを任意の窓へ切り詰めません。

## 3. Route A — phenotype-first（`LLM or agent -> machine`）

`phenotype -> mechanism / pathway / reference gene / family hypothesis -> source / accession / sequence verification -> interval projection`

LLMが仮説を考え、機械が参照identity・配列・区間対応を検査します。初期仮説は区間を見ずにfreezeします。claim未確認でも対応を確認した候補はreviewとして残し、未確認の文章をverified supportにしません。[claim verification](docs/CLAIM-VERIFICATION.ja.md)の`independent_agent_review`も使えますが、`curated_human`へ読み替えません。

## 4. Route B — interval-first（`machine -> LLM or agent -> machine`、任意）

`QTL interval -> all-gene inventory -> family / orthogroup / domain / reference compression -> semantic interpretation -> verification -> interval members`

B0は表現型を見ずにfreezeします。family/subfamily/orthogroupを**中間表現**として資料を**圧縮**し、数百の生の遺伝子をそのまま長いpromptに詰めません。容量を超えたら全unitを保持して分割し、chunkごとのTop-Kや圧縮率で切りません。未処理・失敗・材料不足と、調査済みの空claimを区別します。

通常Bは相同遺伝子の既知機能を調べてから、表現型との整合性を中立的に判断します。Web検索や複数論文の利用を制限せず、取得した資料は新しいbundleへfreezeしてから根拠に使います。詳細は[canonical Route B](docs/CANONICAL-ROUTE-B.ja.md)。

## 5. Direct Route C（必須）

Cは直接報告されたtrait–gene / trait–locus / trait–variantを回収する経路で、一般的な機能から新しい機構を推論するA/Bとは分けます。同種・近縁種・モデル植物の報告を調べられます。

`独立したC担当の調査（または任意のroute-c-direct-search） -> route-c-direct-freeze -> route-c-direct-project -> route-c-direct-members -> Wet`

supportへの昇格には、記名した`curated_human`または`independent_agent_review`が機械検査（凍結sourceの逐語引用・引用の実在・凍結identity・取得範囲）を通した封印済みbindingと、認可された投影が必要です。種をまたぐ対応、orthology/BLASTのみ、矛盾する投影はreviewに留め、区間外・no-hit・未解決は台帳に残します。区間no-hitと`outside_current_qtl`を区別します。詳細は[Direct C](docs/ROUTE-C-DIRECT-RESCUE.ja.md)と[binding producer](docs/ROUTE-C-CLAIM-BINDING-PRODUCER.ja.md)。

## 6. A/B/Cはfreeze前に混ぜない

Aの機構想起、Bの区間情報、Cの既報回収は異なる失敗様式を持ちます。各経路をfreezeしてから比較する**procedural independence**を保ちます。freeze後の一致は来歴であり、独立証拠数やscore bonusではありません。旧来のlexical/LSA/neural retrievalは診断・歴史的なもので、現在のDirect Cとは区別します。

## 7. familyは答えではなく「橋」

family / orthogroup / subfamilyは圧縮、種間の橋渡し、IDへの接地、copy曖昧性の保持のための中間表現です。`family -> one gene`は必須ではありません。orthology・domain・synteny・発現・位置はcopy比較の材料で、未測定はnegativeではありません。区別できないcopyは集合で扱います。

## 8. Wet review budget

補充・採用集合・支持状態の規則は[canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md)を正本にします。

- 通常handoffは全gene台帳を保存し、T1と既存の適格supported・reviewを引き継ぐ。報告根拠のある上位Tier（T1〜T3）があれば、その上位集合と保護されたDirect Cを今回の対象にし、下位は台帳に残す。20未満でも件数だけを理由に追加調査や補充をしない。
- 累積25を初めて超えるTierだけを境界レビューへ渡す。上位は保持し、後続Tierは台帳に残して追加比較しない。25は切断上限ではなく、同Tier内の順位付けや無理な縮小を要求しない。
- 境界の採否は元担当の理由付き判断を機械照合して適用する。保留（defer）は非原因の証明ではなく、明示的に保留した遺伝子を通常の補充で黙って戻さない。
- 固定の設備容量があるときだけ`--budget-policy strict`を使う。

`wet_candidate_trace.tsv`は全geneのunion/review在籍、bound claim数、relevance、第1段の扱いを記録します。これだけでモデルの能力や原因を診断するものではありません。

## 9. LLMに任せること / 任せないこと

再現したいのは「本命は数件だが、自分の仮説が外れたときに備えて約20枠を使う」という研究判断です。

機械は計算・全候補の整理・ID/種/配列/出典/引用箇所/hashの照合・重複除外・件数確認を担当します。LLM/agentは機構推論、文献の意味理解、copyや別仮説の比較、理由付き採用/保留案を担当します。機械照合の成功を生物学的な証明とせず、LLMが強い説明を書けないだけで候補を落としません。

文章だけでfreezeを変更せず、数値のcausal score、未確認の唯一copy、provider障害による除外、架空の引用を作りません。研究者が記録した遺伝子名/IDをDB検索の値で上書きしません。詳細は[役割分担](docs/LLM-BOUNDARY.ja.md)。

## 10. mappingは最後まで独立した証拠

このworkflowはQTL mappingを加速し、区間を実験可能なWet候補へ狭めるためのものです。coding variantに限らず、SV・調節・発現を含む因果経路を保持し、同じ区間に原因遺伝子が複数あればすべて残します。

map-based cloningは置き換えません。候補の生物学的根拠が弱い、copyを分けられない、Wetで否定された場合は、組換体・追加mapping・marker・親配列や発現・annotation rescue・SV/TE/CNV/PAVの調査へ戻れます。`Top candidate != mapping complete`です。一方で、fine mappingを極限まで終えなければWetへ進めない、という条件にもしません。

## 11. repository内の階層

| Layer | 責務 |
|---|---|
| CORE | Route A、Direct C、freeze、supported union、family/copy、Wet handoff |
| OPTIONAL | 区間全体のRoute B、Direct C専用検索、第2段review、追加evidence |
| SCREENING | `discover`、`auto-triage` |
| DIAGNOSTIC | budget portfolio、ablation、旧retrieval |
| HISTORICAL | 過去の凍結結果・旧実装の来歴 |

ファイルがpackageに存在することと、通常経路から使われることを区別します。ベンチマーク用データと評価ツールは公開版に含みません。

公開版の構成：`src/candidate_qtl/`（パッケージと`candidate-qtl` CLI）、`workflow/`（通常案件の入口とテンプレート）、`scripts/`（親Agent用補助とローカルCI）、`docs/`、`templates/`、`config/`、`tests/`。

## 12. 変えてはいけない科学的invariant

`not_evaluated != negative`、annotation absence != sequence absence、interval no-hit != outside-QTL、paralog ambiguity != counterevidenceを維持します。generic homology、種をまたぐ既報、経路の一致、LLMの記憶は、それだけで対象遺伝子の因果的支持になりません。実在の引用やquoteがあるだけで、gene/species/trait/experimentの関係が確認できたとはしません。評価用の既知の答えを、production条件・alias・閾値へ持ち込みません。

## 13. 現状と限界

補充、review込みの採用・保留、適格なDirect C reviewの保持は接続済みです。support確認と実行coverageの表示整理、未見事例での選定性能の確認は残ります。安全に保持することと、有効な候補集合を選べることを同一視しません。

## 14. 何を目指していないか

支持確定件数やTop-1の最大化、単一copyの強制、LLMへの無検査のランキング丸投げ、20への切断、捏造・重複による水増し、provider不在の生物学的解釈、新しいRouteを増やすこと自体は目的ではありません。**正当な補充で約20枠を使うことは目的に含まれます。**

## 15. 文書の読み順と責務

作業前にこの文書、`docs/CORE-CONTRACT.ja.md`、最も近い`AGENTS.md`、必要な実装/runbookの順に読みます。Wet規則は`docs/CANONICAL-WET-MEMBERSHIP.ja.md`へ戻ります。

| 文書 | 責務 |
|---|---|
| このCORE-README | 目的、全体像、階層、未完の境界 |
| README / CLI / runbook | 現在の入口と操作 |
| AGENTS / CORE-CONTRACT | agentの実行契約・科学的な変更境界 |
| CANONICAL-WET-MEMBERSHIP | 補充・採用集合・支持状態の正本 |
| 経路ごとのdocs | source / projection / discoveryの詳細 |

同じ規則を別文書へ書き足すより、担当文書へリンクします。
