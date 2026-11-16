# Canonical Wet membership

初期の `paper_checklist_v1` / `initial_bibliography_v1` 記録では、掲載種別・`source_excerpt`・`source_location` は取得できる場合に残す任意項目です。これらだけの不足では補完再調査を要求しません。値・欠測・明示unknownは表に保全し、`missing_fields_json` と `record_state=unfilled` が残っていても、必須項目が揃えば補完依頼は不要です。論文ID/タイトル、報告遺伝子名・名前空間付きID（取得不可なら理由付きunknown）・種、参照対応、9回答と判断理由は元の研究LLMが記録し、機械取得した値で原記録を上書きしません。必須項目の未記入・具体的矛盾・上位Tierへの具体的手がかりは従来の限定補完対象です。既に確認した取得不可IDを自動追加探索義務やTierペナルティにせず、旧方針・明示 `content_v1` の厳密確認は維持します。

機械照合は意味を推測する曖昧一致を使いません。DOI接頭辞/URL、PMID接頭辞/公式PubMed URL、同じデータベースの既知の名前空間表記、Unicodeの表示幅・空白、明示的な回答/enumの表記だけを正規化します。遺伝子IDの名前空間とアクセッションの版は保持し、共通の遺伝子名があっても異なる明示IDを一致させません。unknown/不明などの欠測語を実遺伝子として扱わず、nullの回答は未記入として差し戻し、明示的不明には確認理由を要求します。タイトルのHTML書式・実体参照は除いて既存の正規化タイトルで照合し、別のタイトルや別の文献IDは同一化しません。Bで同じ凍結済みIDの表示形式だけが異なるときは、既存bundleの正確なmaterial keyを再利用し、新しいIDは従来どおり未凍結として扱います。これらは同じ事実の形式差を吸収する処理で、原文の意味判断や未知の遺伝子対応を機械が推測するものではありません。

通常A/Bの論文根拠は `paper_reports_json` を正本とし、取込時に `paper_evidence_table.tsv` を同じデータから生成します。1提案・1論文・1報告遺伝子ごとに、論文ID/タイトル、原著/レビュー、報告遺伝子名・種・名前空間付きID、参照遺伝子IDと対応関係、観察事実 `source_excerpt` と確認箇所 `source_location`、9つの判定回答、自動Tierを記録します。Aの区間候補IDは盲検投影前には利用できず、Bでは凍結unitのmemberを対応情報として示します。報告遺伝子、配列取得用の参照、区間copyを混同しません。

判定回答の未記入・空文字は `paper_checklist_missing_fields_json` に残し、Tier計算上はunknownのままとします。明示的なunknownは調査結果であり、確認先と不明理由を `unknown_reason` に記録します。論文に遺伝子IDがなければ `reported_gene_id=unknown` と理由を残し、既知の名前・種・対応を保全します。ID欠落や記録不足を自動的なT5確定・候補除外・支持昇格には使いません。表の `computed_tier` は元のチェックから計算した値、`applicable_tier` は参照対応が使える範囲の値、`tier_state` は記録不足・対応未解決・評価済みを分けます。

空欄と、上位Tierを変え得る具体的な手がかりだけを `paper_followup_requests.json` にまとめ、元の調査担当へ原則1回の限定バッチとして返します。単にunknownを全てyesと仮定して上位可能性を作らず、表現型の記載済み事実または元の担当が記録した `followup_question` を使います。確認後は `followup_state=completed` とし、未解決なら理由付きunknownを受け入れます。依頼生成は実調査の実行ではありません。回答は元payload/responseの新しいコピーを修正し、新出力先へ再取込します。新B文献は従来のbundle extension/freezeを経由し、既存freezeや独立A/B入力を変更しません。未修正なら `recording_state=needs_completion` を残し、既存の候補採用条件とは分離します。

通常A/Bでは、元の参照遺伝子または共有できる機構ごとに約20件の異なる文献候補をタイトル・要旨から調べます。各論文の元の観察・9項目の回答から初期Tierを計算し、短い理由とともに記録します。その初期Tierが高い方から有用な5報を選び、次段へ渡す主な根拠論文として記録してから、同じ集合のID・タイトルを書誌照合します。DOI/PMIDなど同じ論文の別IDは1報です。20件は全interval geneごとの全文読解数ではなく、共有調査をcopy間で再利用します。選択時の初期Tierと、後段で参照・区間geneへ適用したTierは区別します。

5報は次段へ渡す根拠論文数の目標です。高Tierの根拠を優先しながら各報の観察・対応・理由を保持し、ある1報に誤読や不適合があっても他の論文を参照して再検討できるようにします。原著とレビューを含め、T1/T2の根拠をT5の注釈より優先します。強い論文が取得不能・不適合なら別の取得先や別文献も調べます。1報見つけただけで調査完了とせず、5報に満たなければ実数、調べたID・検索・取得失敗・不足理由を残します。件数合わせ・論文数の加点・不足による候補の機械的除外はしません。初回の論文Tierがすべて誤っていた場合は、機械的な5報選択だけでは修復できません。

元の全論文記録・観察・反証・選外理由を取込と各経路の凍結監査出力に保全し、gene Tierは既存の候補対応・最強適用根拠の規則を使います。次段向けの選択のために、その規則の適用範囲を狭めたり広げたりしません。各論文記録には`paper_handoff_policy`、`paper_handoff_scope`、`paper_handoff_state`（`primary` / `audit_only`）、`paper_handoff_rank`を記録し、既存の`paper_evidence_table.tsv`にも選択状態を表示します。`audit_only`は次段の主資料に選ばれなかった記録であり、生物学的な否定や候補除外ではありません。第2段requestを作る場合は、その`primary_paper_evidence`へ選択済み5報の根拠を範囲ごとに渡します。候補に対応付いた根拠はWetの`relevance_evidence_json`にも保持しますが、Bの選外claimを自動的に候補の適用根拠へ追加しません。元の全記録は上流の監査出力から再検討できます。

`source-verify`の軽量経路は、この同じ集合に既存のcache-first書誌確認を適用し、申告タイトル、取得タイトル、DOI、照合状態を根拠表へ結合します。選択と選外理由を既存の`bibliography_selection.json`にも記録し、表のための追加取得は行いません。報告遺伝子の名前/ID監査も既存の確認を再利用します。種・報告遺伝子・参照IDが同じ対応を再利用しても、各論文の観察事実と判定回答は個別に保持します。通常CLIは軽量確認を既定とし、明示的な`--verification-policy content_v1`と旧凍結policyの厳密確認は維持します。JSONのtrue/falseとyes/no、同一論文の別IDの配列と区切り文字列は取り込み時に正規化し、意味が同じ回答を形式差だけで失いません。未確認の回答・不明な遺伝子対応は引き続きunknownとして残します。

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

Productionで「最終候補集合」と呼ぶものは **`wet_candidates.tsv`だけ**です。
この文書をWet補充・採用集合・支持状態の正本とし、他の入口文書で別の選抜規則を定義しません。

2026-09-21時点で、実Wet書出しには新relevance属性と`relevance_evidence_json`を含めます。
保存済み原claimと追加入力を再利用し、同一geneの先頭claimだけを採る処理を廃止しました。
支持済みgeneと重なるreview根拠も残し、支持状態や候補数には加点しません。
供給された区間全geneの`wet_relevance_inventory.tsv`は監査台帳であり、実採用集合ではありません。
A/C優先（2026-11-08）でBを供給しない場合も同じ規則です。適格集合はA/Cのsupported・review候補だけで、区間注釈による補充はしません。採用件数が0なら`route_a_rethink.state=required`、`--min-route-candidates`（既定3）未満なら`recommended`をmanifestへ記録し、親がAを新しい世代で一度考え直します。この合図は採用集合を増減させません。
通常の`tiered_interval_fill` + `target` handoffは、#817の`upper_evidence_scoped_v3`を使います。全rosterのrelevance台帳と元の適格supported/review集合を保存し、今回のactive scopeを分けます。適格な報告根拠の新T1〜T3（旧T1/T2も対象）または無条件保護T1があれば、その上位集合とprotected Cだけを今回の採用・比較対象にします。旧`T3_functionally_plausible`は新T3の報告根拠と同一視しません。上位0なら従来の適格scopeを使います。

下位のupstream適格性は消さず、`not_selected_upper_tier_priority`として全台帳へ残します。生物学的陰性、独立比較によるdefer、読了とは異なります。上位が少ない場合はその数で終了し、20への補充・下位追加比較は行いません。累積25を初めて超えるTierだけを既存の境界比較へ渡し、それより上は保持、後続Tierは`not_selected_after_first_boundary`として台帳に残して追加比較しません。25以下への切断、同Tierの順位付け、無理な縮小や追加検索を必須にしません。Tierや表現型を別agentで再審査しません。scopeは元`eligible_gene_ids`と今回の保持prefix・境界・protected Cの`active_gene_ids`、policy版と全inventory hashを封印し、finish/partial後も下位を自動復帰せず追加roundを起動しません。新しい根拠による訂正は新generationで計算します。旧`evidence_scoped_v2`・版なしscopeは元の挙動で再生します。
protected Cは上位優先の発動集合に含めません。Cの保護追加だけで通常A/Bの適格下位候補を押し出さず、通常A/B側に上位がある場合もCは別枠で保持します。
`#750`のboundary requestと`#758`のfinishを再利用し、policy版と適格IDをhash付きimmutable scopeへ保存します。旧v2・版なしscopeはそれぞれ元のscope/continuation規則で再生します。
strict capacity、`supported_only`、明示的refinementは旧診断経路として残ります。
通常target + atlas handoffでは、boundary scopeのlineageを固定するため非空の`case_id`を要求します。

上位Tier優先は、元の調査LLMが報告された表現型・対象種・gene対応を判断し、その記録から決定表でTierを計算するという分担に基づく採用方針です。評価では、実際のTier誤判定、利用可能な根拠の受渡し漏れ、種をまたぐ対応の途切れを確認します。全geneの機能解明や下位Tierの一律再調査を、通常運用の完了条件にはしません。

## 0. 主目標は候補保持であり、支持確定ではない

**支持確定候補を作ることではなく、真の原因候補を失わず、実験可能な約20遺伝子のWet集合へ整理すること**を主目標とします。**支持・検証状態は候補ごとの属性であり、支持確定をWet採用の必須条件にしません。**

主評価は、**実際に採用したWet集合での正解保持と実遺伝子数**です。support/review/fillの内訳、採用段階、完了・未評価、実行負担を別に報告します。reviewやfillとして正解が残っていても保持として数えます。全区間を残して保持率だけを高く見せることは目標の達成ではありません。

### 20は取りこぼし対策の目安

**20件は目安であり、最低件数・切断上限・調査完了の条件ではありません。20件を超える既存候補は数字だけで切らない。** 必要な確認を終えて12件なら12件で終了できます。件数不足だけを理由に追加探索やlower-tier全件比較を要求しません。

既に挙がった弱いreview候補・別機構・区別できないcopyは安価に広く残します。補充に本命と同じ直接証拠・claim確認・具体的な機構説明を要求しません。identity-onlyの全geneを実験候補へ写すこととは区別し、追加調査は採否・identity/copy・作用機構を変える具体的な問いに集中させます。架空ID・重複・根拠の格上げは禁止します。

元の集合を保持した追加は、その集合に既にある正解を失わせません。ただし、20件を埋めただけで未知の正解保持を保証できるわけではありません。元候補の入替や第2段の選定にはこの単調性を適用しません。

| 状況 | 扱い |
|---|---|
| 必要確認済みの8/12件 | 実数で返す。20との差分だけで追加比較を要求しない |
| 上位20件と既存の後続Tier review候補 | v3は上位根拠のpoolを使い、最初の超過Tierだけを比較する。後続は適格性を保って台帳。旧v2は後続も従来scopeへ渡す |
| relevanceで採用した21〜25件 | その採用集合を保持する。20件への切断をしない |
| T1だけで25件を超える | T1全体を保持する。後続の適格候補があれば、その比較scopeは別に作る |
| 適格poolの累積で最初のboundary Tierが生じる | v3は上位を保存し、その境界Tierだけを比較する。後続は台帳に保持。旧v2は境界と後続を同じscopeへ渡す |
| inventoryなし・適格補充候補の不足 | 原因を明示する。架空IDや重複で20件にしない |

全copyを実遺伝子として数え、単一copy確定や集合内順位は成功要件にしません。benchmarkのtruth単位と採用段階は事前固定し、未同定の実案件では保持成功は未判定です。#684のベンチは明示再開までPAUSEDです。

## 1. Supported member universe

```text
A/B verified members + optional Direct-C verified_direct_support
  -> supported_member_union -> family/copy context
  -> supported/review stateを属性として保持（unionとWet membershipは別）
     + T1と既存適格supported/reviewのgene/copyを件数非依存で引き継ぐ
  -> wet_candidates.tsvをfreeze
  -> v3: 最初に25を超えるTierだけboundary比較、後続Tierは台帳に保持
  -> finishで採用/比較待ち/理由付き保留を保存し、実Wetを再freeze
  -> orderingは同じ集合内のみ
```

`supported_member_union.tsv`は支持確認済みの部分集合です。unverified proposal、provider failure、Direct-C review/ledger、区間外、no match、mapping unresolved、generic homologyだけの行はsupported unionへ昇格しません。ただし、unionへの非昇格とWetからの除外は同義ではありません。実在・区間対応を確認できたreview候補は別入力から合流します。

A-only / B-only / C-only / A+B / A+C / B+C / A+B+Cは支持のprovenanceであり、route convergenceはscore bonusではありません。通常はA/C-firstでRoute Bは任意です。Cの探索済み空結果と未実行を区別し、旧A/B凍結世代は元の方針で再生します。

## 2. Membership rule

### #904: 第一段の本候補採用

新bundleは `b_review_membership_policy=initial_bibliography_v1` をfreezeします。
元のA/Bが判断した具体的な機能と表現型への推論に、解決した文献ID・提示されていればタイトル、凍結参照との対応を合わせます。
汎用kinase/TF/GOやIDの存在だけで候補化しません。具体的なcurated function/annotationも従来のprovenanceと対応を保って使えます。
Aはprovider identity・種・配列と区間投影、Bはfrozen memberと参照対応を保ちます。1:1 orthologyは必須でなく、card/lookup/相同性handleの条件付き対応も使います。
一つの参照機能を別の参照に対応するcopyへ暗黙適用せず、確定した候補範囲を読み戻し時にも保持します。
候補採用は既存の非support laneを通り、`supported_member_union`には昇格しません。Bの資料不足はA/Cの独立適格性を拒否しません。
metadataはcache優先・Europe PMC liteを使い、返された実IDと提示titleを照合します。title未提示は照合済みと捏造しません。
同一論文の別IDは `source_identifier_aliases`、別の論文は別source/claimとして扱います。no-hit・取得失敗・未取得・ID/title不一致を区別します。
本文の独立確認状態は第一段の候補表・traceへ追加しません。source provenanceとsupport状態は別に保存します。
通常 `prepare → review → finish` はレビュー回答なしで第一段をfreezeできます。旧prepared/sessionは元のstrict方針で再生します。
明示互換モードはordinary configの `candidate_admission_policy=content_v1` です。

Route A keeps the original hypothesis citations separate from the resolved provider record's full bibliography. Only the identity-checked declared accession (or a uniquely resolved named target) carries those citations; a changed accession stays unresolved. A title and its identifier aliases refer to one paper; other source_ids may cite separate papers without inheriting that title.

For `initial_bibliography_v1`, active Route-A function provenance uses the ID/title-checked original papers, the exact reviewed database-record handle and any explicitly authorized bound citation. Provider bibliography is retained in the verification audit and sealed original inputs, not added as function evidence or copied into active review evidence. Identity/name/species/accession/sequence/copy metadata stays intact. This separation adds no content-review or verified-support admission gate; legacy content/curated/independent-review authorization remains unchanged.

When the parent declares no stable ID and a provider lookup returns multiple candidates, an already recorded unique paper-level cross-reference rebind may carry only that paper's original citation, title and aliases to the identity-checked accession. It does not restore the parent bibliography or override a different declared stable ID. Name-only/shared-ID ambiguity, other candidates and conflicting same-gene species remain unresolved; original reports and retrieval provenance are retained.

### #903: 元のA/B判断で付ける初期Tier

`candidate-relevance-2.5.0` / `reported_function_v2` を使い、旧Tier名を新しい意味へ置換しません。

| 初期Tier | 根拠 |
|---|---|
| T1 | 対象種の同一・同等または非常に近い具体的表現型と遺伝子の既報関係、適格な区間内相同投影（#1002の標準判定。下記「初回A/Bのチェック表」節） |
| T2 | 同種・他植物で同じ/非常に近い細胞・発生現象の既報 |
| T3 | 別器官・時期等の既報表現型と具体的な共通過程 |
| T4 | 類似表現型の既報なし。報告機能・経路からの具体的仮説 |
| T5 | 機能実験は薄いが、出典付きの具体的注釈/domainからの仮説 |
| T6 / U | identityのみ / 未解決 |

`initial_tier`、`initial_tier_reason`、`reported_phenotype/process`、`phenotype_comparison`、既知の器官・時期を元の提案で残します。
Tier専用の追加LLM・検索・独立再読はしません。機能説明だけでT2へ上げません。既報関係の記載が足りない具体的候補もT4で残せ、厳密なTier確定は入場条件ではありません。
既存のrelevance inventory・boundary・Wet表・trace/reportへ引き継ぎます。旧ラベルと旧scopeの順序・再生を保ち、過去freezeを書き換えません。

### 互換性: #901 source_function_v1の出典付き機能条件

旧B bundleは `b_review_membership_policy=source_function_v1` をfreezeする。
参照機能の実際の記述、参照gene/種との対応、区間unitへの対応、今回の表現型への
具体的な条件付き橋渡しを確認したものだけをB単独のWet reviewへ渡す。
対象種の因果論文、同一表現型、機構全体の直接実証は不要。
Type3の `mechanistic_review_seed` は主要なreview候補で、verified supportへ昇格しない。

既存独立claim reviewのB 1.1は、引用span・source/hash・provider gene/species対応を照合し、
参照機能と具体的推論について別のsemantic attestationを記録する。
原著抄録/本文と、provenanceを保った具体的curated function recordを利用できる。
出典ID、論文タイトル、汎用kinase/TF/GO、モデル知識だけでは入場条件を満たさない。
機械によるspan/identity照合もreviewの意味判断も因果の証明ではない。

取得失敗・未探索・未読・対応不明は適格未確定。出典なしやreview不適合も非原因とはしない。
原claim、全gene/copy、理由、request/verification状態を `route_b_review_eligibility.tsv`
および全interval Wet inventoryへ保全する。A/Cの独立適格性に対する拒否権は持たない。
copy ambiguityだけで一つへ切らず、対応のあるgroupの複数copyを保持する。

版のない旧bundleは `legacy_grounded_review` の旧入場規則で再生する。
旧レビューpacket 1.0の読取りも保ち、過去freeze/採点を書き換えない。
新規入場規則は意図的な科学的membership変更であり、機械的parityや性能向上とは主張しない。
1〜2原著のsoft budgetは未採用。検索自由、複数資料、矛盾/copy差の追加調査を禁止しない。

### #817: 有限scopeと終了の運用

開始時に必須確認、対象unit/候補、終了stage、既存のround budgetをmanifest/receiptへ結び付ける。
追加調査には採否/次実験を変え得る具体的な問い、対象ID、入力差分、追加負担を残す。
未解決・資料不足だけで次roundを自動要求しない。既存budget消尽後はpartial/未評価・
保留理由・再開先を報告し、保存した第1段集合を変更せず終了できる。
全gene台帳、現在の適格候補、第1段保持、最終採用、比較待ち、deferを区別する。
新しい終了frameworkや科学的な予算既定値は追加しない。20は終了条件ではない。
停止済みrunはユーザーの再開指示なしに再開しない。資料不足を生物学的陰性に変えない。

### #789: 同じ参照機能判断の共有とcopy固有適用

既存B独立review promptは、同じreference gene/種/source/function/mechanism/phenotype relation/
代替仮説/不確実性についてexact questionをまとめる。superfamilyでgeneを集約しない。
共有source snapshotと凍結された全入力hashを固定し、機能判断は一度、claim/copyへの適用は個別に記録できる。
共有判断のIDはレビュー対象全体のhashに依存し、表現型・source・比較条件の変更時に旧判断を無条件で再利用しない。
copy固有の非適用や未評価はpositive判断を借りられず、区別できないplausible copyはreviewで残せる。

第一段では、出典付き参照機能と具体的な表現型への推論に基づくfamily/相同/近縁paralogの候補保持を許す。
focal copyの参照endpoint完全一致、1:1 orthology、直接実証、精密系統解析は必須にしない。
既存freezeのID/member anchorと実際の参照source identityは照合するが、間接的な機能推論を禁じない。
汎用kinase/TFという分類だけで全memberを候補化する理由にはしない。一方、familyの規模だけでも排除しない。
copy差が採否を左右する具体的な問いがあるときだけ詳細検討し、発現/系統/位置等は利用可能な補助資料として使う。
共有質問ID/理由をreview候補へ渡し、前段判断を後段の比較資料として利用できる。

auditは実際に返された共有機能判断とcopy適用判断の数を別記録する。
保存bytes、実原文読込回数、実tokens/費用は別で、未取得値はnull。
この共有回答形式を実モデルでの読解時間/科学的保持改善の実証としない。

### 人間の候補選定を再現するための役割分担

| 担当 | 責務 |
|---|---|
| 機械 | 配列比較、全候補の整理、ID・種・アセンブリ・区間・引用箇所・hashの照合、重複除外、件数確認、判断記録の保存 |
| 研究者役LLM / agent | 仮説生成、文献の意味理解、コピーや別機構の比較、本命が外れた場合の見落とし検討、理由付き採用/保留案 |
| 既存の検証・handoff | 対応の正確さと出力契約を確認し、明示された選定結果を記録する。照合成功を生物学的な証明と呼ばない |

研究者役LLMはWebで一次資料を調べ、他種の報告された機能・変異表現型と、対象種の相同遺伝子・family・機構を橋渡しすることを主軸にします。文献の観察事実、配列/familyの対応、今回の表現型への推論を分けて記録します。対象種での機能実証や参照との完全な表現型一致は、review候補化の必須条件ではありません。同種の報告は候補の根拠を補強する追加情報として扱い、その既知遺伝子検索を主軸にしません。特定種・形質・ベンチ例専用の選定規則を作りません。未実証の橋渡しは`mechanistic_review_seed`等の明示reviewで扱い、supportとは別にWetへ渡します。

この知見の転用と科学的な意味判断がLLMへ期待する仕事です。モデル植物でも研究の乏しいgeneは、利用可能な資料と不足理由を記録した上で未解決のまま残せます。全区間geneの台帳と依頼済み初期判断の範囲を保ち、調査済みでも根拠不足の状態と、未読・未実行を分けます。全geneへ機能や表現型を付与することは求めません。

境界比較では元研究役が記録済みの根拠を使って実験優先度を判断し、Tierや表現型を別agentで再審査しません。理由付き延期は生物学的陰性と区別し、支持確認未完やcopy ambiguityだけでは延期せず、区別できない適格copyを保持します。v3の封印済み境界だけは、元担当の理由付き`adopt/defer/not_evaluated`をhash・roster・executor identity/modelの機械照合で適用します。`original_assessment_v1`としてmanifest・履歴・保留台帳に記録し、reviewer欄は空、独立reviewは`not_required_original_assessment`です。未評価・欠落は保持して履歴未完とし、全未評価を完了にしません。旧v2/legacyとgeneric consumerの独立合意契約は維持します。

v3では#732のreview込みconsumerを再利用し、coordinatorが指定した既存研究担当の単独判断を使います。identity/modelの照合は申告整合とprovenanceであり、上流A/B担当との同一人物認証ではありません（`executor_identity_basis=declared_identity_model_consistency`）。同じ研究threadの再利用は実行記録で残し、混合A/Bのactorを機械が新たに選びません。`agentic_boundary.py merge --role proposal`へ実executor ID/modelを渡し、`claim_review_session.py finish-wet-review`は同じexecutor ID/modelとproposalを受け取ります。v3ではreviewer引数・独立review artifactを指定しません。

既存の `scripts/agentic_boundary.py prepare --comparison-context` は、記録されたphenotypeとboundary/higher全IDに対応する追加調査材料を任意のhash付きJSONから各shardへ渡します。request policyも各shardへ保持し、mergeは追加材料のhashと回答の `comparison_context_sha256` を照合します。追加材料は比較の根拠であり、support・Tier・higher membershipを変更する入力ではありません。数値順位・固定family枠・Top20を追加しません。

同じprepareはcanonical `shard.json`を保持し、`comparison_index.json`、`material_catalog.json`、shard別の`reading_view.json`を作ります。前段の判断・資料はcomparison contextへ由来付きで供給し、全候補の短い索引から機構の食い違い、参照/copy対応、採否を変える不足を確認します。完全一致する資料の部分木だけを共有し、gene/copyと出典の違いは別bindingへ残します。索引で短縮・省略した値にはcanonical ownerのfile/SHAや原値のhashを残し、必要箇所だけ原資料を開きます。候補ごとの判断や独立reviewは共有資料への参照と区別します。

`--max-material-bytes`（既定262144）は配送の目安です。超過は`chunk_required`等で報告し、資料やIDを切断しません。共通索引自体が超える場合、候補を全て単独shardへ分解して負担を増やさず、親が共通文脈を分けて読みます。catalog/index/viewの物理bytesと、各shardへ展開した読取bytesは別計測です。共有の少ない資料では索引分が増え得るため、token・費用の削減を保証しません。mergeは派生artifactのhashとcanonicalへのprojectionを検査し、旧manifestも互換受理します。

Route Bの凍結promptが供給source IDだけを受け付けるのは引用の取り違えを防ぐ境界です。調査agentのWeb検索を禁止する意味ではありません。新しく見つけた文献は取得範囲・reference対応とともに新しいbundle/promptへ保存してから解釈・importし、既存freezeを上書きしません。

補充では、検証未完の提案、曖昧なcopy、別機構、注釈の薄い候補を見直します。これは検討観点であり、固定比率やfamilyごとの枠、論文数の閾値、新スコアを設ける指示ではありません。20件への不足も到達も調査開始・停止の理由にせず、既存の本命を補充候補で追い出しません。

IDや引用の誤りはその対応・証拠を修正する理由です。一つの引用の失敗や未評価だけで、別の根拠を持つ遺伝子まで非原因としません。捏造引用、誤った種/ID、区間外member、破損freezeを通す変更ではありません。

この役割分担は#732の設計方針です。通常target handoffの適格候補の機械的引継ぎをLLMによる救済比較の完了やcausal rankingと読み替えません。strict、`supported_only`、明示的refinementでは旧tier/ID順fallbackを互換維持し、新しい選定器・Route・providerは追加しません。

## 3. Relevance Tier target と現行runtime vocabulary

### 3.1 target: canonical relevance Tier

**新規A/Bの正本はT1〜T6/U（`reported_function_v2` / `paper_checklist_v1`）**。今回の表現型との関連度であり、verified support・mapping・copyは別の属性とする。

| relevance_tier（新規A/B） | 根拠 | Wetでの位置づけ |
|---|---|---|
| `T1_reported_focal_trait` | 同種の同一・同等または非常に近い具体的表現型の既報と、妥当な区間内相同対応 | 第一段の高優先度本候補。適格な全copyを保持し件数だけで外さない |
| `T2_reported_similar_phenomenon` | 同種・他植物で同じ/非常に近い細胞・発生現象のgene–phenotype既報 | 強い比較候補 |
| `T3_reported_shared_process` | 別の器官・時期等のgene–phenotype既報と具体的な共通過程 | 過程を橋渡しする候補 |
| `T4_reported_function_hypothesis` | 類似表現型の報告がなくても、報告機能・経路からの具体的な推論 | 未知原因の機能仮説 |
| `T5_annotation_hypothesis` | 出典付きの具体的な注釈・domain等による仮説 | 弱い機能候補 |
| `T6_identity_only` | BLAST/相同性等identityのみで具体的trait relevanceは未確立 | 全gene台帳・必要時rescue |
| `U_unresolved` | identity/機能が未解決 | 未評価・未解決台帳 |

**T1〜T3には遺伝子と表現型の既報関係が必要**。機能の一般論だけでT2/T3に上げず、T4/T5も具体的な推論と別の適格条件を満たせば初期Wet候補になり得る。T6/Uは全gene自動Wet採用の対象ではないが、生物学的陰性として捨てない。既存の適格review候補・曖昧copyは、未支持や件数だけで落とさない。

T2/T3は、非常に近い**観察現象**の既報か、異なる観察現象をつなぐ**具体的な過程**かで区別する。T4は報告された**機能**、T5は出典付き**注釈**、T6は**identityのみ**。論文の意味の初回判断はLLM、9項目チェックからのTier計算・ID照合は機械が行う。Tierの高さを支持確定と混同しない。

報告機能の出典は論文だけでなく同定可能なdatabase accessionでもよい。元研究者の `evidence_types=["function"]`・`reported_function`・具体的な形質への橋渡しを保持し、既存のsource/参照対応検査を行う。注釈やdomainだけを機能へ読み替えず、出典IDだけでT4にしない。これは既存T4判断の輸送修正であり、T1〜T3の表現型既報要件やverified supportの基準は変えない。

A/B共通の `paper_reports_json` は従来の配列に加え、`{"bibliography":{"p1":{"source_id":"PMID:...","source_title":"原題"}},"observations":[{"paper_ref":"p1", ...}]}` という省略表現を受理する。LLMは論文ID・タイトルを記録し、各観察には報告遺伝子名・ID・種、形質/実験context、理由・9回答を残す。機械はID/titleだけを展開し、同じレビュー論文の異なるgene/形質を統合しない。不明な参照・矛盾するID/titleは拒否する。既存のsource照合・Tier算出・台帳には展開した同一記録を渡し、別routeの判断は共有しない。

旧凍結世代の`T1_same_species_known / T2_trait_or_mechanism_linked / T3_functionally_plausible / T4_identity_only / U_unresolved`は旧relevance語彙。互換再生のため保持するが、新Tierへ単純に読み替えない。旧凍結policyの候補順・予算処理も新通常運用へ持ち込まない。

### 3.2 Tierとは別に保持する軸

| axis | 既存の流用先 |
|---|---|
| support / verification state | verified / review / context / not_evaluated、A/B/C各routeの現行state |
| mapping / identity quality | `family_candidate_bridge.mapping_tier` |
| family resolution | superfamily -> family -> subfamily/orthogroup -> specific homolog -> focal copy set |
| copy state | `family-to-copy` の preferred/co-preferred/paralog/tandem/unresolved states |
| route provenance | A / B / C |

mappingのunique/ambiguousはrelevance Tierを上下させない。copyが曖昧でも、同じT2 familyに3 plausible genesがあれば3 geneともT2として残し得る。family resolutionを細かくすることと、copyを1 geneへ強制collapseすることは別である。

### 3.3 既存artifactからの流用

- Route B: `direct_functional_bridge` / `mutant_phenotype_bridge` / `mechanistic_process_bridge` / `mechanistic_review_seed` は初回仮説の材料。claim_typeだけからT2/T3を付けず、論文別gene–phenotype報告と観察比較をチェックして新Tierを判断
- Route A: mechanism statement、model-species evidence/mutant phenotype、reference bait、source、interval projectionを同じTierへnormalize
- Direct C: `verified_direct_support / reference_phenotype_review / ledger_only` をTierそのものとして使わず、relevance Tier + support stateへ分解。同種same-trait knownはT1、cross-species direct phenotypeは主にT2
- Route B BLAST: #751のidentity metadataのみなら新T6材料。報告機能があればT4、具体的な注釈ならT5を別に判断し、BLAST hitだけをtrait evidenceへ昇格しない

### 3.4 Wet 20–25との接続

```text
T1 -> T2 -> T3 -> T4 -> T5 -> T6 -> U
```

- T1はknown rescueなので原則全保持
- 通常v3は報告根拠の上位があれば上位とprotected Cだけをactive対象にし、上位0なら従来の適格scopeを使う。20未満でも件数による追加調査なしで終了できる
- T6/Uを含む全gene台帳は保存する。T4/T5の具体的仮説と、identity-onlyのT6/Uを区別し、追加探索は具体的な問いに基づいて別generationへ固定する
- poolをTier単位で加えて25件を初めて超える場合、そのTierだけを#750の既存boundary consumerへ渡す。後続Tierと下位の元適格候補は台帳へ保持し追加比較しない
- higher Tierはfreezeし、再選定しない
- boundary Tierでは原典・mechanism・functionに加え、`superfamily -> family -> subfamily/orthogroup` とidentity解像度を上げる
- copyを安全に区別できなければ複数を保持
- 25は科学的cutoffではないため、安全に縮小不能なら超過を許容

以下のprovisional opt-inは旧whole-tier policyの互換経路です。新通常policyで12件を返すことは、それだけで暫定・未完にはなりません。Bの依頼済み回答欠落・不正claim・必要review未完などの実coverageは引き続きpartial/degradedとして報告します。

通常target handoffのinput coverageは、厳密なverified support件数と分けて報告します。Aの`initial_bibliography_v1`では、hash固定した仮説の書誌照合済み状態、凍結済み配列panelと適格review候補を照合し、supportが0という理由だけで劣化にしません。これは調査・import全体の完了宣言ではありません。producerに依頼範囲と完了statusが結び付いていなければ`execution_state=unknown`を残し、書誌・panelの各状態件数も保持します。結び付いた不完全importや利用可能な候補根拠の欠落は引き続き明示します。

Cは、既存consumerで検証した原記録・投影系譜から、source確認済みの`no_interval_match`、適格な区間内support/review、source未確認、投影未解決を区別します。確認済みの区間一致なしは供給された証拠の評価済み状態で、生物学的陰性や探索全体の完了を意味しません。source pendingや未投影・失敗はpartialとして残し、状態別件数を保全します。必要な系譜がなければunknownです。候補membership・Tier・support・原freezeは変更せず、旧policy、strict、`supported_only`の判定は従来どおりです。

明示的な`--allow-provisional-wet`を指定した場合だけ、全区間遺伝子を網羅するhash固定済みscreen ledgerを検証します。T1–T3と保護対象Cの全候補が25以下で、最初の超過境界が比較根拠のないT4/U identity-only集合なら、それらの具体的候補を暫定採用します。T4/Uは陰性扱いせず、全区間の`interval_unresolved.tsv`へ`not_evaluated_identity_only`と観察状態を残します。結果の`handoff_state`は`provisional_adopted_wet`です。screen欠落、informative境界、具体的候補が0件または25件超では従来のwhole-tier境界reviewに戻ります。この例外は真因保持を保証せず、最終Wet件数と未評価数をともに報告します。

同じ例外はboundary review後のlower-tier continuationにも明示的に指定できます。元の全区間inventoryとscreen ledgerを固定したまま、独立reviewで合意した保留台帳だけから現在の候補集合との差を説明します。保護対象・凍結済み上位候補を除去せず、未完のreviewや不一致がある場合は暫定採用しません。採用・保留・identity-only未評価の全遺伝子を区別して記録し、以前のfreeze、relevance Tier、supportを変更しません。API/`finish-wet-review`は`allow_provisional_wet`/`--allow-provisional-wet`とscreen ledgerのpath・SHAを受け取り、以後はhash固定した元資料と累積保留台帳の系譜を継承します。

### 3.5 現行runtime vocabulary

runtimeは次の旧support/membership labelを互換出力として保持する。

| current label | 現在の意味 |
|---|---|
| `tier1_verified_support` | A/B verificationまたはDirect-C member gateを通ったsupport |
| `unverified_route_review` | 明示route review |
| `tier2_reference_phenotype` | reference phenotype欄等を使うfill |
| `tier3_family_or_domain` | family/domain等のfill |
| `tier4_dark` | 現行fill対象外のdark |

これは**移行前のmembership/support vocabulary**であり、旧T1–T4/Uにも新T1–T6/Uにも該当するrelevance Tierと同一視しない。特に `tier1_verified_support` は新`T1_reported_focal_trait`や旧`T1_same_species_known`を意味しない。

通常のrelevance selectionでは同一Tier内のgene ID順で切らず、適格rosterとpolicy版をimmutable scopeとして保存します。`wet_candidates.tsv`は採用と未完の比較待ちを含み得るため、reportのadopted/pendingを区別します。`wet_candidate_trace.tsv`は全geneのunion/review在籍、bound claim数、relevance、first-stage dispositionを追跡します。traceの`upstream_llm_assessment=not_evaluated_by_handoff`は、handoffだけからLLMの認知性能を診断しないための固定値です。

通常target handoffのfillは`membership_basis=interval_review_coverage`、`rank_contribution=none`です。支持フラグや架空の引用は付けません。明示reviewは`membership_basis=route_review_candidate`として、route・claim・hypothesis・verification・reason・sourceのprovenanceを保持します。review/fillをsupported-unionの支持確認済み集合へ算入しません。strict、`supported_only`、明示的refinementでは旧tier2/tier3 fillを使い、`tier4_dark`はfillしません。

A/Bのreview TSVとproducer manifestは`--review-candidates` / `--review-candidates-manifest`へ同じ順序で渡します。`--wet-fill-mode supported_only`は診断・比較用であり、固定設備容量が明示された場合だけ`--budget-policy strict`を使います。

legacy fillで材料不足・inventoryなしなら`unfilled_budget_slots`が残ります。新通常policyの20未満は実数の記録であり、調査不足の判定ではありません。全区間に20遺伝子あるだけでは常に20件になるとは限りません。`tier4_dark` exclusionはlegacy pathの規則であり、具体的な問いに基づくT4/U rescueを禁止しません。


## 4. Direct-C member statesとの対応

### 第2段とDirect Cの残作業

v3の第2段は最初の超過Tierだけを、元研究役の既存根拠で検討する段階です。20/25への数合わせ、同Tierの順位付け、無理な縮小、追加の原文調査や別agentの妥当性再審査を目的にしません。T1/T2の根拠を十分と判断した理由を明示できるなら、境界のT3を全て実験優先度上の保留とすることも禁止しません。確信度を機械が自動推定したり、Tierを意味の違うラベルへ変更したりはしません。後続台帳と元根拠は残し、保留を非原因の証明にしません。

旧v2/legacyの詳細な残存候補比較・独立selection reviewは元policyで再生します。境界のない25件以下の集合からT3を除く入力は自動追加しません。新しい元判断・根拠への訂正は既存の新generation経路で記録します。

支持状態、今回の採用/保留、生物学的な除外は別の判断です。保留は非原因の証明ではなく、support未確認だけを自動不採用にしません。第2段後の件数不足で全lower-tierを調査せず、明示的な保留を通常fillで黙って復活させません。再採用には理由を記録し、原候補・原根拠は保持します。既存consumerはv3の封印済み境界に限って元担当の判断を消費し、旧policy/genericでは独立selection reviewを要求します。どちらもhash/roster検査を通した新Wetと保留台帳を出力し、元freezeを変更しません。
proposalと独立reviewが不一致の場合は候補単位で扱い、合意した`defer`だけを保留台帳へ移します。不一致候補は`not_evaluated`として保持し、両判定・理由・解決状態を専用decision artifactへ記録します。`defer`は実験優先度の判断であり、生物学的negativeやclaim未検証を意味しません。
独立reviewはproposalの証拠・identity・表現型解釈と候補間の実験優先度を監査するためのもので、同じfavorite 20の再選出や非原因の証明を要求しません。支持不足だけでなく、器官・時期・作用機構・copy coverageなど比較可能な理由を伴う場合に限り、提案への合意deferを記録します。

第1段を上書きせず、要求・実行・完了・採用集合を区別します。`second_stage_state`の省略・要求・完了・round budget exhaustedを混同しません。新policyのfinishは封印した比較範囲の採用・保留・not_evaluatedを反映し、件数不足によるlower-tier全件continuationを作りません。旧scopeのcontinuationは旧policyとして再生します。未完でも比較待ちの有効候補を失わず、要求書や保留台帳だけに残る正解を採用集合のhitにしません。

新policyの未完reviewは`review_incomplete_retry_original_scope`として記録し、元completionと比較scopeのhash付きresume先を保存します。再回答は元のcheckpointから別出力へ`finish-wet-review`を再実行します。未完出力を新たなlower-tier探索の入口にせず、旧回答・freezeを保持したまま、元scopeのpendingを解決します。

最終採用WetでのType3 true causal gene全件の保持を評価し、第1段・expanded set保持とは分けます。採用集合、unique gene数、採用stageを固定し、unique gene数とmappingに対するcost/timeも別々に報告します。Type1/2の回収は先行診断であり、最終Wet保持やcausal proofへ自動昇格しません。複数の原因geneが同一区間にある場合は全件をtruthとして扱います。

Direct Cは`route-c-direct-members`から、`verified_direct_support`をsupported unionへ、適格な`reference_phenotype_review`を`unverified_route_review`へ、`ledger_only`を参照台帳へ渡します。適格reviewは空き枠やgene ID順に依存せず全件保持し、同一geneの複数sourceはprovenanceへ集約します。現行support条件はS0/S1、D0/D1、source確認、区間内、exact/assembly mapping等に限定されます。cross-speciesやBLAST-onlyを無条件にtier1へ昇格しません。

**受理済みDirect-C救出候補は最終Wetへ直接合流し、20/25件を超えても保持する。** relevance Tierが低い場合も、上位Tierで予算が満たされたことを理由に外さない。boundary reviewの比較対象から保護し、finish・lower-tier continuationでも保護を維持する。同じgeneがA/Bにもある場合は一度だけ数えて全route根拠を残す。これはWet所属の保護であり、T1やverified supportへの格上げではない。source・identity・projectionの既存受入条件を通らない行やledger-onlyを、この保護へ混入させない。

#661のsource gene/species/trait/experiment binding、#663のprojection evidence順序不変性、#662の検索系譜を確認し、未完境界を隠しません。記名した`curated_human`または`independent_agent_review`が同じ機械検査を通した完全bindingと認可projectionは、追加人間確認なしでsupportへ進みます。未完・矛盾・未認可はreview-onlyです。C reviewはsource identity/species/citationとprojection ledger/hash/interval rosterを検査した行だけを保護します。cross-species、orthology-only、variant-only、未確認/矛盾sourceやmapping conflictはsupportへ昇格せず、条件を満たさない行はledgerまたは通常fillの対象です。generic fillやA/Bだけで残ったgeneをCの増分救済と呼びません。

#722のreview保護にはfreeze済みのsource bindingが必要です。出典実在の確認状態、種/accession、stable geneまたはlocusのidentity、引用・原文hash・reviewerを確認し、同じsource identityが区間内projection ledgerへ使われたことを照合します。trait perturbationの支持確定は必須にしません。variant-only、source矛盾、mapping conflict、binding/ledger欠損は保護対象になりません。legacy自己申告reviewは従来の容量制限付きtier2 fillに留めます。改変された入力はエラーにし、通常fillへ黙って降格させません。これは既存source bindingの消費であり、#662の原文取得系譜を完成させる変更ではありません。

## 論文単位のTierと遺伝子への集約（strongest_paper_v1）

Tierは、対象の問いに対する各論文の報告内容に付ける。仮説全体の確信度や論文数に付けない。
新しいA/B提案は `paper_reports_json` に論文ID、報告遺伝子・種、表現型・過程、比較関係、
Tierと理由、既知の題名を論文ごとに記録する。同じ論文が複数遺伝子を扱う場合も報告対象を区別する。
PMID/PMCID/DOIの別名は同一論文のものだけを記録し、別論文のTier・題名をコピーしない。

遺伝子の代表Tierは、既存の対応付けを満たす報告の最強Tierとする。
**T1が1本と、T1が1本＋T3が多数は、同じT1である。** 本数、経路間の重複、識別子の別名を
加点・投票・平均・同Tier内の優先理由にしない。代表根拠は `tier_source_id` と
`tier_evidence_json` に記録し、全報告は `relevance_evidence_json` に残す。

新しい論文記録は初回A/B調査のチェック表からTierを決める（paper_checklist_v1）。旧観察形式・旧申告形式は元のポリシーで再生する。Tier再判定の別エージェントは追加しない。
複数論文は根拠をたどるために記録し、T1×2・T2×3やT1原著1本＋T1レビュー3本も保持する。
原著・レビューのどちらも使え、原著であることを必須にせず、論文種別だけでTierを上げ下げしない。
T1を得たあとに無理にT2を探すことも、1本だけに制限することも、必ず追加の裏付けを探すこともしない。
T1の既存の保持ルールに、論文数・種別・独立性の条件を足さない。T3以下の仮説にも利用できる論文根拠を残す。
論文がない場合は明示し、本数合わせの捏造や、既存の適格候補の機械的な除外につなげない。

遺伝子Tierとは別に、`paper_tier_counts_json` にTierごとの論文数、`paper_evidence_summary_json` に
各論文の短い説明を記録する。同じ論文のPMID/DOI別名、同じ論文のA/B重複、同じ論文の複数claimは
重複計上しない。これは申告された論文の内訳であり、独立した実験・再現確認の回数ではない。
レビューが原著を紹介する場合、論文は2本でも報告の根は同じことがある。判明している論文種別、
元の研究のID、裏付けの関係を記録するが、未知の独立性を調べる追加工程は設けない。
論文数は検証可能性のための記録であり、スコア、優先順位の加点、Tier昇格条件にしない。

第2段には代表根拠に加えて、下位Tierも含む各論文の要約を渡す。長い元記録は監査台帳へ残して
必要時に読む。固定の1本上限・最低本数・必須の追加レビューは設けない。追加検索が既に確認済みの
同じ報告を繰り返すだけで、出典・対応関係・機構・既知の矛盾の判断を改善しないときに止める。
Tierの総当たりや、必ずT1に到達するまでの探索は要求しない。既知の反証とコピーの曖昧さは残す。

論文報告がない場合は明示的に `[]` を記録する。注釈だけの具体的仮説はT5、同定だけはT6、
未解決はUであり、論文のT1–T4を借用しない。Tierは元の調査判断であり、独立した内容検証や
verified supportを意味しない。Tier確認だけの追加レビュー工程は設けない。
旧形式でこのフィールドを持たない記録は、元のポリシーで再生し、論文単位へ推定変換しない。

## 5. Ordering boundary

`integrated-narrowing`はfreeze済みWet集合内の順序だけを扱います。集合外geneはpriority番号を消費せず、family/copy contextも同じ集合へ射影します。ordering前後でgene setが変わればfail closedします。採用集合の変更は明示的な別段階で行い、並べ替えに隠しません。

## 6. Artifact classes

| artifact | producer class / role |
|---|---|
| `wet_candidates.tsv` | `canonical_wet_membership`：実際に採用する集合 |
| `supported_member_union.tsv` | 支持確認済みの部分集合。最終Wetの代用ではない |
| `route_c_direct_members.tsv` | 必須Cのsupport / review / ledger adapter |
| `stage_a/` | baseline review budget / screening / diagnostic |
| `stage_a_narrowed` | `ordering_only`：freeze済み集合内の順序 |
| budget-engine Top20 | `diagnostic_portfolio`：本番membershipではない |

「20行のファイルがある」と「正解を含む実験集合を作れた」は同義ではありません。未使用枠を使う方針と、その補充法の実性能を分けて検証します。


## 初回A/Bのチェック表と軽いID照合（paper_checklist_v1）

A/Bには共通して、論文名・PMID/PMCID/DOI、報告遺伝子名・遺伝子ID・植物種、
チェック表と短い根拠を記録させる。`reported_gene` と `reported_gene_id` を分け、
IDのDB名を保存する。分からないIDは空欄・不明とし、穴埋めだけの追加探索を要求しない。
チェック表の意味判断と「その論文がその遺伝子・現象を報告しているか」は初回A/BのLLMに任せ、
既存文献確認以上の検証、本文再読、分類専用エージェントを追加しない。

`paper_checklist` の回答は `yes` / `no` / `unknown`。欠損は `unknown` として保持する。
機械は下の表を上から評価し、最初に成立するTierを使う。自由記入のTierよりチェック表を優先する。

| 条件 | Tier |
|---|---|
| 同種の表現型報告 + 目的と同等または非常に近い具体的な細胞・発生現象 + 妥当な区間内相同投影 | T1 |
| 表現型報告 + 目的と同等、または似た細胞・発生現象 | T2 |
| 別の表現型報告 + 具体的な共通過程 | T3 |
| 報告された機能が仮説を支える | T4 |
| 具体的なアノテーションが仮説を支える | T5 |
| 同定情報のみ | T6 |
| 解釈できる情報が未解決 | U |

表の質問は `phenotype_reported`、`equivalent_target_phenotype`、`focal_species`、
`same_gene_as_focal`、`similar_cellular_or_developmental_phenotype`、`shared_specific_process`、
`function_supports_hypothesis`、`annotation_supports_hypothesis`、`identity_resolved`。
複数論文があれば適用可能な最も強いTierを使い、論文本数・原著/レビューでは加点しない。

#1002を通常の標準判定とする。T1は同種の近い既報表現型に基づく高優先度Wet候補であり、
同一座位・原因変異・各copyの因果性やverified supportの確定ではない。
`same_gene_as_focal=unknown/no`、review、copy ambiguityを残したままT1にできる。
初回の論文側T1は区間投影を条件とし、Wetでは既存の適格な配列投影または凍結gene orthologyと
実member対応を消費する。全長100%、同一annotation、一意1:1、相互一位は要求しない。
一般domain、短い局所一致だけ、理由のないproxy、種・対応の明白な矛盾や区間外hitは新T1にしない。
既存の投影閾値、実参照identity・配列・source/member/hash検証とsupport判定は変更しない。
20/25件はT1の切捨て上限ではない。適格な全copyを保持し、件数と根拠差を報告する。
ユーザー指定により並列のTier方針版は増やさない。既存成果物は保存し、現行mainでの再生は新基準を適用する。

`equivalent_target_phenotype` は入力された表現型の粒度で判断する。入力で未指定の変化方向・
発生段階・alleleまで一致を要求せず、未提示の条件を補って作らない。単なる器官名の一致や一般機能だけで
yesにせず、報告された観察が入力の表現型にどう対応するかを元研究者が説明する。対象caseの原因alleleの
証明とは別である。benchmarkのType番号はTierを決めず、Type2等でも実際の根拠が条件を満たせばT1になる。

### 報告遺伝子・代理参照・対象座位を分ける

`reference_gene_id` は既存の参照/配列のIDを選び、`gene_correspondence` は報告遺伝子から
その参照への関係（same_gene/homolog/different_gene/unknown）を表す。
`same_gene_as_focal` は報告遺伝子と対象座位が同じかという別の質問である。
旧新版の同じ座位だと通常のA/B調査で判断できれば同一としてよいが、広いQTL内にあるだけでは同一にしない。

**MtA → AtB → PsC** は、報告geneの配列を取得できない場合にも使える通常の橋渡しである。
ここでMtA、AtB、PsCは模式名であり、実在する特定遺伝子を指さない。

| 役割 | 模式例 | 保持する内容 |
|---|---|---|
| 原論文で研究された報告遺伝子 | MedicagoのMtA | 元のgene・種・観察された表現型・論文・確認箇所 |
| 配列を使う代理参照 | ArabidopsisのAtB | 実際の参照ID・種・配列と、MtAとの対応理由 |
| 対象区間の候補 | エンドウのPsC | AtBからの投影、凍結member、copy ambiguity |

MtAとAtBの相同対応、MtAの表現型が今回の現象にどうつながるかは、元の調査LLMが資料と生物学的知識から理由付きで判断する。
出典で報告された対応か、今回の推論かを区別し、通常の仮説形成として採用してよい。
MtAの配列取得、AtB自身の同じ表現型の論文、別の生物学reviewerを、この橋渡しの必須条件にはしない。
機械は実際に使うAtBのidentity・配列とAtB→PsCの区間投影・frozen member/hashを既存の仕組みで検査する。

根拠はMtAの原論文のまま保持する。そこで報告された表現型が目的と同等、または似た細胞・発生現象ならT2として使える。
具体的な共通過程による橋渡しならT3など、元の観察を既存決定表へ渡す。
MtAの報告をAtBやPsCで実験された事実へ書き換えず、相同対応だけでT1やverified supportへ昇格させない。
一方、原報が対象種のPsAで同一・実質同等または非常に近い具体的表現型を報告し、
PsA→AtBの相同対応に理由があり、AtB→区間内PsC/PsDが適格に投影される場合は全copyをT1にする。
原報はPsAの論文として残し、AtBの表現型論文・PsA自身の配列・同一座位の厳密な確定を要求しない。
最終paper summary/Tier evidenceの `homolog_projection_json` と `proxy_reference_species` は、
原報gene/speciesと実proxy/対象member、投影basis、全copy群、ambiguity、既存配列根拠を併記する。
`initial_tier_reason` に代理を含む対応と根拠を短く残す。対応不明や複数copyは明示し、勝手に一対一へ縮めない。
DBの明示的な同一レコード対応（例：GenBank/EMBLとUniProt）があれば元論文と参照anchorを引き継ぎ、
元の報告遺伝子IDやチェック表は書き換えない。関係のない別copyに論文を配らない。

### freeze後の同一座位評価を対象geneのTierへ渡す

この任意経路は、同種の報告geneと対象座位の厳密な同一性を後段で追加評価する場合に使う。#1002の第一段T1の入口条件にはしない。上記の異種proxyによる通常の仮説形成の必須工程にもならない。

区間blindの原調査では `same_gene_as_focal=unknown` でもよい。freeze後に親または独立agentが
原著gene、公開登録CDS、実験種、対象modelと全競合を評価できた場合、既存
`wet-candidate-handoff --identity-assessment` または `scripts/claim_review_session.py handoff --identity-assessment`
で `post-freeze-identity-assessment-1.0.0` の明示reviewを新世代へ渡せる。
standaloneのWet CLIでは `--identity-target-proteins` に実対象FASTAを指定し、reviewに拘束した
`candidate_proteins` SHAと一致させる。session handoffは既存 `--proteins` を使う。
元の論文回答・freezeは保存し、既存の同じmemberに結び付いた同じ原論文だけから
gene-boundの派生チェック表を作る。親の区間情報はblind Aへ返さない。

入力はcase/atlas、元paper TSV、配列bridge、実review、原著sourceをpath/SHAで結び、
原著・公式CDS・reviewの正確な引用位置、元proposal/source/reported gene、両taxid、
query SHA、対象geneを検査する。原著geneと配列requestの表示名が異なる場合も
両endpointを明記したreviewが必要で、名前だけのalias推測をしない。
初版の機械防御条件は、同種の実行済みwhole-proteome照合で全長完全一致が1 protein modelに限られ、
部分modelや競合する完全一致modelがないこと。加えて原著→配列、実験種、座位、model、全競合の
明示reviewを全て要する。完全一致だけで自動T1にはしない。
実対象FASTA全件とatlas roster・参照proteomeのgene対応と配列も照合する。
全参照の実配列から完全一致集合を再検査し、alignment台帳の欠落で競合を隠さない。

解決できた `same_gene_as_focal` と原論文の他のチェックを既存決定表へ渡す。
対象形質が同等という元の根拠がなければT1にはならない。未実行、他種、部分対応、曖昧copy、
unknownは根拠不足の記録を保ち、生物学的陰性にしない。元の全paperと派生根拠を最終成果物へ残す。
この評価はrelevanceだけに作用し、候補追加、verified supportやDirect Cの保護資格の昇格、
人間curationの代行は行わない。Direct Cの厳格なsource/perturbation/support条件は維持する。

### 軽い遺伝子ID照合

`reported_gene` / `reported_gene_id` / `reported_species`は元の調査LLMが原資料から記録した値を保持する。
機械がDBから取得した名前・IDで原回答の空欄を埋めたり、原値を上書きしたりしてはならない。
IDから名前も自動入力してその組を照合すると、同じDBレコードの自己確認になり、
LLMが論文・種・遺伝子・表現型を取り違えたことへのハルシネーション対策にならない。
DB取得値は別の監査記録へ保存して原回答と比較する。確認できないIDは確認理由付きunknownを残す。
元担当が既に明示した値の機械転記・再利用、監査用metadataや台帳ID・hashの生成は可能である。
役割と理由は[LLM boundary](LLM-BOUNDARY.ja.md#遺伝子名idを機械で埋めない理由)を参照する。

既存の `initial_bibliography_v1` 文献確認の中で、報告IDのDBメタデータを一意なIDごとに取得・再利用する。
UniProt、NCBI Gene、GenBank/EMBL、RefSeqの直接レコードと、TAIR/Ensembl Plants/Gramene/Phytozomeの
明示的なUniProt相互参照を使う。相互参照がないことは、そのDBのIDが架空だという証明ではない。

照合対象はIDと遺伝子名・登録済み別名とし、種名表記は保存するだけで成否条件にしない。
Aですでに取得した参照メタデータも再利用し、同じIDを再取得しない。
結果は `matched` / `mismatch` / `id_resolved` / `unconfirmed`。大小文字・Unicode・表示上の括弧を正規化し、
DBの正式名・登録済み別名を使う。名前の長さや編集距離で拒否しない。IDが実在して名前が確認できない場合は
`id_resolved`（IDのみ確認済み）。取得失敗、未対応形式、オフライン未取得は `unconfirmed`。明示的なID未発見/返却ID違いは `mismatch` と記録する。
この照合は監査属性であり候補除外やsupport昇格のgateにしない。

`source_verify/gene_identity_checks.json` と出典確認表に結果を残し、manifestにhashと取得件数を記録する。
キャッシュの `gene_records` mapを再入力すれば再利用できる。論文本文と遺伝子の関係は検証しない。
既存の配列・区間投影・凍結hashの整合確認は維持する。

旧 `reported_paper_judgment_v1` / `observation_rules_v1` / declared-Tier記録は互換経路で読める。
凍結済み入力・予測は変更せず、チェック表やメタデータ確認を独立した生物学的正しさの証明とは扱わない。
