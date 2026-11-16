# Codex/OpenCode に通常の QTL 案件を任せる

Direct Cの検索済み空結果は、空membersと既存producer manifestに加え、hash付きstatusを渡す。
native研究のstatusには`research_receipt_path` / `research_receipt_sha256`を記録する。
receiptは`contract=route-c-native-research-receipt-1.0.0`、`case_id`、`worker_id`、`model`、
`interval_blind=true`、`state=completed`、`finding_count=0`、`no_findings_rationale`、
`queries=[{query,status,result_path,result_sha256}]`、`read_sources=[{path,sha256}]`を持つ。
queryは実際の成功検索（`retrieved`または`no_hit`）と保存した実tool payloadに結び付ける。
取得ありなら読んだsourceを保存する。検索hitがあっても原研究者が非該当と判断してよい。
失敗・未実行queryや空artifactは完了の根拠にしない。hashは来歴の検査であり、意味判断や生物学的不存在の証明ではない。
専用searchは`search_manifest_path` / `search_manifest_sha256`から既存のquery/status/provenance generationを検査する。
handoffはnative原receipt/statusとrebound版、実query結果・読了資料を保存し、原hashを残す。

通常A/Bの初期論文記録では、掲載種別・抜粋・確認箇所は取得できる場合に記録し、それらだけの不足で元担当へ補完を返さない。欠測自体は根拠表に残す。論文ID/タイトル、報告遺伝子名・ID/名前空間（取得不可なら確認理由付きunknown）・種と9回答等の必須記録、具体的矛盾・上位Tierへの手がかりに限定して従来のbounded followupを使う。既に確認したIDの取得不可を再探索義務やTier低下にせず、機械取得は元LLMの記録を上書きしない。旧方針と明示 `content_v1` の確認は維持する。詳細は [canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md) を参照。

#1002の標準T1は同種の同一・非常に近い具体的既報表現型と、適格な区間内相同投影に基づく第一段候補である。複数copy・理由付き他種proxyを許容し、適格な全copyを保持する。`same_gene_as_focal=unknown/no` のままでよく、同一座位・全長100%・一意orthologは不要。厳密identity手順は後段で必要な場合の追加調査とし、第一段T1のためだけに起動しない。原報・proxy・実member/配列/hashの対応とsupportは区別して保持する。

## 通常案件の知見の橋渡し

モデル植物・他植物で研究された機能や表現型から、対象区間の相同遺伝子を候補化する。全geneを台帳へ保持し、依頼された初期判断の範囲を記録するが、全geneの機能解明を完了条件にはしない。研究の乏しいgeneは不足理由を付けて未解決のまま残せる。

模式例のMtA→AtB→PsCでは、報告遺伝子MtAと代理参照AtBの対応・表現型の意味を元の調査LLMが判断し、原報の事実と理由を記録する。MtA配列を取得できなくても、AtBのidentity・配列と区間gene PsCへの投影を既存の機械検査に通して候補化できる。AtB自身の表現型論文や別の生物学reviewerを追加要件にしない。詳細は[報告遺伝子・代理参照・対象座位](CANONICAL-WET-MEMBERSHIP.ja.md#報告遺伝子代理参照対象座位を分ける)に従う。

## 同種の報告遺伝子と対象座位の同一性を追加評価する

以下は同種の報告geneと対象座位の同一性を追加評価するときの手順であり、異種proxyを使う通常の候補形成の必須工程ではない。

厳密な同一座位の追加評価が必要な場合、親Agentは公開配列の取得・照合まで実行し、未実行のまま「対応できない」と結論しない。論文調査中に得た原著のaccession/version、CDS protein_id、taxid、品種、確認箇所を記録する。区間blindのAへ対象配列・座標・他routeの結果を返さない。取得できない項目は理由付きnot_evaluatedで残し、既存の候補保持を進める。第一段T1では理由付きproxyと適格な区間内相同対応を使える。

論文・参照選択をfreezeした後、親Agentが既存 `reference-gene-resolve` の別モードを使う。研究者に手作業のTSV作成を依頼するのではなく、親Agentが原調査の出典から `sequence_requests.tsv` を生成する。必須列は `request_id, source_id, reported_gene, reference_accession, reference_taxid, focal_taxid`（TSV）、複数CDSを持つ登録には `cds_protein_id` も明示する。名前からBAC内のCDSを推測しない。`reference_accession` は公開INSDC登録IDで、版指定時は完全一致を要求する。

A の参照が対象種の UniProt entry である論文報告は、凍結後に機械的に request 化できる。`scripts/agentic_identity_requests.py requests --prepared <prepare stage.json> --sha256 <SHA> --out <dir>` は封印済み stage の A 論文表と record cache を読み、報告種が対象種の論文について UniProt record の EMBL 相互参照(INSDC 塩基 accession と CDS protein_id)から `sequence_requests.tsv` を書く。対象外の行は理由付きで manifest の `skipped` に残す。同一性の判定、Tier の変更、外部への問い合わせは行わない。whole proteome の `protein_gene_map.tsv` は、FASTA ヘッダーに明示的な `gene=` がある場合だけ `protein-map --proteome <FASTA> --out <map.tsv>` で作れる(推測はしない)。後段で厳密な同一座位評価が必要な場合に限り、A を凍結し第一段 Wet を作った後で request → 下記 resolver → `source-verify --reference-bridge` → assessment → `finish --identity-assessment` を実行する。A の `same_gene_as_focal=unknown` を追加評価する経路であり、通常T1の入口条件ではない。

```sh
candidate-qtl reference-gene-resolve \
  --sequence-requests sequence_requests.tsv --post-freeze \
  --proteome reference_proteins.fa --protein-gene-map protein_gene_map.tsv \
  --records public_insdc_records --out reference_bridge --runner auto --json
```

`--records` は `{accession}.embl` / `.gb` / `.gbk` の保存レコード。新規取得を許可する場合だけ `--allow-external-api` を付ける。公開accessionだけをENAへ送り、対象FASTAは送信しない。保存cacheはsupplied_record_cacheとして由来を明示し、ライブ取得と混同しない。whole proteomeとその `protein_id, gene_id` 対応表には区間外と未配置scaffoldも含める。供給された全参照を検索するため、区間だけの入力を全ゲノムの一意性証明と呼ばない。

出力manifestはrequest/queryの事前freeze、原レコード、入力SHA、全BLAST HSP、全コピー、両側coverageと使用条件を封印する。末尾stop記号はアミノ酸数から除く。既存のBLAST review閾値は候補投影用であり、同一遺伝子を機械的に認定する数値基準ではない。target coverage不足や、複数の原著遺伝子が一つのモデルの異なる領域に対応する状態を別に記録する。染色体配置を優先してscaffold候補を消さない。

各局所一致は `alignment_scope=single_hsp` として個別保存する。同じqueryがモデル内の別の同族領域にも一致しても、そのcoverageを足して「モデル全長が同一」としない。連続する複数HSPを含む真の全体対応は、親Agentが保存済み全HSPと遺伝子構造から判断する。

**座位対応と注釈モデル全体の同一性は別に評価する。** 親Agentは全長高一致と競合の差、原著・GFFの境界を確認し、記録された同一性の問いに答える。注釈不整合があっても座位対応までunknownに戻さず、一方で同じモデルにhitした別遺伝子をaliasとして統合しない。BLASTは表現型の因果関係を証明せず、元のTier/supportを自動変更しない。

A/Bでは既存 `source-verify` に `--reference-bridge reference_bridge/reference_sequence_bridge_manifest.json` を渡すと、元の論文IDと報告geneに結合した `paper_sequence_correspondence.tsv` ができる。これは **post_freeze_parent_only** のviewで、元のpaper_reportsやTierを上書きしない。原調査の不足を修正する場合は根拠と判断を新generationに記録する。区間結果をblind Aへ戻さない。

Direct Cの報告symbolから配列への関係を同じ凍結論文で明示的にreviewした場合は、[paper identity mapping](ROUTE-C-PAPER-IDENTITY-MAPPING.ja.md)を既存 `--identity-mapping` とreference bridgeで投影・review Wetへ接続できる。名前の共出現だけでは許可しない。原referenceが空欄の場合も原recordを保持し、reviewが認可したtargetを投影用viewだけへ補う。verified supportの条件は変わらない。

そのviewから親または独立agentが同一座位を評価したら、hash-boundの `post-freeze-identity-assessment-1.0.0` を既存 `wet-candidate-handoff --identity-assessment` または `scripts/claim_review_session.py handoff --identity-assessment` へ渡す。元paper、原著source、公式CDS、query/atlas、全競合と実reviewを結び付け、元回答を保存したまま対象memberの派生Tierへ消費する。初版は同種のwhole-proteome全長完全一致1 modelに限定し、これだけで同一性を自動認定せず明示reviewも要する。対象形質の同等性は原論文のチェックを再利用し、support・C保護資格は変更しない。詳細は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md#freeze後の同一座位評価を対象geneのtierへ渡す)を参照する。

standalone Wet CLIの対応先配列は `--identity-target-proteins interval_proteins.fa`、session handoffでは既存 `--proteins`。assessmentのcanonical FASTA SHAと実対象を一致させ、全対象geneの参照対応と完全一致競合を原配列から再検査する。原著・CDS・query/atlas/対象FASTA・bridge・引用位置と実reviewを同じ材料指紋へ拘束し、review後の資料差替えを受理しない。

通常案件の `practical_packet_run.py finish --identity-assessment assessment.json` も同じconsumerへ渡す。
対象配列にはpreparedの `packet_b/candidate_proteins.fa` を使い、assessmentをそのpath/SHAとB0 atlasへ結び付ける。
assessmentは元の場所から読み、内部の相対pathを保つ。指定しない既存runの挙動は変わらない。

親担当の資料準備と配送には `scripts/agentic_identity_review.py` を使える。
`prepare --prepared <stage.json> --sha256 <stage SHA> --bridge <bridge manifest> --members <A review candidates>
--members-manifest <producer manifest> --paper-records <保存済み原著JSON cache> --provider <provider>
--model <実reviewモデル> --run-id <実行ID> --out <新規directory>` は、既存memberに結び付いた
元paper、公式CDS、全競合、atlasと対象配列を封印して `review_template.json` を作る。
これは判断前の雛形であり、T1判定や完了reviewではない。親担当は全資料を確認し、別ファイルへ
`decision/reason/review_axes` の実判断を記録する。部分配列・別copy・原報と異なるproxyはunknown等で残す。
`compile --manifest <prepare_manifest.json> --sha256 <SHA> --review <実回答JSON配列> --out <新規directory>` は
実回答の引用位置と材料を結び、既存identity consumerで検証したassessmentを出力する。
そのassessmentを上記finish/handoffへ渡す。原回答の表現型チェック、候補集合、supportは変更しない。
資料不足と複数の異なる参照配列は `skipped` へ理由付きで残し、全correspondenceから追跡する。
この入口はAの封印済みpaper/memberが対象で、BやCへ別の判定経路を追加しない。
同じ公式レコードを複数要求で保存した場合は内容SHAで同一性を確認し、保存pathの違いだけで拒否しない。
UniProtの旧accessionは保存済み公式secondaryAccessionsだけで対応し、競合するaliasは拒否する。
原報が参照と別遺伝子・homologだと明記した行は同一遺伝子用requestを生成せず理由を残す。

任意Cでは既存 `route-c-direct-project` に同じ `--reference-bridge` を渡す。request_idは凍結Cのrecord_idに対応し、source_id、参照accession、登録CDS上の遺伝子/配列endpointも照合される。報告名とCDS qualifierの表記が異なる場合も、凍結済みの参照accessionが完全一致し、requestの報告名が凍結Cの同じ報告identityと一致し、同種の公式translationとsequence SHAが得られていれば、元の報告identityをreview-onlyの投影provenanceとして保持できる。これは公式名のalias登録や支持昇格ではなく、矛盾する明示的gene ID、別accession、未評価sequenceは許可しない。`--eligible-hits` と同時指定せず、新しい投影generationを使う。既存のappend-only投影→members→Wetへreviewとして届き、最終Wetでもbridgeの祖先SHAを再検証する。取得不能・no-hitは独立したexact支持を否定せず、測定された注釈モデル不整合は支持昇格を保留する。CをA/Bの必須条件にはしない。

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

2026-09-23 / Change class: CORE_INTEGRATION / Issue #758。

通常A/Bの論文報告は [canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md) の根拠表契約に従う。取込が出力する `paper_evidence_table.tsv` とmanifestの `paper_evidence.recording_state` を確認する。`needs_completion` なら `paper_followup_requests.json` の対象項目だけを元の調査担当へ、同じroute入力とともに原則1バッチ返す。担当は元payload/responseの新しいコピーを修正し、`followup_state=completed` と確認した不明理由を記録して新出力先へ再取込する。空欄をunknownに機械変換して完了扱いにせず、実行していない依頼は未実行のまま報告する。未達の記録状態は候補除外・支持昇格に使わない。5報の件数合わせや全候補の再調査、AへのB/区間入力の混入、freezeの書換えは行わない。書誌確認後の表には既存のDOI/タイトル照合結果が追加される。

原調査担当へ渡す作業指示には、9回答・理由・報告geneと参照の対応を先に記録し、取込用writerはその判断を配送することを含める。参照取得の成否、全長参照の選別、`Fragment`表示、alignment coverageから回答を自動生成・変更しない。配列の完全性、名前付き参照の同定、機能注釈の有無、対象座位との同一性は別の問いであり、原担当が資料に基づいて答える。部分配列だからunknown、参照が取得できたからyes、のどちらにも一律変換しない。取込時に同じ担当がrawの9回答・理由・対応と正式表の差分を確認し、不一致があれば既存の限定follow-upで理由を確かめる。単なる配送修正では科学回答を保持し、科学回答自体の訂正が必要ならold/new・根拠を別generationへ記録してから再取込する。追加のTier審査役や全候補の再調査は設けない。

研究者は対象種、表現型、QTL区間、FASTA、GFF3を渡し、「この区間に原因になりうる遺伝子はあるか」と尋ねる。親エージェントが調査・分割・再開・統合・報告を行う。約20遺伝子は実験可能な目安であり、最低件数・切断上限・完了条件にしない。必要確認済みの12件でも終了でき、既存の弱いreview・別機構・曖昧copyは件数だけで落とさない。比較待ち・未回答は完了と区別する。

### 科学担当の起動時に情報境界を確認する

別worker IDや「blind」という指示だけでは、A/Bの独立性を保証しない。親が区間や他routeの調査結果を知っている場合は、親の会話履歴を引き継がない新規contextへ分配する。起動前にcontext mode、実行model、許可する入力artifactを既存plan/receiptへ記録する。Aには表現型側の入力とschemaだけ、Bには自laneの凍結材料と表現型を渡し、親の研究メモや他laneの回答・選択論文を混ぜない。

Codexの `spawn_agent` では `fork_turns="none"` を明示し、そのときに依頼されたmodel/reasoningを設定する。既定のfull-history forkは、子への本文に区間を貼らなくても親履歴を漏らし得る。OpenCodeでは別の新規native sessionとlane専用cwd/入力範囲を使い、実際のread/search/fetch権限を確認する。限定したA研究担当なら、cwd内read/writeとWeb検索/取得を許可し、外部directory・shell・孫agentを拒否する形も使える。親/Bの必要なツールを一律に禁止する規則ではない。

区間を見た親は正式なblind A回答を自作しない。隔離不備が見つかった世代はdiagnosticとして保全し、その発見を渡さず新しい担当で影響laneを再実行する。元freezeを書き換えたり、別session名だけで独立実行済みへ読み替えたりしない。

FASTA/GFF3から候補を抽出する通常案件はWSL/Linuxで実行し、`alignment` extra（`pysam`）を入れる。標準wheelだけのclean installでは抽出に必要な`pysam`が入らない。通常入力には [ordinary config template](../workflow/practical_packet_run.ordinary.template.json) をコピーして使う。既存の `practical_packet_run.template.json` は抽出済みpacket用であり、両形式を混ぜない。採否規則は [canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md) を参照する。

## 共有元ゲノムと区間パケット

通常案件の `prepare` は `reference_storage` を省略しても `shared_readonly` を使う。ordinary templateもこの方式を標準にする。元FASTA/GFFを読み取り専用で直接参照し、案件には抽出された区間遺伝子・配列と由来を保存する。全ゲノムのコピーや管理用参照への初回登録は不要である。
`"reference_version": "使用するassembly/annotation release名"` は必須。明示指定する場合は `"reference_storage": "shared_readonly"` とする。
FASTA/GFFは共有元の絶対パスを指定し、非圧縮FASTAと既存の `<FASTA>.fai` を用意する。共有元の索引は作成・再生成せず、索引欠落や読込失敗は停止する。抽出前に別の一時領域へ独立索引を生成し、全contigの名前・長さ・offset・linebases・linewidthを既存FAIと照合する。不整合は停止し、共有元のFAIを修復しない。この検査は全FASTAの追加読み取りを伴い、所要時間を抽出manifestの `fasta_index_validation` に保存する。一時索引は検査後に破棄する。元ファイルと出力先は別にし、共有元への同時書き込みを避ける。

既存の抽出処理でQTLと重なる全遺伝子、CDS、タンパク質、遺伝子全長の配列を保存する。境界に重なる遺伝子は区間外の部分も保持する。元FASTA/GFFをコピーせず、抽出時に各ファイル全体のSHA-256を一度計算するため、大きい元データの読み取り時間は残る。抽出中の元ファイルの置換・size/mtime/ctime変更は拒否する。これは外部プロセスによる書き込みをロックする機能ではない。

`inputs/reference_asset.json` に元パス、release、SHA、seqid、interval、索引のコピーとSHA、抽出manifestのSHAを記録する。stageのhash対象は凍結区間パケット・索引・由来manifestであり、外部の共有元は由来としてのみ参照する。`review` / `finish` は元ゲノムを再読せず、元データが後で移動・変更されても凍結済み区間を消費する。区間内artifactの改変は従来どおり拒否する。新しい区間の `prepare` は元データを改めて検証・抽出する。既定の `snapshot` モードと既存packet入力は維持する。

GFF3の一次転写産物と成熟RNA産物は `Parent` / `Derives_from` の階層を保持する。成熟産物や既知のtranscript IDをplaceholder geneへ置き換えない。正規のgene-less transcript／CDS-onlyモデルの補完は維持する。循環・曖昧な祖先・異なるfeature type間のID衝突・未解決産物は `gff_dialect_qc` に記録し、生物学的陰性とは扱わない。

### 原本を保護して同じ参照・区間を再利用する

元データから独立した全体コピーと共通区間cacheが必要な場合だけ、ordinary config に次を追加して `managed_cache` を明示指定する。これは任意の方式であり、ordinary templateの初期設定ではない。元の `reference_fasta` / `annotation_gff3` はそのまま入力し、参照release名と、案件出力とは別の共通保存先を指定する。

```json
{
  "reference_storage": "managed_cache",
  "reference_version": "assembly-and-annotation-release",
  "reference_cache_dir": "/data/qtl-reference-cache"
}
```

初回は原本から独立した解析用コピーを作り、そのコピー側にFASTA索引を生成する。原本のFAIは不要で、原本への索引作成や書換えは行わない。FASTAは非圧縮を使う。以後、同じ参照内容とreleaseを指定すれば共通コピーを再利用する。管理ファイルは読取専用にするが、同じOS権限を持つ別プロセスの書込みを禁止するロックではない。複数の実行者で共有する場合は保存先のアクセス権や読取専用mountも運用側で設定する。

区間抽出は同じ保存先の `intervals` 以下に再利用可能な形で保存する。必要なら `interval_cache_dir` で別の共通保存先を指定できる。FASTA/GFF/FAIの全文SHA、区間、seqid mapping、抽出コードと依存版で対応を確認するため、名前やmtimeが同じでも内容が変われば以前の抽出を流用しない。変更・欠落したcacheは明示的なエラーにし、未完成のものを完成済みとして扱わない。

案件側に保存するのは従来と同じ7つの抽出ファイル、QC・由来manifest、案件固有のA/B材料である。QTLと重なる全gene、境界外まで続くgene/CDS、no-CDS・noncoding・pseudogeneの台帳を現行の抽出規則どおり保持する。case ID、phenotype、A/B回答、独立review、B0の解釈結果は区間cacheへ混ぜない。案件の抽出ファイルも独立したコピーであり、編集が共通cacheへ伝わるhardlinkは作らない。

`inputs/reference_asset.json` の `managed_asset.cache_hit` と `interval_cache.hit` から、共通コピーと区間抽出それぞれの再利用を確認できる。毎回、原本と管理コピーの全文hashによる確認は残る。初回のコピー・索引作成も必要な時間なので、初回登録と再利用時を分けて評価する。主に削減するのは案件ごとの巨大ファイル書込み、同一区間のGFF解析・配列抽出、後段の全ゲノム再読である。

### 明示したprotein DBを共有する

通常のBLAST protein DBを各案件へ複製せずに使う場合は、既存の `protein_db` に加え、`"protein_db_storage": "shared_readonly"` と `"protein_db_version": "database-release"` を指定する。この設定はFASTA/GFFの保存方式とは独立して選べる。省略時は従来どおりDBを案件へコピーする。

検索対象は指定DB全体のままで、区間候補のうち材料が未解決のタンパク質を照会する既存B0処理を使う。BLASTの閾値やhit数は変えない。prepareは指定prefixの構成ファイルの集合・全文SHA・releaseを確認し、B0前後の変更を拒否する。`inputs/protein_db_asset.json` に由来を保存し、検索結果とqueryは案件側へ凍結する。外部prefixを参照する `.pal` alias DBはこの共有方式では受理しない。単一DBまたは通常のvolume構成を直接指定する。prefix外の共通taxonomy DB等はこのmanifestの固定対象外だが、取得した分類情報を含む実検索出力は案件のstageへ保存する。review/finishで元DBの再読は不要である。

## WSLとclean wheelの準備

以下はWSL内のbash例。`repo`は利用するcheckout、`wheel`はそのcommitからbuild済みのwheelの絶対パスに置き換える。Pythonは3.12以上を使い、空のvenvへ導入する。FASTA/GFFやconfig内のpathもWSLから読めるLinux pathにする。

```bash
repo=/mnt/c/cqtl-repo
wheel=/mnt/c/cqtl-dist/candidate_qtl-VERSION-py3-none-any.whl
python3 -m venv /tmp/cqtl-ordinary-env
source /tmp/cqtl-ordinary-env/bin/activate
unset PYTHONPATH
python -m pip install -c "$repo/constraints.txt" "${wheel}[alignment]"
python -c 'import sys, candidate_qtl, pysam; print(sys.executable); print(candidate_qtl.__file__); print(pysam.__version__)'
candidate-qtl doctor --json
cd "$repo"
```

`candidate_qtl.__file__`がこのvenvの`site-packages`を指すことを保存する。`pip install -e .`や`PYTHONPATH=src`でsourceへfallbackするとclean wheelの到達性確認にならない。wheelはproduction packageを提供し、下記の親Agent用scriptとpacket入口はcheckoutから実行する。必要なローカルBLAST等は`doctor`の結果と実際に使う資源に合わせて確認する。

### 既存成果物の読み取り専用E2E監査

通常案件の実行能力・根拠配送・第一段/最終Wet状態・再開・実測負担の受入は E2E監査runbookへまとめる。既存 `practical_packet_run.py audit --finished STAGE --sha256 SHA` は封印済み成果物を読むのみで、検索能力・token・時間・truthの推定や新しい採否は行わない。現実のWeb検索ツールと原文配送は別途、実行workerと取得receiptで確認する。

## 実行モデルの原則

この手順は「1つのLLMに全部を解かせる」ためのものではない。親エージェントがcase全体を保持し、Route A/B/C、Route B shard、boundary review等をsubagentへ渡す。役割分担やモデルは実行環境に応じて変えてよい。

Sol / CodexとOpenCodeの枠管理、資料配送、独立レビュー、export-first回復は
協業手順を参照する。運用toy試験は科学性能や通常QTLのE2E確認とは区別する。

途中のTSV/JSON/manifestは残してよい。これらは **Agent間の型付き通信 + checkpoint + provenance** であり、削除することを自動化の目標にしない。親Agentが生成・結合してよいが、以下はmachine側で必ず検査する。

- `case_id` / source・input SHA-256
- assigned `prompt_id` / `route_b_unit_id` / `candidate_gene_id` のroster完全性
- unknown / duplicate ID
- executor / model / reviewer provenance
- missing / failed / not-run を `not_evaluated` として保存すること

研究者に求めるのは通常、対象種・表現型・QTL・FASTA・GFF3等の案件入力であり、per-gene表の作成、1000 geneの手分割、A/B/C結果の手結合ではない。

## 親エージェントの仕事

**通常はA/C優先（2026-11-08）。** A経路とC経路を必ず実行し、区間の全遺伝子を調べるB経路（手順4）は、ユーザーが求めたとき、または全遺伝子の確認が必要な具体的な問いがあるときだけ実行する任意経路とする。Bを実行しない場合、Wetは A/C のsupported・review候補だけで構成され、区間注釈による補充は行わない。`finish` / `wet-candidate-handoff` の `route_a_rethink` を必ず確認する。

- `required`（Wetが空）または `recommended`（`--min-route-candidates` 未満、既定3）なら、Aを新しい世代で一度考え直す。新しいA担当に同じ表現型側の入力だけを渡し、出力先を新しくし、元のA freezeは保持する。区間遺伝子の一覧や前回のA/Cの候補名は渡さない。
- 考え直しは既定で1回とする。その後も足りなければ、件数不足をそのまま報告し、任意のBを提案する。
- 閾値は作業量の合図であり科学的な切断値ではない。Wetの構成を増減させない。

1. 入力と候補数を確認し、case IDと作業ディレクトリを決める。FASTA/GFF/区間からの候補抽出は既存 `practical_packet_run.py prepare` が行う。中間候補TSVを手書きしない。既存のBLAST DB、source/record cache、reference function等があればconfigで指定し、同一世代の資源を再利用する。
   通常案件では config に `"external_retrieval": true` と実際の `"a_payload_source"`(Agent回答なら `agent`)を指定する(#986/#987、ordinary template の既定)。`prepare` は自分が取り込んだ A 仮説表から参照 record cache を取得し、A の引用 PMID/DOI の ID・タイトル書誌を1回取得して stage に封印する。prepared がこの設定なら、`review` も B claim の引用書誌を同様に取得して封印する。その後の照合は封印済み cache で offline に行う。cache を親が手で組み立てる必要はない。指定しない場合(固定資料の再生など)、A の引用は `not_attempted_offline` のまま A の区間投影が 0 になり得る。これは陰性ではなく未照合であり、report の degraded coverage を確認する。provider 障害は取得 ledger に残り、健全な record は使える。
2. A担当には表現型と対象種を渡し、機構、関連する既知遺伝子・相同遺伝子、原著を調べさせる。区間遺伝子やB/Cの回答をA凍結前に渡さない。`a_payload`は既存 `mechanism-import` 契約のJSONで保存する。候補が見つからなければ空rosterと探索状態を記録し、創作しない。
3. C担当は**A/Bと別の新規context**で、対象種と表現型だけを受け取り「この形質に関連する既報gene/locus/variantは何か」を自由に調査する。ハーネスで実際に利用できる標準Web検索・取得を優先し、`route-c-direct-search` は任意の追加手段とする。専用の事前gene/family一覧、A/B回答、区間・truth、親の調査履歴は渡さない。Aと重複した報告や後で区間外と判明する報告も、pre-freezeで選別しない。Cは直接の既報事実・関連を回収する経路で、一般機能だけの新仮説はA/Bの仕事とする。研究後に既存17列proposal（`route_c_direct_proposals.empty.tsv`）へ整形し、`route-c-direct-freeze`へ渡す。通常調査プロンプトは`render_harness_native_direct_prompt()`を使う。既存`route_c_practical_prompt.txt`は凍結ベンチ用として変更しない。検索失敗と検索済みの該当なしを分ける。
4. （任意のB経路）B0は表現型を見ずに全区間遺伝子を圧縮する。`prepare` の `bundle.prompts` を `scripts/agentic_route_b.py prepare` で分割し、各shardを担当エージェントへ渡す。親は全prompt ID/unit IDを統合する。shardごとのTop-Kや未回答の黙殺は禁止。Bの担当には各promptに含まれる表現型と材料を使って、マクロ/ミクロな記述差、機構段階の違い、複数の原因・コピー、SV/発現経路を検討させる。根拠の無いclaimは作らない。
5. 新policyでは独立レビュー回答なしでfinishする。以下の原文レビュー手順は明示した旧content policy用である。A/Bのclaimに取得済み原文があれば、提案者とは別contextの担当が原文、種、遺伝子、実験、形質の結合をレビューする。`review`はrouteごとの凍結`source_review_records.records`が厳密に空の場合だけ、既存no-source helperによる機械未実行回答をstageへ保存する。これは独立レビューの代用や生物学的否定ではない。Aの出典・参照recordは `prepare` 前に取得する。Bの提案後に新しい出典・参照recordを取得した場合は `review --source-cache <file|directory> --record-cache <records directory>` へ明示して、新しいreview stageを作る。追加資料はsnapshotとSHA・取得元pathをstageに保存し、source cache指定時は既に解決済みの同じA hypothesis rowsも再検証する。Aのraw proposal・解決ID・旧preparedは変更しない。供給cacheは各軸の使用cacheを置き換えるため、既存の必要な資料も含める。未指定の軸はprepare時のcacheを使い、暗黙の取得・旧cacheへのfallback・支持昇格は行わない。
6. `review` が封印したA/B checkpointを `finish` へ、Cの検証済みmembersとproducer manifestも同じhandoffへ渡す。通常の `initial_bibliography_v1` では `--reviews` を渡さない。旧content policyの独立レビュー回答は、その世代の契約に従って渡す。Cに有効な結果が無い場合は理由を報告し、A/Bを捨てない。Cの同種direct supportと異種・未確定reviewを混同しない。現行のC source binding/projectionの暫定条件は `docs/ROUTE-C-DIRECT-RESCUE.ja.md` に従う。
7. `boundary_scope_manifest` が出たら `scripts/agentic_boundary.py` で境界候補を分割し、元担当の理由付き採用/保留提案をまとめ、`finish-wet-review` を実行する。通常の `upper_evidence_scoped_v3` は累積25を初めて超えるTierだけを比較し、上位は保持、後続Tierは台帳へ残す。`original_assessment_v1` としてhash・roster・executor ID/modelを機械照合し、別のTier/本文reviewerや下位追加比較を要求しない。旧v2・版なしscopeは旧policyで再生する。欠落は `not_evaluated` で保持し、支持未確認だけを `defer` 理由にしない。最終実件数、採用/比較待ち/保留、A/B/C由来、実施範囲・未確認事項を返す。

## 固定入口とcheckpoint

### 通常Bの検索と資料の引渡し

B担当は固定資料を読むだけの役ではなく、必要に応じて通常のWeb検索・関連論文検索・追加資料取得を行う。初期判断から考察途中まで、注釈・相同遺伝子の機能・内部知識・検索・論文確認を行き来してよい。一次判断が終わるまで検索を禁じる段階ゲート、1遺伝子1論文の上限や一律検索ノルマは設けない。弱いが具体的な機能的つながりも検討し、未知という理由だけでは採用しない。論文で確認した参照機能と、対象表現型への条件付き推論を分けて記録する。

通常Bは、まず区間gene/相同先の既知機能・作用方向・器官・時期と対応の不確実性を確認し、その後に今回の表現型との整合性を中立的に判断する。調査した機能から具体的bridgeが見つからない場合は理由付きの空claimを許し、未探索・取得不能と区別する。旧 `source_function_v1` bundleのB-only Wet reviewには、実際の外部参照機能、区間対応、具体的な条件付き推論の独立reviewが必要である。対象種の因果論文・同一表現型は必須でなく、汎用class/GOや出典IDの存在だけでは足りない。資料不足は全gene台帳へ残し、A/Cの独立適格性を拒否しない。詳細は [canonical Route B](CANONICAL-ROUTE-B.ja.md) と [Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md)。

開始時に必須確認・対象ID・終了stage・既存round budgetを記録し、追加調査は採否を変える具体的な問いへ結び付ける。未解決だからという理由だけで次roundへ進めず、予算内に終わらない分はpartial/未評価・再開先を残す。原著1〜2報のsoft budgetは未採用であり、初期Web検索や矛盾・copy差に必要な複数資料を禁止しない。同じ参照機能質問の共有とcopy固有適用は [review session](CLAIM-REVIEW-SESSION.ja.md) を使い、保存bytes削減と実読解・判断削減を区別する。

再開時に元回答・資料を再利用しても、それだけで追加検索を禁止したり、1roundへ制限したりしない。固定資料のみの再生を指定された場合はその指示・policyを既存planへ結び付ける。親が独自に範囲を狭めた場合も親判断と記録し、フローの指定やprovider障害へ付け替えない。roundの予算と実施数を分け、未消化で終えた理由と、既存の具体的な手掛かりに対する未実施作業を残す。全資料不足geneへの一律検索や停止済みrunの再開は意味しない。

Tierの説明では全paper reportsを調べる。代表Tierが他種論文でも、同じ論文・別レコードに同種の観察が残ることがある。同種原報がない、未検索、原報はあるが対象座位との対応が未確定、計算が誤っている、を区別する。原著accessionと対象モデルの対応・競合copyを既存の同一性評価手順で確認し、部分配列の完全一致を一意座位の証明にしない。Typeラベルや期待した件数にTierを合わせない。

親はBを割り当てる時点で、その担当に実際のsearch/fetchツールが提供されていることを確認する。CodexではホストのWeb検索・取得ツールを使う。`agentic_route_b.py`や`practical_packet_run.py`は資料と回答を運ぶ処理であり、モデルのツールを追加するランナーではない。OpenCodeではmanagerのB権限と使用providerの利用可能性を確認する。使えなければ担当・クエリ・取得不能理由を残し、外部障害や公開制限はそのまま報告する。利用可能な通常の取得方法を使えばよく、新しい認証や有料providerを前提にしない。

検索語・取得URL/PMID/DOI・取得時刻・読んだ本文/抄録の範囲と参照遺伝子/種の結び付きを作業資料へ保存する。新しい文献を提案に使う場合は、既存の`literature` TSV入力（`interval_gene_ids`, `reference_gene`, `source_id`, `source_url`, `pmid`/`doi`等）へ結び付け、元のB0・全gene/copyを保持した新しい`prepare`出力でbundleをfreezeする。旧bundle・旧回答を編集せず、新promptに対して回答する。新旧で変わらない回答の再利用は既存reuse plannerのhash条件に従う。本文/抄録は既存`source_cache`へ保存し、`review --source-cache`へ既存の必要資料も含めて渡す。参照同定recordの追加には`--record-cache`を使う。検索結果や内部知識だけを検証済み根拠へ昇格させない。

検証器のoffline設定は検索禁止ではない。独立レビューの再現可能な入力を固定するための設定であり、検索担当が取得した資料を受け取る。明示的に固定資料だけを読む比較ベンチは別契約で、途中から検索を追加して同じ比較と呼ばない。通常調査へその比較専用の禁止指示を流用しない。

OpenCodeの`opencode_manager.opencode_config`は通常のB担当に`websearch`/`webfetch`を許可するが、`science_adapter`では固定回答用に`websearch`を除去し、`webfetch`も`allow_public_webfetch`時だけ許可する。このadapterを通常Bの探索担当として起動しない。通常の検索可能な担当が調査と資料の引渡しを行い、adapterを使う場合はその後の凍結済み回答・検証に用途を限定する。設定生成はツールの利用成功を保証しないため、実際の担当環境で検索/取得が利用可能かを確認し、API認証エラー・通信失敗・取得成功を分けて記録する。

上記venvを有効にした、通常の `initial_bibliography_v1` のWSL bash例。`$run`以下の各出力は別directoryとし、既存stageの中身を書き換えない。各段のJSONに記録された `manifest` と `sha256` を次段へ渡す。

```bash
h=workflow/practical_packet_run.py
run=/mnt/c/cqtl-run/my-case
python $h prepare --config "$run/config.json" --out-dir "$run/01-prepared"
# prepareのJSONにある bundle.prompts を指定。親がshardを各担当へ配る。
python scripts/agentic_route_b.py prepare --prompts $PROMPTS --out "$run/b-shards" --max-prompts 8 --max-material-bytes 262144
python scripts/agentic_route_b.py preflight --prompts $PROMPTS --response $B_RESPONSE_1 --response $B_RESPONSE_2 --executor $B_EXECUTOR --model $B_MODEL --reading-manifest "$run/b-shards/prepare_manifest.json"
python scripts/agentic_route_b.py merge --prompts $PROMPTS --response $B_RESPONSE_1 --response $B_RESPONSE_2 --out "$run/b-responses.json" --executor $B_EXECUTOR --model $B_MODEL --reading-manifest "$run/b-shards/prepare_manifest.json"
python $h review --prepared $PREPARED_MANIFEST --sha256 $PREPARED_SHA --responses "$run/b-responses.json" --b-execution-manifest "$run/b-responses_manifest.json" --out-dir "$run/02-reviewed"
# B提案後の取得資料がある場合は、同じprepared/responseから別stageへ渡す。
# cache取得は外部取得を許可した既存入口等で明示的に行い、review自体はoffline。
python $h review --prepared $PREPARED_MANIFEST --sha256 $PREPARED_SHA --responses "$run/b-responses.json" --b-execution-manifest "$run/b-responses_manifest.json" --source-cache "$run/post-proposal-sources" --record-cache "$run/post-proposal-records" --out-dir "$run/02-reviewed-with-materials"
# 通常policyのreview結果は ready_for_initial_candidates。使用するreview stageのmanifest/SHAを指定する。
# 独立原文レビュー回答は第一段に要求せず、--reviewsを渡さない。
python $h finish --reviewed $REVIEWED_MANIFEST --sha256 $REVIEWED_SHA --route-c-members $C_MEMBERS --route-c-manifest $C_MANIFEST --max-adaptive-rounds 2 --out-dir "$run/03-first-wet"
```

任意のBを実行しない場合（A/C優先）は、`scripts/agentic_route_b.py merge --prompts $PROMPTS --out "$run/b-not-run.json"` で全promptを`not_evaluated` にした応答を作り、`review ... --responses "$run/b-not-run.json" --route-b-not-run "<理由>"` を使う。完了済みのB回答が1件でもあれば拒否し、`--b-execution-manifest` とは併用できない。reportとauditではBを `not_run_optional` と表示し、実行済みや陰性として数えない。

通常の `initial_bibliography_v1` では、`review` は `state=ready_for_initial_candidates` を封印し、独立原文レビューを要求しない。`finish` は非空の `--reviews` 回答を受け取ると停止するため、原文の有無によらずこの引数を省く。この状態は初期候補の採用へ進めるcheckpointであり、独立した意味レビューの実施やWet境界比較の完了を表さない。Bの機能解釈、転移検証、必要なWet境界の採用・保留比較は既存の手順で行う。

#### 旧content policyのcheckpointを再開する場合

`state=awaiting_independent_reviews` の旧content policyには元の契約を使う。A/Bのreview promptはreview結果JSONの `sessions.a.prompt` / `sessions.b.prompt` にあり、`independent_review_routes` にあるrouteの実reviewを `reviews.json` へ記録して `finish --reviews "$run/reviews.json"` に渡す。

機械未実行の補完に対応したcontent-policy stageは `machine_reviews` と `independent_review_routes` をsealする。後者が空でも実レビュー完了を意味しない。`finish` は省略された空原文routeだけを機械記録で補い、全source不在なら `--reviews` なしで進める。原文ありrouteが未供給なら出力作成前に停止する。1件でも原文があるrouteは、原文不足claimが混在していても全体を機械未実行へ置き換えない。`finish` はstage hashと空原文条件、回答の全ID・内容、`no-source-machine-*` / `not_invoked` を再検証する。明示的に供給した実reviewは優先し、機械記録のない旧stageには自動補完を遡及適用しない。その場合は既存の `scripts/agentic_no_source_review.py` を使う。記録は `not_executed_no_source`、全claimは `not_evaluated` であり、意味レビュー実施数に数えない。

取込前の形式検査には、Aは `candidate-qtl mechanism-import --payload "$A_PAYLOAD" --preflight --json`、Bは上記 `preflight` を使う。元回答を読み取り専用で検査し、問題のrecord/claim位置をJSONで返す。Aのfreeze manifestを指定した場合はそのhashも照合する。Bは原回答の形式・provenance・prompt/unit対応を検査する入口であり、checkpointのmanifest/hash検査は既存の正式mergeに残る。どちらもモデル実行・科学的支持の検証・調査完了の認定ではない。不正な論文を削除して機能claimだけに置き換えず、理由付きで元回答を保持し、新しいattemptで修正する。

各 `paper_reports_json` 要素では、科学的な論文項目と調査メモを分ける。`paper_checklist` は従来の9問の回答、`source_excerpt` 等は従来の文字列のままとし、読了範囲・反証・回答理由等は任意の `research_annotations` JSON objectへ保存する。例えば `"research_annotations":{"full_paper_read":false,"paper_checklist_reasons":{"phenotype_reported":"抄録の観察記載"},"counterevidence":["copyは未分離"]}` とする。この欄は原型を保持する監査情報であり、Tier・source/gene binding・支持状態を書き換えない。旧来の誤配置に限り `unit_id`, `counterevidence`, `stage`, `read_scope`, `conditional_focal_reasoning`, `source_url`, `reference_species`, `formal_claim_state`, `retrieved_at`, `full_paper_read`, `paper_checklist_reasons` を専用欄へ移す。同名の値が競合すれば拒否し、その他の未知項目や科学項目の不正な型も従来どおり拒否する。注記のない旧recordには新しい空欄を追加しない。注記がある世代の論文表には `research_annotations_json`、論文summaryには全ての異なる注記を保持し、旧凍結成果物は変更しない。

### 通常Bの共通台帳・軽量screenから追加調査へ

通常のB0→原担当による全区間geneの初期機能判断は維持したまま、**追加の深い文献探索**だけを条件付きで省くための入口。新しいTier決定表、Wet selector、別のscientific reviewerではない。共通の参照遺伝子・speciesについて軽く調べた論文検索結果を一つのscreen記録へ集約し、区間内の各copyは**それぞれ**暫定Tier・出典と元の初期判断を保持する。T4/T5は現時点の根拠を示すだけで、上位論文の不在や既存のT1〜T3を否定しない。正式なTierは従来のpaper checklist（9項目）と既存の区間投影・Wet consumerだけで決める。

1. `prepare --research-bundle ...` で生成される `research/initial-material-*.json` 等の共通資料を使って、元B担当が機能・GO・pathway・mutant phenotypeを含む初期判断を行い、既存形式の `assessments.json` を作る。最低限 `bundle_sha256` と各geneの `unit_id/gene_id/state/reason/used_material_pointers/assessed_fields` を記録する。単に資料を一覧表示しただけでは判断完了にならない。
2. 同じ元担当が、特定のreference gene/speciesの同種・他種の表現型既報を軽くscreenし、検査したquery、取得状態、論文ID、確認範囲、具体的な上位手掛かりを記録する。`no_concrete_upper_lead` は**調べた範囲での具体的手掛かりなし**であり、文献不存在の証明ではない。未実施・検索失敗・未読は `not_assessed` / `retrieval_failed` として区別。単なるfamily/domain名や別speciesの同名遺伝子からscreenを流用しない。
3. 以下の小さいplanを生成し、`gene_actions` / `unit_actions` と全screen台帳を親が読む。**`light_stop_eligible`は広い追加深掘りを優先しなくてよいという助言に限る**。全原geneの正式回答、異なるcopyの対応、9回答、資料/claimのfreeze、既存支持確認と最終Wetを飛ばす許可ではない。

~~~bash
python scripts/agentic_route_b.py summarize-initial \
  --bundle "$BUNDLE" --bundle-sha256 "$BUNDLE_SHA" \
  --assessments "$run/b-initial-assessments.json" \
  --out "$run/b-initial-summary"

python scripts/agentic_route_b.py plan-light-followup \
  --summary "$run/b-initial-summary/initial_assessment_summary.json" \
  --screening "$run/b-light-screening.json" \
  --out "$run/b-light-followup"
~~~

通常の `prepare --research-bundle` は `research/manifest.json` に `phase=initial_screen` と原担当用taskを出す。親はこのtaskから初期判断・共通screenを受け取り、計画を付けた同じ担当laneの継続へ進む。この継続は会話履歴の続きではなく、保存した台帳・planのpath/SHA-256から新しい文脈で再開してよい。これは既存prepareの配送・調査指示で、モデルを自動起動する新runnerではない。

通常は原判断・暫定判断・screenを**同じ** `initial_assessments.json` に保存する。完全に同じ機能contextに対し、元担当が明示した適用memberだけを既存集計が展開する。異なるGO/pathway/表現型・copy固有の差は別recordへ分け、同じ理由をgeneごとに再生成しない。

~~~json
{
  "bundle_sha256": "<pinned bundle SHA>",
  "shared_assessments": [{
    "assessment_id": "J1", "material_id": "<initial-material material_id>",
    "member_gene_ids": ["<copy1>", "<copy2>"],
    "state": "assessed_claim", "reason": "<original shared functional judgment>",
    "assessed_fields": ["reference_function"],
    "research_planning": {"provisional_tier": "T4", "screen_ids": ["S1"]}
  }],
  "research_screening": {"reference_screens": [{
    "screen_id": "S1", "reference_gene": "<frozen endpoint>",
    "reference_species": "Arabidopsis thaliana", "member_gene_ids": ["<copy1>", "<copy2>"],
    "state": "no_concrete_upper_lead", "retrieval_state": "completed",
    "queries": ["<actual query>"], "checked_species_scopes": ["focal_species", "other_plants"],
    "screened_source_ids": [], "reason": "<actually inspected scope and remaining limits>"
  }]}
}
~~~

`material_id` は完全な機能contextのhashであり、単なるfamily名やbest hitではない。`member_gene_ids` は原担当の明示的な適用判断で、機械が全copyへ勝手に配るものではない。機能材料が本当にないmemberだけは `material_id=""` と空の `assessed_fields` で `assessed_insufficient_material` / `not_assessed` を共有できる。これは外部情報の不存在や陰性の判断ではない。従来のgene別 `assessments` とcopy固有の例外も併用できる。

`research_planning` は元判断と同じrecordに置く。`tier_reason` を省略した場合はその原判断の `reason` を再利用する。複数recordの計画が競合すれば未解決とし、必要なgeneだけ `research_screening.gene_assessments` で元担当が明示する。既存の上位判断は下位の計画で消さず、正式Tierは従来のconsumerだけで決める。

~~~bash
# 原判断内のresearch_screeningを使い、同じ1000geneの説明を別fileへ書き直さない。
python scripts/agentic_route_b.py plan-light-followup \
  --summary "$run/b-initial-summary/initial_assessment_summary.json" \
  --out "$run/b-light-followup"
# PLAN_SHAは上の出力plan.sha256。同じlaneの原担当を新しい文脈で再開してよい（会話の続きは条件ではない）。
python scripts/agentic_route_b.py prepare --prompts "$PROMPTS" \
  --research-bundle "$BUNDLE" --research-bundle-sha256 "$BUNDLE_SHA" \
  --followup-plan "$run/b-light-followup/followup_plan.json" --followup-plan-sha256 "$PLAN_SHA" \
  --out "$run/b-followup-view"
~~~

継続の `research/manifest.json` は `phase=followup`、共通の原台帳への `reuse_initial_ledger`、全geneの `research_action_groups` を持つ。`pages` は追加判断対象unitの索引、`initial_material.pages` は読了記録がない機能contextの差分だけになる。既読の全bundle・全理由・全論文を再表示せず、同じlaneの元担当が新しい文脈で未解決の問いに集中する（会話の続きを条件にしない、#1004）。元bundleと全canonical prompt/read receiptは省略・変更しない。全geneを再研究済みにしたり未評価をno-claimへ変えたりせず、正式回答は既存mergeへ返す。新資料が得られたら既存extend/reuse-planで再凍結する。旧bundleのplanを新bundleへ流用するとhash検査で止まるため、別世代の計画か既存の正式回答差分経路へ進む。

以下の別file `b-light-screening.json` は従来形式・外部producer向けの互換例。通常は同じ内容を上記 `research_screening` に保存し、`--screening` を省略する（`<...>` は実際の凍結済みIDや元調査回答で置き換える。架空論文を入れない）。

~~~json
{
  "contract": "route-b-light-tier-screen-1.0.0",
  "case_id": "<case_id>",
  "initial_summary_sha256": "<summary file SHA-256>",
  "reference_screens": [
    {
      "screen_id": "S1",
      "reference_gene": "<frozen reference endpoint>",
      "reference_species": "Arabidopsis thaliana",
      "member_gene_ids": ["<interval gene ID>"],
      "state": "no_concrete_upper_lead",
      "retrieval_state": "completed",
      "queries": ["<actual search query>"],
      "checked_species_scopes": ["focal_species", "other_plants"],
      "screened_source_ids": [],
      "reason": "<what was actually checked; not a biological negative>"
    }
  ],
  "gene_assessments": [
    {
      "gene_id": "<interval gene ID>",
      "unit_id": "<frozen unit ID>",
      "provisional_tier": "T4",
      "tier_reason": "<specific function-linked hypothesis>",
      "screen_ids": ["S1"]
    }
  ]
}
~~~

1個のscreenは、同一参照gene/speciesと明示的に対応する**複数copy**へ共有できる。全gene台帳にないID、別reference/speciesのscreen借用、source/summary改変、検索未実施なのに `no_concrete_upper_lead` を使う入力は拒否する。暫定T4は元のclaimに紐づく実際のfunction/pathway材料、T5はannotation/GO/domain材料を必要とする。元の判断に具体的上位Tier手掛かりがある、他に未screenの明示参照がある、既存bundleの`literature_link_candidates`中の当該参照文献ID/URLを`screened_source_ids`で実際に確認した証跡がない、`literature_query_status`に未解決の取得状態がある、発生/GO等の重要な材料が未読である、解釈不足・検索失敗が残る場合は、終了対象にしない。全memberが該当しないmixed unitもunit全体を終了しない。

このplanは**研究作業の配分**だけであり、`agentic_route_b.py prepare/merge`の全prompt rosterを自動削除・短縮しない。元資料・会話/workerの実トークン・時間が減ったかは別の実走比較が必要。保存済みrunの初期判断や9回答を推測で補完しない。旧世代の再生は従来どおり可能。詳細調査を止めた範囲と未解決の問いは次世代で再開可能な状態として保持する。

`checked_species_scopes` は元担当が実際にscreenした対象を記録する。`focal_species` と `other_plants` の両方が記録されるまでは `needs_light_screen` とし、query文字列から調査済みと推測しない。同じ論文IDでも、あるreference/speciesでの確認を別referenceの読了記録に流用しない。geneに結び付いた別screenの具体的上位lead・検索失敗・未実施は、`screen_ids` から省略されていても保持し、`unlinked_screen_ids` と元screen台帳から再開できる。

### 文脈の引継ぎと出力の規律

長い会話を1本のスレッドとして続けると、累積履歴の再送が新規材料を圧迫し、途中省略（truncated）された資料を「読んだ」ものとして扱う危険が生じる。#1004のa101診断（Route B、実行役はCodex gpt-6.1-sol）では、1本の会話で128 request、input 23.55Mのうち約81%が累積履歴の再送、新規材料は約2.4%だった。約40k字のtool出力が以後最大115回再送され、14件のtool出力は `…tokens truncated…` で中間省略されていた。研究をしない形式化だけの2ターンが研究履歴を丸ごと抱えてinputの13.7%を占めた。削減効果は未実測であり、#1004 step 3の同条件比較で判断する。

Route Bの保存回答の一括機械処理は [採用範囲と手順](ROUTE-B-SESSION-ECONOMY.ja.md) に従う。#1020/#1023は全工程で平均1.46%増となりtoken削減を確認できなかったが、理由配送と利便性を理由に採用継続する。通常batch20/originalを維持し、追加の小分割・reset/compaction・起動profile変更・新resume-planは採用しない。旧5案の記録 に失敗した組合せと未測定範囲を残す。組み込みの検索・取得と調査量は節約の対象にしない。

段階の境目（初期判断checkpoint、追加調査、正式回答）では、判断・使用材料pointer・hash・未解決の問い・予算消化を既存artifact（`initial_assessments.json`、discovery、論文表など）へ保存し、会話履歴にしかない判断を残さない。既存の明示batchや再開時はそのartifactを使うが、節約のための追加の強制resetは行わない。「原担当」の同一性は会話の継続ではなく、lane・台帳・provenance（path/SHA-256）で定義する。長い出力（目安4,000字超）はartifactへ保存し、path・SHA-256・件数・行数だけを示す。読むときは連続した行範囲で1回あたり目安8,000字以内とし、読んだ範囲を記録する。表示が途中省略された部分は未読であり、範囲を指定して読み直すまで読了・no-claim・判断済みにしない。これらはローカル/クラウドやproviderに依存しない規則で、同じartifactから別環境が再開・比較できる。

追加調査のtask（`phase=followup`）も、`research/manifest.json`・元台帳・plan・discoveryのpath/SHA-256から新しい文脈で再開してよい。会話の続きを同じ担当の条件にしない。Light契約で原科学回答が保存済みなら、`light-postprocess`等の既存Python処理へ渡し、転記だけのために新LLMを起動しない。既存の`--research-ledger`付き固定回答worker経路は互換として残すが、利用するときも台帳を明示し、節約目的の新規文脈を強制しない。未記録の科学判断が必要な場合は原担当へ戻し、fresh/continuedは親が既存契約に従って選ぶ。

~~~bash
# prepare時に研究台帳のSHA-256を検証し、manifestへpath/SHAを記録する。
LEDGER_SHA=$(sha256sum "$run/b-initial-assessments.json" | cut -d' ' -f1)
python scripts/agentic_route_b.py prepare --prompts "$PROMPTS" --out "$run/b-transcribe-view" \
  --research-ledger "$run/b-initial-assessments.json" --research-ledger-sha256 "$LEDGER_SHA"
~~~

`--research-ledger` を指定した時だけ、検証済みの台帳bytesをprepare出力内の `research-ledger/` へ複製し、各 `worker_task.txt` へその相対pathと転記規則を追記する。workerのread rootはprepare出力のまま変わらず、task文面はホストの絶対pathに依存しない。台帳にある原判断・九回答・source IDs・member対応だけを転記し、新しい研究・検索・判断を追加しない。判断のないunitは理由「台帳に判断なし」で `not_evaluated` とし、陰性やno-claimにしない。台帳は表示資料・read_receiptの代用にならず、資料の全読み要件は維持する。未指定時の `worker_task.txt` と `prepare` 出力bytesは従来と同一で、既存凍結物のhash互換を保つ。SHA不一致・path/SHAの片方だけの指定は出力生成前に拒否する。

親のreceiptには、その段階を新しい文脈で開始したか会話を継続したかを `context_mode=fresh|continued` として記録する。`fresh` でもlane・台帳・provenanceの対応は同じrunのものを使う。

#### 研究段の分割実行

初期研究は、opt-inの `--research-batch-units` で複数のfresh contextへ分割できる。`prepare --research-bundle` に `--research-batch-units 40` を渡すと、requested unitを共有reference endpoint（`reference_material_cards` の同一 `reference_gene`、なければ `reference_gene_id`）でまとめ、`batch_units` 以下のbatchへ詰める。単独で上限を超えるgroupだけを分割し、その分割を索引へ記録する。出力は `research/batches/<bid>/task.md`・`manifest.json`・batch用の索引/初期資料で、グローバルの `research/task.md` とページは従来どおり残る。`--research-batch-units 0`（既定）の出力bytesは変えない。

~~~bash
python scripts/agentic_route_b.py prepare --prompts "$PROMPTS" --out "$run/b-research-view" \
  --research-bundle "$BUNDLE" --research-bundle-sha256 "$BUNDLE_SHA" --research-batch-units 40
python scripts/agentic_route_b.py merge-initial --bundle "$BUNDLE" --bundle-sha256 "$BUNDLE_SHA" \
  --batches "$run/b-research-view/research/batches/index.json" --batches-sha256 "$INDEX_SHA" \
  --out "$run/b-initial-merge"
~~~

batchごとに `research/batches/<bid>/task.md` から新しい文脈のセッションを1つ起動する。並行実行してよく、このflowはexecutorに依存しない。各sessionは担当batchのgeneだけを判断し、原台帳を `research/batches/<bid>/initial_assessments.json` へ既存と同じschemaで書き、`assessment_id`/`screen_id` は `<bid>:` で始める。親は `merge-initial` で台帳をbatch順に連結し（欠落batchは `missing` として記録し、判断を作らない）、以後は既存の `summarize-initial` → `plan-light-followup` → `prepare --followup-plan` → 台帳を明示した既存の `--research-ledger` 経路へ進む。fresh/continuedは親が既存契約に従って選ぶ。結合は原記録の連結だけで、判断の作成・変更・推測をしない。`merge-initial` は各batch台帳の記録を、後段（`summarize-initial`・`plan-light-followup`）と同じ語彙・binding・pointerで1件ずつ検査する。通らない記録と、外したscreenを参照する記録だけを理由付きで `dropped_records` に残して結合から外し、同じbatchの他の記録は使う（1件の誤りでbatch全体の判断を捨てない）。bundleに結び付かない台帳や表の形が壊れた台帳だけはbatchごと `invalid` とする。記録がすべて外れたgeneは未評価のまま残り、親は必要なbatchだけを再実行して結合し直し、新規文脈で始めるか継続するかは既存契約に従って選ぶ。`assessed_fields` は列挙した資料のいずれかにある欄であればよい。TASKには台帳の語彙（暫定Tier・screen state・species scope・欄名・必須欄）を明記している。これは契約の明確化であり、特定の実行役向けの調整ではない（OpenCode deepseek-v4.1-flashを実行役にした試行で、batchごとに書式が分かれた例があった）。

文献探索の停止条件は分割の有無で変えない。Tierが付き得るgeneでは上限約20報までできるだけ論文を集め（原著に限らずレビューも可）、集めた中から高Tier順に有用な最大5報を選ぶ。もともと文献が少ない場合はそれでよく、尽きたら「文献が尽きた」と記録するだけで詳しい理由は不要。件数補充や全論文精読の義務にはしない。

参考：a101（1000 gene、559 requested unit）で、OpenCode deepseek-v4.1-flashを実行役にした1例・1回の比較では、1本の会話はinput約16M・文献screen 9件、≤40 unitのbatch分割はinput約111M・62 screen・193論文でgeneごとの判断が行われ、真のgeneはどちらもT1だった。この数字と振る舞いは実行役・ハーネスに依存する観測であり、flowの性質や既定値の根拠ではない。

### 軽量Route B

#1013の新形式は `light-prepare --output-contract science-v2` で明示的に選ぶ。既定の `legacy-v1` と旧freezeは維持する。新形式でも論文ID・正確なタイトル、報告遺伝子名・ID（namespace付き、未取得は明示unknown）・種は原研究LLMが記録し、取得した書誌で上書きしない。Tier番号とgeneごとの書誌管理文字列は不要で、元の9回答・理由・報告遺伝子との対応と既存の投影検査から機械導出する。

`results.jsonl` に1個のJSON objectを保存し、`contract: agentic-route-b-science-2.0.0`、`paper_observations`、`gene_decisions`を持たせる。観察は `observation_id` で明示参照し、論文が同じでも報告gene・形質・実験・反証が違えば別観察を残す。copyごとの `paper_assessments` に対応・9回答・理由を記録する。`basis=paper_evidence` は独立DB根拠を推定しない。原研究者が明示する `source_bound_function` / `annotation` は、正しい凍結DB出典がある場合だけ独立の機能・注釈判断として論文と併存できる。

`--packing reference-affinity` は同じ正式参照・種・資料snapshotを持つunitを近くへ配置するopt-inで、既定は `original`。batch上限・全roster・copy差分は保持し、Tierや採否は共有しない。convertのCLI返却 `receipt_sha256` を保存し、最終case reportでは `agentic_route_b_light.verify_light_completion(receipt_path, receipt_sha256)` で入力・回答・coverageのhashを検証した上で、返された `state`・`member_counts`・`partial_units` を必ず表示する。unitのtransport `completed` を全memberの科学調査完了と解釈しない。正常claimと同じunit内の未評価・不正回答・取得失敗は別々に保持する。新形式の実トークン削減率は、同条件の実モデル比較を行うまで未測定とする。

上の研究台帳＋正式転記の代わりに、opt-inの軽量経路を使える。従来の `legacy-v1` では元研究者は `material.md` の各geneについてTier・理由・採用論文の9回答を書き、機械は幻覚検査（batch内gene ID、frozen reference endpoint、Tier語彙、論文ID解決とtitle一致・abstract取得）を行う。既存の台帳経路は互換のため残し、この経路は指定したときだけ動く。

~~~bash
# 1. requested unitを最大N unitのbatchへ詰め、material/task/rosterを書く（model呼び出しなし）
python scripts/agentic_route_b.py light-prepare --bundle "$BUNDLE" --bundle-sha256 "$BUNDLE_SHA" \
  --out "$run/b-light" --batch-units 20 [--units-file units.txt]
# 2. 各batchで results.jsonl を書く（この工程だけが研究者の仕事）
# 3. 幻覚検査・discovery・llm-run用agent file・source-verify用cacheへ変換
python scripts/agentic_route_b.py light-convert --bundle "$BUNDLE" --bundle-sha256 "$BUNDLE_SHA" \
  --light-index "$run/b-light/light/index.json" --light-index-sha256 "$INDEX_SHA" \
  --out "$run/b-light-convert" --cache "$CACHE" [--offline]
# 4. 以降は既存canonical経路（receiptのnext_commandsに実コマンドを記録）
#    route-b-semantic-bundle --extend-bundle -> llm-run --provider codex-agent --agent-file ->
#    route-b-semantic-import -> route-b-semantic-source-verify --cache -> transfer-verify -> candidates
~~~

`light-prepare` はpost-B0 bundleのSHAを検査し、requested unitをbundle順のまま `light/<bid>/material.md`・`material.full.json`・`task.md`・`roster.json` へ書き、`light/index.json` にbundle path/SHA、batch_units、各fileのSHAを記録する。`--units-file` は小規模試行用にunit IDを限定する。taskは担当material・完全資料・resultsの絶対pathとgene一覧を明記する。見つからなければそのpathの失敗を報告し、別batchを探索・代用しない。配送先を変える場合はそこでprepareし直し、新しいhashを使う。

`material.md` はreference cardの機能・GO・組織・発生段階・domain等の生物学的欄を文字数で切らずに表示する。表示する生物学的欄・参照・出典が完全一致したcardだけをbatch内で一度表示し、各geneから参照する。異なるGO・組織・発生段階などは別cardとし、gene固有contextと `reference_gene` / `orthology_reference_gene` も保持する。元の `family_id` / `family_label` / `reference_family_id` / `evidence_completeness` / `notes` は各geneに、`unit.identity` はunit IDに結び付けて表示し、注釈の矛盾や不明点をそのまま渡す。取得時刻等の管理欄と生の投影・source情報を含む全unit原本は、hash検査付き `material.full.json` に保持する。必要なunit/fieldだけを開き、巨大JSONを会話へ再表示しない。長い資料は分割して読み、ツールの返却省略を読了としない。旧sidecarなしの凍結世代も再生できる。`roster.json` のendpointはsemantic importがgroundingに使うのと同じ欄（gene reference_gene、orthology_reference_gene、family_index_subject_endpoints、sequence_similarity_hits、reference_material_cards / literature_link_candidates の reference_gene）から作る。

軽量回答の各観察には、LLMが `paper_id`・`title`、`reported_gene`・`reported_gene_id`・`reported_species`、所見・対応・9回答を書く。原報のIDが不明なら明示的unknownとし、proxy/reference IDを代入しない。書誌情報だけは、`results.jsonl` の先頭に `{"bibliography":{"p1":{"source_id":"PMID:...","source_title":"LLMが原資料から写した題名"}}}` を一度書き、各観察で `paper_ref: "p1"` としてよい。gene・形質・実験・組織/時期・所見・対応・9回答は共有しない。同じレビュー論文の異なる観察を別々に保持し、同じ論文IDというだけで2件目を捨てない。完全一致の観察だけは重複を除く。5報は異なる論文の数であり観察行の上限ではない。超過した原記録も保持し、既存canonical handoffの選択へ渡す。書誌参照の不明・不一致は記録し、別geneの正常回答は残す。

`light-convert` の機械検査は上記の幻覚検査だけで、意味判断を追加しない。`interval_gene_id` がbatch外、`reference_gene` がそのgeneのfrozen endpoint外、Tier語彙外、解決しない・題名不一致・abstractなしの論文は、gene/paper単位で `light_check.json` に記録して落とす。落とした論文でTierを支えられないgeneはclaimを残さず、そのunitは `not_evaluated` として記録する（no-claimsにしない）。`not_assessed`・未回答geneを含みclaimのないunit、light run外のunitも同じく `not_evaluated` とする。この表現はcanonical `llm-run --agent-file` の `status: not_evaluated` 回答（理由付き）として渡し、semantic importではclaimなし・`no_claims` ではないstatus行になる。

機械検査は研究者の判断を作らず、直さず、格上げしない。論文は、機械が `Tier` を計算し直すためではなく、実在・題名一致・checklist完全性（9キーすべてがyes/no/unknown。欠けた・空の・語彙外の回答を `unknown` で埋めない）・新規sourceならdiscoveryを組めること（abstract、https URL、retrieved_at、frozen reference species、reported species）を満たすときだけ残す。落とした論文は `light_check.json` の `paper_issues` に記録する。

採用論文を落とした結果、残った論文のchecklistが示すTier（`checklist_tier`）が宣言Tierと異なるgeneは、claimを残した上で `light_check.json` の `tier_comparison` に `declared_tier`・`computed_tier`・`dropped_papers` を記録する（宣言Tierを主張する文言はclaimに書かない）。T1–T3で残る論文が無ければclaimを残さない。論文なしのT4/T5は、そのunitにfrozenされたdatabase accessionとreference gene単位のprovenanceがあるときだけ、`evidence_types`・`reported_function`/`reported_annotation` で出す。元研究者が報告機能と具体的な形質への橋渡しを記録したT4は、論文ではなくDBを引用しているという輸送上の理由でT6へ落とさない。注釈だけはT5、identityだけはT6のまま。A/B共通のimporterも明示的な空の `paper_reports_json` によってDB機能・注釈を消さない。機能の創作、出典なしの昇格、verified supportへの昇格は行わず、既存のsource/identity/member検査を維持する。accessionが無ければclaimを残さず記録する。

論文IDは `PMID:n` / `PMCID:PMCn` / `DOI:x` の正規形で出し、既にfrozenなIDはそのエイリアス（PMID/PMCID/DOI）すべてで照合してfrozen側の表記を使う。Europe PMC問い合わせは429/5xx/通信エラーで最大3回、指数backoffで再試行し、取得失敗（`provider_failure`）は未発見（`not_resolved`）と別に数える。Europe PMCでDOIを解決できない場合はCrossrefの同一DOI recordへfallbackする。返ったrecordが要求したID自身を持たなければ採用しない。題名比較はHTMLタグ除去・実体参照展開後に行い、保存する題名も除去後のものとする。Crossrefにabstractがない場合も創作せず、取得した出版社資料を`abstract_provenance`付きcacheで渡すか、未取得のまま記録する。出版社403を直接取得成功としない。取得元・URL・hash等の`metadata_provenance`と`abstract_provenance`は両cacheに保持する。`--offline`では追加取得しない。

取得したmetadataは `metadata_cache.json` に残り、再実行では `--cache <前回の出力>/metadata_cache.json` で再利用できる。同一geneの重複行は、同一内容なら1行にまとめて記録し、内容が食い違えばそのgeneを `not_assessed` として記録する。receiptにはbatchごとの `results.jsonl` のSHA-256も残す。`material.md` などハッシュするfileはWindowsでもLF固定のbyteで書く。

canonical B sessionは、出力に隣接する既存のimport statusとmanifestがある場合、そのstatusを凍結済み回答・bundleから再導出して照合し、sourceのhashとともにsnapshotへ渡す。missing/failed/skippedは完了へ変換しない。statusのない旧sessionは未確認のまま再生でき、supportの有無と実行coverageを混同しない。

新しい論文は既存schemaの `discoveries.json` にまとめ、`route-b-semantic-bundle --extend-bundle` で新generationへ凍結してからimportする。既存 `stable_source_ids` にある論文はdiscoveryにしない。`agent_responses.json` は全requested prompt分を作り、T1–T5のgeneだけをclaimにする。T6/Uはclaimを作らない。`metadata_cache.json` と、source-verifyが読む形式の `source_verify_cache.json` を同じmetadataから書き、`--offline` ではcacheにある論文だけを使う（未取得は `not_attempted`）。

研究台帳・共有判断・詳細screen・`summarize-initial`/`plan-light-followup`/`--research-ledger` 転記の経路は削除せず、指定した場合だけ軽量経路を使う。参考：q001（Type2）とd001（Type3）の真因周辺約100 geneで、OpenCode deepseek-v4.1-flashを実行役にした各1回の試行では、軽量経路はextend→import→source-verify→transfer→freezeまで下流のコード変更なしで通り、真のgeneはいずれもT2で `adopted_first_stage`、論文IDと題名はすべて検査を通った。同じ実行役でのinput tokensはbatched ledger経路の約1/6だった。数字・Tierは実行役・ハーネスに依存する1例・1回の観測で、flowの性質、一般性能や保持率の保証ではない。flow自体はハーネス・モデルを問わない。

### B担当への資料配送

Bの`prepare`はcanonical `b*/prompts.jsonl`を保持し、追加で共有資料付きの読取りビューを作る。親は`prepare_manifest.json`の`reading_view.shards`からassignmentを配る。`reading/b*/assignment.json`には全assigned prompt/unit ID、unit本体、各geneの対応、参照gene・species・source handleがある。`materials`にあるpathはprepare出力directoryからの相対pathで、`shared-material/*.json`を必要な範囲だけ開く。system、phenotypeを含むuser前置き、参照機能・出典の完全なrecordは内容と種類が同一のときだけ1ファイルへ共有される。unit/geneごとのbindingと同一行の反復は保持され、同じfamily名でも機能・citation等が違えば別資料になる。既読の資料は再利用し、各unitのclaim判断と回答は別々に行う。

`prepare --prefer-shared-material`は、同一内容の参照資料を共有するunitを同じreading assignmentへ寄せる明示opt-in。既定はOFFで、canonical prompt shardsと各unitの科学promptは変えない。元の順序と配置案を比較し、総読取りbytesが減り、shard数・最大assignment読取りbytes・oversized数のいずれも増えない場合だけ配置案を使う。それ以外は元の順序へ戻し、manifestの`packing_strategy=canonical_order_no_safe_affinity_gain`に記録する。別資料の共有を崩してかえって遅くする配置は採用しない。これは配送量の検査であり、token・実時間・候補保持性能の改善を保証しない。

このopt-inでは`reading_view.shards`のID集合がworkerの担当範囲になる。canonical `shards`の連続番号へ戻して担当を分けると、共有資料の再配置効果が失われる。既存のbounded readerを利用する場合は、`workflow/bounded_b_reader.py init --reading-shard reading/bNNNN/assignment.json`でmanifest掲載のassignmentを指定し、通常の`--view-root`・`--state`・`--expected-source-sha256`・実executor/modelも渡す。`--start/--stop`とは併用しない。元のsource indexを保存したままassignment順に読み、同じassignmentのpartial checkpointから未完IDだけを再開できる。従来のrange指定と既存cursorはそのまま再生する。完全promptへ再展開するnative packet経路には共有資料の読取り削減を自動適用しないため、その経路の効果は別に計測する。

新しいprepareは各reading shardへ`worker_task.txt`を添え、path/SHA/bytesを`worker_task`欄に記録する（CORE_FIX、#757/#789）。親はこのtaskを起動指示として渡し、prepare出力directoryをworking directory/read rootにする。assignment配下だけを許可して共有資料を読めなくする起動文は使わない。隔離したjobへ配送する場合は、そのassignmentと参照する共有資料の相対配置・bytesを保ち、起動文内のrootをjob rootへ対応させる。隣のassignmentやA/truthを担当へ渡さない。taskはlaunch指示であり科学promptや新たな支持根拠ではない。`--max-material-bytes`は従来どおりassignment＋参照資料の読取り量を数え、launch指示のbytesは別欄に記録する。既存のtaskなし世代も再生できる。

行単位readの長行切捨てを避けるため、新規prepareは`display/`配下へassignmentと全参照shared-materialの表示用JSONも生成する。元JSON・material ID・prompt/hash・packing予算は変えない。`reading_view.shards[].display_view`に元/表示path・SHA・bytes・行数・最大行長を記録し、生成・復元時に元canonical bytesへの完全一致を検査する。表示は最大1000文字/行で、worker_taskも同じ行長内に収める。担当はtaskに従い`display/reading/b*/assignment.json`を連続する行範囲で読み、その`value.materials[].file`ごとに`display/<元path>`へ進む。

表示envelopeの`value`は通常のフィールド構造を保ち、`chunked_string_paths`に載った位置だけが長文字列の`__cqtl_display_string_chunks__`表示になる。配列を順に文字追加なしで結合する。元データの同名keyを自動的にmarkerと解釈しない。長いkeyの置換は`renamed_keys`に元keyのchunksと位置を保持する。表示注記は科学回答schemaではなく、元system/user-prefixの全内容と元の回答形式を使う。表示分の保存bytesは増えるが資料を要約・削除せず、共有資料は重複排除後に一度だけ生成する。旧sidecarなし世代はそのまま復元できる。表示生成や復元成功を実際の読了・研究完了とは扱わない。

Aの起動指示も既存`mechanism-import`の形式へ揃える。回答は余分な括弧・fence・説明文のないJSONとし、`source_ids`には認識可能なPMID/PMCID/DOI、または型付きdatabase accessionを入れる。単なるURLを入れる指示は使わない。識別できない出典はresearch log/uncertaintyへ残し、IDを推測しない。親はraw/native exportを保存して、既存`packet_to_canonical.import_a_payload`を新しい検証出力へ実行し、rejected rowsも確認する。import成功はsourceの意味検証やsupport確認ではない。形式エラー時はrawを自動修復・上書きせず、実行が許可されたrunの予算内で当該回答のみ最大1回再入力させ、別attemptと変更履歴を残す。直らなければ理由付き未評価checkpointで止める。Bも既存`merge`で形式・ID/hashを検証し、healthy checkpointを保持したまま失敗unitだけを同じfrozen promptsへ再入力する。準備・検証だけではモデルを呼ばず、性能ベンチのPAUSED状態も解除しない。

読取りassignmentは`--max-prompts`と、assignment本文＋そのassignmentが参照する全unique資料のUTF-8 bytesで分割する。`--max-material-bytes`の既定値は262144で、科学的な選抜基準ではない。単一unitが超える場合は専用assignmentに`reading_state=oversized_full_material`、超過bytes、完全な資料参照を保存し、他のassignmentは継続する。親はこのunitを省略せず資料を分けて読む。canonical外のuser書式は推測せず全文inlineで保持する。

生成時に全field・user文字列をcanonical promptへ完全復元し、ID/順序・case・source SHA・各artifact hashを照合する。再読込みの機械検査には`from scripts.agentic_route_b_material import restore_reading_shard`を使い、`restore_reading_shard("<prepare出力directory>", "reading/b0000/assignment.json")`を呼べる。`exact_roundtrip_verified=true`は輸送の完全性だけを表し、資料の科学的検証やLLM実行を意味しない。既存`merge --prompts`には引き続き元のfrozen `$PROMPTS`を渡す。prompt数、Bのrequested条件、全interval roster、未依頼/未回答の状態、候補選抜や支持属性は変更しない。共有ビューは全geneを索引から検討するための配送改善で、familyの科学的一括判断や詳細調査の省略は実装していない。

`config.json`には `case_id`, `focal_species`, `phenotype`, `interval`, `reference_fasta`, `annotation_gff3`, `a_payload`, `a_proposer_id`, `a_model`, `b_proposer_id`, `b_model` が必要。pathはconfigのdirectoryを基準に解決される。`interval`はFASTA/GFFと同じassemblyの`seqid:start-end`を指定し、seqid対応が必要な場合だけ既存`seqid_map`を使う。既存の `record_cache`, `record_cache_manifest`, `source_cache`, `reference_corpus`, `family_index`, `domains`, `orthology`, `annotations`, `literature`, `protein_db` は存在する場合に指定する。`protein_db`はFASTAファイルではなくformatted BLAST DB prefix。普通のFASTA/GFF入力では`reference_corpus`は必須ではない。各レビューの `reviews.json` は既存入口の形式で、`a`/`b`ごとに `{ "response": "path-relative-to-reviews-json", "reviewer_id": "...", "reviewer_model": "..." }`。実行済みreviewとmachine未実行ではID/モデルを区別する（後者は `no-source-machine-a|b` と `not_invoked`）。`scripts/agentic_route_b.py merge` の `not_evaluated` はそのunitの未処理であり「claimなし」ではない。partial merge outputを同じfrozen promptsに再入力し、retry shardのraw responseとともに新attemptへ渡すことで、生成済みcheckpointを再利用できる。未応答promptだけでなく一度`failed/not_evaluated`になったpromptも別attemptでcompletedへ更新できるが、同一attempt内duplicateと別case/unknown IDは拒否する。healthy responseへの後発failed/not_evaluatedは無視して既存の正常回答を保持し、downgradeしない。mergeが返す execution manifest は `review --b-execution-manifest` で同じfrozen prompt/response hash、prompt→unit、parse status、raw response hash、executor/modelとretry historyまで検証される。common handoffもB/C state、件数、case、source artifact hashを再検査する。manifestを省略した場合は `not_recorded` として残す。

たとえばBのpartial checkpointを保存した後は、同じ`$PROMPTS`に対して次の新attemptを作る。checkpointに入ったhealthy rawを再度`--response`へ重複投入しない。retry回答にも元prompt ID/unit IDと実executor/modelを保持する。

材料追加などでprompt generationが変わる場合は、既存checkpointと原rawを読み取り専用で照合できる。

```bash
python scripts/agentic_route_b.py reuse-plan --old-prompts "$OLD_PROMPTS" --new-prompts "$NEW_PROMPTS" --old-response "$OLD_B_RESPONSE" --import-manifest "$OLD_B_IMPORT_MANIFEST" --out "$run/b-reuse-plan.json"
```

この計画は全field（system/user、case、unit/member、source、追加metadataを含む）のdigestで`same/changed/new/removed`を記録し、旧execution manifestに束ねられた原rawのpath/SHA、実executor/model、raw本文とcheckpointの一致を確認する。`same`でも未完・invalid・既知import-rejectedは再利用適格にしない。import manifestのprompt/response SHAが無い旧形式は、raw claimと一致したrejectの除外にだけ使い、健全な回答の適格性は`import_audit_input_not_bound`とする。import監査を省略した場合も`import_audit_not_recorded`で未確定にする。必要なら旧bundle/prompt/responseを現行importerへ別出力で再投入し、元の実行由来を保ったhash-bound監査を作る。原本の書換え、scope補正、回答の自動merge、モデル実行はしない。`new_model_execution_count=0`であり、再利用を新実行数や科学的支持へ読み替えない。system説明だけを変えたpromptも`changed`なので、旧回答の全文同一再利用にはならない。

名前・種の複合参照は、対応が明示された範囲でatomic pairに分解し、全てのprovider一致を保持する。`reference-gene-cache`の一部queryが失敗しても健全なrecordは使える。`cache_complete=false`とquery ledgerの未評価は別に報告し、同じmechanismsを`--resume-manifest <旧manifest>`と新しい`--out`へ渡して失敗queryだけ再試行する。`prepare`はrecord・mechanisms・query ledgerをhash検証し、移設後もledgerと元manifestの由来を保存する。

独立source reviewで `decision=supported` を返す際、`reported_organism` はauthoritative recordの学名またはcommon nameと文字通り一致するtokenを使い、同じliteral `organism_quote` にもそのtokenが含まれることを確認する。原文がcommon nameだけなら、引用に無い学名へ置き換えて支持扱いにしない。凍結済みreview promptのinstructionsやquoteを変更して通さない。

B bundle schema 1.9は、semantic call前に`source_material_class`と決定的な入力量見積もりを記録する。裸のID・hit・記事linkだけではcallを要求せず、全unit/memberを未評価台帳に保持する。見積もりは実測tokenや課金額ではない。表示を短くする場合も、具体的bridgeの材料や切れた文章を全文参照へ戻し、全文を読めなければその範囲を未評価として記録する。JSONが構文またはclaims envelope/member形式に違反した応答はunit単位の`not_evaluated`とし、rawと失敗理由を保存する。unknown/duplicate ID、case/hashの不整合、正しいraw JSONと`model_json`の矛盾は引き続き拒否する。

`agentic_route_b.py merge`の`response_state_counts`は元回答から生成し、通常のreview/finishと最終報告へ渡す。完了claimあり／空claim、担当者の明示未完／明示失敗、形式違反、一部claim隔離、未回答、旧形式で分類不能を分ける。理由文から資料不足・未読・provider障害・生物学的陰性を推測しない。新しい内訳はconsumerで回答と再照合し、内訳のない旧manifestは元回答から再計算して旧artifactを保持する。これは回答の実行状態であり、読了・文献探索・科学判断の正しさや完了を証明するものではない。

新しい `agentic_route_b.py prepare` は表示ファイルの末尾に `read_receipt` を置き、`reading_view.shards[].read_receipts_by_prompt` に各unitが読むべき表示ファイルのreceiptを封印する。担当は完了回答ごとに `read_receipts`(assignmentとそのunitが参照する全表示ファイルのreceipt)を返す。`merge --reading-manifest <prepare_manifest.json>` は照合し、receipt不足の空claim回答を `not_evaluated`(`read_receipt_incomplete`)へ戻し、claim付き回答は保持して `read_receipt_state=incomplete` と記録する。manifestの`read_receipts`に件数が残る。receiptはファイルを末尾まで開いた配送証跡で、理解や判断の正しさの証明ではない。grep・抽出スクリプトによる拾い読みの空claimを「調査済みclaimなし」と数えないための検査であり、生物学的陰性ではない。旧prepare世代(receiptなし)は `--reading-manifest` を渡さず従来どおり再生する。`preflight` も同じ `--executor/--model/--reading-manifest` を受け取り、merge前にreceipt不足を報告する。

```bash
python scripts/agentic_route_b.py merge --prompts "$PROMPTS" --response "$run/b-responses.json" --response "$B_RETRY_RAW" --out "$run/b-attempt-02.json" --executor "$B_EXECUTOR" --model "$B_MODEL" --reading-manifest "$run/b-shards/prepare_manifest.json"
python "$h" review --prepared "$PREPARED_MANIFEST" --sha256 "$PREPARED_SHA" --responses "$run/b-attempt-02.json" --b-execution-manifest "$run/b-attempt-02_manifest.json" --out-dir "$run/02-reviewed-retry"
```

通常Cの順番は `独立したハーネス標準Web検索/取得 → proposal TSV → route-c-direct-freeze → claim-binding prepare/seal（独立agent） → route-c-direct-project → route-c-direct-members`。recordが報告記号だけを持ち配列accessionが無い場合は、bound manifestの前後どちらでもよいので `route-c-provider-bait → blast-panel-freeze → blast-panel-search` でbaitとhitsを作り、`route-c-direct-project --eligible-hits <hits> --identity-mapping <provider_identity_edges.tsv>` に渡す。provider edgeが無いBLAST hitは、問い合わせ配列が報告記号に結び付かないため、記号でbindingしたrecordでもWetのreview集合に入らない。`route-c-direct-search`は専用検索の任意fallbackであり、正常なC調査の必須CLIではない。`project` は `prepare` が返す `b0.interval_gene_atlas` を使い、`members` が返すmembers TSVとmanifestの**元パス**を `finish` に渡す。詳しいフラグと提案17列は `docs/ROUTE-C-DIRECT-RESCUE.ja.md` と `candidate-qtl route-c-direct-* --help` を参照。Cが検索済みで空なら有効な空artifactを通す。provider障害時は `{"case_id":"<case_id>","state":"provider_failure","reason":"..."}`、その他の未評価は `case_id` を含むJSONに `state=not_evaluated` と理由を保存して `finish --route-c-status <path>` に渡す。statusはhash付きartifactとしてhandoff内にsealされ、case/controlled stateをbindする。C未供給でstatusも無い場合、reportは `not_evaluated_no_input` と表示する。A/BのWetが生成されても三経路の完了とは呼ばない。

C memberの生成件数とWetへ渡ったreview件数は別に点検する。review保護には凍結source identity bindingが必要で、binding欠損のmemberは保護対象にならない。既存 `python -m candidate_qtl.route_c_claim_binding_producer prepare|seal` を使い、独立agentが原sourceの種・identityと投影ledgerのsource identity対応を確認したreview TSVを新しいmanifestへsealできる。agent reviewは `reviewed_by=agent:...`、`curation_method=independent_agent_review` とし、確認していない軸は未評価のままにする。全軸を確認し、引用・逐語引用・凍結identity・取得範囲の機械検査を通ったbindingは `verified_complete` になり、projectionが一致すればsupportへ昇格する（人間の別確認は不要）。軸が残る `complete=False` / `review_required` のままでも、既存source identity・projection条件を満たすreview候補はWetへ渡せる。品種・系統はgene identityではないので、review候補では `source_accession_or_cultivar_state=not_evaluated` を許し、`contradicted` だけを除外する。報告アレルは従来どおり検証済みか該当なしを要する。supportは全軸の検証を要する。agentが `curated_human` を名乗らない。旧freeze/manifestを変更せず、sealが返す新manifestでproject/members/finishを再実行する。

検索障害の再開では、失敗cassetteもimmutableな実行証跡として元stageに残す。現在の検索は一致するcassetteをstatusに関係なくreplayするため、失敗を含むcacheをそのまま`--cache`へ渡しても再取得にはならない。親Agentがcassette JSONの`retrieval_status`を確認し、`retrieved`または`no_hit`のものだけを**bytesを変えず**別のhealthy cache directoryへコピーする。同じ検索条件で新しいoutを指定し、外部取得を許可した検索を行う。健全なqueryはreplay、cacheに含めなかった障害queryは再取得される。旧cassetteを上書き・削除したり、障害を該当なしへ置換したりしない。

```bash
candidate-qtl route-c-direct-search --case-id "$CASE_ID" --phenotype "$PHENOTYPE" --focal-species "$FOCAL_SPECIES" --species-scope "$SPECIES_SCOPE" --limit 10 --cache "$run/c-healthy-cache" --allow-external-api --out-dir "$run/c-search-attempt-02"
```

同じ条件には`--limit`とspecies scopeも含める。検索のstatus/provenanceとraw cassetteを各attemptに保存し、新しいC proposal/freeze/project/membersを作る。A/Bのsource cache追加は元prepared stageを変更せず、同じproposalとprepared anchorを `review --source-cache` へ渡して新しいreview stageから再開する。reference resolution自体を更新する必要がある場合は、別の `prepare` 世代にする。

### 通常v3のWet境界と再試行

新しいreading projectionは`structured_json_catalog_refs_v1`として、`candidate_input`の各行、
上位context、主論文panel、request policyも既存catalogの同じinline/refで共有する。
妥当なobject/list形式の`*_json`だけを読み用に構造化し、invalid JSONやscalarの原文字列は保つ。
`expand_reading_view(view, catalog)`で全reading fieldを展開できる。元のJSON文字列の表記、
全raw audit、科学候補行はcanonical shardにbyte/hash不変で保持し、final consumerは原shardから
同じprojectionを再計算して照合する。旧projection世代は元の方式で読む。
共有panelを一度読んでもcandidate/copyの各bindingや固有差の比較は省略しない。
inline/refの費用比較はunique候補数ではなく各shardへの実配置回数を数える。全shardに渡すread-only上位候補の材料もその回数で共有し、候補ごとのbindingと原shardは保持する。
新しいpage/task hashは配送の検査であり、旧世代の既読記録や科学的な読了へ読み替えない。

新しい`initial_bibliography_v1`の通常A/B handoffは、凍結B bundleから境界候補とread-only上位候補の機能・対応情報を`comparison_context.json`へ保存し、境界requestへpath/SHAを結び付ける。`agentic_boundary.py prepare`はその材料を自動で読み、共有catalogと候補ごとのviewへ渡す。対象外の下位geneは追加せず、同じ参照の資料は共有する。これは既存材料の配送であり、モデル実行・機能確認・支持昇格ではない。Bに材料が無いgeneはその状態を残す。生のhit metadataは元bundleのpath/SHAとsource pointerから必要時に参照できる。

追加の`--comparison-context`を指定した場合も、requestに束ねた既存材料を空欄で上書きしない。同じ記録表現型・scopeを確認して両方を保存する。旧requestは従来の明示contextとその検証を継続する。

元担当は論文のTierを変更せず、コピーごとの機能・aligned regions・利用可能なdomain・別identityや既知の矛盾を比較する。`defer`は今回の実験を後回しにする理由付き判断で、非原因の証明を要しない。比較した理由と再検討に必要な問いを残す。ただし部分一致・支持未確認・材料不足だけで自動保留したり、20/25件へ強制的に合わせたりしない。未完の比較は`not_evaluated`とし、追加のTier reviewerや固定の追加検索を要求しない。

境界がある場合は `03-first-wet/handoff/wet` の境界requestパスとpolicy版を確認する。以下は `upper_evidence_scoped_v3` の手順である。各shard回答には `shard_id`, `case_id`, `source_wet_sha256`, `responses`（`candidate_gene_id`, `decision`, `reason`）を含め、元担当のexecutor ID/モデルを付ける。`agentic_boundary.py` は未回答候補を明示的に `not_evaluated` として全rosterを出す。v3のfinishには同じexecutor ID/modelを渡し、独立review artifactやreviewer引数を指定しない。

```bash
python scripts/agentic_boundary.py prepare --request "$BOUNDARY_REQUEST" --out "$run/boundary-shards" --max-candidates 8
python scripts/agentic_boundary.py merge --prepared "$run/boundary-shards" --response "$PROPOSAL_REPLY" --role proposal --out "$run/proposals.tsv" --executor-id "$EXECUTOR" --executor-model "$MODEL"
python scripts/claim_review_session.py finish-wet-review --handoff "$HANDOFF" --completion-sha256 "$HANDOFF_SHA" --out "$run/04-adopted" --proposals "$run/proposals.tsv" --executor-id "$EXECUTOR" --executor-model "$MODEL"
```

boundary担当へは`manifest.json`と担当shardの`worker_task.md`を渡し、case/source SHA、comparison-context SHA、shard IDと既存回答schemaを確認させます。新しいprepareは`reading_pages`に、完全な`reading_view.json`・共有`comparison_index.json`・`material_catalog.json`のlossless表示ページを出力します。taskの親directoryの`..`がreading rootで、task内の全pathはそのrootを基準に解決します。nested workerの現在directoryに依存させず、共有資料を含むprepare出力全体をread許可範囲にします。共有catalog/indexは一度だけ生成・閲読し、各候補の`candidate_input`と`primary_paper_evidence`、資料binding、read-only上位contextを含む全viewをたどります。`comparison_context`だけでは元論文の要約・選定報告を含む全比較資料にはなりません。

各表示ページはUTF-8で12 KiB・160行以内、各行1000文字以内です。これは観測された約56 KiBの返却上限と行番号等の追加分に余裕を置いた運用値で、全harnessで切断されない保証ではありません。番号順に全ページを読み、実返却最終行・EOFを確認し、途中で切れたら同じページの次の実返却行から続けます。細分化でファイル数・read回数が増えるため、担当viewから読み始め、共有catalog/indexは一度読んだ範囲を再利用し、未取得ページと再開位置を記録します。ページ数は調査件数や新しい判断単位ではありません。ページを余分な文字なしで連結すると表示用envelopeに戻り、元fieldは`value`にあります。Route Bと同じlossless表示形式を再利用し、長文字列は明示されたpathだけchunksを連結します。元データのmarker風keyや長keyも区別して復元します。これは回答schemaではありません。

完全資料を再現する機械検査は`from scripts.agentic_boundary_material import restore_reading_shard`の`restore_reading_shard(prepared_path, manifest, shard_id)`です。source/page/taskのhash・順序・全field復元を検査し、旧ページなし世代も読めます。新prepareの`reading_delivery_contract` marker、またはpage/taskの残存がある場合、ページmetadataの欠落/nullを旧世代として受理しません。canonical shard・科学判断・回答schemaは変えません。ページの存在・復元一致・回答roster完了は、実際の資料受領や理解、科学比較の完了を証明しません。

`prepare --comparison-context "$CONTEXT" --max-material-bytes 262144`で前段資料・判断を由来付きで供給できます。bytesは目安であり、`chunk_required`または共通文脈の超過があれば完全資料を分けて読み、担当geneを省略しません。共通索引は1ファイルで保存されますが、全shardで読み直すと展開量は増えます。実際のtoken/モデル費用は別に測る必要があります。詳細は[canonical contract](CANONICAL-WET-MEMBERSHIP.ja.md)を参照してください。

boundary shardが未回答/失敗なら、元requestとprepared shardsを保存し、そのshardだけ別attemptで再実行する。Bのaggregate checkpointとは異なり、boundaryのmergeは同一candidate IDの重複を拒否する。再mergeにはhealthy shardのraw回答とretry shardのraw回答を渡し、置き換え対象のfailed rawや前attemptのpartial TSVは渡さない。v3では更新されたproposalと同じexecutor ID/modelをfinishへ渡す。独立review回答を新proposal hashに対応させて再取得するのは、独立reviewを要する旧policyの場合である。

再試行の`finish-wet-review`は、**未完の子stageのcompletionではなく、元のfrozen handoffとそのcompletion SHA**をconsumeし、新しいout directoryへ出す。子stageの`review_history_state=incomplete`は後続stageにも累積されるため、未完の子を入力にして同じ境界のretryを完了扱いにしない。新通常policyは未完出力に元checkpointのresume先を保存する。旧policyで完了した境界から次のlower-tier境界へ進む場合は、その完了stageを次のhandoffにする。

```bash
python scripts/agentic_boundary.py merge --prepared "$run/boundary-shards" --response "$BOUNDARY_HEALTHY_RAW" --response "$BOUNDARY_RETRY_RAW" --role proposal --out "$run/proposals-attempt-02.tsv" --executor-id "$EXECUTOR" --executor-model "$MODEL"
python scripts/claim_review_session.py finish-wet-review --handoff "$HANDOFF" --completion-sha256 "$HANDOFF_SHA" --out "$run/04-adopted-retry" --proposals "$run/proposals-attempt-02.tsv" --executor-id "$EXECUTOR" --executor-model "$MODEL"
```

#### 旧v2・legacyの独立review

凍結済みrequestが旧policyなら、その世代のscopeと独立review契約で再生する。上記のv3 finishと混ぜず、proposal作成後に別担当のreviewを取得して次を使う。review回答には `proposal_sha256` と別のreviewer ID/modelも必要で、retry後は新proposalへ結び直す。通常v3へこの手順を追加しない。

```bash
python scripts/agentic_boundary.py merge --prepared "$run/boundary-shards" --response "$REVIEW_REPLY" --role review --proposal "$run/proposals.tsv" --out "$run/independent-review.tsv" --reviewer-id "$REVIEWER" --reviewer-model "$REVIEW_MODEL"
python scripts/claim_review_session.py finish-wet-review --handoff "$HANDOFF" --completion-sha256 "$HANDOFF_SHA" --out "$run/04-adopted-legacy" --proposals "$run/proposals.tsv" --independent-review "$run/independent-review.tsv" --executor-id "$EXECUTOR" --executor-model "$MODEL" --reviewer-id "$REVIEWER" --reviewer-model "$REVIEW_MODEL"
```

`03-first-wet/report.md` は第一段の説明。boundary処理後は `finish-wet-review` の `adopted_wet_candidates` が正本で、第一段レポートを最終結果として流用しない。最終 `wet_review_report.md` にはWet SHA-256/adopted stage、A/B/C execution state、decision reason、deferral、未評価、review history、next boundary/actionが統合される。`review_history_state=incomplete` が残る限り `next_boundary_review_required` より `boundary_review_incomplete_retained` を優先し、要求生成や未応答完了と表示しない。support、review候補、汎用fill、境界未評価を区別し、Wet保持を因果確認と呼ばない。

### 候補一覧と任意のLLMによる整形

`wet_review_report.md` の「候補一覧」は、同じWetの全行を次の6列で表示する。第一段とboundary後の両方に出力し、元の行順・ID・Tierを保持する。保留中の行は一言メモに明示し、採用件数と保留件数はレポート冒頭で確認できる。

| 列 | 表示内容 |
| --- | --- |
| 遺伝子名・ID | 対象遺伝子の安定IDと、Wetに保存された対象遺伝子名があれば併記。名前がない場合はIDのみ。 |
| 位置 | 同じWet行の染色体・開始・終了座標。明示的な参照buildが保存されていれば併記。座標系やbuildを推定しない。 |
| 相同遺伝子（種・ID） | 元claimに保持された参照と種。元論文の報告遺伝子が別にある場合は「報告 → 参照」と区別する。集合化されたgene列とspecies列を順番で対応付けない。 |
| Tier | 正本Wetの `relevance_tier`。 |
| 論文 | 保存済みPMID・PMCID・DOIへのリンク。主なTier根拠を先に置き、他の保持論文も表示する。リンクの作成では追加検索しない。 |
| 一言メモ | 保存済みの判断理由の冒頭と、記録済みのcopyの曖昧さ・保留状態。理由の全文と来歴は後続の詳細表に残す。 |

「未記録」は不在の証拠ではない。欄を埋めるために全区間geneの機能を調べ直したり、原報・参照の名前を対象gene名へ流用したりしない。

読みやすい日本語の一言メモが必要なら、**整理担当の別LLMを使ってよい**。同じstageのWet・report・対応する保存済み論文表を渡し、次の要領で表示用の `wet_candidate_summary.md` を別に作る。これは任意の文章整理で、追加の科学reviewや候補採用の条件ではない。

> 入力Wetの全候補を、遺伝子名・ID／位置／相同遺伝子（種・ID）／Tier／論文／一言メモの6列に整理する。元の候補集合・行順・ID・位置・Tier・論文IDを保持し、MtAの報告事実、AtBの配列参照、対象PsCを区別する。一言メモは保存済みの理由から「なぜ今回の表現型の候補か」を一文で示し、記録にある場合だけcopyの曖昧さや保留を添える。未記録の事実・引用・実験提案を補作しない。表の冒頭に元WetのSHA-256とstageを記録する。

整理後は、IDの集合と行順、Tier・位置、論文IDが同じ入力に対応することを確認する。正本Wetとhashが記録された既存reportは上書きしない。科学的な修正が必要と分かった場合は、表示文の修正と混ぜず、既存手順の新generationに理由と判断を記録する。

## 未採点の接続確認例

エンジニアリング上の接続確認には、実QTLとの一致を要求せず、公開FASTA/GFFから原因情報を使わずに固定した遺伝子区間を使える。表現型・機構は一次論文から準備し、truthは与えず採点しない。区間座標、全gene ID、入力由来とSHAを凍結し、候補調査・資料不足・source review・最終Wetを別々に記録する。

2026-10-01の確認では、公開[Cameor V2](https://doi.org/10.57745/5MDE99)の `chr7:10000001-62712084`（1000遺伝子）と節間短縮表現型を使用した。抽出タンパク質1000件は公開注釈の対応isoformと全件一致し、clean wheelから実BLASTを実行した。Bの816単位は442依頼単位と374資料不足単位に分かれ、3担当の実回答442件を統合した。147未回答からの再開、および別attemptの1件失敗からの再開で、healthy raw回答を保持した。Cの5件はprojection未確定のledger-onlyであり、C provider failureを注入した別stageでもA/B候補ID・支持属性は保持された。

この例の第一段120候補は、上位14候補とT3境界106候補に分かれる。第一段の件数やdriver内部の `prediction_calls=0` は科学性能やagent実行数を表さない。実際のagent回答・原典レビュー有無・境界採用stageを併記して、上記のfrozen checkpointから完了させる。

## B Web discovery return and re-freeze

複数担当のdiscoveryを再配布へつなぐ場合は、既存処理をまとめた `agentic_route_b.py extend` を使える。モデルや検索は起動せず、元のfreezeを保存したまま新bundle、読取りassignment、reuse-plan、再利用回答、再配布一覧を作る。

```sh
python scripts/agentic_route_b.py extend \
  --bundle "$OLD_BUNDLE" --bundle-sha256 "$OLD_BUNDLE_SHA" \
  --old-prompts "$OLD_PROMPTS" --old-response "$OLD_MERGED_RESPONSE" \
  --discovery discovery-01.json --discovery discovery-02.json --out "$run/b-extended"
```

`--old-response` はmerge済み回答とそのsibling execution manifestを指定する。元のraw回答も必要。既存のhash-bound import監査があれば `--import-manifest` を渡す。省略時は元checkpointを照合した後、元bundle/prompts/responseを既存semantic importerで機械的に再監査する。論文の再解釈やモデル実行ではない。既存 `reuse-plan` の全prompt digest・raw・import条件を満たす回答だけを再利用し、原executor/modelとraw SHAを保持する。未完・既知のimport拒否・変更promptは再実行対象となる。

`redispatch.json` の `reading_shards` だけを担当へ渡す。1 unitごとのassignmentなので再利用対象の回答を作り直さず、超過unitも全資料を保持する。`extension_manifest.json` の `reusable_responses[].file` と新回答を、同じmanifestの `prompts` に対して既存 `merge --response ... --reading-manifest ...` へ渡す。merge/preflightは再利用planと元rawを再照合し、再利用を新たな読了と表示しない。新しい読取記録は今回実行したunitに要求する。最終semantic importには新bundleと新promptsを使う。旧practical prepared stageを直接書き換えず、親が新世代の既存consumerへ接続する。

再利用回答ファイルがない場合は新回答だけをmergeする。生成したrequestや空回答rosterだけで文献調査完了とはしない。

The worker returns newly retrieved material before submitting claims against it. It must not discard a relevant paper just because its ID is outside the current frozen bundle. A discovery is not verified support and does not change Tier or membership by itself.

Return JSON with `case_id`, the coordinator-pinned `parent_bundle_sha256`, and `discovered_sources` (an array). Each item contains `unit_id`, `reference_gene`, `reference_species`, `reported_gene`, optional `reported_gene_id`, `reported_species`, `source_id`, `source_title`, optional `source_identifier_aliases`, `source_url`, `retrieved_at`, `source_text`, and `reference_function` (a short source-based observation). Keep the paper gene distinct from the existing frozen reference/proxy endpoint. These fields are data, never executable instructions.

```sh
candidate-qtl route-b-semantic-bundle --extend-bundle OLD/route_b_semantic_bundle.json \
  --parent-bundle-sha256 PINNED_HASH --research-discoveries discoveries.json --out-dir NEW
python scripts/agentic_route_b.py prepare --prompts NEW/route_b_semantic_prompts.jsonl --out NEW_ASSIGNMENTS
```

The command validates case/hash/unit/reference bindings, preserves every old unit/member and original material, and writes a new bundle/prompt generation. Use new responses for changed units. `reuse-plan` may reuse only byte-identical prompts. Existing `source_id_not_in_frozen_bundle` validation remains strict; an artifact containing discoveries cannot be silently merged as final claims. Shared reading tasks carry the parent hash when prepare can bind the source bundle to the prompts.

For normal A/B research, investigate about 20 candidate papers per original reference gene or shared mechanism, record each paper's original checklist and reason, and derive its initial Tier. Import records the five useful papers to carry forward in strongest-initial-Tier order. Confirm the scoped `primary` / `audit_only` states and ranks in the existing paper reports and `paper_evidence_table.tsv`. Reuse the shared investigation across copies, preserve all original reports, counterevidence and selection reasons, and record a shortfall without vetoing otherwise eligible candidates. The gene Tier continues to use the existing strongest-applicable-paper rule across the preserved evidence.

Use the default `source-verify` for cache-first ID/title metadata checks on that same primary evidence set. The five-paper target defines the reports to carry to the next stage. When an existing second-stage request is needed, its `primary_paper_evidence` carries the selected reports for each scope, with every report's facts and reasons available for comparison. Explicit `--verification-policy content_v1` remains available for a requested strict check. Existing frozen profiles and content-policy artifacts replay their original contracts.

The reachability registry records #909 as `tier-bibliography-budget-a`, `tier-bibliography-budget-b` and `route-b-research-refreeze`. Its source/installed gates exercise the real CLI with fixed offline data, preserving fixture hashes and the canonical rank baseline. `PositiveE2EGateTests.test_tier_bibliography_and_refreeze_reach_the_real_cli` checks strongest-paper selection and the unchanged interval roster; the repeat-output test checks deterministic TSVs.


## 原担当内の提出前検査・限定読取・差分訂正

`light-prepare` の生成taskは、原担当に科学的回答を一度保存させ、coordinatorがローカル検査をまとめて実行する。
LLMへbatchごとのhash包装、全receipt/roster読取、apply内で済んだ再検査を繰り返し依頼しない。
通常は `light-review-batch` / `light-apply-batch` を使い、下記単発コマンドは限定診断用とする。
coordinatorが生成receiptからbundle/indexの絶対pathとSHA256を渡す。indexがtask自体をhashするため、task本文に最終index SHAを埋め込まない。
次の `$bundle` / `$index` / SHA / `$batch` はその固定値、`$results` は保存した原回答を指す。

```powershell
python scripts/agentic_route_b.py light-preflight --bundle $bundle --bundle-sha256 $bundleSha --light-index $index --light-index-sha256 $indexSha --batch-id $batch --results $results --out preflight.json
python scripts/agentic_route_b.py light-read-material --bundle $bundle --bundle-sha256 $bundleSha --light-index $index --light-index-sha256 $indexSha --batch-id $batch --gene-id GENE_ID --field /genes/0/reference_function
python scripts/agentic_route_b.py light-read-source --source saved-source.json --source-sha256 $sourceSha --pointer /abstract
python scripts/agentic_route_b.py light-read-source --source saved-source.txt --source-sha256 $sourceSha --line-start 10 --line-end 30
python scripts/agentic_route_b.py light-apply --bundle $bundle --bundle-sha256 $bundleSha --light-index $index --light-index-sha256 $indexSha --batch-id $batch --results $results --patch correction.json --patch-sha256 $patchSha --out-results corrected.json --receipt correction-receipt.json
```

`--field` はunit内、`--pointer` は保存JSON内のRFC6901 pointer。selector未指定は索引だけを返し、読了とはしない。
全文は改変せず、path/SHA/pointerまたは行範囲を返す。`--max-chars`（既定12000）を超える選択は `narrow_selection_required` となる。400文字等で切り捨てず、選択範囲を狭めて読み直す。source読取はA/Bどちらの保存資料にも使え、新規取得・要約・生物学判断は行わない。

preflightはJSONの `state` と `findings` を確認する（正常起動の終了codeだけでは判定しない）。
ID/title・報告gene名/ID/species・9回答・roster/bindingの欠落を一括で返す。
`ready` は提出形式が整った意味であり、論文確認済み・全gene考察済みとは異なる。
`not_assessed` / `retrieval_failed` やunknownを埋めるための推測は禁止し、科学的完了は既存member coverageとcompletion gateで報告する。

patchは次の形式で、preflightの `base` をそのまま使う。`expected` は元値、存在しないkeyだけ `{"$missing":true}` とする。

```json
{
  "contract": "agentic-route-b-light-patch-1.0.0",
  "base": {"results_sha256":"...", "bundle_sha256":"...", "light_index_sha256":"...", "batch_id":"b01"},
  "operations": [{"op":"set", "record":{"kind":"gene_decision","id":"GENE_ID"},
    "path":["paper_assessments",0,"checklist","focal_species"],
    "expected":{"$missing":true}, "value":"unknown"}]
}
```

例のunknownも原担当が確認した回答である場合だけ記録する。kindはv2の `gene_decision` / `paper_observation`、v1の `legacy_gene`。
共有observationの訂正には全利用geneの `affected_gene_ids` が必要で、同一論文でも別gene/形質/観察を自動統合しない。
旧hash・元値の不一致、重複/重なった訂正、record追加削除やID変更を拒否する。原回答と非対象recordを保持し、新回答・訂正receipt・再preflightを残す。
訂正後の回答を既存 `light-convert` に渡し、source確認とcanonical consumerは従来どおり実行する。新規科学判断が必要な訂正は単なる転記修正として扱わない。

### 一括処理とLLMへ返す範囲

```powershell
python scripts/agentic_route_b.py light-review-batch --request review-request.json --request-sha256 $requestSha --out review-generation
python scripts/agentic_route_b.py light-apply-batch --review review-generation/review.json --review-sha256 $reviewSha --decisions decisions.json --decisions-sha256 $decisionsSha --out corrected-generation
```

requestは `contract=agentic-route-b-light-review-request-1.0.0`、共通 `binding={bundle,bundle_sha256,light_index,light_index_sha256}`、`batches=[{batch_id,results,sources:[{source,source_sha256,pointers:[...]}]}]`。sourceは任意で、明示した保存資料・選択範囲だけを返す。相対pathはrequestの所在から解決する。
reviewは詳細preflightをファイルに保存し、stdoutには集計・不足項目・選択した原資料値・reviewのpath/hashだけを返す。完全なmember状態やreceiptを原担当に再読させない。

decisionsは `contract=agentic-route-b-light-batch-decisions-1.0.0`、`batches=[{batch_id,operations:[...]}]`。operationsは前記set操作と同じで、変更しないbatchは `{batch_id,unchanged:true}` を明示する。全review対象batchが一度ずつ必要。Pythonが既存のhash-bound patchを構築し、全操作を検査してから新成果物を作る。原回答・非対象記録・source binding・共有観察の全copy確認を維持する。

原LLMがすでに書いた明示的な訂正は、その保存ファイル/hashを使ってPythonだけで配送できる。転記のためだけに新しいLLMを起動しない。一方、新しい科学判断や未記録の回答は機械で補わず、必要な原資料と不足項目だけを原担当へ返す。これは研究判断の自動化ではない。

計測は既知訂正の機械配送、判断を依頼するモデル処理、元の文献調査を分ける。保存JSON bytesの短縮や機械配送のモデル呼出し0を、全研究フローのtoken削減率と呼ばない。source/入力とモデル条件を揃えた実counter比較と科学的記録の一致を確認してから効果を報告する。

## Prompt / transport failure attribution

Use the current generated prompt and exact schema/validator from the execution
checkout. Preserve their hashes with the first raw response and actual tool/input
receipt; a prompt requesting Web search does not prove a tool call executed.
Run the existing read-only payload preflight on a fixed response before import.
A failed preflight is a transport result, not a biological exclusion. Repair only
unambiguous representation in a new generation, or return a bounded clarification
request; do not infer missing scientific facts or silently normalize ambiguous enums.

For Direct C, `render_direct_proposal_prompt()` in
`src/candidate_qtl/route_c_direct_pipeline.py` renders exact proposal columns and
full enum literals from the validator constants. The harness
`route_c_practical_prompt.txt` is its generated copy. Use full values such as
`S1_same_species` and `D0_causal_perturbation_or_rescue`, not display shorthand.
For A, `research_annotations` belongs within each `paper_reports_json` element;
it is not an arbitrary hypothesis-level extension.

Distinguish flow/transport defects, contradictory or underspecified prompts,
model instruction-following failures, research reasoning, and user/executor input
or tool failures. Attribute a clear model format violation only after checking the
actual supplied instructions against the rejected response. Request generation,
request dispatch, completed model execution and validated import are separate states.
An unread source or failed search is not biological absence. Ordinary A permits
focal-species primary references; only explicit `benchmark_mode` excludes them.
Focal-species and cross-species research remain mechanism-driven, with recorded
queries/coverage and reported gene, sequence proxy and focal locus kept separate;
there is no T1 quota or Tier-only research stage.

通常Bの研究taskは、資料の一覧化・閲覧、原研究者による初期機能判断、文献調査を区別する。
全区間geneの初期機能判断を依頼された場合、上位候補発見後も共有判断とunit/member対応で全geneを覆う。追加文献は候補・Tier・重要な対応を変え得る具体的な手掛かりに限定する。
この依頼は、利用可能な機能と今回の表現型の関係を各geneで少なくとも一度考えることを含む。同じ機能・同じ関連contextは考察と適用gene IDsを共有し、重要なcopy差だけを別記してよい。索引の閲覧・状態行の一括生成・一般的な不足理由の貼付は考察の代わりにしない。
元の理由・読了範囲では、保存入力の機能材料なし、材料あり未読/未考察、考察済みだが機能/対応が未解決、具体的仮説なし、具体的leadありを分ける。未読/未考察はnot_assessedとし、assessed_insufficient_materialへ一括変換しない。入力の材料欠如は外部情報の不存在ではない。短いfunction prefixだけを読んだ場合、未読の補助材料を明示する。
親は共有考察の全gene対応と情報不足・未調査・考察済み未解決を別々に確認し、未調査が残れば部分実施と報告する。Tier・採否規則・追加agent数は変えず、全件文献検索や未知の機能解明を新たな必須条件にはしない。過去の曖昧な不足記録は後付けで考察済み/情報なしへ再分類しない。
ここで全geneを覆うとは、共有判断・情報不足・未評価を含む各geneの状態を記録することであり、全geneの機能解明や表現型説明を成功させる義務ではない。未評価を台帳へ記録しても、依頼された初期判断の実施完了には数えない。利用できる材料を届けて判断に使うことと、モデル植物でも不足する知識を新たに得ることを区別する。
実際に参照機能と表現型を検討して具体的claimなしと判断したunitは、読了範囲と理由を研究台帳へ残し、既存のcompleted/空claimsで記録できる。全unitの文献検索は一律条件ではない。
一般名だけで解釈不能、未読、取得失敗、未処理の具体的手掛かりはnot_evaluatedを維持し、no-claimを生物学的陰性にしない。
共有資料は対応・例外を保持して再利用し、原文や巨大JSONを会話へ反復出力しない。通常initial_bibliography_v1では別の意味reviewerを増やさず、原九回答を既存機械照合へ渡す。
同種既報のreference対応は取得済みID/assembly情報を優先し、複数copyの候補保持と同一copyの確定を区別する。これらはpromptの明確化であり、候補保持性能やtoken削減の実測保証ではない。

`reference-gene-resolve` は、実際に取得した同種 provider record の ID と明示名で
確認できた参照対応だけを、派生行の `paper_reports_json.reference_gene_id` へ渡す。
原入力の単一参照表示名（括弧内の明示名を含む）はその atomic pair 内に限定し、
任意の表示名を新しい alias として登録しない。名前だけで複数 copy を選ばず、
競合 ID・別種・複合参照は未解決のまま保持する。受渡しの原 anchor・新 ID・根拠は
`retrieval_provenance.paper_reference_rebindings` に記録する。原回答ファイル、
reported gene、9回答、`same_gene_as_focal=unknown`、理由、source 検証方針は変更しない。

照合の種は、記録の学名が宣言した二名法の亜種・変種（`subsp.`/`var.`/`cv.` 等）であれば同種とみなす
（例: 宣言 `Oryza sativa` と記録 `Oryza sativa subsp. japonica`）。宣言側が亜種で記録が種のときは確認にならない。
名前は宣言種の略号接頭辞（`SlGLK2` → `GLK2`）を外した形も照合する。
A が自分で宣言した UniProt ID の記録が実在し種も合うが名前だけ記録に無いときは、留保
（`declared_accession_species_name_unconfirmed`）として採用する。留保は review-only の投影に参加するが、
verified support にも Tier の確定にも使わない。同じ名前・種で名前と種が完全に合う別記録があれば、宣言 ID は採用しない。
留保の行は `identity_review_requests.tsv` に出る。親は独立の判定役（原 A 担当以外）に、論文と記録を見て同一遺伝子かを判定させ、
`reference_gene_name, reference_species, accession, record_gene_name, verdict, evidence_source_ids, reviewer, rationale`
の TSV を作る。`verdict=same_gene` で引用と判定者があり、`record_gene_name` が実際に記録上の名前である行だけが
`reviewed_name_alias` として確定する。ledger は `reference-gene-resolve --identity-review`、
ハーネスでは config の `identity_review` で渡す。判定が無い・否定なら留保のまま残る。

### 同じ親担当による中間・終了時の振り返り

実行前に想定した機構、Tier別候補数とその幅、二段目が必要になる条件を記録する。A/Bの凍結・合流・必要な境界比較が終わるたびに、親担当は保存済み成果物で予測との差を確かめる。benchmarkでは予測・最終候補のfreeze後に真値を開き、真値に合わせた修正を元世代へ混ぜない。振り返りは結果の説明と運用改善であり、追加のTier審査段階や候補採否gateではない。

- 正解保持、実遺伝子数、根拠Tier、区間遺伝子への対応の強さ、未評価範囲を別々に確認する。少数候補、A/B一致、既報の直接回収だけでTierの正確さや未知形質への性能を証明しない。
- 予測との差を、原LLMの発想・解釈、挙げたが未追跡のlead、資料取得、記録/配送、投影・統合のどこで生じたか、raw・実行receipt・正式表から切り分ける。原因不明なら断定せず、問いに対応する最小限の再生・解析を行う。
- 記録/配送の不整合は、同じ原担当・同じ保存済み資料の限定修正を別generationで比較する。新検索・科学判断の変更があれば区別し、修正後の成功を独立のblind試行に数えない。既存指示の運用逸脱ならdispatchや意味保存の確認を直し、既存契約を重複させる仕組みやcase別救済を足さない。

振り返りの結論には、確認できた改善、検証していない一般化、残った未解決点を記す。次に実行するのはその不確実性を解く条件とし、すでに確認できた工程を一律に全再実行しない。
