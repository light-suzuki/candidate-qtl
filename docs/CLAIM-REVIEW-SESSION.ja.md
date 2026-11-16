# 固定入力からA/Bレビューを再開する

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

Change class: `CORE_INTEGRATION`。Refs #708。既存の独立claimレビュー、A投影、
B transfer、family/copy、Wetを接続する。新しい科学的判定・順位・Routeは導入しない。

`PYTHONPATH=src python scripts/claim_review_session.py --help`で操作を確認する。
各操作は別の新規出力ディレクトリへ書く。既存の科学結果を上書きしない。

## 新しい第一段の操作

通常のordinary prepare/review/finishは新policyを既定で保存し、finishに独立レビュー回答を渡しません。
個別CLIでは `source-verify --verification-policy initial_bibliography_v1`、sessionの
`prepare --review-policy initial_bibliography_v1` を明示します。`complete --session ... --request-sha256 ...` に回答やreviewerを捏造して渡しません。
Tierは元のA/B提案の属性を使い、別のTierレビューを起動しません。第一段を本候補としてfreezeした後、必要な本文確認だけを第二段で扱います。

通常のWet境界は`upper_evidence_scoped_v3`です。累積25を初めて超えるTierだけを比較し、上位は保持、後続は台帳に残します。元の調査担当者の理由付き判断を`original_assessment_v1`として機械照合し、別の生物学reviewerを要求しません。25で強制切断せず、後続Tierの追加比較も起動しません。旧`evidence_scoped_v2`・版なしscopeは下記の旧policyで再生します。詳細は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)に従います。

## 互換性用の独立contentレビュー手順

1. Aの`source-verify`出力、Bのproposal/source-verify出力、参照providerレコードを固定する。
   A/Bの提案はそれぞれ別文脈で生成・凍結してから扱う。
2. `prepare --route a|b`へ`--hypotheses --sources --proposal --records`と実際の
   `--proposer-id --proposer-model --case-id`を渡す。Aは`--focal-species`も指定する。
   Bは`--bundle --source-verification`を追加する。
   ファイルはsessionへコピーされ、正確な入力を含む`prompt.json`とhash付きrequestができる。
3. 通常のCodex/OpenCode等の実行担当が、提案者とは別の実行文脈に`prompt.json`を渡す。
   helper自身はモデルを起動しない。実回答・実executor/modelを保存し、推定で補わない。
4. `complete --session ... --request-sha256 ... --response ... --reviewer-id ...
   --reviewer-model ...`で、同じ固定入力から既存の検証器とcandidate freezeを実行する。
   Aはbait freeze、Bはtransferとinterval candidate freezeまで進む。
   再検索・原文の再取得は行わない。支持がないという回答も正常なレビュー結果である。
5. `handoff --a-session ... --b-session ... --genes ... --proteins ... --atlas ...
   --out ...`でAのローカルBLAST投影、A/B union、family/copy、Wetを実行する。
   A/Bのcase ID、全回答の存在、固定済みファイルhashを確認する。
   Aの既存interval-family状態表も生成してbridgeへ渡す。既存のreference_genome
   family indexは`--family-index`で明示指定でき、`--resources`、環境変数
   `CANDIDATE_QTL_FAMILY_INDEX`、標準キャッシュからも検証して解決する。
   family indexと隣接manifestのパス・SHA・解決経路をhandoffの`inputs.json`へ保存する。
   欠損・不適格な自動発見indexは未評価として記録し、未実行を生物学的陰性としない。
   明示指定した不適格scopeはfail-closedで拒否する。
   `completed_partial_review`を完全実行へ読み替えない。未評価と明記した回答は回答済みとして扱う。
6. 旧`evidence_scoped_v2`の再生では、最初のoverflow Tierと後続の既存適格候補を固定したboundary requestを独立実行担当に渡す。
   回答は既存のproposal TSV schemaで保存し、boundary requestの`source_wet_sha256`を各行へ保持する。
   `finish-wet-review --handoff ... --completion-sha256 ... --proposals ... --independent-review ...
   --reviewer-id ... --reviewer-model ... --executor-id ... --executor-model ... --out ...`で取り込む。
   理由付き採用/保留案と別担当の独立確認を既存boundary consumerへ渡し、新しいWet世代を保存する。検証済みsupportへの昇格とは別であり、legacyな`--verification`でboundary reviewを迂回できない。
   `next_boundary_scope_manifest`がある場合だけ、その世代のcompletion hashと次scopeを使う。後続scopeの有無と継続は、保存された旧policyに従う。旧unversioned scopeは元のwhole-tier/continuationで再生する。
   不完全な比較はpending候補を保持し、`resume_from_completion`が指定する元の第一段completionから別出力へ再開する。不完全な選択出力自体を再開元にしない。既存round budget消尽はpartial/未評価・理由・再開先を記録して終了する。

## Bの参照機能reviewと同じ判断の共有

旧B bundleの `source_function_v1` は、出典IDの存在だけでなく参照機能の実内容、参照gene/種と区間unitの対応、具体的な条件付き表現型推論を確認する。対象種の因果論文や同一表現型は不要で、Type3はreviewのままWetへ渡せる。旧bundleのlegacy方針と旧review packetは変更せず読む。原著1〜2報のsoft budgetは既定化していない。採否条件の正本は [Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md)。

既存B review promptには、凍結入力の完全一致する `reference_function_questions` と対応claim IDsが入る。任意の `shared_reference_function_reviews` はquestionごとに1つの参照機能判断を返し、各claimの `reviews` は `reference_function_question_id`、`copy_application_state`、`copy_application_reason` を別々に記録する。資料を共有してもcopy差を共有判断で消さない。適用不能・未読はpositive判断を借りられず、適用未確定のplausible copyはreviewに保てる。

question IDはhypotheses/source全体のhashに依存し、参照gene・種・source・機能・機構・phenotype relation・代替仮説・不確実性の違いを区別する。変更後の資料や別世代へ旧判断を無条件流用しない。auditは共有機能判断とclaim/copy適用の件数を別に保存するが、実読解回数・token・費用削減は別途測定が必要である。

`complete`の同一入力による再呼び出しはhashを再確認して保存済み結果を返す。
別回答・同一executor・改変された入力/出力は拒否する。処理中のクラッシュや失敗は
そのディレクトリを保持し、同じ場所で暗黙に再実行しない。失敗原因を確認し、同じ固定入力から
別の出力session/attemptを作成する。モデル回答の品質改善のための再試行はこの機構の役割ではない。

各completionの`files`が成果物の固定境界であり、工程全体の完了と支持の有無は別に読む。
一段目Wetのfreezeは二段目で変更しない。ベンチの採点はこのhelperには含めず、予測freeze後に行う。

少数ケースの接続確認はmethod validationであり、合成入力や既知のfixtureから自律的な科学性能を主張しない。
