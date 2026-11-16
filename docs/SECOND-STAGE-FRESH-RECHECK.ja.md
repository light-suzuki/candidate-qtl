# 第2段：残った候補を新しく調べ直す

CORE_INTEGRATION / #732。第1段の説明や全証拠を再構成するのではなく、残った候補IDと配列から新しく調査する。

## 動作

既存の `guided-profile --profile triple-route --execute` は、requestのstageにより次の処理を分ける。現行 `boundary_review_request` はsealed scope・hash・全candidate IDを検査し、`boundary_review` に未実行状態と既存の `candidate_qtl.boundary_tier_review.consume_boundary_selection` を記録する。候補は切り詰めず、要求作成時に新規BLASTやモデル実行を追加しない。採用の手順は後述の `finish-wet-review` と[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)を使う。

明示的な旧 `adaptive_synthesis_request` だけは、既に指定されている `--proteins` と明示指定の `--protein-db` によるローカル再調査へ進む。以下の配列再調査はこの旧contractか、単独入口を明示的に使う場合の説明である。dry run・screening・要求なしの実行に追加処理は行わない。

1. 第2段requestの全候補IDを受け取る。前段の回答や説明は再調査promptへ移植しない。
2. 対象候補の蛋白質配列を抽出し、明示された参照FASTAまたはローカルBLAST DBに対して `blastp` を新しく実行する。
3. 残存候補の蛋白質同士も比較する。同じgeneの複数isoformと、別geneを区別できるID対応表を保存する。
4. 全HSP原本と、query/subjectの両coverage・残基集約identity・最大HSP identity・配列長・E-value・参照ID/titleを保存する。返された複数参照をbest-hit 1件に縮約しない。
5. 各候補について、文献・機能・部分一致・類似copyを新しく調べるためのprovider-neutral JSONL promptを出す。

`second_stage_recheck/second_stage_recheck_manifest.json` から原本、結果、promptを辿れる。外側のguided manifestにもその参照を記録する。

外側の `guided_profile_manifest.json` は `canonical_execution` にcanonical manifestとWetの内容hash、現在のmembership/input coverage状態、必要作業の残りと再開先を保存する。semantic入力やboundary reviewが未完なら `overall_status=partial`、`final=false` とする。旧profileの段階別完了情報は `legacy_profile` に残し、canonicalの完了と混同しない。`completed_stage=first_stage` は第1段の完了範囲を示し、科学的な因果確定や第2段レビュー完了を意味しない。

## 単独で再調査する

既に第1段が終わっている場合や、参照FASTA・新しいDBを使って再調査する場合は、**新しい出力directory**を指定する。第1段の再実行は不要。

```powershell
python -m candidate_qtl.second_stage_recheck `
  --request 'run/second_stage_request/adaptive_synthesis_request.json' `
  --out-dir 'run/recheck-r2' `
  --phenotype '観察した表現型' `
  --species '対象種' `
  --proteins 'inputs/candidate_proteins.fa' `
  --reference-proteins 'inputs/reference_proteins.fa'
```

実際のrequest位置はWet handoff manifestの `adaptive_synthesis_request` を使う。FASTAではなくBLAST DBを使う場合は `--reference-proteins` の代わりに `--protein-db PREFIX` を指定する。両方の同時指定はしない。

FASTAの蛋白質IDがgene IDと違う場合は `--protein-gene-map FILE` で `protein_id` / `gene_id` 列を持つTSVを渡す。`.t1`等を推測で削って対応させない。明示対応した複数isoformは全て残す。

## 調査結果の意味

- BLASTは**相同性の調査材料**。一致率・bit scoreだけでtrait-specific support、唯一のortholog、最終候補順位を作らない。
- HSP集約はPR #725/#726の `aggregate_hsps()` を再利用する。両方向のcoverageは同じalignmentから算出し、subject側への表の転置を「逆方向BLAST」と呼ばない。
- 候補同士の比較は残存蛋白質の範囲だけ。全ゲノムのparalog探索、系統樹、synteny検証、reciprocal-best-hit証明を実装した意味ではない。
- 機能・種・系統・実験対象・表現型の対応は再調査/レビューで確認する。旧A/B/Cの一致を新しい独立証拠に数えない。
- 複数copy、配列欠損、実行不能、正常なヒットなしを区別する。いずれも元のgene IDを消さずpromptを用意する。取得情報がないことだけで負の生物学的判断をしない。

## 数値は運用設定

約20は目安。現行boundaryの対象と開始条件は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)を正本とする。旧adaptive要求の25件超triggerとの互換性は保持し、候補数の上限・強制padding・Top-Kは追加しない。
参照FASTAと候補同士の比較は、パネル内の全targetが報告可能な件数を指定する。DB検索の既定取得件数100、探索E-value 10、DB検索2 threads、各呼出600秒は検索/実行の設定であり、候補の生物学的採否基準ではない。単独呼出では `--max-targets` / `--evalue` / `--threads` / `--timeout` で変更できる。
`-max_hsps`は付けず、取得したHSPを後処理で1件へ制限しない。BLAST自体の検索・report設定による未取得の一致まで回収したとは主張しない。

## 実行・再開と失敗

ローカルBLASTのみ実行し、外部モデル/APIや文献検索を勝手に起動しない。蛋白質配列そのものはpromptへ入れない。候補ID・参照title・数値等は入るので、外部サービスへ送る際は通常の入力確認を行う。

同じguided出力の再開では、入力hash/設定と作成済み結果のhashが一致する再調査記録を再利用する。新しい検索をしたとは表示せず `reused_previous_report=true` とする。入力/結果が変わった場合は旧結果を上書きしない。BLAST導入後の再試行や新しい参照での実行には、単独コマンドで新しい出力先を使う。

FASTA入力と結果は内容SHA-256を記録する。BLAST DB prefixはlocatorのみで、alias先を含む全DBの内容freezeではない。DBの厳密な再現が必要な実行では、運用側でDB snapshotを固定する。過去記録の再利用は現在のDB内容を再検索したという意味ではない。

## ここまでで実装したもの／していないもの

実装：残存配列抽出、実 `subprocess` による再BLAST、候補同士の比較、既存HSP集約の再利用、状態保存、新規調査用prompt、通常guided経路からの接続、同一出力の検査付き再利用。

**LLMの文献再調査・科学的レビュー実行は既存agent/runner側の仕事。** JSONLは既存 `llm-run` の `prompt_id/gene_id/system/user` 形式。rawなレビューを既存検証へ持ち込む判断は別です。review-inclusive selectionを行う場合は、proposal TSVの全行に`source_wet_sha256`、明示的な`decision`、非空の`executor_id`・`executor_model`・`reason`を付ける。独立review TSVを使う場合は全行に同じWet hash・`proposal_sha256`・reviewer identity/model・reasonを付け、`finish_wet_review`へproposalとは別のexecutor/reviewer IDを渡す。proposalがWetの一部だけでも未レビュー候補は保持され、manifestの`review_coverage=partial`として記録される。独立review未供給時はdeferできず、selection stateは`not_evaluated_no_independent_review`になる。

この再調査結果自体はproposalです。提案単独では`wet_candidates.tsv`・supported union・support flagsを変更しません。独立reviewをハッシュ・ID roster・理由つきで検査したreview-inclusive consumerだけが、元Wet行をコピーした新しい`wet_candidates.tsv`と保留台帳を別directoryへ出力します。保留は非原因の証明ではなく、元のfreezeは変更されません。

入口は `scripts/claim_review_session.py finish-wet-review`（内部APIは `candidate_qtl.claim_review_session.finish_wet_review`）です。LLM自動実行は行いません。boundary requestの`source_wet_sha256`をproposal各行へコピーし、`--completion-sha256`にはfirst-stage `completion.json`のSHA-256を渡します。`decision=adopt`は採用、`defer`は採用集合から保留台帳への移動、`not_evaluated`は未評価のまま元候補を保持します。意味が曖昧な`hold`は拒否します。支持未確認だけを保留理由にせず、実際の比較判断と理由を残してください。

proposalと独立reviewの不一致は候補単位で解決します。両者が同意した`defer`だけを保留台帳へ移し、不一致候補は`not_evaluated`として採用集合に残します。全候補を一括で未評価に戻さず、両判定・両理由・解決状態は`review_decisions.tsv`へ記録します。`defer`は実験上の優先度判断であり、生物学的negativeやclaim未検証の証明ではありません。
独立reviewはproposalの証拠・identity・表現型解釈と相対的な実験優先度を監査するもので、同じfavorite 20を独自に再選出したり非原因を証明したりする必要はありません。別の有力候補があるだけでは不一致理由にならず、器官・時期・作用機構・copy coverageなど候補間の比較を伴う理由を残します。

```python
from candidate_qtl.claim_review_session import finish_wet_review, sha256

result = finish_wet_review(
    "first_stage", sha256("first_stage/completion.json"), "second_stage",
    proposals_arg="selection_proposals.tsv",
    executor_id="luna-proposer-run", executor_model="gpt-5.6-luna",
    independent_review_arg="independent_review.tsv",
    reviewer_id="luna-reviewer-run", reviewer_model="gpt-5.6-luna",
)
adopted_wet = result["adopted_wet_candidates"]
```

CLIでは次のように実行します（`REQUEST`はboundary request、`PROPOSALS`と`REVIEW`はそのrequestの`candidate_gene_ids`だけを対象にしたTSVです）。

```powershell
$first = 'C:/cqtl-run/case/first_stage'
$completionSha = (Get-FileHash "$first/completion.json" -Algorithm SHA256).Hash.ToLower()
python scripts/claim_review_session.py finish-wet-review `
  --handoff $first `
  --completion-sha256 $completionSha `
  --out 'C:/cqtl-run/case/boundary_round_1' `
  --proposals $PROPOSALS `
  --independent-review $REVIEW `
  --reviewer-id 'independent-reviewer-run' `
  --reviewer-model 'gpt-5.6-luna' `
  --executor-id 'proposal-executor-run' `
  --executor-model 'gpt-5.6-luna'
```

結果JSONの`adopted_wet_candidates`、`wet.wet_candidates`、`files`を保存し、`wet_candidates.tsv`の実パスと採用件数を報告します。`next_boundary_scope_manifest`があり、`adaptive_round < max_adaptive_rounds`なら、次のscopeの`source_wet_sha256`を次のproposal/review各行へ束縛し、出力directoryの`completion.json`を新しい`--handoff`として同じコマンドを再実行します。

```powershell
$next = 'C:/cqtl-run/case/boundary_round_1'
$nextSha = (Get-FileHash "$next/completion.json" -Algorithm SHA256).Hash.ToLower()
python scripts/claim_review_session.py finish-wet-review `
  --handoff $next `
  --completion-sha256 $nextSha `
  --out 'C:/cqtl-run/case/boundary_round_2' `
  --proposals $NEXT_PROPOSALS `
  --independent-review $NEXT_REVIEW `
  --reviewer-id 'independent-reviewer-run-2' `
  --reviewer-model 'gpt-5.6-luna' `
  --executor-id 'proposal-executor-run-2' `
  --executor-model 'gpt-5.6-luna'
```

`max_adaptive_rounds`到達時、requestが無い場合、またはreviewが未完の場合は、first-stageまたは直前世代の採用集合を保持して終了状態を記録します。要求生成だけを実行完了や因果性能と読み替えません。

独立reviewと提案の原本コピーは第2段directoryへ保存し、completionのhash一覧に含めます。採用集合の正本は`adopted_wet_candidates`、元集合はfirst_stageのWetです。第2段の出力を新たなgeneric fillへ戻して保留候補を再採用しません。

## 検証

`python -m unittest discover -s tests -p 'test_second_stage*.py' -v`

ローカルでは初版44件と追補28件、合計72件の機能/合成テストを実行。BLAST部分はfake executableによる実子processを使用し、BLASTの数値精度・実データ性能の検証ではない。実guided parser/handlerの接続試験では、無関係のupstream producerをstubに置換し、再調査本体とファイルI/Oを実行する。全repo suite・clean wheel・実NCBI BLAST・実LLM/実データE2Eは未実施。

NCBI公式BLAST+ manual Table C1（outfmt / subject / max_target_seqs / max_hsps）:
https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.options_common_to_all_blast/

## 追補：閾値・family/subfamily・copyも詳しく再点検する

再BLASTだけでなく、同じ測定結果を複数の条件で再評価し、nativeなfamily資料も
新しい調査promptへ渡す。元の候補数・参照hit・isoform・HSP原本は減らさない。

### 相同性の条件を変えて比較

既定で次の2つの**調査用プリセット**を並べる。第1段の実採否条件の再現ではなく、
全植物・全familyに妥当性を検証したcutoffでもない。

| 調査用条件 | 両方向の残基一致率 | query coverage | subject coverage | E-value |
|---|---:|---:|---:|---:|
| broad_alignment_review | 30%以上 | 0.5以上 | 0.5以上 | 1e-5以下 |
| strict_full_length_review | 50%以上 | 0.8以上 | 0.8以上 | 1e-10以下 |

すべて調整可能な初期値。`illustrative_review_heuristic_not_calibrated`と記録する。
50%や80%に科学的な絶対性を与えず、短いペプチド・遠縁・部分gene modelなどは
未達でも調査対象から除かない。探索時のE-valueと、取得後のこの確認条件は別。

各hitについて `meets_review_conditions` / `below_review_conditions` と未達の軸を出す。
厳しい条件を満たすhitだけをpromptに渡すのではなく、**未達のhitも並べて見せる**。
一致率を満たしてもsubject coverageが低ければ、全長での対応とは混同しない。
検索上限/E-valueのため未取得の参照まで評価したとは扱わない。

単独module CLIでは `--review-profiles profiles.json` で比較条件を明示できる。
JSONは上表のprofileに対応する `profile_id`, `basis`, `min_identity`,
`min_query_coverage`, `min_subject_coverage`, `max_evalue` を持つobjectのlist。
`basis`にはその調査条件を使う理由/検証範囲を記す。未知field、重複profile ID、
非有限値、範囲不正は入力エラー。gene rank/scoreは受け付けない。

### familyを粗い名前の一致で済ませない

既存形式を使い、`--interval-atlas interval_gene_atlas.tsv` と
`--family-copy family_copy_discrimination.tsv` を任意に渡せる。

- discovery bucket/legacy familyと、evolutionary familyを分離する。
- subfamilyが不明・orthogroup proxy・symbol stem由来なら、詳しい確認を要求する。
- 既存のtarget相対のsubfamily/domain一致・不一致を、対象targetごとに報告する。
  targetの検証状態/IDが不足すれば `not_evaluated_target_not_bound`。
- familyの多数派ドメイン構成との差はQCであり、target相対の機能差とは区別する。
- one-to-many/paralog、部分gene model、競合するfamily/targetは消さず可視化する。
- 対応はnative gene IDの完全一致。同じgeneに複数family/target行があれば全て保持し、
  一方のsource/target判定を他方へ付け替えない。順位・score列は読み込まない。

通常guided実行では、同じ実行のB0 atlasと、canonical結果が明示する
`effective_family_copy`または`family_copy`を自動で渡す。未知のhomeや資源を探索しない。
資料不足・対象geneの行なしは未評価のまま新しい調査へ渡す。6軸を全部そろえる必要はない。

**今回のfamily確認は既存/明示供給された分類の再点検であり、新規HMM/InterProScanや
系統樹計算を実行したものではない。** 新しく作成したnative形式の分類を単独入口へ
渡して比較できる。報告されたtarget検証状態も独立に再検証済みとは呼ばない。

### 保存と実装境界

追加出力は `alignment_review_checks.tsv` と `second_stage_review_checks.json`。
候補別promptにも条件比較とfamily資料/確認事項が入り、hashは再開時に検査する。
`candidate_membership_effect=none`、`review_state=not_executed`を維持する。
数値未達は `support_falsified`に変換せず、機能supportや最終採用集合も自動では作らない。

新contractは `second-stage-fresh-homology-1.1`。旧1.0の報告を新しい確認まで済んだと
読み替えないため、旧出力の再利用時は新しい出力revisionを要求する。元の第1段は不変。
通常guidedのCLI引数は増やさず、数値条件を変更する場合は単独入口と別revisionを使う。

資料の意味の参照（機能同一性の閾値を規定する根拠ではない）:
- EMBL-EBI InterPro entry types: family/superfamily/domainの区別。
  https://www.ebi.ac.uk/training/online/courses/interpro-functional-and-structural-analysis/what-is-an-interpro-entry/interpro-entry-types/
