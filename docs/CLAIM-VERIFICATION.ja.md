# 引用の実在確認から、主張の裏付け確認へ

引用listの分割は共通`split_citations`を使い、DOI suffix内のsemicolon/commaを保持する。
次が明示的なPMID/PMCID/DOI開始なら区切り、pipeもlist区切りとして扱う。
曖昧な`DOI;裸数字`を勝手に分割せず、複数citationは`PMID:`等のprefixで明示する。
原citationの文字列・source hashと、実在照合・引用範囲・bindingの検査は引き続き別々に保持する。

`source-verify` は「捏造引用」の穴を塞いだ。解決しない PMID/DOI は `rejected` で、
モデルの記憶は引用元にならない。しかし**引用が実在することは、主張が支持されている
ことより弱い**。実在する論文が、引用された因果的主張には一切触れていないことがある。

Route A の `source_ids` が UniProt/TAIR ID だけの場合、それ自体は論文引用ではない。
`reference-gene-resolve` は、単一の矛盾しない参照IDから取得したレコードの論文IDを
元のIDに追加し、取得元を記録する。名前と種を指定しないaccession-only入力では、既存の PMID/DOI や不明な文字列は変更しない。
複数ID・宣言したIDとの不一致・論文のないレコードから支持を作らず、取得失敗も
生物学的な否定として扱わない。追加した論文は引き続き本文と主張の検証を必要とする。

名前と種を指定した仮説では、モデルのaccession/引用を根拠に名前解決を省略しない。
機械で解決した別identityへ元の引用を移さず、入力履歴として保持する。
claim未確認の候補も区間対応を確認できればWetへreview候補として残せるが、
本書のverified support判定とは区別する。

通常の`initial_bibliography_v1`はcache優先のID・タイトル照合だけを行う。Europe PMC liteでPMCIDが未取得の場合に限り、公式NLM PMC ID Converterの単一対応からPMID/DOIを確認し、既存cacheまたは同じlite照合でタイトルを取得する。対応の矛盾・複数対応・取得失敗は成功にせず、宣言済み全aliasの照合を維持する。取得URL・時刻・応答hash・ID対応はmetadata provenanceへ保存し、抄録・本文の取得や追加の意味判断は行わない。

区別しなければならないもの:

| 実際に論文がしていること | evidence_relation | claim_support_status |
|---|---|---|
| ノックアウト＋相補性検定で表現型が再現/回復 | `direct_genetic_test` | `supported_direct_genetic_test` |
| 過剰発現・RNAi・VIGS 等の摂動表現型 | `perturbation_phenotype` | `supported_perturbation_phenotype` |
| 発現量が相関している | `expression_association` | `associative_support_only` |
| GWAS/QTL/連鎖の関連のみ | `genetic_association` | `associative_support_only` |
| 総説・背景記述 | `review_or_background` | `background_or_review_only` |
| 遺伝子名は出るが機構の記述が無い | `mention_without_mechanism` | `review_required` |
| 引用対象の遺伝子/family/pathwayに触れていない | `no_claim_term_in_source` | `unsupported_mention` |
| 取得できるテキストが無い | `not_evaluated_no_text` | `not_evaluated` |

## 2軸を分ける

`source_validation.tsv` に列が2本ある。

- `citation_verified` — **識別子が実在の記録に解決したか**（引用の軸）
- `claim_support_status` — **取得できた本文が、その主張を何によって支持しているか**（証拠の軸）

実在する論文が遺伝子名を挙げているだけの場合は
`citation_verified=True` かつ `claim_support_status=unsupported_mention` になる。
以前はこれが `verification_status=review_required` に潰れていて、
「引用が怪しい」のか「引用は本物だが主張を支えていない」のかが読めなかった。

`hypothesis_source_status.tsv` には `claim_support_summary`
（`mechanistic_support` / `associative_support_only` / `review_required` /
`unsupported_mention_only` / `background_or_review_only` / `not_evaluated`）と
機構的/関連的/言及のみの引用数が付く。

## 意図的な限界

- **cue照合は判読ではなくトリアージ**。「相補性検定の記述がある」と「相関の報告である」
  は区別できるが、その実験が*この*機構を支持するかは判定していない。cueが立たない場合は
  常に `review_required` で、passにもrejectにも倒さない。
- `retrieval_scope`（`title_abstract` / `title_only` / `full_text_cached` / `none`）を
  記録する。abstractしか見ていない判定を全文判定と同じ顔で出さないため。
- **スコアは一切書かない**。分類が決めるのは「人間が何を読み直すべきか」であって
  候補の価値ではない。

## agent抽出のゲート

`source-verify` は `source_review_records.json` に取得済みのtitle、abstract、
cached full textを元の値・空白・markupのまま保存する。出典ID、本文hash、
同じ論文を引用する全hypothesis IDを保持し、manifestにartifact全体のSHA-256を記録する。
既存の語句照合で支持されなかった論文も、取得できていればレビュー対象として残る。

`available_text_scope` は実際に保存した範囲、`classification_scope` は既存の分類が
記録した範囲であり、両者を混同しない。本文保存によって分類や候補昇格を変更しない。
このartifactの `semantic_binding_state=not_evaluated`、`membership_effect=none` は
明示的な境界であり、**本文があるだけではbound claimは完成しない**。

`--claim-extraction` に agent/人間の抽出（`citation` / `claim_type` /
`evidence_relation` / `quoted_support` / `locator`）を渡せる。ゲートは決定的で、
**引用した文字列が実際に取得済み記録の本文に出現しなければ `rejected`**。

任意の `hypothesis_id` を指定すると、その仮説と `citation` の組だけへ抽出を適用する。
Direct C では `record_id` がこの仮説IDに対応する。同じ論文の異なる引用を別の仮説へ渡す場合は、
各抽出に対象IDを明記する。未知のID、対象が引用していないcitation、同じ組への異なる抽出、
同一citationでの対象指定あり・なしの混在はエラーとして拒否する。同一内容の重複は一つにまとめる。
対象未指定の単一抽出は従来どおり、そのcitationを引用する仮説へ適用する。
引用本文から対象を推測したり、引用を連結したりはしない。対象指定は既存のsource・identity・support gateを迂回しない。

- `rejected_quote_not_in_retrieved_text` — 位置できない引用は、誰が作ったものでも証拠ではない
- `rejected_unknown_citation` / `rejected_unknown_relation` / `rejected_missing_field`

通過した抽出は行に注釈を付けるだけで、**解決しなかった引用を支持済みへ昇格させることは
できない**し、スコアも作らない。

## Route A: 独立agentレビューの取り込み

`hypothesis-verify --agent-claim-review packet.json --source-review-records source_review_records.json`
は、提案とは別実行のagentによる意味レビューを取り込む。`--bound-claims`とは排他的。
既存の`--hypotheses`、`--record-cache`、`--source-records`、`--focal-species`、`--out-dir`も従来通り指定する。
新しいRouteやscoreは作らない。`curation_method=independent_agent_review`をvalidationとbaitへ保持する。

1. `source-verify`のannotated hypothesesと`source_review_records.json`を固定する。
2. `candidate_qtl.agent_claim_review.prepare_review_prompt(hypotheses_path, sources_path, prompt_path)`
   でJSON promptを出力する。このhelperは文章判断もmodel実行もしない。
3. 提案者と**別の実行コンテキスト**でprompt全体を読むagentへ渡し、JSON回答をそのまま保存する。
   AのreviewerにもB/区間候補/truthを見せない。実際のexecutor ID/modelとraw回答を保存する。
4. 次のpacketを作成し、既存の`hypothesis-verify`へ渡す。

```json
{
  "contract": "independent-agent-claim-review-1.0",
  "hypotheses_sha256": "<exact annotated TSV sha256>",
  "source_review_records_sha256": "<source snapshot sha256>",
  "proposal_run": {
    "executor_id": "<actual proposal task/run id>",
    "model": "<actual model>",
    "raw_response": {"path": "proposal-raw.json", "sha256": "<sha256>"}
  },
  "review_run": {
    "executor_id": "<different actual review task/run id>",
    "model": "<actual model>",
    "prompt": {"path": "prompt.json", "sha256": "<sha256>"},
    "raw_response": {"path": "review-raw.json", "sha256": "<sha256>"}
  }
}
```

prompt内に回答schemaと全review軸が含まれる。`reviews`にはhypothesisごとに最大1件を置き、
未読を機械的な空回答で埋めない。未回答は`not_evaluated`。positive quotaは無い。
複数出典から1件を選ぶ判断はreviewerが行い、選んだ出典/箇所/理由を記録する。

機械照合は以下を要求する。

- prompt内のhypothesis rosterとsource snapshotがfreezeされた実入力そのものと一致する。
- citationが取得recordのPMID/PMCID/DOIと一致し、このhypothesisに紐付く。
- abstractまたはfull textの指定文字位置にverbatim quoteが存在する。titleのみでは昇格しない。
- gene token、種、実験摂動、表現型が同じ引用内に存在し、agentが同じ実験の関係と判断している。
- gene tokenが同じprovider recordのaccession/gene name/synonym/locus/xrefに完全一致する。
  quoteにaccession自体が無いsymbol対応は、そのprovider recordから当該論文へのcitation linkも要求する。
- organismが同じprovider recordのscientific/common nameと一致する。
- 既存のsupporting predicateと全意味review軸が満たされる。不一致/不足をbiological missにしない。

`agent_claim_review_audit.json`には元review、provider record、入力/回答hash、方式、
executor provenanceを保存する。意味判断はagentの明示的なattestation（レビュー結果）であり、
語句一致が関係を証明したとは主張しない。executor IDはorchestratorが記録するもので、
model identityや隔離の暗号学的証明ではない。proposal rawは履歴であり、意味的正しさの証明ではない。
評価結果に合わせたreviewの書き換えをしない。

checked claimは実行中だけtyped receiptを持つ。method名や`verified_complete`を記したTSVだけを
`--bound-claims`へ渡しても受理しない。再利用は元packet/source/hypothesesを再検証する。
`bait-panel-freeze`はagent方式の場合、隣接verification manifestとaudit、hypothesis/validation hashも確認する。
旧human/deterministic方式やDirect Cのゲートは変更しない。

`guided-profile --profile triple-route --agent-claim-review ...`もA04へ接続する。
これは**正確に同じsource/hypothesis入力に対して事前実施済みのreview**を供給する方式で、
profile内でagentを自動起動するものではない。取得時刻等で入力bytesが変われば再利用を拒否するため、
通常は固定した`source-verify`出力から既存stageを再開する。`audit-predict`等への無視された入力にはしない。

## Route B: 同じ独立レビュー方式の接続

原文・入力を固定して実レビューを受け取り、既存のA/B候補生成とWetへ渡す操作は
[固定入力からのレビュー再開](CLAIM-REVIEW-SESSION.ja.md)を参照。

`prepare_review_prompt(proposals, sources, prompt, route="b")` は凍結済み
`semantic_claim_proposals.tsv` と B source-verify の `source_review_records.json` を渡す。
packet contractは`independent-agent-route-b-claim-review-1.0`、回答の主キーは`claim_id`。
提案者とは別contextで実際にレビューし、raw回答とproviderレコードを保持する。

既存`route-b-semantic-transfer-verify`へ`--agent-claim-review`、
`--source-review-records`、`--record-cache`を指定する。確認方式はAと同じで、
citation linkはモデルの申告ではなくproviderの実レコードから取得する。
種名が欠ける場合は同一reference endpointの凍結metadata、または一致するtaxon IDでのみ補完する。
既存のreference-to-interval transfer gateと全copy保持は変更しない。
意味レビューは順位やスコアを作らず、最終Wet集合の構成員となる根拠を確認するだけである。

guided-profileではAのpacketと混同しないよう`--route-b-agent-claim-review`と
`--route-b-source-review-records`を使用する。profileはagentを自動起動しない。
再利用する場合も元proposal/source/prompt/raw回答を再検証し、
`agent_claim_review_audit.json`とtransfer manifestへhash・provider・根拠を記録する。

## ベンチマーク

```bash
candidate-qtl claim-verification-benchmark \
  --cases cases.tsv --cache <europe pmc cache> --cache-is-complete --out-dir <dir>
```

`cases.tsv` は `case_id` / `citation` / `expected_claim_support_status`
（＋任意で `reference_gene_name` / `gene_family` / `pathway_or_process` /
`expected_citation_verified`）。真の機構的支持・関連のみ・無関係な言及・捏造引用の
4クラスを別々に採点し、クラス別accuracyと混同表を出す。

- 記録を取得できなかったケースは `not_evaluated` で**分母から外れる**。
  providerが無いことは誤分類ではない。
- 捏造引用ケースをオフラインで採点するには `--cache-is-complete`
  （「このcacheはこれらの引用について完全である」という運用者の明示的宣言）が要る。
  宣言なしのcache missは `not_evaluated` に倒す。

## 非目標

実在する論文に遺伝子名が共起していることを機構的支持として扱わない。
