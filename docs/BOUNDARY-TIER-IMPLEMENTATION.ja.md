# Boundary Tier review: implementation and integration status

採否規則の正本は [CANONICAL-WET-MEMBERSHIP.ja.md](CANONICAL-WET-MEMBERSHIP.ja.md) です。
この文書はAPIの責務と未接続箇所の説明であり、別のTop20 producerを正本にしません。

## 実装したAPI

`boundary_tier_review.partition_relevance_tiers` は正規化済みの全gene rosterを
frozen / boundary / lowerの3範囲に分けます。同Tierの一部をgene ID順で採用する処理はありません。
IDによるsortはserializationだけです。既にレビューを終えた上位prefixを渡すと下位Tierへ継続できます。

`freeze_boundary_scope` は元inventoryのSHA-256とcaseを必須にし、範囲表とmanifestを新しいdirectoryへ保存します。
元inventoryや過去世代は上書きしません。`_load_scope` はhashだけでなく元inventoryから再計算した全行との一致を確認します。

`consume_boundary_selection` は既存の
`adaptive_synthesis.consume_review_inclusive_selection`に**境界Tierだけ**を渡します。
独立レビュー・executor/reviewer identity・proposal/source hash・deferの検証は既存関数に委譲します。
上位/下位のgeneをreview行へ混ぜた入力はconsumerを呼ぶ前に拒否します。
戻った境界行は元の全列と比較し、上位行は変更せず合流します。

レビュー後にtarget未満なら下位Tierを丸ごと進め、次の境界があれば`next_scope_manifest`として明示します。
未解決・不一致は既存契約に従って保持し、数だけを理由に切り捨てません。
これはLLM実行ではなく、既存のレビュー結果を消費するAPIです。

## 通常handoffへの接続

通常の `wet_candidate_handoff` は、全区間rosterを入力にしたrelevance selectionでこのAPIを呼びます。
第1段の `wet_candidates.tsv` と `wet_candidate_trace.tsv`、boundary scope、review requestを別々にfreezeし、
`claim_review_session.finish_wet_review` がhashを検証して境界選定を消費します。
finish後は採用結果と未採用lower Tierを次のscopeへ渡し、round budget内で再handoffできます。
これはモデル実行を意味せず、実際のproposal/reviewer結果が無い状態は未完・not_evaluatedとして残ります。

統合後も次の境界は残ります。

1. strict capacity、`supported_only`、明示的refinementは互換性のためlegacy diagnosticとして維持する。
2. `wet_candidate_trace.tsv`はhandoff observableの監査記録であり、LLM性能やcausal truthを単独で判定しない。
3. Type3の最終Wet20 retentionとactual unique gene countは、first-stage retentionやexpanded-set retentionと別に評価する。

## 実行した検証

新規17テストは、3+8+35、21〜25保持、T1=31保持、T4/U救済、全copy保持、入力順/ID置換不変、
上位行改変検出、hash/case不一致、独立レビューなしの除外拒否、T3再検討後にT4へ進む2世代を検査します。
監査時にsource抽出とpath shimを除去し、既存の実selection consumerを直接呼ぶ試験へ変更しました。
通常session/CLIからこのadapterへの接続は別の受入条件です。

通常handoffへの接続を含むCORE_INTEGRATIONです。実モデル性能の主張やbenchmark再開は行いません。benchmark pauseは維持します。
