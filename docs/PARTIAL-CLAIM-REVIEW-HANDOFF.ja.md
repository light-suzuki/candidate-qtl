# A/B独立レビューの部分完了からWetへ渡す

`claim_review_session.handoff_sessions()`（`scripts/claim_review_session.py handoff`）は、凍結済みの`completed`と`completed_partial_review`を受け取る。未回答のclaimを支持に変更せず、既存のverification/member gateが生成したsupported候補とreview候補を、そのまま既存Wet consumerへ渡す。

## 出力

`inputs.json`、第一段の`completion.json`、`wet_candidate_handoff_manifest.json`には`claim_review`を保存する。

```json
{
  "state": "partial",
  "routes": {
    "a": {"state": "completed_partial_review", "requested": 10, "returned": 9, "missing": 1},
    "b": {"state": "completed", "requested": 5, "returned": 5, "missing": 0}
  }
}
```

部分完了で通常の`wet_reviewable`が返ったときは`wet_reviewable_degraded_coverage`へ変更する。既にreview-required等の状態なら上書きしない。`first_stage_frozen`は候補集合を凍結したという意味であり、全claimのレビュー完了や原因遺伝子同定を意味しない。

件数は非負整数で`returned + missing == requested`が必要。stateとmissing数の矛盾、failed/awaiting状態、request/出力hashの不整合は拒否する。providerや後段consumer自体の障害を成功へ変換する変更ではない。

## 保持する境界

元のrequest、completion、raw responseは変更しない。回答の差替えや、既存completionへの追加レビューをこの操作では行わない。不足分を別attemptとして追記するAPIは未実装。

候補のidentity/source/projection/member gate、support policy、20件の運用目標、25件超の第二段要求は変更しない。未回答をverified supportにせず、候補数で切断しない。候補TSVを編集せず、handoffの実行coverageだけを追記してから第一段のhashを固定する。

関連: Issue #714。回帰テスト: `python -m unittest tests.test_partial_claim_review_handoff -v`。このテストはdownstreamのBLAST/biology consumerをmockしたorchestration検証であり、実BLASTや生物学的性能の評価ではない。
