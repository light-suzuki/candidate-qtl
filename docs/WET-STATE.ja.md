# Wet assay state と陰性結果の管理

`wet-state` は、assayの進捗とcandidate/mechanismの解釈を分離して正規化する決定論的な段階です。LLMは使いません。`completed`、非増幅、部分的な配列確認を原因性や真の陰性へ読み替えません。

## 入力と出力

canonical入力templateは次の3表です。

- `templates/wet_assay_registry.tsv`: gene/transcript/hypothesisとplanned region、assay state
- `templates/wet_assay_target_coverage.tsv`: planned regionに対する実測covered region、coverage、quality、read direction
- `templates/wet_results.tsv`: parent A/B、F1、F2/F2:3、recombinantのindividual-level result

出力は `assay_registry.tsv`、`assay_target_coverage.tsv`、`wet_results.tsv`、`candidate_wet_state.tsv`、`wet_state_history.tsv`、`wet_review_required.tsv` の6表です。`wet_state_history.tsv`は前回履歴を`--previous-history`で読み、許可された遷移を追記します。禁止遷移も`rejected_forbidden_transition`として消さずに残します。

```powershell
candidate-qtl wet-state `
  --assays templates\wet_assay_registry.tsv `
  --coverage templates\wet_assay_target_coverage.tsv `
  --results templates\wet_results.tsv `
  --out-dir results\wet_state
```

## 安全な解釈規則

- PrimerBlastとexon tilingのcoverageは設計上の予測であり、Wet coverage gateを通しません。
- `coding_region_cleared`には、宣言した全CDS/splice assayの完了、実測coverage、quality pass、既定ではbidirectional read、parent A/B双方のreview済み陰性sourceが必要です。
- coding regionをclearしてもpromoter、UTR、regulatory、SV/PAV hypothesisは残します。gene全体のrejectにはなりません。
- assay failureや片親非増幅は陰性ではありません。primer mismatch、SV、PAVを再検討するreview itemになります。
- variant confirmationと遺伝学的validationは、reviewed、quality pass、source付きでなければ状態を昇格しません。
- F2/F2:3/recombinant resultには`individual_id`が必須です。単純な共分離支持は`causal`にせず`requires_independent_validation`に留めます。
- `fails_cosegregation`と`excluded_by_recombinant`は別状態です。

`case-run`ではcanonical Wet resourceまたはPrimer/exon resultが指定された場合にこの段階を実行し、正規化済み`candidate_wet_state.tsv`を`next-experiment`のSV/PAV branch判定へ渡します。
