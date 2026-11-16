# Locus causal scenario ledger

Issue #362 のための、gene rankingとは独立してlocusの因果形式を並列管理する `locus_causal_scenarios.tsv` contractです。

A-only / B-only / neither を単純なcandidate missへ潰さず、coding、regulatory、SV/PAV、reference-parent mismatch、mapping/liftover、multi-gene、environment-specific といった別scenarioを同一case内で共存させます。

このledgerは **ranking tableを上書きしません**。全rowで `rank_contribution=none` です。

## Scenario type

許容する最小scenario type:

- `coding_gene`
- `regulatory_element`
- `sv_te_cnv_pav`
- `reference_parent_mismatch`
- `mapping_or_liftover_error`
- `multi_gene_or_epistasis`
- `environment_specific`
- `novel_or_unannotated_mechanism`
- `phenotype_context_review`
- `mapping_adequacy_review`
- `provider_missing`

未知typeはfail-closedです。

### #363 由来のreview category (schema拡張の理由)

`novel_or_unannotated_mechanism`、`phenotype_context_review`、`mapping_adequacy_review` の3種はIssue #363のA/B discordance adapterがreview候補として出すcategoryです。既存type (`multi_gene_or_epistasis` / `environment_specific` / `mapping_or_liftover_error`) へ**潰さず**独立typeにしています。理由:

- `novel_or_unannotated_mechanism` は既知のmulti-gene仮説ではなく「未知・未注釈機構のレビュー」であり、`multi_gene_or_epistasis` (既知の複数遺伝子仮説) とは情報量が異なる。
- `phenotype_context_review` は表現型測定の妥当性レビューであり、環境要因の仮説 (`environment_specific`) ではない。
- `mapping_adequacy_review` はマッピング妥当性のレビューであり、mapping/liftoverエラーの断定 (`mapping_or_liftover_error`) ではない。

`provider_missing` はroute provider欠落でdiscordanceの四状態を判定できなかったrow専用です。どの因果形式も名指ししないため、いずれのreview categoryにも丸めません (reasonに「route_X provider missing; discordance cannot be classified」を持つ `not_evaluated` rowになります)。

## Canonical columns

`locus_causal_scenarios.tsv` は次の列を持ちます。

- `case_id`
- `scenario_id`
- `scenario_type`
- `supporting_evidence`
- `contradicting_evidence`
- `assessment_state`
- `recommended_next_experiment`
- `reason`
- `source_ids`
- `rank_contribution`

`scenario_id` が未指定なら、case/type/evidence/next experiment/reasonからdeterministicに生成します。同一caseに複数scenarioを保持できます。

## Assessment state

- `unassessed`
- `supported`
- `contradicted`
- `mixed`
- `review_required`
- `not_evaluated`

supporting evidenceもcontradicting evidenceも無いrowは必ず `unassessed` または `not_evaluated` のままです。証拠が無いことをcontradictionへ変換しません。

### `not_evaluated` と `unassessed` の区別

- `unassessed`: 評価の試行自体がまだ無い状態。
- `not_evaluated`: provider欠落・意味論的に評価不能のため、評価の結果を作ることができなかった状態 (A/B discordance adapterのprovider gap row)。

両者は別stateです。`not_evaluated` は証拠なしで成立しうるため、`normalize_scenario` のguardは「証拠なし+非{`unassessed`, `not_evaluated`} → ValueError」です。`scenario_summary()` は `not_evaluated_rows` を `unassessed_rows` と分けて集計します。

## A/B discordance adapter

`discordance_scenarios.discordance_review_rows()` は `dual_route_hypothesis_comparison.tsv` / `triple_route_hypothesis_comparison.tsv` の各rowを四状態に分類し、review候補scenario rowへ決定論的に変換します。出力はそのまま `build_scenario_ledger()` へ渡せます (ledger側で `rank_contribution=none`、sort、dedup)。

| 四状態 | review候補scenario_type |
|---|---|
| A∩B | `coding_gene`, `regulatory_element` |
| A-only | `mapping_or_liftover_error`, `reference_parent_mismatch`, `sv_te_cnv_pav` |
| B-only | `coding_gene`, `novel_or_unannotated_mechanism`, `multi_gene_or_epistasis` |
| neither | `phenotype_context_review`, `mapping_adequacy_review`, `environment_specific` |
| route provider欠落 | `provider_missing` 1本 (`not_evaluated`) |

四状態のrowは全て `assessment_state=review_required` で止まり、reasonに判定理由 (`four-state=<state>; judgment: ...`) を出力します。candidateの自動除外・自動昇格はしません。same-root非独立収束 (shared source ids) は分類を変えず、reasonとsource_idsにsource dependenceを明記します。provider欠落 (`route_a_evaluated=False` / `route_b_evaluated=False`) のrowはneither bucketへ混ぜず、`not_evaluated` 1本 (reason: `route_X provider missing; discordance cannot be classified`、recommended_next_experiment: `provide route X provider or record route X as intentionally absent`) として分離します。`input_fingerprint` は呼び出し側が計算して渡し、source_idsに `comparison=<id>;fingerprint=<fp>;route_a_source_ids=...;route_b_source_ids=...` として保存します。

## 作成

`build_scenario_ledger()` はmappingの列からcanonical DataFrameを作ります。`write_scenario_ledger()` はそのまま `locus_causal_scenarios.tsv` を書き出します。

最小例:

```python
rows = [
    {
        "scenario_type": "coding_gene",
        "supporting_evidence": "candidate coding variant",
        "assessment_state": "supported",
        "recommended_next_experiment": "sequence candidate CDS",
        "reason": "coding hypothesis remains plausible",
    },
    {
        "scenario_type": "sv_te_cnv_pav",
        "assessment_state": "unassessed",
        "recommended_next_experiment": "inspect long-read/SV evidence",
        "reason": "SV/PAV space has not been evaluated",
    },
]
```

## 表示

`scenario_summary()` はscenario typeごとに:

- scenario数
- supported rows
- contradicted rows
- unassessed rows
- supporting evidenceが実際に記録されたrow数
- contradicting evidenceが実際に記録されたrow数

を分けて集計します。

`render_scenario_report()` はsupport / contradiction / assessment / next experiment / reasonをscenarioごとに表示しますが、scenario間のwinnerを選びません。

## 現在の配線範囲

このPRではschema、deterministic作成、TSV writer、summary、Markdown表示に加え、`case-run` への自動配線を追加します。

- 新step `causal_scenarios` は next_experiment の後・report の前に実行され、常に `<out-dir>/locus_causal_scenarios.tsv` を出力します。scenario rowリソース (case.yaml の `inputs.scenario_rows_resource`、TSVまたはYAML list-of-mapping) が無い旧caseでも、ヘッダのみのledgerと report の「No causal scenarios have been registered; status is unassessed.」で未評価を誠実に示します。
- このstepはranking tableを読み書きしません。全row `rank_contribution=none` のまま、deterministic rank・Stage-A・scoreへは一切干渉しません。
- #360 qtl_effect_contract は未mergeのため、このstepからは該当モジュールをimportしません。任意入力 `inputs.qtl_effect_contract_resource` (YAML `qtl_effect_contract:` キー配下またはトップレベルlist、あるいはTSV) を防御的に読み、該当caseの `sign_concordance` が `incompatible` なら各scenario rowの `contradicting_evidence` に、`concordant` なら `supporting_evidence` に `sign discordance/concordance: favorable_parent/effect_direction` を追記します (source_idsにも追記)。ファイル欠損・parse失敗・`unknown` はcrashさせず、`run_manifest.json` の `scenario_ledger_status.sign_contract_status` (`present|absent|invalid|unknown`) に記録するだけです。
- `run_manifest.json` に `scenario_ledger_status` (rows_input / sign_contract_status / scenario_count)、`case_state.json` のfinalに `scenario_ledger` / `scenario_count` を追加します。stepは `cached()`/`complete()` の通常パターンでresumeキャッシュされ、digestはscenario rows・contractファイルのfingerprintとcase_idを含みます。

### A/B discordance adapter のcase-run配線

- 任意入力 `inputs.dual_route_comparison_resource` (A/B比較TSV) が登録されている場合のみ、`causal_scenarios` stepは `discordance_review_rows()` の出力rowをledgerの構築rowsへ追加します。追加rowの `source_ids` には `comparison=<comparison_id>;fingerprint=<SHA-256 hex>` が入ります。fingerprintは `discordance_scenarios.comparison_fingerprint(path)` が計算します: ファイル内容bytesのSHA-256 (64桁hex) で、データ行は辞書順に正規化 (ヘッダ行は先頭固定) してからhashするため、同じ比較行の並び替えではfingerprintもledgerも変わりません。
- `inputs.dual_route_route_a_manifest_resource` / `inputs.dual_route_route_b_manifest_resource` を登録すると、route status は明示artifactから**fail-closed**に読みます。Route A は `phenotype_first_status.json` (`route=phenotype_first`, `verified_bait_count >= 1`) を、Route B は `interval_first_status.json` (`route=interval_first`, `semantic_status=evaluable`) を評価条件にします。`not_evaluated*`、route不一致、壊れたJSON/TSV、必要キー欠落は該当routeを評価不能として、全比較rowを `provider_missing` / `not_evaluated` の1本に潰します。従って provider欠落は `B-only` / `neither` という生物学的主張へ誤分類されません。旧precomputed dual comparisonでstatusも直接route tableも無い場合だけは既存互換のlegacy inferenceを維持します。
- canonical `guided-profile --profile triple-route` は `triple_route_hypothesis_comparison.tsv` と両routeのstatus JSONを `case-run` のoverrideとして渡します。このためtriple-route comparisonもcausal scenario ledgerへ到達します。比較・statusはreview queueの入力に限られ、rank / score / Stage-Aを変更しません。
- 比較TSVの読み込み・分類に失敗した場合はcrashせず、`discordance_input=invalid` を記録し、reportに `discordance input invalid` の警告行を追加します (scenario_rows invalidの流儀)。
- step digestはscenario rows・contractに加え、比較TSV (存在時) とroute A/B manifest (存在時) のpath/fingerprintを含みます。resume時に同じstep hashなら `cached` でbyte-identicalのledgerが再利用され、**ledgerは毎run、全入力 (scenario_rows + sign contract + discordance) からゼロから再構築**されるため、duplicate scenarioは構造的に増殖しません。
- `complete()` のextraと `scenario_ledger_status` に `discordance_input` (`present|absent|invalid`)、`discordance_rows` (追加row数)、`route_a_status` / `route_b_status` (`absent|evaluated|not_evaluated|unknown`) を記録します。reportの `Ledger summary` 行に discordance 由来のscenario数 (`discordance-derived: N (present|absent)`) を明記します。比較リソースが無い旧caseの挙動 (ledger・manifest・report) は変わりません。
