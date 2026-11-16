# 次実験選択（action priority）

`next-experiment`は、生物学的なcandidate順位を変更せず、「次にどの実験を行うと候補を効率よく減らせるか」を別artifactへ出すrule-based decision supportである。原因遺伝子確率、Bayesian information gain、外部発注の自動決定ではない。primary recommendationを実行する前に人が承認する。

同じ選択器は`case-run`のdeterministic rankとprimer planningの後にも自動実行され、`next_experiment/`へ保存される。LLM-assisted rankはaction選択へ混ぜない。

## 比較する経路

- additional grid markers / recombinant genotyping
- Top-N parental Sanger
- multiplex amplicon sequencing / target capture
- SV・TE・CNV・PAV branch
- expression assay
- annotation rescue（比較annotation投影でreference-missing locusを補完）
- non-gene follow-up（linked QTL・低組換え領域・association/haplotype・マルチ遺伝子座機構）

細粒度markerと組換え個体で区間を先に狭める設計は、候補遺伝子を逐次読むより少ないgenotypingで候補集合を減らせる場合がある。QIRによるfine mappingの実例ではrecombinant-focused genotypingが作業量を抑え、rice `qLTG3a`の同定でも追加markerとrecombinant breakpointによる区間圧縮後に候補解析・発現比較へ進んでいる。

- [Quantitative interval and recombinant-focused mapping](https://pmc.ncbi.nlm.nih.gov/articles/PMC1456828/)
- [Rice qLTG3a fine mapping and functional analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8674402/)

## 新規 action class

- **annotation_rescue**: hypothesis ledger に reviewed かつ manual/workflow origin の
  `annotation_missing` 仮説がある場合（LLM仮説はpending_reviewのまま除外）、または
  仮説の `required_experiment_ids` が annotation rescue を要求する場合に優先される。
  比較annotation資源（reference proteome等）が必要
- **non_gene_follow_up**: mapping gate の status が `linked_qtl_suspected` /
  `low_recombination_region`、または reviewed かつ manual/workflow origin の
  `unknown_locus_mechanism` / `linked_qtl` / `multi_gene_combination` 仮説がある
  場合に優先される

mapping refinement の token マップには `fine_mapping_required` と
`linked_qtl_suspected` を追加した（どちらも localization不足 → 区間refinement）。
`low_recombination_region` は追加しない。marker追加が解決にならないためで、
non-gene follow-up 枝のreview信号として残す。

## 入力と制約

case manifestの`decision_support`へ`budget_class`、`allowed_technologies`、`available_samples`、必要なら実施設定に基づく`local_cost_classes`を記載する。価格情報がなければ`unknown`のまま保持し、数値やlow/medium/highを推測しない。

既存実験はassay outcome表で次の2軸を分ける。

- state: `not_started` / `partial` / `failed` / `complete`
- outcome: `technical_failure` / `informative_negative` / `positive` / `inconclusive` / `unknown`

`technical_failure`は仮説の反証ではない。`informative_negative`も、測定範囲外のregulatory/SV原因や未測定candidateまで除外しない。

### candidate-elimination 後の remaining set

`--remaining-candidates` に candidate-elimination の canonical
`remaining_candidate_set.tsv` を渡すと、selector は candidate count を
**elimination 後の生存者数**に切り替えます（`candidate_count` /
`top_or_reserve_count` は remaining set と ranked table の intersection から再計算）。
`input_facts` には `candidate_count_before_elimination` と
`candidate_count_after_elimination` の両方を残し、決定JSONの入力指紋にも
remaining set が載ります。

生存状態のもう1つの消費者として、`integrated-narrowing` も同 canonical
`remaining_candidate_set.tsv` を `--remaining` で受けます（設計・tier規則・
現在の配線と残る限界は [INTEGRATED-NARROWING.ja.md](INTEGRATED-NARROWING.ja.md)）。

- 渡さない場合は**従来挙動のまま**（byte-stable）。
- remaining set が空の場合は `remaining_candidates_empty` の review を出し、
  勝手に候補を復活させません（「原因遺伝子が無い」とは断定しない）。
- remaining set に ranked table に無い gene_id があれば fail-closed で
  ValueError（silent drop しない）。
- 生物学的 rank / `dry_score` は一切書き換えません。

## 実行例

```powershell
candidate-qtl next-experiment `
  --ranked results\case01\rank\ranked_candidates.tsv `
  --case templates\case_manifest.yaml `
  --out-dir results\case01\next_experiment `
  --mapping-gate results\case01\mapping_gate.json `
  --assay-outcomes results\case01\assay_outcomes.tsv
```

## 出力

- `next_experiment_options.tsv`: candidate rankとは独立した`action_rank`、feasibility、ordinal cost/information、rules、入力facts
- `next_experiment_decision.json`: primary、alternatives、policy、入力指紋、`human_approval=pending`
- `experiment_dependencies.tsv`: 各optionのblocking dependency
- `experiment_review_required.tsv`: 不明なsample、cost、technology、coverage等

順位はhard dependency、feasibility、ordinal information、cost/risk/runtime/reuseのPareto関係、stable tie-breakの順で決める。同じ入力はbyte-stableである。candidate scoreや`dry_group`には費用・機器可否・LLM出力を混ぜない。

## ベンチマーク

`next-experiment-benchmark` は frozen な per-case artifact に対して
`select_next_experiment` を実行し、推奨の escalation を採点する。

```bash
PYTHONPATH=src python3 -m candidate_qtl.cli next-experiment-benchmark \
  --manifest benchmarks/next_experiment_v1/cases.tsv --out-dir bench_out
```

- `exact_action_accuracy`: primary recommendation が expected action class と一致
- `appropriate_escalation_accuracy`: expected または manifest の
  `acceptable_alternatives`（|区切り）のいずれか
- `expected_blocked_technology` のcaseは exact/appropriate の分母から除外し、
  `blocked_detection_rate` として別報告する（blocked≠inappropriate）
- manifest で指定した入力ファイルが存在しないcaseは `not_evaluated` として
  全分母から除外・明示報告する（必須sample/provider欠落は experiment class への
  反証にしない）
- 特別指標: `metric_group=repeated_informative_negatives` での coding-gene
  continuation 回避率、`metric_group=weak_mapping_override` での
  `additional_grid_markers` 選択率

出力は `next_experiment_benchmark_cases.tsv` / `_summary.tsv` /
`_manifest.json`（入力fingerprint・語彙・policy文）/ `_report.md`。policy は
「推奨はordinalでhuman-reviewableのまま、causal probability・正確なcostを捏造
しない」ことを manifest に記録する。コミット済みfixture
（`benchmarks/next_experiment_v1/`）はsyntheticなcanaryであり、checksum検証は
`tests/test_next_experiment_benchmark.py` が行う。
