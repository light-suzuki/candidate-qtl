# QTL effect context contract

Issue #360 のための、QTLで観測されたアレル効果と候補仮説の機構方向を分離して記録する監査用contractです。

このcontractは **ranking scoreではありません**。`incompatible` であっても候補を自動除外・自動減点せず、QTL effect と機構仮説の矛盾を review 可能な状態として残します。

## 入力contract

YAMLではcase manifest内の `qtl_effect_contract` として、TSVではcaseごとの1行として次の列を使用できます。

| field | 意味 | 未指定時 |
|---|---|---|
| `case_id` | case識別子 | 呼び出し側のcase_id、なければ `unknown` |
| `favorable_parent` | 有利アレルを持つ親 | `unknown` |
| `effect_direction` | QTL効果方向（例: `increase`, `decrease`, `gain`, `loss`） | `unknown` |
| `additive_effect` | 加算効果 | `unknown` |
| `dominance_class` | 優性様式 | `unknown` |
| `effect_ci` | effectの信頼区間 | `unknown` |
| `environment_id` | 環境・試験識別子 | `unknown` |
| `phenotyping_protocol` | 表現型測定法 | `unknown` |
| `stable_across_environments` | `true` / `false` / `unknown` | `unknown` |
| `source_id` | 出所の安定ID | `unknown` |
| `source_detail` | table/figure/local artifact等の補足 | `unknown` |

既存caseにこのcontractが無くても有効です。欠損はnegative evidenceに変換せず、明示的な `unknown` になります。

### YAML例

```yaml
case_id: pea_height
qtl_effect_contract:
  favorable_parent: ParentA
  effect_direction: increase
  additive_effect: "3.2"
  dominance_class: additive
  effect_ci: "2.1..4.3"
  environment_id: field_2026
  phenotyping_protocol: mature stem height
  stable_across_environments: true
  source_id: PMID:12345678
  source_detail: Table 2
```

### TSV例

```text
case_id	favorable_parent	effect_direction	additive_effect	dominance_class	effect_ci	environment_id	phenotyping_protocol	stable_across_environments	source_id	source_detail
pea_height	ParentA	increase	3.2	additive	2.1..4.3	field_2026	mature stem height	true	PMID:12345678	Table 2
```

複数caseを1つのTSVに入れる場合、loaderには `case_id` を渡して1行へ解決します。該当caseが無い、または1caseに複数行ある場合はfail-closedです。

## Sign concordance

`candidate_qtl.qtl_effect_contract.annotate_sign_concordance()` は、hypothesis側の `expected_mechanism_sign` とcontractの `effect_direction` を比較して次を付加します。

- `sign_concordance`: `concordant` / `incompatible` / `unknown`
- `sign_concordance_reason`: 正規化した両signを記録
- `qtl_effect_source_id`
- `qtl_effect_source_detail`
- `qtl_effect_environment_id`

`increase/gain/positive/+` はpositive sign、`decrease/loss/negative/-` はnegative signとして比較します。どちらかをsignへ正規化できない場合は `unknown` です。

この処理は入力DataFrameのコピーへ監査列を追加するだけで、`dry_score`、rank、candidate集合には触れません。

## case-run への配線

contractは `case-run` のoptional inputです。case manifestに次のように登録すると、`resources.yaml` 経由で解決されます。

```yaml
inputs:
  qtl_effect_contract_resource: qtl_contract
```

- inputが未宣言の場合: 全列 `unknown` のcontractとして扱い、既存caseはそのまま有効です（backward compatible）。
- inputが宣言済みだが解決できない、またはYAMLのtop-level/nested `case_id` がcaseと一致しない、header-only TSVでrequested caseの行が無い場合は fail-closed（`ValueError` でrunが失敗）です。他caseのQTL効果を静かに流用しません。
- contractデータは `prioritize_evidence` やrank/evidenceのdigestに渡りません。`rank_contribution: none` です。

### 生成されるaudit artifact

| パス | 内容 |
|---|---|
| `audit/qtl_effect_contract.tsv` | canonicalな1行contract（contract提供時のみ） |
| `hypotheses/locus_hypothesis_sign_audit.tsv` | `locus_hypotheses.tsv` のコピーにsign判定列を付加した監査表（ledger生成時は常に出力） |

`locus_hypotheses.tsv` 自体は書き換えません（schemaと `_stable_id` をbyte-stableに保つ）。auditの `expected_mechanism_sign` はhypothesis proposal入力の同名列から `hypothesis_id` で突合します。列順はcanonical（ledger列 + 監査6列）で、`hypothesis_id` のstable sortです。

### manifest / case state への記録

`run_manifest.json` と `case_state.json` に `qtl_effect_contract` ブロックが入ります。hypothesis stepのstateにも同ブロックが記録されます。

```json
{
  "status": "used",
  "artifact": ".../audit/qtl_effect_contract.tsv",
  "fingerprint": "<contractのsha256>",
  "sign_audit": ".../hypotheses/locus_hypothesis_sign_audit.tsv",
  "sign_concordance_counts": {"concordant": 0, "incompatible": 1, "unknown": 0},
  "sign_concordance_review_required": true,
  "rank_contribution": "none"
}
```

- `status`: contract使用時 `used`、未提供時 `missing`（`artifact`/`fingerprint` はnull）。
- `sign_concordance_review_required` は `incompatible` が1件以上あるとき `true` です。**review可視のみ**で、候補の自動除外・減点・`reviewed` 化はしません。`incompatible` なhypothesisは `locus_hypotheses.tsv` に完全なまま残ります。
- contractのsha256は `hypothesis_hash` のdigest configに含まれ、contract artifactとsign auditは `final_hash` のartifact listに含まれるため、contract編集でhypothesis/ledger/reportのstale cacheは失効します。`rank_hash` / `evidence_hash` には一切含めません。

## 現在の配線範囲

`case-run` への配線は上記のaudit専用の形で完了しています。contractが無い旧caseでは `expected_mechanism_sign=unknown` / `sign_concordance=unknown` / provenance unknown となり、有効な状態として扱われます（エラーではありません）。
