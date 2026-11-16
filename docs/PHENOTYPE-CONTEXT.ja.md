# Structured phenotype context

Issue #361 のための、自由文 phenotype を置き換えずに発達段階・組織・環境・測定法などを機械可読に保持する監査用contractです。

このcontractは **ranking scoreではありません**。構造化情報の有無を候補の加点・減点に使わず、Route A/C がどの文脈の phenotype statement を扱っているかを追跡するための情報として保持します。

## 後方互換

既存caseが `phenotype_text` しか持たなくても有効です。その場合、自由文はそのまま保存し、構造化項目は `unknown`、`term_usage=raw_text`、`structured_context_status=structured_context_missing` として扱えます。欠損をnegative evidenceへ変換しません。

## `phenotype_context.yaml`

```yaml
case_id: flowering_case
trait:
  phenotype_text: delayed flowering
phenotype_context:
  trait_ontology_id: TO:0002616
  trait_ontology_term: flowering time trait
  unmapped_term:
  developmental_stage: reproductive transition
  tissue: shoot apical meristem
  environment_or_treatment: long day
  measurement_method: days to first flower
  comparison_or_control: wild type
  effect_direction: increase
  replicate_info: 3 biological replicates
  source_id: PMID:12345678
  source_detail: Figure 2
```

canonical fields:

| field | 意味 | 未指定時 |
|---|---|---|
| `case_id` | case識別子 | 呼び出し側case_id、なければ `unknown` |
| `phenotype_text` | 元の自由文 | `unknown` |
| `trait_ontology_id` | ontology ID | `unknown` |
| `trait_ontology_term` | normalized ontology term | `unknown` |
| `unmapped_term` | ontologyへ安全に写せない原語 | `unknown` |
| `developmental_stage` | 発達段階 | `unknown` |
| `tissue` | 組織 | `unknown` |
| `environment_or_treatment` | 環境・処理 | `unknown` |
| `measurement_method` | 測定法 | `unknown` |
| `comparison_or_control` | 比較対象・control | `unknown` |
| `effect_direction` | 効果方向 | `unknown` |
| `replicate_info` | replicate情報 | `unknown` |
| `source_id` / `source_detail` | provenance | `unknown` |

## Term usage

loaderはontologyを検索・推論しません。与えられた情報だけから、どの表現を使用しているかを次のように記録します。

- `normalized_term`: ontology ID または normalized term が入力済み
- `unmapped_term`: ontologyに強制変換せず原語を保持
- `raw_text`: 自由文のみ
- `unknown`: phenotype表現自体が未指定

normalized term が無いからといって `unmapped_term` を自動生成したり、自由文からontology IDを推測したりしません。

## Context fingerprint

同じ `phenotype_text` でも、developmental stage、tissue、environment/treatment、measurement method、comparison/control、effect direction、replicate情報が異なれば別の biological context として識別できるよう、これらの意味的フィールドから deterministic SHA-256 fingerprint を生成します。

`case_id` と `source_id/source_detail` はfingerprintから除外します。同じ生物学的文脈を別case・別出所で記録した場合に、provenance差だけで別contextへ分裂させないためです。

## Structured-context status

`structured_context_status()` はreport/manifestへ載せやすい形で以下を返します。

- `status`: `structured_context_used` / `structured_context_missing`
- `term_usage`
- `present_fields`
- `context_fingerprint`
- `rank_contribution=none`

## 配線（CLI・manifest）

`phenotype_context.yaml` は Route A/C の実フローに「監査記録」として配線されています。自由文は常に保持され、構造化フィールドが自由文を置き換えたり、クエリ・プロンプト・retrieval入力に注入されたりすることはありません。

### CLI

- `guided-profile` / `retrieval-run` に任意オプション `--phenotype-context <path>` を追加（存在しないパス・不正なYAML root・case_id不整合は `ValueError` で fail closed）。
- `retrieval-run` では既定の lexical provider のみがブロックを記録します。`--provider lsa/neural` と併用すると fail closed します（黙って記録を落とさないため）。
- コンテキストファイル未指定でもフローは通常通り実行され、manifest には `structured_context_missing` が記録されます。

### manifest の `structured_context` ブロック

以下のmanifestに同一形式の `structured_context` ブロックが書かれます:

- `phenotype_investigation_manifest.json`（Route A: a01 phenotype-investigate）
- `phenotype_first_status.json`（Route A: a06 phenotype-first-status）
- `route_c_local_retrieval_status.json`（Route C: c01 retrieval-run）
- `route_c_retrieval_status.json`（Route C: c02 retrieval-import）
- `guided_profile_manifest.json`（profile レベル）

ブロックのフィールド:

| field | 値 |
|---|---|
| `status` | `structured_context_used` / `structured_context_missing` |
| `term_usage` | `normalized_term` / `unmapped_term` / `raw_text` / `unknown` |
| `phenotype_text` | 正規化後の自由文（未指定は `unknown`） |
| `unmapped_term` | 原語のまま保持（未指定は `unknown`） |
| `present_fields` | 構造化フィールド名の `;` 区切り、無ければ `none` |
| `context_fingerprint` | SHA-256（deterministic、case/source provenance は除外） |
| `rank_contribution` | 常に `none` |
| `query_usage` | そのステージの誠実な利用記録（下記） |

### query_usage の誠実性契約

`query_usage.usage` は「そのステージの下流クエリが実際に構造化フィールドを消費したか」を記録します。現在の実装では構造化フィールドをクエリへ流し込む消費ステージは存在しないため、全ステージが `"usage": "stored_only"` と理由を記録します（例: 「TF-IDF retrieval query は自由文のまま、構造化フィールドは追加されない」）。`usage` が `stored_only` 以外になるのは、実際にフィールドを下流クエリへ渡した時だけです。記録だけして利用を主張することはしません。

`phenotype_used: bool`（Route C status 既存フィールド）は従来どおり自由文の有無を真に反映します。構造化フィールドの有無では変わりません。

### fingerprint と cache

- `context_fingerprint` は各 manifest に記録されるため、同じ自由文でも stage/environment が異なる run は manifest 上で判別でき、古い cache の再利用による混同は検出可能です。
- クエリ自体が無変更のため、retrieval の index cache キーや reverse-literature のクエリキャッシュキーは変更しません（コンテキストはそれらのクエリをfeedしない）。
- fingerprint は `case_run` の rank/evidence ダイジェストには一切含まれません。rank/score/weight 出力は byte 同一のままです。

### 後方互換（実行レベル）

- コンテキストファイル無しの既存実行: 全出力（`ranked_candidates.tsv` 等）は従来と byte 同一。manifest に `structured_context.status=structured_context_missing` が追加されるのみ。
- 任意オプションのため、既存のCLI呼び出しは変更不要です。
