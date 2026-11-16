# Canonical schema と検証

`candidate_qtl.schemas` は、表(table)の列契約・enum語彙・列間チェックを**一元管理する単一の真実源**です。以前は `discover` / `evidence-build` / `case-run` / 各 importer がそれぞれ個別に列・型・enum を検査しており、同じ列名でも nullable 条件や意味がコマンド間でずれる余地がありました。定義をここへ集約することで、契約の一貫性と、契約が変わったときに `SCHEMA_VERSION` を1か所で上げられる可搬性を確保します。

## 共有 enum

| 名前 | 値 | 用途 |
|---|---|---|
| `DIRECTIONS` | `support` / `against` / `rescue` / `context` | evidence の方向。`context` は決して加点されない |
| `CONFIDENCE_LEVELS` | `high` / `medium` / `low` | evidence の定性的確度 |
| `DRY_GROUPS` | `top` / `reserve` / `rescue` / `not_selected_yet` | Dry優先順位づけのバケット |
| `ASSAY_STATES` | `not_planned`から`completed`/`failed`/`needs_redesign`等 | Wet assayの実施状態。candidate解釈とは分離 |
| `CANDIDATE_WET_STATES` | `not_tested`、`partially_covered`、`coding_region_cleared`、genetic validation等 | coverage/quality gateを通したcandidate/mechanism状態。`causal`は含めない |
| `HYPOTHESIS_TYPES` | coding/splice/promoter/TE/SV/annotation-missing/multi-gene/linked-QTL/intergenic等 | gene順位とは独立したmechanism/locus hypothesis |
| `HYPOTHESIS_STATUSES` | proposed/pending-review/testing/weakened/rejected等 | hypothesis状態。gene dispositionとは分離 |
| `PROPOSAL_ORIGINS` | `manual` / `workflow` / `llm` | proposal経路。LLMはpending-review gateを通る |
| `CONFIDENCE_MULTIPLIERS` | high=1.0 / medium=0.75 / low=0.5 | 確度→スコア係数(enum の隣に置き乖離を防ぐ) |

`evidence.py` の `EVIDENCE_COLUMNS` / `CONFIDENCE_MULTIPLIERS` / `VALID_DIRECTIONS` は、この `schemas` からの再エクスポートです。既存の `from .evidence import EVIDENCE_COLUMNS` 等はそのまま動きます。

## canonical table 契約

- `EVIDENCE_SCHEMA` — canonical 8列(gene_id / evidence_code / direction / confidence / value / source_type / source_id / description)。`direction`・`confidence` は enum、`gene_id`・`evidence_code`・`direction`・`confidence` は非空必須。
- `CANDIDATE_SCHEMA` — `extract` が出力する candidates 表。`gene_id` 非空かつ一意、座標 `chrom`/`start`/`end` は `start <= end`・正の1-based・seqid 非空。
- `WET_ASSAY_SCHEMA` / `WET_HISTORY_SCHEMA` — assay IDと追記型state historyのcanonical契約。列順は`WET_*_COLUMNS`に集約する。詳細な遷移とcoverage gateは[Wet state](WET-STATE.ja.md)を参照。
- `LOCUS_HYPOTHESIS_SCHEMA` — hypothesis type/status/origin、gene/interval membership、parent/supersedes relation、gene dispositionとrejection scope。詳細はLocus-hypothesis ledgerを参照。
- `ANNOTATION_RESCUE_SCHEMA` — query provenance、projection quality/provider state、native overlap、review/promotion state。詳細は[QTL区間限定 annotation rescue](ANNOTATION-RESCUE.ja.md)を参照。
- `REFERENCE_RESOURCE_PANEL_SCHEMA` — 比較ゲノムresourceのrole/type/taxon/version/source、availabilityとskip理由。LLM由来のknown-gene panelとは別物。詳細は[比較 reference-resource panel](REFERENCE-RESOURCE-PANEL.ja.md)を参照。
- `EXTERNAL_RANK_INPUT_SCHEMA` / `EXTERNAL_RANK_LEDGER_SCHEMA` — 外部rankerのrank/coverage/family/source入力と、coverage-aware ledger(consensus_top/native_only/external_rescue/discordant/uncovered)。生スコアは`dry_score`へ加算しない。詳細は外部比較rank ledgerを参照。
- `PROPAGATED_TERM_SCHEMA` — 参照ゲノムのオーソログ機能語を区間候補へ伝播した語ledger。provenance(source genome/role/ontology/relationship)、`match_type=ortholog_propagated`、`review_status=pending`、trait_relevance。詳細は[trait検索語伝播](TERM-PROPAGATION.ja.md)を参照。
- `QUERY_ELIGIBILITY_SCHEMA` — homolog query の species/genus/assembly/selection basis/sequence origin と provenance監査結果、`benchmark_eligibility`(eligible/excluded_pisum/excluded_answer_derived/excluded_unknown_provenance)。Pisum由来queryを禁止。詳細は[query eligibility gate](QUERY-ELIGIBILITY.ja.md)を参照。

## 検証の使い方

```python
from candidate_qtl import schemas

errors = schemas.collect_errors(frame, schemas.EVIDENCE_SCHEMA)   # 非raise、行/列単位のリスト
schemas.validate(frame, schemas.EVIDENCE_SCHEMA)                  # 失敗時 SchemaValidationError（.errors に全件）
```

`SchemaError` は `(table, code, message, column, row)` を持ち、pandera の lazy 集約レポートと同様に**行・列単位で失敗を保存**します。

## どこで強制されるか(gate)

- **evidence**: `rank`(`hardening.rank`)が読み込んだ evidence 表を `EVIDENCE_SCHEMA` で検証します。`discover` を含む全ての evidence 生成物は最終的にここを通るため、consumer 側の単一 gate で enum ・非空ドリフトを捕捉します。
- **candidates**: `extract` が `candidates.tsv` を書く直前に `CANDIDATE_SCHEMA` で検証します。壊れた GFF3 が gene ID を重複させた場合(`gff_dialect_qc.duplicate_parent_ids` に報告)、gene_id-keyed な evidence の二重結合を防ぐため先頭行を残して重複を落とします。

## pandera(optional)

`pandera` は**任意の追加依存**です。既定の検証は依存ゼロで、local-first ランタイムが pandera を要求することはありません。

```bash
pip install candidate-qtl[schema]   # pandera を導入
```

導入時、`schemas.pandera_schema(EVIDENCE_SCHEMA)` は**同じ `TableSchema` descriptor から** `pandera.DataFrameSchema` を生成します。pandera ビューと組み込み検証は同一 descriptor 由来なので乖離しません。`schemas.pandera_available()` で有無を判定できます。

## SCHEMA_VERSION

table 契約を後方非互換に変更したら `schemas.SCHEMA_VERSION` を上げてください。これは canonical *table* 契約のバージョンであり、各コマンドの run manifest が持つ `schema_version`(manifest フォーマットのバージョン)とは別の名前空間です。
