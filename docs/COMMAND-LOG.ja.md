# command_log.tsv 形式仕様と監査

`command_log.tsv`（benchmark run の実行コマンド台帳）の canonical 形式、legacy 12 variant の扱い、
機械監査コマンド `command-log-validate` の仕様を定義する。

## 背景

T11 run（`t11-130-fresh-20260812T045741Z`）の `cases/*/command_log.tsv` 98ファイルには形式仕様が無く、
12種のヘッダ variant・3種の時刻表記・2種の欠落マーカー（`-` と空文字）が混在していた。
さらに extract 8 case・dryrun 15-16 case の行欠落（計23-24行）、時刻/wall 欠落40行が確認された。
本仕様は (1) 新しい run が書くべき canonical 形式、(2) 既存 legacy ファイルの機械監査、の2つを提供する。
**既存の T11/T8 frozen run・`benchmarks/**` ・`freeze/**` の実データは書き換えない**（監査は読み取り専用）。

## 1. canonical 形式（新規 run が記録する形式）

- ファイル: UTF-8、タブ区切り、行末 `\n`（`csv.writer(delimiter="\t", lineterminator="\n")` による正規エスケープ）
- ヘッダ: 以下13列を**この順序で完全一致**させる（`CANONICAL_COLUMNS`）
- 欠落値は空セル（欠落マーカー）。空や `0` で補完しない
- 時刻は UTC aware（`YYYY-MM-DDTHH:MM:SSZ`）のみ。naive 値は書けない（writer が拒否）

| 列 | 型 | 必須 | 例 |
|---|---|---|---|
| `schema_version` | 文字列 `"1.0"` | 必須 | `1.0` |
| `case_id` | 文字列 | 必須 | `pea_brc1` |
| `stage` | `extract` / `dry_run` / `full_execute` | 必須 | `dry_run` |
| `command` | 文字列（人間向け要約） | 必須 | `guided-profile --profile triple-route` |
| `argv_json` | JSON配列（秘密は redact 済み） | 必須 | `["python","-m","candidate_qtl.cli","extract",...]` |
| `status` | `success` / `failure` / `skipped` | 必須 | `success` |
| `exit_code` | 整数 / 空 | 任意 | `0` |
| `started_at_utc` | `YYYY-MM-DDTHH:MM:SSZ` | 必須 | `2026-08-12T06:55:00Z` |
| `completed_at_utc` | `YYYY-MM-DDTHH:MM:SSZ` | 必須 | `2026-08-12T06:56:05Z` |
| `wall_time_seconds` | 非負浮動小数 / 空 | 任意 | `65.0` |
| `working_directory` | 文字列 | 必須 | `/work/run-root` |
| `tool_version` | 文字列 | 必須 | `0.22.0` |
| `notes` | 自由文 | 任意 | `16/16 stages completed` |

- `exit_code` / `wall_time_seconds` / `notes` は記録時点で不明でも空で書ける（任意列）。他の列の欠落は strict 検証で error
- 実行しなかったことが分かっている stage は `status=skipped` の行として正直に記録する（行を消さない）
- `wall_time_seconds` は `timing()` ヘルパー（monotonic clock、3桁小数）で計測する

### 秘密値の禁止

`argv` は構造化値（`argv_json` 配列）として記録し、**秘密値をそのまま記録しない**。
`--api-key` / `--token` / `--password` / `--secret` 等のフラグ値、`key=value` 形式の代入は
`<redacted:名前>` に置換される（組み込み規則 + 呼び出し側の追加パターン）。

### writer 利用例

```python
from datetime import datetime, timezone
from candidate_qtl.command_log import CommandLogWriter, timing

with CommandLogWriter("cases/pea_x/command_log.tsv", tool_version="0.22.0") as writer:
    with timing() as elapsed:
        result = run_extract(...)
    writer.append(
        case_id="pea_x", stage="extract", command="candidate-qtl extract",
        argv=["candidate-qtl", "extract", "--interval", "chr1:1-2", "--api-key", api_key],
        status="success" if result.ok else "failure",
        exit_code=result.returncode,
        started_at_utc=started, completed_at_utc=datetime.now(timezone.utc),
        wall_time_seconds=elapsed(), working_directory=str(workdir),
    )
```

- ヘッダはファイル作成時のみ書く。既存ファイルのヘッダが canonical と一致しなければ追記を**拒否**する（legacy 混在ファイルを canonical 行で汚さない）
- 書き込み後に `flush()` + `os.fsync()` し、writer 障害は例外として伝播させる（成功に見せない）

### 期待 stage manifest（`--manifest`）

```json
{"schema_version": "1.0", "expected_stages": ["extract", "dry_run", "full_execute"]}
```

- stage 名は canonical 名のみ。未知の stage 名は即 error（タイポを黙殺しない）
- **manifest が無い限り「この stage が欠落している」とは断定しない**。無い場合は
  `not_evaluated_expected_stage_set_missing` として report に載せ、欠落判定をしない

## 2. legacy 12 variant（T11 実測ヘッダ）

variant 判定はヘッダ列タプルの**完全一致**（`LEGACY_VARIANTS`）。集約ファイルの
`case_id` 前置（`case_id,case_id` 重複含む）と、既知 subset ヘッダ（T8 の V1 から
`wall_seconds` が無い形など）は、前置除去・順序部分列で判定する。

| ID | case数 | ヘッダ（タブ区切り） |
|---|---|---|
| V1 | 10 | `case_id,stage,command,exit_code,start_utc,end_utc,wall_seconds,notes` |
| V2 | 8 | `argv,exit_code,log,started_utc,ended_utc,wall_seconds`（case_id列なし → 親ディレクトリ名） |
| V3 | 8 | `case_id,command,argv,exit_code,log,start_utc,end_utc,wall_sec` |
| V4 | 8 | `case_id,phase,argv,exit_code,stdout_log,stderr_log,started_utc,finished_utc,wall_seconds` |
| V5 | 8 | `case_id,step,command,exit_code,log_path,run_id,started_utc,ended_utc,wall_seconds,notes` |
| V6 | 8 | `command,argv,exit_code,log_stdout,log_stderr,started_at_utc,finished_at_utc,wall_seconds` |
| V7 | 8 | `phase,command,argv,exit_code,stdout_log,stderr_log,started_at_utc,ended_at_utc,wall_seconds` |
| V8 | 8 | `run_type,argv,exit_code,log_file,started_utc,ended_utc,wall_seconds` |
| V9 | 8 | `started_at_utc,finished_at_utc,wall_seconds,argv,exit_code,log_path,status_summary` |
| V10 | 8 | `started_at_utc,finished_at_utc,wall_seconds,exit_code,argv,log_path` |
| V11 | 8 | `step,argv,exit_code,stdout_log,stderr_log,started_utc,ended_utc,wall_seconds` |
| V12 | 8 | `step,command_argv,exit_code,stdout_log,stderr_log,start_utc,end_utc,wall_seconds` |

### 時刻表記3種（実測）

| 分類 | 実例 | 備考 |
|---|---|---|
| `iso_z` | `2026-08-12T06:58:41Z` | Z と任意の `±HH:MM` オフセットを同分類（aware として扱い、UTC 換算はしない） |
| `run_id_style` | `20260812T070045411092Z` | yyyyMMddTHHmmssffffffZ |
| `datetime_no_zone` | `08/12/2026T15:52:15` | MM/dd/yyyyT HH:mm:ss、zone 無し → naive |

欠落マーカーは `-` と空文字の両方（V11/V12 は `-`、V5 は空文字、V9 は wall のみ `-`）。
V2 の `(mtime-based)` は「非数値」として数え、数値として扱わない。

### stage 名の正規化（実測）

- `extract` → `extract`
- `dryrun`, `dry-run`, `dry_run`, `dryrun-guided-profile` → `dry_run`
- `full-execute`, `execute`, `full_execute`, `guided-profile-execute`, `guided-profile-full`, `guided-profile-triple-route` → `full_execute`
- 未知の stage 名は正規化せず `unknown` として報告（`unnormalized_stage_values`）
- stage 列を持たない variant（V2, V6, V7, V9, V10, V11, V12）は argv から推測する
  （`(dry-run)`/`/dryrun` パス → dry_run、`--execute` → full_execute、先頭 `extract` → extract）。
  推測由来の行は `stage_source.inferred` として区別し、記録された事実と混同しない
- T8 の実測 stage 値 `guided-profile` は曖昧（`--execute` 付き=full_execute、無し=dry_run）なので
  argv（V1 形では argv 列が無いため command 列）から推測する。推測できない場合は
  `unnormalized_stage_values` に入り、manifest 付きでも欠落断定はしない
- 非空だが正規化できない stage 値は `stage_rows.unknown` / `stage_source.unknown` として
  報告する（欠落 `missing` と区別。欠落は stage 列自体が無いか空/`-`）

### 列マッピング（legacy → canonical）

`argv_json`←`{argv,command_argv}`、`command`←`{command,command_argv}`、`started_at_utc`←
`{start_utc,started_utc,started_at_utc}`、`completed_at_utc`←`{end_utc,ended_utc,ended_at_utc,finished_at_utc}`、
`wall_time_seconds`←`{wall_seconds,wall_sec}`、`status`←`{status_summary}`、`stage`←`{stage,phase,step,run_type}`。
`log`/`log_path`/`run_id`/`stdout_log` 等の付随列は情報列として無視（canonical 項目に捏造しない）。

## 3. `command-log-validate` コマンド

```
candidate-qtl-bench command-log-validate --log PATH [--mode strict|legacy] [--manifest MANIFEST] [--json]
```

- `--log`: ファイル、またはディレクトリ（`**/command_log.tsv` を再帰収集。`cases/*/command_log.tsv` を包含）
- `--mode`: `legacy`（既定）/ `strict`
- 戻り値: `{"files": [...], "totals": {...}, "exit_code": 1 if any errors else 0}`（複数ファイル集計）
- exit code: **error があれば 1**、finding（warning）のみなら 0。**legacy 監査は historical log の欠落で失敗しない**

### legacy モード

`audit_legacy_command_log` の finding を warning 化。error になるのは構造的破損のみ:
ファイルが読めない / UTF-8 でない / ヘッダがどの variant にも一致しない / 行の列数がヘッダと不一致 /
`case_id,case_id` 重複ヘッダ / データ行がヘッダより前にある。

report の主キー:

- `variant_id`（複数ヘッダ混在ファイルは `mixed` + `blocks[]` 単位の variant 分類。全ブロックが
  同一ヘッダ形状の場合はその variant を返す。例: T8 集約は 10 ブロックすべて V1 形 → `V1`）
- `coverage.argv_complete`: `full_argv` / `medium_argv` / `short_argv` / `no_argv`
  （full=プログラム呼び出しを含む、medium=シェル包囲・中略 `...` を含む、short=ラベル・オプション羅列のみ、no=欠落。
  判断不能を完全扱いしない）
- `coverage.exit_code` / `coverage.time` / `coverage.wall`（欠落・非数値を正直に数える）
- `missing_fields`: canonical 13項目ごとの欠落数と行番号（列が variant に無い場合も欠落として数える。
  `rows` は先頭100行で打ち切り、`count` は全行分）
- `timestamp_formats`: `iso_z` / `run_id_style` / `datetime_no_zone` / `missing` / `unparseable`
- `stage_rows` / `stage_source`（`column`/`inferred`/`unknown`/`missing`）/ `unnormalized_stage_values`
- `missing_expected_stages`（manifest 指定時のみ。無ければ `not_evaluated_expected_stage_set_missing`。
  stage 値が正規化できない場合は `not_evaluated_unnormalized_stage_values`、複数ヘッダ混在ファイルは
  `not_evaluated_mixed_blocks` となり、いずれも欠落断定をせず warning で報告する）
- `parse_errors` / `duplicate_stage_rows` / `end_before_start_rows` / `naive_timestamp_rows`

### strict モード（canonical）

ヘッダが13列と完全一致しないと error。行ごとに:
- required 列（`exit_code`/`wall_time_seconds`/`notes` 以外）の欠落 → error。
  欠落は空セルと `-` の両方（`MISSING_MARKERS`）で判定し、`-` を「値あり」と誤認しない
- `schema_version` 不一致 / `stage` 非 canonical / `status` 非許容値 → error
- `argv_json` が valid JSON 配列でない → error
- 時刻が naive または不正 → error、`completed_at_utc < started_at_utc` → error。
  必須である `started_at_utc`/`completed_at_utc` の `-` は required 欠落 error
- `wall_time_seconds < 0` または非有限値（NaN/inf） → error。`|wall - (completed-started)| > 30秒` → warning
- 任意列（`exit_code`/`wall_time_seconds`/`notes`）の欠落（空または `-`）→ warning

## 4. 既存データの扱い

- T11 run root・`benchmarks/t11_130_from_scratch_20260812`・`benchmarks/t8_130_from_scratch_20260811`
  の `executor/command_log.tsv`（集約・ヘッダ混在・`case_id,case_id` 重複を含む）は**読み取り専用**。
  監査は sha256 を変えない（テストで検証）
- 集約ファイルはブロック分割され、ブロックごとに variant 分類・集計。重複ヘッダは `parse_errors`
- log 品質は candidate score / rank / benchmark hit-miss には一切影響しない（スコア計算コードに触れない）

## 5. 不変条件

1. 欠落フィールドの補完（空文字/0/プレースホルダ）による「完全」扱いはしない。missingness と coverage を出す
2. 時刻は UTC aware へ正規化。naive / 不正 / 終了<開始は strict で error、legacy で finding
3. 期待 stage の欠落断定は manifest が与えられた場合のみ
4. 秘密値の記録禁止（redact 済み argv_json のみ）
5. 既存 frozen run・benchmarks/**・freeze/** への書き込み・変更・削除・移動禁止
