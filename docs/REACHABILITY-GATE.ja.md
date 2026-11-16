# 統合到達可能性ゲート（Canonical Feature-Reachability Gate, issue #373）

## 目的と適用範囲

このワークフローのコマンド（CLI）は「配線（wiring）」が途中までしかなくても存在できてしまう。CI のグリーンと個々のテストだけでは「その機能がパイプラインに統合されている」ことは保証されない。本ゲートは、**統合済みと宣言された機能を、リポジトリ所有の機械可読レジストリに基づいて実 CLI エントリポイントを subprocess で実行し、成果物・マニフェストキー・決定性・ランク非関与・fixture 不変性を機械的に検証する**。

宣言は2種類しかない。

- `integrated_closing_slice` — 正準配線が完了している機能。ゲートが合格して初めて「統合済み」と数えられる。統合配線を担う PR は、この宣言と registry entry を追加し、クローズする issue を閉じることができる。
- `core_only_non_closing_slice` — 配線が未完の機能。`remaining_wiring`（残りの正準配線）と `doc_reference` だけを記録し、**統合カバレッジには数えない**。

**グリーンテストだけでは統合受容にならない。** 統合済み宣言の技術的事実（正準エントリポイントが存在し、fixture で実行でき、宣言された成果物とマニフェストキーを出し、ランク成果物を書かず、fixture を変更しない）は registry + gate で検証される。チェックボックスや説明文だけでは不十分である。

## registry の場所と schema

- 場所: `tests/fixtures/reachability/registry.json`
- 検証実装: `src/candidate_qtl/reachability.py`（標準ライブラリのみ。`jsonschema` 非依存で手書き検証）
- ランナー: `scripts/reachability_gate.py`
- ゲートの CI 実行: `.github/workflows/ci.yml` の `test` ジョブ（source モード）と `build` ジョブ（installed モード）
- 負の fixture: `tests/fixtures/reachability/negative/orphan_registry.json`（正準配線のない合成孤児。検証が必ず落ちることの見本）

トップレベル: `schema_version`（== 1）、`entries`（feature_id 昇順ソート必須）。

エントリフィールド:

| フィールド | integrated 必須 | core_only 必須 | 内容 |
|---|---|---|---|
| `feature_id` | 必須 | 必須 | `^[a-z0-9][a-z0-9-]*$`、全体で一意 |
| `integration_status` | 必須 | 必須 | 上記2値のいずれか |
| `entrypoint` | 必須 | 禁止 | `command`（`candidate_qtl.commands.COMMAND_NAMES` に存在必須）+ `argv`（非空の文字列リスト） |
| `fixtures` | 必須 | — | 非空リスト。リポジトリ相対パスのみ（絶対パス・`..` 拒否） |
| `expected_artifacts` | 必須 | 禁止 | 非空リスト。`out-dir` 相対の `/` 区切りだけ（絶対パス、空/`.`/`..` 成分、`\\` を拒否）。fnmatch グロブは通常成果物だけで可 |
| `expected_manifest_keys` | 必須 | 禁止 | `{"artifact": "declared_manifest.json", "keys": [...]}` のみ。JSON manifest は exact の宣言済み成果物へ必ず束縛する |
| `expected_exit_codes` | 任意 | — | 既定 `[0]` |
| `rank_contribution` | 必須 | — | `none` または `policy:<name>` |
| `rank_assertion` | 必須 | — | `type` ∈ `{deterministic_and_no_rank_artifacts, emits_canonical_rank_unchanged}`、非空 `reason` |
| `clean_install_smoke` / `cache_fingerprint_assertion` | 必須 | — | 真偽値 |
| `owning_issue` / `doc_reference` | 必須 | 必須(後者のみ) | 非空文字列 |
| `remaining_wiring` | — | 必須 | 非空文字列 |
| `assert_values` | 任意 | — | `{"artifact": "declared_manifest.json", "values": {dot-path: 期待リテラル}}`。値検査も同じ JSON manifest へ束縛 |

argv プレースホルダ: `{out}`（そのエントリの out-dir）と `{fixture:<相対パス>}`（宣言済み fixture の絶対パスに解決）。シェル解釈は一切ない（argv は構造化リスト）。

fail-closed 規則: duplicate feature_id / fixture または成果物のパス走査（`..`）/ 絶対パス / 未知 status / 必須フィールド欠落 / 未知コマンド / 未宣言 fixture 参照 / 未知キーはすべてエラー。CLI に存在するだけの helper は統合済みではない。feature_id ごとに `reachability.py` の canonical orchestration と一致し、そこが要求する JSON manifest を registry が露出しなければ orphan としてエラーになる。core_only が `entrypoint`・`expected_artifacts`・`expected_manifest_keys` を持つことも schema エラー。統合エントリが0件の registry はゲートがハードフェイルする。

## 新規 entry 追加手順

1. **feature_id 命名**: `^[a-z0-9][a-z0-9-]*$`。CLI コマンド名と同一にすると追跡しやすい。
2. **fixture 追加**: `tests/fixtures/reachability/features/<feature_id>/` に最小かつ決定的な入力（synthetic で可、ネットワーク・外部バイナリ不使用）。絶対パス禁止、パスは常に `/` 区切り。
3. **マニフェストキーの調べ方**: 実際にコマンドを走らせて出力 JSON を読み、その manifest 自身を `artifact` として明記する。
   ```
   python -m candidate_qtl <command> <...fixture 引数...> --out-dir /tmp/probe
   cat /tmp/probe/*_manifest.json
   ```
   安定キー（schema_version、カウント類）だけを `expected_manifest_keys` / `assert_values.values` に書く。`wall_seconds`・`computed_at_utc`・`mtime_ns`・絶対パスなど volatile な値は絶対に assert しない。manifest JSON に glob を使ったり、他の JSON を探索して偶然同名キーを見つけることはできない。
4. **決定性確認**: 同じ入力で2回実行し、宣言する TSV が byte 一致すること（`cmp` または diff）を確認してから commit する。JSON マニフェストは volatile キーのみの差異を許容（ゲートは TSV の byte 一致のみを強制するわけではなく、テストが `run_entry` の再実行間の TSV byte 一致を検証する）。さらに全ゲート実行の前後で canonical `tests/fixtures/case_run` を実行し、`candidates.tsv`、`rank/ranked_candidates.tsv`、`rank/score_contributions.tsv`、`stage_a/deterministic_shortlist.tsv` の SHA-256 が同一であることを必須にする。
5. **rank 方針の理由明記**: `rank_contribution=none` なら、その機能がスコア・ランク成果物を書かない理由を `rank_assertion.reason` に書く。ゲートは実行後、out-dir に `rank/ranked_candidates.tsv` / `score_contributions.tsv` / `candidates.tsv` / `flagged_candidates.tsv` / `stage_a/**` / `dry_rank*` / `ranked_candidates*` が無いことも検証する。`policy:<name>`（例: `policy:deterministic_v1`）は、その機能が正準ランク段を変更せずに実行する場合だけ使う（`type: emits_canonical_rank_unchanged`）。この宣言の entry は no-rank-artifact 走査の対象外になるため、canonical ランク成果物が**変更されていないこと**の保証（gate の canonical-baseline before/after SHA-256 比較、および既存 E2E テストの byte 同一性）とセットで宣言する。
6. **PR**: `.github/PULL_REQUEST_TEMPLATE.md` の統合宣言セクションに feature_id とテスト名を書く。`tests/test_reachability_gate.py` のコミット済み registry 自己検査（`CommittedRegistrySelfCheckTests`）が CI で回る。

## 失敗診断

| エラー | 意味 | 対処 |
|---|---|---|
| `unknown command 'x'` | コマンドが `COMMAND_NAMES` に無い（typo、未配線、古い wheel） | `python -m candidate_qtl --help` で実在を確認し、entry を修正 |
| `fixture does not exist` / `repo-relative` / `.. component` | fixture パスが不正・未 commit | fixture を追加するかパスを修正 |
| `expected artifact not produced: X` | コマンドは成功したが宣言成果物が出ない | 実際の出力名を確認し registry を直す（あるいはコマンドの欠陥） |
| `manifest key not found in 'X.json': 'a.b'` | 宣言した JSON manifest にキーが無い | 実マニフェストのキーと artifact 束縛を確認 |
| `not wired into canonical orchestration` | 実在 CLI/helper だが統合の正準経路と manifest 露出が未登録 | feature を `core_only` に留めるか、正準配線と `CANONICAL_INTEGRATIONS` を同一PRで追加 |
| `assert_values mismatch` | 実値が期待リテラルと異なる | 実値を確認して registry を直す（機能変更なら直す方） |
| `rank-impacting artifact(s) written` | rank_contribution=none なのにランク成果物を出力 | コマンド側の配線を確認 |
| `fixture(s) changed during the run` | 実行が入力 fixture を書き換えた | コマンドの入出力分離を修正 |
| `command exited N; expected [0]` | 実行失敗（exit 1=入力エラー、2=内部エラー） | stderr tail を見て入力かコマンドかを判断 |
| installed モード: `import falls back to the source checkout` / `not from site-packages` | 検証対象 python が editable インストール・PYTHONPATH 経由でソース解決している | クリーン venv の非 editable インストール（wheel）を検証対象にする |
| `candidate_qtl is not importable` | 検証対象 python にパッケージが無い | インストールして再実行 |

## CI での実行方法

- source モード（`test` ジョブ、editable install）:
  `python scripts/reachability_gate.py --mode source --repo-root "$GITHUB_WORKSPACE"`
- installed モード（`build` ジョブ、クリーン wheel インストール）:
  `/tmp/cleanenv/bin/python scripts/reachability_gate.py --mode installed --python /tmp/cleanenv/bin/python --repo-root "$GITHUB_WORKSPACE"`
  （ゲート自体をクリーン環境の python で動かすことで、スクリプト自身の import も site-packages から行われ、ソースフォールバックの余地がない。`candidate_qtl.commands` の import 連鎖が pandas 依存のため、パッケージ未導入の素の python では検証できないためである。）

ローカル再現:

```
python scripts/reachability_gate.py --mode source --repo-root <worktree>
python scripts/reachability_gate.py --mode installed --python <cleanenv>/bin/python --repo-root <worktree>
```

installed モードは最初に `assert_no_source_fallback` を実行し、検証対象 python が `candidate_qtl` を site-packages から import すること（ソース checkout にフォールバックしないこと）を証明してから、各エントリをリポジトリ外の cwd で実行する。この検査の失敗は即ハードフェイル。

注意:

1. **installed モードの venv はリポジトリの外に作る**。リポジトリ内（`repo_root` 配下）に作った venv は、たとえ正規の wheel を非 editable インストールしていても、import パスが `repo_root` 配下に解決されるため source-fallback ガードに引っかかる。CI は `/tmp/cleanenv`（リポジトリ外）を使う。
2. **source モードは cwd=repo_root で実行される**ため、子プロセスの import がソース checkout を解決しても当たり前であり、検証できるのは「宣言された fixture の不変性」までである。完全に密閉された（ソース非依存の）証明は installed モードである。どちらのモードも PYTHONPATH を子プロセスから除去して実行する（呼び出し側の PYTHONPATH による import 迂回を防ぐ）。
3. **`clean_install_smoke` / `cache_fingerprint_assertion` は宣言的メタデータ**であり、人間が読むための記録である。実行カバレッジは両モードの全 integrated エントリで担保され、この2フラグの真偽で実行内容は変わらない。

## 現在の登録状況（2026-08-30 時点）

- 統合済み: `budget-run`（issue #375 の決定論的バジェットエンジン）、`discovery-chain-benchmark`、`family-candidate-bridge`（#374 の canonical bridge 配線）、`family-index`、`mechanism-import`（#380 の payload adapter 配線）、`mechanism-import-coverage`（#658 のphenotype-only causal-layer coverage監査）、`route-a-interval-match`（`discovery-chain-benchmark` / `family-index` / `route-a-interval-match` は issue #373 の見本配線）、`route-family-metrics-not-in-index`、`triple-route-e2e`（#374 の guided-profile triple-route 全鎖 clean-install E2E）。
- **#360 / #361 / #362 / #364 / #367 は未統合。** 2026-08-12 時点で main に未 merge であり、registry には登録されていない。将来の各 PR が、統合配線と共に registry entry を追加する契約である（core_only 宣言はその中間状態の記録手段）。

### family-candidate-bridge

- entrypoint: `family-candidate-bridge --route-b-atlas <atlas fixture> --out-dir <out> --case-id demo`
- fixture: `tests/fixtures/reachability/features/family_candidate_bridge/atlas.tsv`（route B atlas 単独の最小セット。curated 軸 1 行 = tier1、naming bucket 1 行 = tier5。atlas 単独で bridge が通ることは `tests/test_family_candidate_bridge.py::test_route_b_atlas_is_a_sufficient_standalone_input` で担保）
- 検証するテスト: `tests/test_reachability_gate.py::test_committed_registry_passes_source_mode_gate`（7 エントリ）と `::test_family_candidate_bridge_output_carries_the_mapping_tiers`、`tests/test_family_candidate_bridge.py`（22 列 schema、fingerprint 決定性、not_evaluated 語彙）。
- rank 方針: `rank_contribution=none`。bridge は 3 route の family 証拠を 1 表に正規化する report 専用 artifact で、score 列・rank artifact を書かない。`guided_profiles.py` の triple-route は compare と family-to-copy の間で bridge を実行し、family-to-copy は bridge の mapping tier を shortlist 昇格ゲートとして使う（`family_to_copy.py` の `bridge_promotion_gate`）。
- dependency clusters 配線: triple-route profile は `route-b-prepare` が生成する `evidence_dependency_clusters.tsv`（`<b01>/evidence_dependency_clusters.tsv`）を bridge の `--dependency-clusters` 入力へ自動 chain する。同表は Route B の orthology 表と凍結 atlas に現れた orthogroup/HOG ごとに 1 cluster を持ち（同一 orthogroup/HOG 由来の行 = 同一 source/homology 依存）、schema は `dependency_cluster` / `dependency_status` / `root_source_id` / `source_ids`。bridge の Route B 行は orthology の `source_id` が無い場合 `orthogroup:<id>` トークンを source として持つため、この表を通じてのみ cluster へ解決される（表が無ければ `dependency_unknown` のまま、自動昇格なし）。**cluster は証拠依存のラベルであり加点しない**: 全行 `rank_contribution=none` のまま。
- dependency clusters が無い場合: 供給 orthology 表も atlas も orthogroup/HOG 値を持たない場合（`route-b-prepare` 自身の local-blastp 由来 orthology 表は orthogroup 列が空）、cluster は 0 件で表は header-only になり、bridge は `dependency_clusters` 入力を空として受け取り全 Route B 行が `dependency_unknown` のまま。捏造はしない（`route_b_prepare.py` の `dependency_clusters.no_data_rule`）。

### mechanism-import

- entrypoint: `mechanism-import --payload <legacy 7列 payload fixture> --out-dir <out>`
- fixture: `tests/fixtures/reachability/features/mechanism_import/mechanism_payload.tsv`（凍結 T11 形式の 7 列、3 行）
- 検証するテスト: `tests/test_reachability_gate.py::test_committed_registry_passes_source_mode_gate`（同値 7 エントリ）と `::test_mechanism_import_output_preserves_legacy_values`（mechanism_statement / selection_rationale / adapter map / accepted_payload_schema の実値）。
- rank 方針: `rank_contribution=none`。hypothesis 行は score を持たず、rank artifact を書かない。

### mechanism-import-coverage

- entrypoint: `mechanism-coverage-audit --mechanisms <12行 fixture> --mechanism-import-manifest <import manifest> --localization-freeze-manifest <manifest> --investigation <phenotype> --out-dir <out>`
- fixture: `tests/fixtures/reachability/features/mechanism_coverage/`。高confidence signaling局在、隣接regulation 1件、alternative 2件を含むため、偏りを均等化せず`status=pass`になる。
- assert: `interval_blind=true`、`rank_contribution=none`、`freeze_ready=true`、merged mechanisms / auditのSHA-256がmanifestへ固定される。
- rank 方針: causal-layer labelはcoverage監査専用。score・rank artifactは書かず、下流のsource/accession/sequence verificationを迂回しない。

### triple-route-e2e（#374 Acceptance Criterion 12）

- feature_id: `triple-route-e2e`（canonical orchestration は `guided-profile` と `guided_profile_manifest.json` の露出を要求する）。
- entrypoint: `python -m candidate_qtl guided-profile --profile triple-route --execute <fixtures> --case/--resources <case fixture>`。1 回の実行で **bridge → family-to-copy → case-run（integrated narrowing）→ stage_a_narrowed → next-experiment** の全鎖を生成する。合成 fixture のみ（record cache / source cache / retrieval index すべて synthetic、private 配列・外部 API 不使用）。
- blast 層の範囲外の明記: **blast 層（`blast-panel-freeze` / `blast-panel-search`）は `--proteins` 必須のため本 E2E には含まれない**（`--proteins` 非指定のため profile は blast 段を計画しない）。blast 層の実装検証は `tests/test_blast_panel_integration.py`（実 BLAST+ の freeze→search E2E）と `tests/test_triple_route_profile_bridge.py`（MockBlastTools 経由の guided 配線、`--proteins` 指定時に blast 段が計画・実行されること）がカバーする。
- fixture: `tests/fixtures/reachability/features/triple_route_e2e/`。`records/`（UniProt-shaped 3 件）と `pmc/`（Europe PMC-shaped 3 件）と `retrieval_index/plant_ref_test/`（`retrieval-index-build` で生成したローカル index）は**ディレクトリ fixture** で再帰ハッシュ検証される。
- 成果物: `guided_profile_manifest.json` / `z01b_family_candidate_bridge/{family_candidate_bridge.tsv, _manifest.json}` / `z02_family_to_copy/{family_copy_shortlist.tsv, family_copy_discrimination.tsv, family_copy_manifest.json}` / `z03_case_run/integrated_narrowing/{integrated_narrowing.tsv, narrowing_manifest.json}` / `z03_case_run/stage_a_narrowed/deterministic_shortlist.tsv` / `z03_case_run/next_experiment/{next_experiment_options.tsv, next_experiment_decision.json}` / `z03_case_run/rank/ranked_candidates.tsv`。
- assert: `guided_profile_manifest.json` の `legacy_profile.overall_status=completed`・`legacy_profile.final=true`・`profile=triple-route`、`family_copy_manifest.json` の `bridge.supplied=true`（bridge が family-to-copy へ届いた証明）、`family_candidate_bridge_manifest.json` の `rank_contribution=none`。
- #859: このfixtureはcanonical semantic reviewが未実行のため、外側は `overall_status=partial`・`final=false`。`canonical_execution.state=semantic_review_required`、manifestの内容hashと再開先を検査する。旧profileのbridge/copy配線完了を、canonical Wet完成と読み替えない。
- rank 方針: `rank_contribution=policy:deterministic_v1` + `type: emits_canonical_rank_unchanged`。case-run 段が正準 deterministic rank 表（`ranked_candidates.tsv` 等）を `z03_case_run/` 配下へ変更なしで出力するため。byte 同一性は `tests/test_triple_route_profile_bridge.py`（`test_legacy_stage_a_and_rank_unchanged` / `test_two_runs_produce_identical_artifacts`）が担保し、gate は canonical-baseline before/after の SHA-256 比較でランク配管が不変であることを追加検証する。
- 検証テスト: `tests/test_reachability_gate.py::test_committed_registry_passes_source_mode_gate`・`::test_triple_route_e2e_bridge_reaches_family_copy`・`::test_triple_route_e2e_rank_outputs_are_confined_to_the_case_run_subtree`、および決定性テスト（宣言 TSV の 3 実行 byte 一致）。

## budget-run（決定論的バジェットエンジン, issue #375）

- feature_id: `budget-run`（CLI コマンド `budget-run` と同一）。canonical orchestration は `budget_manifest.json` の露出を要求する（`src/candidate_qtl/reachability.py` の `CANONICAL_INTEGRATIONS`）。
- entrypoint: `python -m candidate_qtl budget-run --extract-dir {fixture:.../extract} --out-dir {out} --budget-profile zero`。
- fixture: `tests/fixtures/reachability/features/budget_run/extract/`。`--extract-dir` は抽出ディレクトリ（ディレクトリ fixture）を指すため、**ディレクトリ fixture の不変性は再帰ハッシュで検証**される（`_fixture_sha256` が相対パス順にツリー全体をハッシュ。`test_declared_fixtures_are_unchanged_after_the_gate` も同じ関数を使う）。
- 成果物: `tier0_inventory.tsv` / `tier1_clusters.tsv` / `tier2_plan.tsv` / `budget_handoff.tsv` / `budget_manifest.json` / `benchmark_resource_usage.json`。ランク成果物は書かない（`rank_contribution=none`）。
- assert: `budget_manifest.json` の `schema_version=1.0` と `status=completed`（zero プロファイルはクエリを1件も実行しないため `budget_exhausted` にならない。ステータス規則は docs/BUDGET-ENGINE.ja.md 参照）。
- 検証テスト: `tests/test_reachability_gate.py` の `PositiveE2EGateTests.test_budget_run_manifest_reports_completed_zero_profile` ほか、committed registry の全エントリを駆動するゲートテストが自動で実行する。

## PR template の使い方

`.github/PULL_REQUEST_TEMPLATE.md` の「Integration declaration」は必須セクションであり、どちらか一方だけにチェックする。

- 統合を主張する PR: `integrated_closing_slice` を選び、registry の feature_id と、それを検証する `tests/test_reachability_gate.py` のテスト名を明記する。ゲート（source + installed 両モード）がグリーンであることが統合の前提であり、チェックボックスはその記録である。
- 配線途上の PR: `core_only_non_closing_slice` を選び、残りの正準配線を書く。この状態ではクローズする issue を閉じてはならない。
