# CI と開発ワークフロー

## 開発環境

~~~bash
python -m pip install -c constraints.txt -e ".[dev,alignment]"   # 固定版 + ruff, coverage, build, pysam
python -m unittest discover -s tests
ruff check src scripts tests
~~~

## 依存の固定と再現性

`constraints.txt` に主要依存(biopython / pandas / PyYAML)とoptional alignment依存(pysam)の固定版を置く。default installはpysamを要求せず、配列I/Oを使うLinux/WSL環境では`pip install -c constraints.txt -e ".[alignment]"`を使う。

- **更新手順**: 依存を上げる場合は `constraints.txt` の版を編集 → 新環境で `python -m pip install -c constraints.txt -e ".[dev,alignment]"` → `python -m unittest discover -s tests` が緑になることを確認 → コミット。version は release tag ごとに固定される。
- `candidate-qtl doctor --json` は `python_version` と `package_versions`(biopython/pandas/pysam/PyYAML/candidate-qtl)、および利用可能な外部ツール(BLAST+/Primer3/minimap2/HMMER/bcftools)の `tool_versions` を記録するので、run のたびに実際に使われた版を残せる。CI の `test` job でも毎回出力する。

WSL/Ubuntu(または Linux)を前提とする。pysam は`alignment` extraで、**配列 I/O を行う command はextra導入済みLinux/WSLで実行する**。package import、表計算command、default clean-wheel検証はpysamなしで動く。どちらに属するかは`doctor`の`pysam_available`と`pysam_dependent_commands`で判定できる。

## CI ジョブ(`.github/workflows/ci.yml`)

| job | 内容 |
|---|---|
| `test` | Ubuntu / Python 3.12 でproduction/core partition、core coverage閾値、`compileall`、helper、CLI reference、reachabilityを常時検証 |
| `schema-extra` | optional な `schema` extra(pandera)を install し、pinned `test` job では skip される `PanderaBackendTests` を実際に走らせる |
| `lint` | `ruff check src scripts tests`。実バグ中心(`select = ["F", "E9"]`)で、既存の密な1行スタイルへ広範なstyle強制はしない |
| `build` | `python -m build` で sdist/wheel を作り、clean venv へ wheel を install して `candidate-qtl --version` と import を検証 |
| `integration-blast` | `ncbi-blast+` を導入し、実 `makeblastdb`+`blastp` で `discover` の local_blastp 経路を通す(`tests/test_integration_blast.py`) |
| `windows-wrapper` | `windows-latest` で PowerShell ラッパー `scripts/candidate-qtl.ps1` を構文検証。pysam の制約により、Windows ではパッケージ install ではなく Windows 向け入口(.ps1)を検証する |

benchmark・frozen fixture・audit harnessのtestは`.github/workflows/benchmark-ci.yml`へ分離した。関連path変更時、手動、週次scheduleで実行し、通常のproduction変更で重いhistorical suiteを毎回起動しない。

全範囲はローカルでも再現できる(下記「ローカルで CI と同等を回す」)。

`tests/test_integration_blast.py` は BLAST+ が PATH に無ければ自動 skip するため、通常の unit test は外部ツール無しで hermetic に走る。

## ローカルで CI と同等を回す(`scripts/local_ci.sh` / `scripts/local_ci.ps1`)

private リポジトリの GitHub ホストランナーは Actions 分数を消費する(Windows ランナーは分単価2倍)。支払い停止や spending limit でランナーが起動しなくなると、ジョブは**ステップ0件・数秒で失敗**し、アノテーションに `The job was not started because recent account payments have failed or your spending limit needs to be increased` が入る。この状態は待っても回復しないので、ローカルで同等の検証ができる必要がある。

~~~bash
scripts/local_ci.sh                 # Linux 側のfull validation（core + benchmarkを含む全test）
python scripts/run_test_partition.py --partition core
python scripts/run_test_partition.py --partition benchmark
scripts/local_ci.sh lint build      # job を指定
scripts/local_ci.sh --list
pwsh -File scripts/local_ci.ps1     # windows-wrapper job (PowerShell 専用)
~~~

`scripts/local_ci.sh`は`.github/workflows/ci.yml`と`benchmark-ci.yml`を包含する非課金full gateである。BLAST+ が PATH に無い場合`integration-blast`はSKIPと表示され、PASSを騙らない。

### venv を job ごと・repo root ごとに作る理由(重要な罠)

reachability gate は `python -m candidate_qtl` を**サブプロセスで、PYTHONPATH を意図的に剥がして**起動する(呼び出し側が別のソースツリーを注入できないようにするため)。したがってサブプロセスが import するのは**インストール済みパッケージ**であり、venv がどのチェックアウトを指しているかが結果を決める。

複数の git worktree で1つの venv(= 別チェックアウトの editable install)を共有すると、**gate は「別チェックアウトのコード」を「このチェックアウトの fixture」に対して実行する**。その結果、canonical baseline が `stage_a/deterministic_shortlist.tsv` を produce できず、`tests/test_reachability_gate.py` が丸ごと(13件)失敗する。コードの問題に見えるが原因は環境であり、しかも「環境が壊れている」のではなく**venv の向き先が違うだけ**である。

`scripts/local_ci.sh` は repo root のパスをキーにして `~/.cache/candidate-qtl-local-ci/<hash>/` 配下に venv を作るため、worktree 同士が環境を共有しない。手で回す場合も worktree ごとに venv を作ること。

~~~bash
cd <worktree>
python3 -m venv ~/venv-<name>
~/venv-<name>/bin/pip install -c constraints.txt -e .
~/venv-<name>/bin/python -m unittest discover -s tests
~~~

## CI ゲート

計測して表示するだけで落ちない検査は、実質的に検査ではない。次の2つは閾値/差分で job を失敗させる。

### CLI reference の drift

`docs/CLI-REFERENCE.md`はdefault production surface、`docs/CLI-BENCH-REFERENCE.md`は`candidate-qtl-bench`の全research/benchmark surfaceをargparse registryから生成する。旧`BENCHMARK-CLI-REFERENCE.md`は重複本文を持たずcanonical版への案内だけを生成する。CLIを意図的に変更したらgeneratorを実行してコミットする。command数は生成物から取得し、手書き文書へ固定値を転記しない。

~~~bash
python scripts/generate_cli_reference.py          # 再生成
python scripts/generate_cli_reference.py --check  # CI と同じ drift 検査
~~~

`--check` は `test` job の専用ステップで走る。`tests/test_cli_reference.py` も同じ一致を検証しているが、専用ステップにしておくと「CLI を変えて再生成を忘れた」ことが test failure ではなく drift として直接分かる。

### coverage 閾値

`pyproject.toml` の `[tool.coverage.report] fail_under` にゲートを置く。CI のフラグではなく設定ファイルに置いてあるので、ローカルの `python -m coverage report` も同じ閾値で落ちる。

~~~bash
python -m coverage run scripts/run_test_partition.py --partition core
python -m coverage report --omit="src/candidate_qtl/*benchmark*.py,src/candidate_qtl/audit_freeze.py,src/candidate_qtl/commands/benchmark.py,src/candidate_qtl/challenge_aggregation.py,src/candidate_qtl/command_log_validation.py,src/candidate_qtl/corpus_masking.py,src/candidate_qtl/knowledge_coverage.py,src/candidate_qtl/lane_coverage.py,src/candidate_qtl/query_eligibility.py"
~~~

現在値は 80(導入時の実測 branch coverage は 84%)。これは目標値ではなく ratchet で、**PR を通すために下げない**。上げるときは意図的に上げる。

## リリース / バージョンガバナンス

### 権威バージョンソース

`src/candidate_qtl/__init__.py` の `__version__` が唯一の version リテラルである。`pyproject.toml` は静的 `version` を持たず `[tool.setuptools.dynamic]` でこの属性から解決するため、`candidate-qtl --version`、wheel の package metadata、benchmark manifest の `tool_version` はすべてここから導出され、互いに食い違う経路がない。README の `Version: \`...\`` 行と `CITATION.cff` の `version` は文芸的な重複面として意図的に保持しているが、`tests/test_version_single_source.py` が機械的に一致を検証する。

### benchmark provenance の決定(0.22.0)

`benchmarks/supervised_58case_20260810/run_manifest.json` が記録する `tool_version: "0.22.0"`(2026-08-10 実行)は、**未 release の dev 状態**で走ったものである。実行当時の `__init__.py` リテラルは 0.22.0 だが、それはブランチ `agent/issues-225-230-231-followup` 上のローカル値で、git tag も GitHub release も存在しない。さらに実行当時 `pyproject.toml` は静的 `version = "0.21.0"` だったため、同じツリーから build した wheel は 0.21.0 を刻印し得た(package-version レベルの provenance は曖昧)。このため同 manifest が保証するのは artifact レベル、すなわち同 manifest の `file_sha256` が固定する配布物の同一性であり、`tool_version` は「その時のソースの `__version__`」と読む。実行時の commit SHA は manifest に記録されておらず、branch 名を除いて復元不能である。`#266`(commit 2192b3e)により version は一元化されたため、以後の実行ではこの曖昧さは解消される。fixture は書き換えずそのまま保持する(既存テストが checksum を固定している)。

### release checklist

1. `src/candidate_qtl/__init__.py` の `__version__` を bump(唯一の変更点)
2. `python scripts/generate_cli_reference.py` で `docs/CLI-REFERENCE.md` を再生成
3. フルテストを実行し、`python -m coverage report` が `fail_under` を下回らないことを確認
4. `python -m build` → clean venv へ wheel install して `candidate-qtl --version` / import / weight registry 同梱を確認(ci.yml の `build` job と同じ手順)
5. 凍結 benchmark を再実行する場合のみ `run_manifest.json` の `tool_version` と commit SHA を更新
6. git tag と GitHub release を作成
7. `tests/test_version_single_source.py` が README の `Version:` 行と `CITATION.cff` の `version` に一致することを確認

## Lint 方針

`ruff` は現状 pyflakes 系(未使用 import/変数、未定義名など)と構文エラーのみを対象とする。style ルールを増やす場合は `pyproject.toml` の `[tool.ruff.lint] select` を段階的に広げ、既存コードを併せて整える。導入時点で、重複 dict キーや dead 変数などの実バグを検出・修正した。
