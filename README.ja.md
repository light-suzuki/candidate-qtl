# candidate-qtl

[English](README.md) | **日本語**

`candidate-qtl` は、植物のQTL区間と表現型から、Wet実験で検証できる規模の候補遺伝子集合を作る研究用workflowとCLIです。モデル植物や他植物の文献・相同遺伝子の知見を使って原因候補を失わずに残し、各判断の根拠を保存します。約20遺伝子は目安であり、最低件数・切断上限・完了条件ではありません。mappingは独立した証拠として保持し、map-based cloningを置き換えたり原因遺伝子を自動で断定したりはしません。

## 全体の流れ

```text
表現型 + QTL区間 + 参照FASTA/GFF3
  ├─ A経路（必須）  phenotype-first：機構・family・参照遺伝子の仮説
  │                 → 出典・配列の検証 → 区間への決定的な投影
  ├─ C経路（必須）  Direct C：区間を見ずに形質と遺伝子/座/変異の既報を回収
  │                 → freeze → 区間への投影
  └─ B経路（任意）  interval-first：区間の全遺伝子 → family/orthogroup単位へ圧縮
                    → 意味解釈 → 検証
        各経路は独立にfreeze
  → claim review → 論文ごとのチェック表（9回答） → Tier T1〜T6/U
  → canonical Wet membership → wet_candidates.tsv
       review予算を最初に超えたTierだけを境界レビューへ
```

- **A経路**は表現型から機構と参照遺伝子を考え、出典・accession・配列を確認したうえで区間へ機械的に投影します。
- **C経路（Direct C）**は区間や他経路の結果を見ずに、形質と遺伝子・座・変異の直接の既報を集めて投影します。検索した結果の空と、未実行を区別して記録します。
- **B経路**は任意です。区間の全遺伝子をfamily/orthogroup単位に整理し（中間表現であり「1 familyにつき1遺伝子」へ縮約するものではありません）、表現型との関係を解釈します。
- **独立性**：A・B・Cは各々freezeするまで手続き上独立に保ちます。freeze後の一致は来歴であり、証拠数の加算やscore bonusではありません。
- **未実施は報告し、作らない**：未実行・失敗は`not_evaluated`であり、生物学的陰性ではありません。Bを使わない場合、WetにはA/Cの候補だけが入り、注釈で数を埋めません。候補が無い・少ないときはhandoffが`route_a_rethink`を出し、親がAを新しい世代で一度考え直します。

## Tier判定

原調査担当（人間またはLLM agent）が論文ごとに、論文ID/タイトル、報告遺伝子の名前・ID・種、9項目のyes/no/unknown回答と短い理由を記録します。決定表が回答からTierを計算し、遺伝子のTierは適用できる最も強い論文で決まります。論文数による加点はありません。

| Tier | 意味 |
|---|---|
| T1 | 同種で同一または非常に近い表現型の既報 + 妥当な区間内相同投影 |
| T2 | 同じ・非常に近い細胞/発生現象についての遺伝子–表現型の既報 |
| T3 | 別の表現型の既報を、具体的な共通過程で橋渡し |
| T4 | 報告された機能に基づく仮説 |
| T5 | 出典付きの具体的な注釈に基づく仮説 |
| T6 / U | identityのみ / 未解決。全遺伝子台帳に残し、陰性とは扱わない |

Tierは今回の表現型との関連度であり、因果確率ではありません。関連度、verified support、Wet採用、調査の実施範囲は別の軸です。機械はID・配列・引用・引用箇所・区間対応・hashを照合しますが、研究者が記録した遺伝子名/IDをDB検索の値で置き換えません。詳細は[canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md)と[役割分担](docs/LLM-BOUNDARY.ja.md)を正本とします。

## 出力

集合と一緒に、未解決・保留の状態と実際の調査範囲を報告してください。

- `wet_candidates.tsv` — 採用したWet集合（唯一の「最終候補集合」）
- `wet_relevance_inventory.tsv` — 区間の全遺伝子と関連度・評価状態の台帳
- `wet_candidate_trace.tsv` — 遺伝子ごとのhandoff観測値

## リポジトリ構成

| パス | 内容 |
|---|---|
| `src/candidate_qtl/` | Pythonパッケージと`candidate-qtl` CLI |
| `workflow/` | 通常案件の入口`practical_packet_run.py`と設定テンプレート`practical_packet_run.ordinary.template.json` |
| `scripts/` | 親Agent用の補助（Route Bの分割、claim review session、境界レビュー）とローカルCI |
| `docs/` | 契約・実行手順・CLI文書（主に日本語） |
| `templates/`、`config/` | 入力テンプレートと設定例 |
| `tests/` | 合成fixtureを使うunit/integrationテスト |

## インストール

- Version: `0.24.0`
- Python 3.12+。FASTA/GFFを扱う解析はLinuxまたはWSLと`alignment` extraを使います。

```bash
git clone https://github.com/light-suzuki/candidate-qtl.git
cd candidate-qtl
python3 -m venv .venv
source .venv/bin/activate
pip install -c constraints.txt -e ".[alignment]"
candidate-qtl doctor --json
```

Windowsでは`./scripts/install-local.ps1 -AddToUserPath`も使えます。

## 通常案件の実行

親Agent（または人）が案件を持ち、A・C・任意のBを、親の会話履歴を引き継がない別の担当へ分配します。下のwrapperは保存済み回答の配送と検査を行うもので、調査担当の起動や検索ツールの提供はしません。

```bash
h=workflow/practical_packet_run.py
run=/path/to/run
cp workflow/practical_packet_run.ordinary.template.json "$run/config.json"   # 編集して使う
python $h prepare --config "$run/config.json" --out-dir "$run/01-prepared"
# Bを使わないA/C優先：全Bプロンプトをnot_evaluatedにしてからreview
python scripts/agentic_route_b.py merge --prompts "$PROMPTS" --out "$run/b-not-run.json"
python $h review --prepared "$PREPARED_MANIFEST" --sha256 "$PREPARED_SHA" \
  --responses "$run/b-not-run.json" --route-b-not-run "<理由>" --out-dir "$run/02-reviewed"
python $h finish --reviewed "$REVIEWED_MANIFEST" --sha256 "$REVIEWED_SHA" \
  --route-c-members "$C_MEMBERS" --route-c-manifest "$C_MANIFEST" --out-dir "$run/03-first-wet"
```

各段は次段に渡す`manifest`と`sha256`を出力します。Bを実行する場合は実際のshard回答をmergeし、`--b-execution-manifest`を渡します。

A経路の回答形式、Direct C（`route-c-direct-freeze` → `route-c-direct-project` → `route-c-direct-members`）、`scripts/agentic_route_b.py`によるBの分割、境界レビューは[通常案件の手順](docs/AGENTIC-QTL-CASE.ja.md)を参照してください。コマンドの詳細は[CLI案内](docs/CLI.ja.md)と[CLIリファレンス](docs/CLI-REFERENCE.md)にあります。

## テスト

```bash
python scripts/run_test_partition.py --partition core
scripts/local_ci.sh test lint
```

## 文書

| 知りたいこと | 文書 |
|---|---|
| 目的と全体構成 | [CORE-README](CORE-README.ja.md) / [core contract](docs/CORE-CONTRACT.ja.md) |
| Tierと採否の規則 | [canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md) |
| 案件の実行 | [通常案件の手順](docs/AGENTIC-QTL-CASE.ja.md) / [三経路の手順](docs/TRIPLE-ROUTE-RUNBOOK.ja.md) |
| Direct C | [Direct Cの契約](docs/ROUTE-C-DIRECT-RESCUE.ja.md) |
| coding/research agent向けの規則 | [AGENTS.md](AGENTS.md) |

## 現状と限界

開発中の事例で既知の原因遺伝子が残ったことは診断的な証拠にとどまり、未見事例での精度は確立していません。20枠を埋めても未知の原因遺伝子の保持は保証されません。ベンチマーク用データと評価ツールはこの公開版に含みません。

## ライセンスと引用

[MIT License](LICENSE)。引用情報は[CITATION.cff](CITATION.cff)、貢献者は[CONTRIBUTORS.md](CONTRIBUTORS.md)を参照してください。
