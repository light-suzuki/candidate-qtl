# 調査LLMと機械処理の役割分担

元の研究者が出典の意味を読み、機械はその記録を保全・照合する。科学的な判断規則は[canonical Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md)、実行順序は[親Agentの手順](AGENTIC-QTL-CASE.ja.md)を正本とする。

| 工程 | 元の調査LLM・親Agentの判断 | 機械が担当すること |
|---|---|---|
| A | 表現型から機構・既知gene・参照を調べ、原報を解釈する | 記録形式、参照identity/配列、freeze後の区間投影を照合する |
| B0/B | B担当は機能・文献と表現型の関係を判断する | B0は表現型なしで全gene/copy台帳を整理。資料の圧縮・配送、source/member/hashを管理する |
| 論文記録 | 論文ID/タイトル、報告gene名/ID/種、9回答、理由・unknown・反証を記録する | 原値を転記し、書誌・IDを別途監査し、既存決定表でTierを算出する |
| Direct C | 原sourceのgene・種・形質・実験の対応を確認する | strict source/identity binding、freeze、projection、memberの条件を検査する |
| Wet境界 | 封印済みv3では元担当が既存候補の採用・保留を理由付きで判断する | policy・roster・hashを検査し、新Wetと全gene台帳を生成する |

## 遺伝子名・IDを機械で埋めない理由

研究者が「論文Pはgene名A・ID Xを扱う」と記録し、DBから取得したXの正式名・別名と**別の値として**比較することで、取り違えを監査できる。

DBのID Xから名前Aも自動取得し、その二つを研究者の回答欄へ入力すると、照合は同じDBレコードの自己確認になる。実在する名前・IDの組を作れても、元のLLMが論文・種・gene・表現型を取り違えたことは検出できない。**この自動入力をハルシネーション対策とは扱わない。**

機械取得した正式名・別名・IDは監査側へ保存し、`reported_gene` / `reported_gene_id` / `reported_species`を上書きしない。原調査でIDを確認できなければ、確認範囲と理由付きのunknownを残す。名前が分かるのにIDを取得できないことだけで候補を下げたり、探索を際限なく追加したりしない。

元担当が明示した値の機械転記・保存・再利用は可能。台帳IDやhashの生成も機械の担当である。禁止しているのは、科学的な原回答を別の推測やDB値で補って「元研究者が確認した」と扱うこと。詳細は[軽い遺伝子ID照合](CANONICAL-WET-MEMBERSHIP.ja.md#軽い遺伝子id照合)を参照する。

## Tier・支持・採否を混ぜない

Tierは形質との関連度、verified supportは確認済み証拠の属性、Wetは今回実験する候補の採用集合である。source照合の成功やBLAST類似性から、表現型の因果関係を自動生成しない。報告gene・配列proxy・対象座位を分け、妥当な複数copyと未解決状態を保持する。

元LLMの不明点や具体的な矛盾は既存のbounded follow-upへ戻す。writerが9回答を作り直したり、別のTier専用reviewerで全件を再調査したりしない。

## 現在のWet境界

`upper_evidence_scoped_v3`では、累積25を初めて超えるTierだけを元担当の比較へ渡す。20/25は機械的な切断上限ではない。理由付きの`defer`は今回の実験を後回しにする判断で、非原因の証明ではない。未回答は保持・未完として残す。

この採否判断は「検証済み反証があるときだけ候補を除く」という古いrefinement説明と区別する。生物学的排除、support昇格、原freezeの改変にはそれぞれの検証が必要。旧policy・generic consumerのreview要件は旧世代のまま維持する。
