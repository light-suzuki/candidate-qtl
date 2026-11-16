# このPC・既存ツールへの接続

## 確認済みのローカル環境

2026-07-10に確認した時点で、WindowsネイティブPATHには blastn、makeblastdb、primer3_core はない。一方、WSLには blastn、makeblastdb、primer3_core があり、minimap2 も利用可能である。未公開ゲノムの解析・primer特異性確認はWSLを実行環境にする。

## QTLift

ローカル配置：`/path/to/QTLift`

QTLiftは、旧assemblyまたは別accession上の `chr:start-end` を、anchor gene/synteny、任意のmarker配列、必要ならminimap2による全ゲノムevidenceで別assemblyへ移す用途に使う。候補先行ワークフローでは **Step 0: 区間座標の整合性確認** に置く。

- source/targetそれぞれのFASTA + GFF3が必要。
- chromosome/scaffold名を揃える。
- 大きなゲノムでは `.fai` とBLAST DBを推奨する。
- QTLiftの結果がMedium以下、split、Manual-checkなら、候補順位付けを始めずにmarkerまたはassemblyを再確認する。
- 元研究の座標は上書きせず、`templates/qtl_sources.tsv` にsourceとtargetの両方を保存する。

QTLiftは因果遺伝子を同定するツールではない。候補抽出に入れる座標の版違いを防ぎ、異なる研究のQTL・候補遺伝子位置を同じtarget assembly上で比較するツールである。

### ローカルAPI連携

QTLiftを起動するとローカルREST APIを利用できる。`scripts/qtlift_client.py` はAPIへjobを送り、入力と結果をJSONとして固定する。

~~~powershell
python .\scripts\qtlift_client.py health
python .\scripts\qtlift_client.py scan --genome-root C:\private\genomes --out C:\private\qtl\genomes.json
python .\scripts\qtlift_client.py submit `
  --genome-root C:\private\genomes `
  --target-ref RefA `
  --source-ref RefB `
  --contig Chr1 --start 100000 --end 500000 `
  --peak 300000 --name QTL_example `
  --wait --out C:\private\qtl\QTL_example.qtlift.json
~~~

marker配列を追加する場合はJSON objectとして渡す。

~~~powershell
python .\scripts\qtlift_client.py submit `
  --genome-root C:\private\genomes --target-ref RefA --source-ref RefB `
  --contig Chr1 --start 100000 --end 500000 `
  --markers-json '{"left":"ACTG...","right":"TGCA..."}' `
  --wait --out C:\private\qtl\QTL_example.qtlift.json
~~~

API clientはゲノム配列を外部へ送らない。QTLift serverが同じPC上のローカルファイルを読む。

### 変換済みQTLの研究間比較

QTLiftで同じtarget assemblyへ移した後、`scripts/compare_qtl_sources.py` で区間重複、reciprocal overlap、Jaccard係数、target peak間距離を計算する。

~~~powershell
python .\scripts\compare_qtl_sources.py `
  --input C:\private\qtl\qtl_sources.tsv `
  --out C:\private\qtl\qtl_pairwise_comparison.tsv
~~~

同じ形質名の記録だけを比較する場合は `--same-trait-only` を付ける。ただし文字列の完全一致だけなので、trait ontology正規化前の補助フィルターである。

出力の `high_reciprocal`、`nested`、`partial` は区間形状の分類にすぎない。区間が重なっても同じ原因座とは断定せず、形質定義、効果方向、交配親、環境、marker、anchor gene、QTLift confidenceを併記する。

## AI非依存の候補圧縮

標準経路は候補proteinだけをqueryとして、次を組み合わせる。

- BLASTP/DIAMOND：モデル植物・近縁種への相同性
- HMMER/Pfam：domain
- eggNOG-mapper：orthologyとGO
- Europe PMC：候補別の文献検索
- 任意のOrthoFinder：厳密なorthogroup、paralog関係
- 任意のMCScanX/jcvi：シンテニー
- 対応作物がある場合のKnetMiner：knowledge graph探索

詳細は `docs/KNOWLEDGE-PRIORITIZATION.ja.md` を参照する。KnetMiner、Ensembl Plants、Grameneなどの外部サービスは、対応種・assembly・API安定性が異なるため、現段階では固定依存にせずspecies adapter候補として扱う。

## primerblast-oss

GitHub: https://github.com/light-suzuki/primerblast-oss

primerblast-oss はPrimer3とローカルBLAST+を使い、未公開ゲノムを外部サービスへ送らずに以下を行える。

- candidate geneまたはQTL intervalのprimer設計とoff-target確認
- in-silico PCR
- intervalに沿ったgrid marker設計
- GFF3/VCF/intervalを用いるbreeding assay
- CAPS/dCAPS候補の探索

| 段階 | 入力 | 使う目的 |
|---|---|---|
| 組換え体の再分解 | interval + genome + BLAST DB | grid markerを作る |
| 親間候補変異の検証 | gene/amplicon + genome | Sanger/PCR primerの特異性を確認 |
| 変異marker化 | SNP/Indel + VCF + genome | CAPS/dCAPSまたはIndel assay候補を出す |
| アッセイ採用前 | candidate primer pair + 各親ゲノムDB | multi-reference off-targetを確認 |

ツール出力のA–D/riskは設計優先順位であり、実験PCRを置き換えない。両親・F1・代表F2で実測してから全F2へ展開する。

## 最小データレイアウト（Git管理外）

~~~text
project-data/
  reference/                        # FASTA, GFF3, fai, BLAST/DIAMOND/Pfam DB
  qtl/
    qtl_sources.tsv                 # source/target assemblyとQTLift confidence
    qtlift_jobs/                    # QTLift JSON
    qtl_pairwise_comparison.tsv     # 研究間の区間比較
  annotation/                       # ortholog/domain/GO/synteny evidence
  parents/                          # parental reads and Sanger
  f2/                               # phenotype and DNA manifest
  results/                          # ledger and marker calls
~~~

Gitに残すのは、入力ファイルのパスを伏せたmanifest、ハッシュ、実行コマンド、結果要約、primer配列、解析版だけにする。
