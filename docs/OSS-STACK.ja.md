# Win/WSL OSSスタック

## このPCでMVPに使うもの

| 役割 | OSS | このPCの状態 | 組込み位置 |
|---|---|---|---|
| 局所FASTA random access | pysam | 利用可 | candidate-qtl extract |
| 配列操作 | Biopython | 利用可 | CDS/逆鎖の抽出 |
| 表形式・要約 | pandas | 利用可 | ledger/F2集計 |
| primer設計 | Primer3 | WSLで利用可 | primerblast-oss経由 |
| 特異性確認 | BLAST+ | WSLで利用可 | primerblast-oss、QTLift |
| assembly比較 | minimap2 | WSLで利用可 | QTLiftの独立evidence |
| ドメイン検索 | HMMER | WSLで利用可 | Pfam DBを配置後の任意段階 |
| 座標移送 | QTLift | Windows配置済み | Step 0 |

## 追加すると価値が高いが、MVPの必須にしないもの

| OSS | 追加する時 | 理由 |
|---|---|---|
| primerblast-oss | marker設計を実行する時 | interval/grid marker、CAPS/dCAPS、multi-reference specificity |
| bcftools | 親VCFやtargeted-seq VCFが得られた時 | 変異抽出・正規化・フィルタ |
| SnpEff または VEP | VCFを多数扱う時 | coding/splice等の注釈を機械化 |
| R/qtl2 | QTL検出・再解析まで一体化する時 | F2遺伝地図・QTLモデルを再計算 |
| Snakemake | 複数QTL・複数交配を繰り返す時 | 入力ハッシュ・依存関係・再実行の管理 |
| OrthoFinder | 近縁種protein setを揃えられる時 | orthogroupをBLAST単独より堅くする |
| eggNOG-mapper / InterProScan | 大量候補の機能注釈を揃える時 | GO/ドメインを一括付与 |

## 実行順

1. QTLiftでassemblyと座標を固定する。
2. candidate-qtl extractでQTL区間だけを抽出する。
3. BLAST/orthology/HMMER/文献をcandidate ledgerへ根拠URLまたはDOI付きで入れる。
4. candidate-qtl rankで、変更可能で説明可能な実験順にする。
5. primerblast-ossで親間Sanger/marker assayを設計する。
6. candidate-qtl cosegregationでF2またはF2:3を集計し、ledgerを更新する。
7. 複数QTLを反復する段階でSnakemakeを追加する。

## データ境界

Gitへ入れるのはコード、TSVヘッダー、ハッシュ、コマンド、匿名化した結果要約だけである。未公開ゲノム、VCF、表現型、Sanger chromatogram、BLAST DBはWindows/WSLのprivate data rootに置く。
