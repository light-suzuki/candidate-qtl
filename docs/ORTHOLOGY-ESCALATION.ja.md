# QTL区間限定解析とorthology精度の両立

## 結論

このワークフローで限定するのは **query** であり、検索先までQTL区間へ限定しない。

- query: QTL区間内のcandidate proteinだけ
- focal search DB: 対象種の全proteome。1回作成し、複数QTLで再利用する
- reference search DB: 少数のモデル植物・近縁種proteome。1回作成し、再利用する
- gene-family tree: 全候補には作らず、上位かつorthologyが曖昧な候補だけ

対象種全proteomeを検索先に含める理由は、QTL区間外にある近縁paralogを見落とさないためである。区間内配列だけで相互検索すると、実際には多コピーの遺伝子を誤って1:1 orthologと判定し得る。

計算量削減の単位は「全ゲノムを解析しない」ではなく、**全proteome同士のall-vs-allをQTLごとに繰り返さず、少数のcandidate queryだけを既製・再利用DBへ当てること**とする。

## Core i7・RAM 16GB向け三段階運用

### Tier 1: 全候補に行う軽量解析

対象: QTL区間内の全候補

目的: 厳密なortholog確定ではなく、機能クラスと既知情報への入口を得る。

推奨処理:

1. representative transcriptを1本選ぶ
2. DIAMONDまたはBLASTPで4〜8種程度のreference panelへ検索
3. 対象種全proteomeへ自己検索し、近縁paralog数とbitscore差を記録
4. HMMER/PfamまたはInterProScanでdomain architectureを取得
5. eggNOG-mapperまたはOMAmerで既存ortholog group/HOGへ配置
6. Europe PMC、作物DB、モデル植物DBの検索語を作成

出力ラベル:

- `homology_only`: 類似配列はあるがorthology未確定
- `orthology_probable`: RBH、HOG、domain、syntenyの複数根拠が整合
- `orthology_ambiguous`: 多コピー、上位hitが近い、domain差、系統特異的重複など
- `no_reliable_assignment`: 十分な根拠なし

Tier 1では、BLAST best hitをorthologと断定しない。

### Tier 2: 上位候補に行う相互検索と軽量な遺伝子ファミリー確認

`orthology-plan`は一次順位のTop/Rescueから指定数だけを抽出し、近縁種proteomeへの順検索、上位hit抽出、対象種全proteomeへの逆検索を行うscriptを生成する。`orthology-import`はRBHと非RBHをcanonical evidenceへ戻す。シンテニーは `synteny-import`（MCScanX/jcvi）で別途取り込む。

対象:

- Dry順位のTop/Reserve候補
- `orthology_ambiguous`
- paralog冗長性が原因遺伝子解釈に影響する候補

処理:

1. 各reference speciesから上位hitを複数取得する
2. reference hitを対象種全proteomeへ逆検索する
3. RBH、一対多、多対多を記録する
4. 必要ならMCScanX/jcviまたはQTLiftのanchor情報で局所シンテニーを確認する
5. OMAmer HOG、eggNOG OG、RBH、syntenyを別列で保存する

Tier 2で根拠が一致すれば、通常は系統樹を作らず親間direct sequencingへ進めてよい。

### Tier 3: 曖昧なTop候補だけにtargeted gene-family tree

対象:

- Top候補である
- 複数paralogが同程度に近い
- モデル植物機能をどのコピーへ移すかでWet優先順位が変わる
- duplication前後の機能分化が疑われる

処理:

1. 近縁種・モデル種5〜10種から、各種上位3〜10配列を取得
2. domain architectureと配列長で明らかな別familyを除く
3. MAFFTでアラインメント
4. ClipKITまたはtrimAlで必要に応じてtrim
5. FastTreeで初期tree、必要な場合だけIQ-TREEで再解析
6. species tree、duplication、syntenyを併用してortholog/paralogを解釈

16GB運用では、一度に1 familyを処理し、原則として数十〜200配列程度へ抑える。NLR、RLK、F-box、巨大転写因子familyなどで配列数が膨らむ場合は、domain構成、近縁種、syntenic locus、HOG subfamilyを使って先に絞る。

## エスカレーション条件

以下のいずれかを満たす場合だけTier 3へ進む。

- 対象種内の第1・第2hitのbitscore差が小さい
- reference geneから対象種への逆検索で複数コピーが同程度
- RBHとHOG assignmentが一致しない
- domain architectureがbest hitと異なる
- 1:manyまたはmany:many orthology
- 対象系統でgene family expansionがある
- シンテニーとsequence similarityが異なる候補を示す
- LLM/文献調査でparalog固有機能が重要と判断された

これらがなく、RBH・HOG・domain・機能が整合する場合は、`orthology_probable` として次へ進む。

## reference panel

全植物を毎回入れない。対象種に応じて4〜8種程度を選ぶ。

- 最も近い利用可能な注釈済み種: 1〜2種
- 同科または同目の高品質ゲノム: 1〜2種
- 機能情報が豊富なモデル植物: 1〜2種
- 単子葉/双子葉など必要な外側比較: 1種
- tree用外群: 必要時のみ1種

reference panelとDB版はmanifestに固定する。

## 16GB RAM運用ルール

- CPU threadは原則4、最大でも物理コア数以下
- representative transcriptのみを標準入力にする
- query FASTAはQTL候補だけ
- reference/focal proteome DBは作成後に再利用する
- family treeは逐次処理し、並列に複数familyを走らせない
- DIAMOND/MMseqs2の一時ファイルを保存する十分なdisk容量を確保する
- local eggNOG/InterPro DBが大きすぎる場合はweb serviceまたは必要なclade DBを使う
- 中断再開可能な出力とcacheを優先する
- full OrthoFinderをQTLごとの標準処理にしない

## OrthoFinderの位置づけ

OrthoFinderはcomplete proteomeを種ごとに入力し、orthogroup、gene tree、species tree、duplication、orthologを包括的に推定するため、精度確認には有用である。一方、QTLごとに複数種の全proteomeで実行すると、本ワークフローの軽量性と競合する。

使用する場合は次のどちらかとする。

1. 研究室で一度だけ小さなclade coreを作り、複数QTLで再利用する
2. OrthoFinder v3のcore/assign機能を利用し、対象種を既存coreへ追加する

coreを作れない場合、標準経路はOMAmer/eggNOG placement + targeted RBH + 必要時のfamily treeとする。

## 類似OSS・サービスの役割

| ツール | 得意な作業 | このワークフローでの位置 |
|---|---|---|
| QTG-Finder | QTL内原因遺伝子の機械学習順位付け | 対応種では比較用baseline。マイナー種の中核にはしない |
| KnetMiner | region/gene listとtrait keywordからknowledge graph探索・ranking | 対応speciesがある場合の強力な外部prior |
| eggNOG-mapper | 既存ortholog groupとphylogenyを使う機能注釈 | Tier 1標準候補 |
| OMAmer | precomputed HOGへの高速・tree-informed配置 | マイナー植物・低計算環境で有力なTier 1標準候補 |
| OrthoFinder | 複数complete proteomeから包括的orthology・gene tree推定 | 一回限りのcore構築または難しい候補の検証 |
| fDOG | feature/domain architectureを考慮したdirected ortholog search | Top候補の厳密検索、未注釈assemblyへの探索 |
| DIAMOND / MMseqs2 | 高速protein similarity search | 全候補の基本検索とfamily sequence回収 |
| HMMER / InterProScan | domain・family・functional site | 全候補の機能分類とparalog誤認防止 |
| MAFFT + FastTree/IQ-TREE | targeted family alignment/tree | Tier 3のみ |
| MCScanX / jcvi | collinearity・synteny | sequenceだけで曖昧な候補の補強 |

## このワークフローの独自範囲

既存OSSは各作業を個別に実行できるが、次を一つの低予算ゲートへまとめる部分が本リポジトリの役割である。

- QTL区間からqueryを限定する
- 軽量解析から厳密解析へのエスカレーション条件を固定する
- 根拠を同じcandidate ledgerへ統合する
- Top/Reserve/Rescueへ分類する
- 親間direct sequencing用の対象遺伝子・領域へ接続する
- Wet結果をDry順位とは分離して記録する
