# QTL区間限定 annotation rescue

`annotation-rescue` は、source-backed homolog proteinをQTL区間配列だけへ投影し、新しいreference GFF3に無い・壊れた・部分的なmodelを別レーンで監査する。既存候補の順位、score、fixed primary 20は変更しない。

## 実行

~~~powershell
candidate-qtl annotation-rescue --fasta genome.fa --candidates candidates.tsv `
  --query-proteins homologs.faa --query-metadata annotation_rescue_query_metadata.tsv `
  --interval chr1:100000-200000 --out-dir rescue --execute --json
~~~

miniprotが無い場合もcommandは完走し、`run_local_annotation_rescue.sh` と `not_evaluated` reviewを出す。これは`unmapped`ではない。既計算GFFをimportする場合、全queryを含む成功runだと検証できる時だけ`--projection-complete`を付ける。付けなければ、GFFに無いqueryを不存在扱いしない。

## 出力と判定境界

- `annotation_rescue.tsv`: 全projectionと全queryのledger。`native`、`native_discordant`、`reference_rescue`、`ambiguous_overlap`、`unmapped`、`not_evaluated`を分離する。
- `annotation_rescue_review.tsv`: 上限付きreview appendix。全ledgerは削らない。
- `annotation_rescue_review_overflow.tsv`: appendix上限外の行。
- `annotation_rescue_manifest.json`: query/FASTA/candidate/projection fingerprint、command、tool状態、外部送信、閾値、科学的限界。

全行の`promotion_status`は`not_promoted`である。near-complete projectionは「query coverageが高く、miniprotがframeshift/stopを報告しなかった」という計算結果であり、機能するORFや因果性の証明ではない。no-hitもgene absenceの証明ではなく、配列分岐、assembly gap、閾値、query選択が代替説明として残る。

`high` qualityはunique rank-1、query coverage 0.8以上、identity 0.5以上、CDSあり、reported disruptionなし、非canonical splice報告なし、区間端非接触を満たすscreening policyである。この閾値は生物学的真理ではない。

miniprotはprotein-to-genome alignmentでsplicingとframeshiftを扱い、GFF3へ`Target`、`Rank`、`Identity`、`Positive`、`Frameshift`、`StopCodon`等を出す。[公式repository](https://github.com/lh3/miniprot)、[manual](https://lh3.github.io/miniprot/miniprot.html)、[原著](https://pmc.ncbi.nlm.nih.gov/articles/PMC9869432/)を参照。
