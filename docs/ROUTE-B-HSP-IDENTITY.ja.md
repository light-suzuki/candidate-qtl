# Route B のHSP集約一致率

ローカルBLASTの入口は `query_residue_consensus_v1` を使用する。
`percent_identity` と `query_coverage` は同じ**重複を除いたquery残基集合**に基づく。
BLAST標準のalignment-column `pident`とは区別する。

- query残基が同じsubject残基・文字へ重複して対応する場合は1回だけ数える。
- 対応位置または文字が重複HSP間で食い違う残基は、保守的に非一致とする。
- subject側gapはquery残基の非一致。query側gapにはquery残基が無く、分母から除く。
- 同一アラインメントの分割・重複・行順変更で、一致率とcoverageの判定を変えない。
- `max_hsp_identity`は別列の測定値として保存するが、採否には使用しない。

200残基×40%一致と800残基×25%一致が非重複・gapなしである場合、
1,000残基の280残基が一致するため28%となる。最大値40%ではない。

`route_b_local_blastp_alignment.tsv`にはqseq/sseqを含む原本を保存する。
`route_b_local_blastp_pair_evidence.tsv`には分母・一致数・conflict数・policy・採否を保存する。
既存family-index向けの13列は別の`route_b_local_blastp.tsv`に保持する。
各ファイルのhashをpreparation manifestのhomology providerへ記録する。

新しい測定政策であり、過去benchmarkの値は書き換えない。新generationで比較する。
この変更は局所配列相同性の判定であり、orthologyや形質特異的supportを認定しない。
既存外部orthology表は変更しない。旧13列だけから重複残基のidentityを推測しない。

## 検証

`python -m unittest tests.test_route_b_blast_evidence tests.test_issue597_merged_hsp_reference_identity -v`

合成アラインメントとmock BLAST subprocessを使う回帰テストである。
実データのWet保持率・実BLAST全体の精度検証とは区別する。

出力fieldの定義: https://www.ncbi.nlm.nih.gov/books/NBK279684/
