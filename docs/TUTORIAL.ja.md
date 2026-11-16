# 再現可能チュートリアルと既知QTGベンチマーク

`docs/tutorial/` に再配布可能な小型fixture（人工ゲノム・GFF3・protein・注釈・発現・QTL source）を置く。CIはこのチュートリアルをオフラインで最後まで実行し（`tests/test_tutorial.py`）、既知の原因遺伝子が上位へ回収されること、および回収できないケースを**失敗として記録する**ことを検証する。

## fixture

| ファイル | 内容 |
|---|---|
| `genome.fa` / `genes.gff3` / `proteins.fa` | chr1(900 bp)に3遺伝子。`GENE_CAUSAL` が原因遺伝子、`DECOY1`/`DECOY2` はダミー |
| `candidate_eggnog.tsv` | `GENE_CAUSAL` のorthogroup/GO/機能記述 |
| `expression.tsv` | 対象組織(meristem)特異発現。原因遺伝子のみ濃縮 |
| `manual_evidence.tsv` | 近縁種の既知QTG（`known_qtg`、source ID付き） |
| `qtl_sources.tsv` | 2つのassemblyから同じ模擬QTLを参照へliftした比較例 |
| `known_qtg_benchmark.tsv` | ベンチマーク定義（原因遺伝子・変異クラス・期待回収） |

## オフライン解析チェーン

外部DB・ネットワーク・BLAST/primerblastを必要としない範囲でチュートリアルを完走できる。Windowsでは各コマンドを `.\scripts\candidate-qtl.ps1` に置き換える。

~~~bash
# 1. QTL区間から候補を抽出
candidate-qtl extract --fasta docs/tutorial/genome.fa --gff3 docs/tutorial/genes.gff3 \
  --interval chr1:1-900 --out-dir out/extract

# 2. 既知知識を1つの証拠台帳へ統合（発現・orthology・既知QTG）
candidate-qtl evidence-build --candidates out/extract/candidates.tsv \
  --eggnog docs/tutorial/candidate_eggnog.tsv \
  --expression docs/tutorial/expression.tsv --target-tissue meristem \
  --manual-evidence docs/tutorial/manual_evidence.tsv --out-dir out/evidence

# 3. source-backed evidenceだけでWet前順位付け
candidate-qtl evidence-rank --candidates out/extract/candidates.tsv \
  --evidence out/evidence/evidence.tsv --out-dir out/rank

# 4. 別assembly間のQTL座標整合
python scripts/compare_qtl_sources.py --input docs/tutorial/qtl_sources.tsv --out out/qtl_compare.tsv
~~~

`annotate-plan`・`research`・`primer-plan`・`orthology-tier1-plan` などはBLAST/DIAMOND/eggNOG/PrimerBlast OSSやAPIを要するため、実データでの追加ステップとして `docs/` の各章に従う。

## 期待される出力

- `out/extract/candidates.tsv`：`GENE_CAUSAL` / `DECOY1` / `DECOY2` の3行、いずれも `complete_orf`。
- `out/rank/ranked_candidates.tsv`：`GENE_CAUSAL` が `dry_group=top`・`dry_rank=1`・`top_eligibility=strong_genetic_evidence`。decoyは `evidence_insufficient`。
- `out/rank/weights_used.tsv`：使用重みが固定され、指標が再現可能。証拠のsource ID（DOI・GO・発現matrix版）は台帳に残る。
- `out/qtl_compare.tsv`：2 assembly由来の同一QTLが1件の比較として整合。

## 失敗ケースの記録

ベンチマークは原因遺伝子を外したケースも隠さず記録する。

- `tut_coding`（`GENE_CAUSAL`, coding）：機能・遺伝学的証拠があり **top-1で回収**。
- `tut_regulatory`（`DECOY2`, regulatory）：coding証拠が無いため **Topに入らず miss として記録**。

これは「AIや知識解析があれば必ず当たる」ことを主張しないためで、`top-1/top-5 recall`・`median rank`・`Wet対象数`・変異クラス別（coding/regulatory/SV）の失敗率を、既知遺伝子バイアスと分けて評価する枠組みを与える。LLM補助あり/なしは別集計とし、LLM出力はsource ID付きで人が確認した根拠だけを台帳へ戻す。
