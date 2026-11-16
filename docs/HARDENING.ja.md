# 実データ利用前の安全化

この文書は、候補順位付け、親間配列比較、PrimerBlast OSS、LLM連携を実データへ使う際の安全条件をまとめる。

## 1. 親間・assembly間配列比較

`allele-compare` は minimap2 の最良PAF alignmentについて、次を記録する。

- reference coverage
- alternate coverage
- alignment identity
- mapping quality
- alignment strand
- alignment本数

referenceまたはalternateのcoverageが90%未満の場合、結果は `partial_alignment` とする。整列部分のSNV/InDelは出力するが、「変異なし」を配列全体の一致とは解釈せず、genomic marker候補の出力を停止する。

`--reference-candidates`を指定した場合は、入力をgene-oriented genomic FASTAとして扱う。CDS FASTAのoffsetを単純にgene startへ加算するとintronを無視した誤座標になるため、CDS比較から直接genomic marker座標を生成しない。

最良alignmentが逆向きの場合も、alternate alleleのforward-genome表現を自動確定せずmarker出力を停止する。複雑な逆位、重複、複数alignment、PAV、SVは別Issueで扱う。

## 2. trait keyword evidence

`annotation-keywords`は、best-hitの説明文にkeywordがあるだけでは採点しない。既定では次をすべて満たすhitだけを `trait_function_match` に使う。

- identity 30%以上
- query coverage 50%以上
- subject coverage 30%以上
- E-value 1e-5以下

全best hitとgate判定はaudit TSVへ残す。keyword一致は機能仮説であり、原因遺伝子証拠ではない。

## 3. 支持と反証の競合

同じ `evidence_code` にsupportとagainstがある場合、どちらか一方を行順で捨てない。各方向の最も強い1件を別々にcapして加減点し、`evidence_conflicts`へ記録する。

競合が残る候補は自動Topにしない。sourceを人が確認し、誤対応、異なるparalog、異なる器官・発達段階、gain-of-functionとloss-of-functionの違いなどを解消する。

## 4. Rescue候補

Rescueは「残り全部」ではなく、最初にWetへ送る少数の見落とし防止枠とする。既定は3遺伝子で、次を優先する。

1. CDS欠落・不完全annotation
2. target tissue/stage特異的発現
3. QTL peak近傍
4. 保存されるが機能未知
5. synteny支持

Rescue枠を超えた候補は `evidence_insufficient` に残す。除外ではなく、Top/Reserve/Rescueが外れた場合に再検討する。

## 5. exon-aware親間Sanger計画

一般的な `assay --interval` は、指定区間内の任意位置にprimerを置くため、CDS全域を覆う保証がない。

`scripts/candidate_qtl_extensions.py exon-plan` は、代表transcriptの各CDSを最大500 bpの重複windowへ分け、左右flankを含む実ゲノム配列をFASTAとして切り出す。その各templateをPrimerBlast OSSの `tile` へ渡し、coding baseとsplice-boundary近傍を親間PCR/Sangerで覆う計画を作る。

~~~powershell
python .\scripts\candidate_qtl_extensions.py exon-plan `
  --candidates C:\private\qtl\ranked_candidates.tsv `
  --gff3 C:\private\reference\genes.gff3 `
  --genome C:\private\reference\genome.fa `
  --db C:\private\blastdb\genome `
  --primerblast-root C:\tools\primerblast-oss `
  --out-dir C:\private\qtl\exon_plan `
  --top-n 7
~~~

GFFとFASTAのseqidが異なる場合は `--seqid-map GFF_ID=FASTA_ID` を使う。生成primerは親A、親B、F1でPCRと波形を確認してからF2へ展開する。

## 6. OSSからLLMへの往復

LLMへ生の候補一覧だけを渡さず、ranked candidates、canonical `evidence.tsv`、`score_contributions.tsv`をJSON packetへ埋め込む。

~~~powershell
python .\scripts\candidate_qtl_extensions.py research-packet `
  --candidates C:\private\qtl\ranked_candidates.tsv `
  --evidence C:\private\qtl\evidence.tsv `
  --contributions C:\private\qtl\score_contributions.tsv `
  --species "Pisum sativum" `
  --trait "seed shape" `
  --out-dir C:\private\qtl\research_packet
~~~

LLMの返答はそのまま採点しない。人がsourceを確認し、`reviewed=true`、レビュー者名、stable source ID、canonical evidence codeを持つclaimだけをimportする。

~~~json
{
  "gene_id": "Psat.cameor.v2.1g000010",
  "claims": [
    {
      "evidence_code": "trait_function_match",
      "direction": "support",
      "confidence": "medium",
      "value": "ortholog mutant is fasciated",
      "source_id": "UniProt:Q9LZS3",
      "description": "...",
      "reviewed": true,
      "reviewed_by": "Reviewer A",
      "reviewed_at": "2026-08-08"
    }
  ]
}
~~~

~~~powershell
python .\scripts\candidate_qtl_extensions.py llm-evidence-import `
  --candidates C:\private\qtl\ranked_candidates.tsv `
  --input C:\private\qtl\reviewed_llm_output.json `
  --out C:\private\qtl\llm_reviewed_evidence.tsv `
  --resolve-against D:\ref\arabidopsis_uniprot.fasta D:\ref\citations.tsv
~~~

### source_idの実在性

`classify_source()` は source_id の**書式**を判定する。書式が正しいことは、そのaccessionが**実在する**ことを意味しない。`UniProt:Q00000` は形は正しく、存在しない。この区別が無いと、幻覚accessionを2件並べるだけで独立グループが2つでき、`independent_source_count >= 2` を満たして `top_eligibility` が `independent_experimental_support` へ昇格してしまう。

`--resolve-against` は、source_id が実際に見つかるべきローカル資料(reference proteome FASTA、annotation TSV、引用リスト等)を指定する。オフラインで完結し、形式は問わない(tokenを収集するので、`UniProt:Q9LZS3` はFASTA headerの `Q9LZS3` にも一致する)。

| `--resolve-against` | 解決したclaim | 解決しないclaim |
|---|---|---|
| 指定あり | 従来どおり accession が独自の独立グループを持つ | `unresolvable_source_id` として拒否 |
| 指定なし | — | 受理するが `independence_group=llm_curation_source_unverified` を付与 |

指定なしの場合、未検証rowは**全件が単一の独立グループを共有する**ため、何件並べても単独では `independent_experimental_support` へ到達できない。検証を行ったか否かは結果の `source_verification`(`performed` / `not_performed`)に必ず記録され、沈黙して通ることはない。

### `reviewed=true` の限界

`reviewed=true` は**JSON内の自己申告**であり、人間がレビューした証明ではない。モデルがJSON全体を書けば自分で `true` と書ける。そこで `reviewed_by` を必須にし、誰の主張なのかを追える形にした。これは帰属であって検証ではないため、結果には常に `human_review_attestation: self_declared_not_verified` を記録し、`<out>_review_audit.tsv` に claim 別のレビュー者・レビュー日・source検証状態を残す。

## 7. SBEIケーススタディの位置づけ

ZW6–Cameor V2 SBEI例は、既知座周辺の擬似QTLと明示的trait keywordを使ったend-to-end sanity checkである。

示せること:

- transcript-aware抽出から候補順位まで実データで動く
- 既知機能候補を早期に回収できる
- assembly間coding差がないnegative resultを扱える
- 未確認assembly SNPをmarker設計へ流さない

まだ示せないこと:

- 未知QTLに対する一般的なTop-1精度
- 多様な遺伝子クラスに対するrecall
- LLM追加による性能向上
- 必要Sanger本数の平均削減率

一般性能はIssue #5の複数既知QTG benchmarkで評価する。
