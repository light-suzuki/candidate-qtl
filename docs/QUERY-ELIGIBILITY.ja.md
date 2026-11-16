# Query eligibility gate（Pisum由来query禁止）

Cameor V2 の既知遺伝子ベンチが主張するのは、**「他のPisum reference配列を使わず、非Pisum植物の既知機能遺伝子と V2 native annotation だけを用いた、literature-informed な map-based candidate prioritization」** に限定した性能。V2 の gene ID が旧referenceと違っても、Cameor v1a・ZW6・JI系統等のPisum reference にある同一原因遺伝子配列を V2 区間へ BLAST すれば、実質的に正解を直接検索したことになる。

`query-eligibility` は homolog query manifest を **ベンチ実行前に**検証し、**非Pisumの機能ホモログだけ**を通す gate。`discovery-benchmark` の real-answer-source leave-out は別の追加監査で、この gate の代替にはしない（責務分担）。

## 禁止 / 許可

**禁止（除外）**:
- `Pisum sativum` / Pisum 属由来配列（genus単位で安全側に禁止）
- Cameor v1a・ZW6・JI系統等、V2以外のPisum assembly由来
- Pisum配列から作った consensus / projection / renamed record / translated CDS
- 正解ID（V2/旧reference）から逆算した配列（answer-backmapped）

**許可（eligible）**:
- Medicago / Lotus / Glycine / Arabidopsis 等の非Pisum植物の機能ホモログ
- 表現型・機構・pathway・domain 仮説から選んだ既知機能遺伝子
- V2固有の FASTA / GFF3 / proteome / eggNOG・GO・domain annotation

非Pisum query は「正解copyを保証する因果証拠」ではなく、通常のマップベースクローニングにおける候補圧縮用の comparative-function evidence として扱う。

## 入力

`--query-manifest` TSV。必須列: `case_id`・`homolog_query_id`・`query_species`・`sequence_origin`・`selection_basis`。任意列: `query_genus`（無ければspeciesから導出）・`query_taxon_id`・`source_assembly`・`homolog_source_id`・`homolog_source_type`。テンプレートは [`templates/homolog_query_manifest.tsv`](../templates/homolog_query_manifest.tsv)。

`--query-fasta`（任意）: query配列FASTA。header の `OS=`/`OX=` または `[Species]` から taxon を監査し、SHA256 を記録。**self-申告TSVだけに依存せず**、renamed record が Pisum provenance を隠せないようにする（FASTA が Pisum を示せば申告より優先して除外）。

`--known-answer-sources`（任意）: 正解論文ID（PMID/DOI）、正解アクセッション、正解エイリアスが記録されたTSVやファイルパス。これらに一致する情報源から選抜されたhomolog queryを`excluded_answer_derived`として遮断する。

## 判定（順序）

```
Pisum genus / Pisum assembly / Pisum-derived origin / FASTA taxon=Pisum -> excluded_pisum
selection_basis in {answer_backmapped, pisum_sequence_based}          -> excluded_answer_derived
known_answer_sources との一致 (PMID, DOI, alias, accession)            -> excluded_answer_derived
taxon provenance を確認できない / sequence_origin=unknown             -> excluded_unknown_provenance
申告species と FASTA監査species の不一致                                 -> excluded_unknown_provenance
strictモードでquery_fasta指定時にSHA256が取得できない                    -> excluded_unknown_provenance
非Pisum + provenance確認可                                             -> eligible
```

## 実行

~~~powershell
# audit mode（除外して続行、理由を記録）
candidate-qtl-bench query-eligibility --query-manifest queries.tsv --query-fasta queries.faa --known-answer-sources truth.tsv --out-dir gate --json

# strict mode（禁止queryが1件でもあれば実行前にfail）
candidate-qtl-bench query-eligibility --query-manifest queries.tsv --query-fasta queries.faa --known-answer-sources truth.tsv --out-dir gate --strict

# reference-benchmark との連携（事前監査manifestを渡し、fail-closedでquery SHA256をバインド）
candidate-qtl benchmark-reference --config benchmark.tsv --fasta Cameor_v2.fa --gff3 Cameor_v2.gff3 --protein-fasta Cameor_v2_prot.fa --eggnog Cameor_v2_eggnog.tsv --out-dir out --homolog-proteins queries.faa --query-eligibility-manifest gate/query_eligibility_manifest.json
~~~

## 出力

- `eligible_queries.tsv` / `excluded_queries.tsv`
- `query_provenance_audit.tsv`（全query: 申告 vs 監査 species/genus、SHA256、provenance_method、eligibility、理由）
- `eligibility_summary.tsv`（eligibility×理由の件数）
- `query_eligibility_manifest.json`（`strict` / `pisum_queries_present` / `pisum_queries_executed=false` / `answer_derived_queries` / `fingerprints_verified` / `queries` 台帳 / unknown数 / eligible・excluded数 / 理由 / taxon検証方法 / `benchmark_claim`）

## 責務分担とガードレール

- gate は **実行前**の query フィルタ。`discovery-benchmark` の source-based leave-out は **評価後**の追加監査。両方を保持し、gate の代替にはしない。
- gene ID が違うことを非リークの根拠にしない。Pisum別reference・renamed record・projection・answer-backmapped を一律禁止。provenance不明は主ベンチから除外。
- `benchmark-reference` 実行時に `--query-eligibility-manifest` を指定した場合、manifestが `strict: true` かつ `excluded_count == 0`, `pisum_queries_present == 0`, `answer_derived_queries == 0` であることを検証し、実行される query 配列の SHA-256 が監査済み SHA-256 と完全一致することを強制する（不一致または未監査queryは即時fail-closed停止）。
- これは**他作物一般化のde novo discoveryではなく**、Cameor V2 上の literature-informed map-based benchmark である。
