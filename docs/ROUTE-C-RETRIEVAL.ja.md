# Legacy Route C retrieval — phenotype/protein similarity lane

> Historical origin: issue #193
> Current status: **DIAGNOSTIC / HISTORICAL / EXPERIMENTAL**
> Current Direct Route C: [`ROUTE-C-DIRECT-RESCUE.ja.md`](ROUTE-C-DIRECT-RESCUE.ja.md)

この文書は、lexical / LSA / neural phenotype retrieval系を説明します。

このlaneは現在も再現、比較、語彙miss診断に利用できますが、**current Direct Route Cのprimary meaningではありません**。

```text
Legacy retrieval C
  = phenotype textとreference recordのsimilarity retrieval

Direct C
  = source-backed direct trait–gene/locus/variant factのrescue
```

legacy retrieval score、rank、interval matchをDirect-C direct evidenceやcanonical supportへ読み替えません。

---

## 1. 役割

A/Bとは異なるfailure modeを観察するためのphenotype-only retrieval laneです。

| lane | 起点 | 処理 | current role |
|---|---|---|---|
| Route A | phenotype | mechanism/family/referenceを推論しmachine verification | CORE |
| Route B | QTL interval | all genesをfamily/orthogroupへ圧縮しsemantic interpretation | CORE |
| Legacy retrieval C | phenotype text | reference recordとのlexical/latent/neural similarity | DIAGNOSTIC / HISTORICAL |
| Direct C | phenotype + species/cultivar | directly reported factをsource verification/freeze/projection | OPTIONAL |

legacy retrievalは「類似記述を探す」ものであり、「sourceで直接検証されたtrait–gene relation」を作るものではありません。

---

## 2. Interval-blind contract

retrievalはintervalを見る前にfreezeします。

```text
source-backed reference corpus
  -> retrieval-index-build
  -> retrieval-run --phenotype ...
  -> retrieval-import / freeze
  |  interval-blind up to here
  v
retrieval-interval-match
  -> historical/diagnostic comparison
```

構造的な境界:

1. `retrieval-run` / `retrieval-import`にQTL intervalやcandidate IDを渡さない。
2. input contract外の列を拒否する。
3. `retrieval-import`でfreezeした後だけ`retrieval-interval-match`する。
4. interval matchはretrieval rowやscoreを書き換えない。
5. convergenceは`rank_contribution=none`。

---

## 3. Retrieval contract

`templates/route_c_retrieval.tsv`の主な列:

```text
reference_gene_id
reference_species
family_id
orthogroup_id
retrieval_score
matched_modality
matched_description_or_phenotype
model_name
model_version
database_version
rank
```

- family_id / orthogroup_idは少なくともどちらか必要。
- retrieval_scoreはsimilarity measureでありcausal probabilityではない。
- model/database provenanceを空にしない。
- contract外のinterval/candidate/truth列を拒否する。

---

## 4. Provider modes

### `lexical`

決定論的なlocal TF-IDF系。第三者providerに依存せず、historical baselineやreplayに向きます。

### `lsa`

latent-semantic retrieval。neural embeddingではありません。

### `neural`

optional embedding adapter。model/version/database provenanceを必ず保存します。

どのproviderでもsimilarityはmechanism evidenceやdirect genetic evidenceではありません。

---

## 5. 主経路

### 5.1 Reference corpus

```bash
candidate-qtl reference-corpus-build \
  --from-uniprot refcorpus/uniprot_records.tsv \
  --corpus-version legume_v1 \
  --out-dir corpus
```

corpusはphenotype-blind、versioned、fingerprintedにします。

### 5.2 Build index

```bash
candidate-qtl retrieval-index-build \
  --corpus corpus/legume_v1/reference_corpus.tsv \
  --database-version legume_v1 \
  --out-dir index
```

### 5.3 Phenotype-only retrieval

```bash
candidate-qtl retrieval-run \
  --index index/legume_v1/index_manifest.json \
  --phenotype "non-nodulating" \
  --top-k 50 \
  --out-dir out/legacy_retrieval \
  --case-id pea_nin
```

### 5.4 Validate and freeze

```bash
candidate-qtl retrieval-import \
  --retrieval out/legacy_retrieval/route_c_local_retrieval.tsv \
  --out-dir out/legacy_retrieval \
  --phenotype "non-nodulating" \
  --case-id pea_nin
```

### 5.5 Post-freeze interval match

```bash
candidate-qtl retrieval-interval-match \
  --retrieval out/legacy_retrieval/route_c_retrieval_frozen.tsv \
  --atlas out/interval_atlas/interval_gene_atlas.tsv \
  --out-dir out/legacy_retrieval_match \
  --case-id pea_nin
```

---

## 6. External retrieval compatibility

外部embedding/search toolの出力も11列contractへnormalizeして`retrieval-import`できます。

複数DBをmergeする例:

```bash
candidate-qtl retrieval-merge \
  --retrieval arabidopsis_hits.tsv \
  --retrieval legume_hits.tsv \
  --top-k 50 \
  --out-dir out/legacy_retrieval_merge \
  --case-id pea_fas
```

外部toolを使ってもinterval列、candidate ID、truthをpre-freeze inputへ入れません。

---

## 7. Interval matchの意味

matching basisの例:

- same family ID
- same orthogroup ID
- exact `orthology_reference_gene` match

symbol alias、version、species identityを推測して補完しません。

0 matchの場合:

- namespace mismatch
- family/orthogroup coverage不足
- reference corpus不足
- phenotype vocabulary mismatch
- provider/model coverage不足

などを診断します。0 matchはbiology不存在やcandidate negativeではありません。

---

## 8. `triple-route`との関係

`guided-profile --profile triple-route`にはlegacy retrieval stagesが残っています。

```text
c01 retrieval-run
c02 retrieval-import
c03 retrieval-interval-match
```

これはhistorical/diagnostic comparisonです。

一方、current Direct Cは別contractです。

```text
pre-frozen Direct-C record
  -> post-B0 route-c-direct-project
  -> route-c-direct-members
  -> optional supported-member-union / Wet
```

legacy retrieval artifactをDirect-C proposal/frozen recordへ偽装しません。similarity convergenceをverified supportにしません。

---

## 9. Outputs

- `route_c_local_retrieval.tsv`
- `route_c_local_retrieval_status.json`
- `route_c_retrieval_frozen.tsv`
- `route_c_retrieval_status.json`
- `route_c_interval_match.tsv`
- `route_c_interval_match_status.json`
- optional historical comparison artifacts

これらはcurrent Direct-C `route_c_direct_rescue_frozen.tsv`、`route_c_direct_projection.tsv`、`route_c_direct_members.tsv`とは別schemaです。

---

## 10. Safety boundary

- retrieval scoreはcausal probabilityではない。
- high similarityはdirect trait–gene evidenceではない。
- interval matchはsupport promotionではない。
- route convergenceはscore bonusではない。
- no hit / provider failure / not evaluatedはnegative evidenceではない。
- legacy retrievalからsingle focal copyを強制しない。
- benchmark truthを見てsynonymやpromptを同generationで調整しない。
- historical benchmark metricをcurrent Direct-C performanceとして報告しない。

---

## 11. 使うべき場面

- vocabulary mismatchの診断
- deterministic lexical baseline
- embedding provider比較
- historical benchmark再現
- Direct CやRoute Aが見落としたreference descriptionのreview候補探索

Direct factを探してcanonical rescue evidenceへ進めたい場合は、[`ROUTE-C-DIRECT-RESCUE.ja.md`](ROUTE-C-DIRECT-RESCUE.ja.md)を使います。

---

## 12. 一文での定義

> **Legacy Route C retrievalは、phenotype textとfrozen reference corpusの類似度をinterval-blindに測り、freeze後にQTL intervalと照合するdiagnostic/historical laneである。Direct trait–gene factを検証するcurrent Direct Cではなく、similarity scoreやconvergenceをcanonical supportへ昇格しない。**
