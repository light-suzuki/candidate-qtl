# Interval atlas（phenotype非依存の区間俯瞰、issue #168）

## 何のためのものか

このworkflowは、正しいpathway / gene family / homolog bait が与えられた後の区間内候補回収には強い。最大の不確実性はその手前にある。

```text
phenotype -> どのpathwayか -> どのgene familyか -> どのreference homologをbaitにするか
```

ここは研究者またはLLMが**あらかじめ正解を思いつけるか**に依存している。`interval-atlas` は順序を逆にしてこの依存を減らす。

```text
QTL interval
  -> 区間内の全遺伝子
  -> family / domain / orthogroup / broad function
  -> conserved / relative-enriched / no-detectable-homology / unresolved か
  -> reference種との対応
  -> 人が読める interval atlas
```

正解のbaitを先に当てる必要を減らし、**区間内に実在する有限の選択肢を見てから仮説を作る**routeをfirst-classにする。

## phenotypeを使わない

このstageには**phenotypeが一切入力されない**。設計上そうしてある。

- 仮説を持つ前の区間の見え方を固定し、後付け仮説で見え方が変わっていないことを監査できるようにするため
- manifestに `phenotype_used: false` と `atlas_sha256` を残すので、後段のoverlayがatlasを黙って書き換えていないことを検証できる

phenotypeのoverlayはatlas生成後の別stage（[#169のRoute B](PHENOTYPE-PATHWAY-PROPOSAL.ja.md)）が担当する。

## 分類レイヤー

| 列 | 内容 | 由来 |
| --- | --- | --- |
| `family_id` / `family_label` | 家系（best-hit reference geneのstem） | `family-index` を再利用（再計算しない） |
| `discovery_bucket_id` / `discovery_bucket_label` / `discovery_bucket_method` | **探索用naming bucket**（`symbol_stem` / `curated_keyword` 等）。仮説生成のrecall-first index | `family-index` |
| `evolutionary_family_id` / `evolutionary_family_source` / `evolutionary_family_status` | **homologous family**（`supported` / `provisional` / `unresolved`）。sourceは `curated_family` / orthogroup / domain+sequence等 | `family-index`のcurated family map等 |
| `subfamily_id` / `subfamily_source` | homologous family内の機能的/進化的枝。atlasはcurated subfamily mapを持たないため、**orthogroupを最も正直な代理**とし、無ければ `unresolved` | orthology |
| `domain_ids` / `domain_architecture` | domain集合と**順序** | `--domains`（Pfam/InterPro相当） |
| `orthogroup` / `orthology_relationship` / `orthology_providers` | orthogroup、1:1 / 1:many / many:many、provider | `--orthology`（orthofinder/omamer/synteny import等） |
| `copy_set_id` / `copy_set_status` | **copy集合の予告**。homology layer id（curated evolutionary family / orthogroup / domain architecture）があれば `validated`、discovery bucketのみなら `provisional`、family無しなら `not_applicable`、dark/未割当なら `unresolved`。**最終決定は family-to-copy** | 上記の合成 |
| `broad_function_bins` / `primary_broad_function` | 粗い生物プロセスbin | 統制resource（下記） |
| `conservation_class` | `conserved_cross_species` / `relative_enriched` / `no_detectable_homology` / `unresolved`（`lineage_restricted_supported` は reserved） | orthology + family + domain |
| `focal_paralog_count` / `tandem_expansion` | 区間内の同family数、tandem集積 | 座標 + family |
| `evidence_completeness` | 分類根拠の充足度 | 下記 |

domain **順序** を保持しているのは、後段（#170 family-to-copy）で似たfamily memberを識別するのに効くからです。

## 「family」という語のレイヤー

workflow中の `family` は1語で2つの別物を指していました。これを**別列・別概念として分離**します。

```text
1. discovery bucket（探索bucket）
   命名 / 広い機能検索bucket。例: CLV, LRR-RLK, WOX-like
   symbol stemやcurated keywordから作る。recall-firstの仮説生成には使ってよいが、
   単独では「同じhomologous family」を保証しない

2. homologous family（進化的family）
   sequence / domain / gene-tree / orthogroup等でhomologyが支持されたfamily

3. orthogroup / HOG（provider定義の進化的group）

4. subfamily（homologous family内の機能的/進化的枝）

5. paralog / copy set（同一focal genome内で、同じhomologous family/subfamily由来として
   比較可能なcopy集合）
```

- 同じsymbol stem（CLV1/CLV2/CLV3 → CLV）は**同じbucket**には入りますが、**同じhomologous family・copy setを意味しません**
- 既存の `family_id` / `family_label` は**移行期間のalias**として残ります（値はbucketとcurated familyのどちらになり得るか入力依存）。意味を読み取るには `discovery_bucket_*` / `evolutionary_family_*` / `orthogroup` / `copy_set_*` を使ってください
- manifestの `family_semantics` に各レイヤーの定義を記録します
- dark / unknown / 未解決familyは**除外されません**。`no evolutionary family assignment` はnegative evidenceにしません（short peptide / small ORF / fast-evolving geneは `family_resolution_required` / `unresolved` として残る）

### broad function bin は統制resource

binのルールはコード中のkeyword listではなく、版管理された `templates/broad_function_bins.tsv` にあります（wheelにも `candidate_qtl/resources/` として同梱、drift testで同期を保証）。

- `match_type=go`: GO IDの完全一致
- `match_type=keyword`: description中の部分一致
- GO由来の割り当てはkeyword由来より優先して `primary_broad_function` になる
- 各遺伝子の `broad_function_evidence` に **どのパターンが当たったか** が残る（例 `go_term:GO:0010073`）

`--function-bins` で差し替え可能です。

### evidence_completeness

「分類不能」を消さず独立bucketとして残します。

- `high_confidence_family_assignment` — family + orthogroup
- `family_assignment` — familyのみ
- `domain_only_assignment` — domainのみ
- `orthogroup_only_assignment` — orthologyのみ
- `conflicting_annotation` — family indexとorthology referenceが**別の家系名を指す**（自動解決せず記録）
- `dark_or_short_orf`
- `no_informative_annotation`

後ろ2つと `conflicting_annotation` は `unclassified_dark_candidates.tsv` にも出て、`annotation_rescue_needed` と `required_next_action` が付きます。

### conservation_class と lineage claim

**検索で見つからなかったことと、進化的に存在しないことは分離する。** family assignmentもreference orthologyも得られないことは、focal lineageに特異的なgeneであることを意味しません（reference proteome不足、remote homology、short peptide、partial gene model、fast-evolving、domain-only conservation、provider/threshold failure等が原因になり得るため）。

| 値 | 意味 |
| --- | --- |
| `conserved_cross_species` | one-to-one orthology row（orthology入力由来） |
| `relative_enriched` | orthologyのreference-gene match（1:many / many:many等の明示的なreference遺伝子対応） |
| `unresolved` | family / orthogroup / domain の部分的なsignalはあるが、reference-gene matchがなくenrichmentを断定できない |
| `no_detectable_homology` | family / orthology / domain のいずれもsignalなし（**safe default**。scoreを下げず、gene lossやnovel geneの証拠にしない） |
| `lineage_restricted_supported` | **reserved**。複数reference/outgroupの十分なcoverageを持つnegative search等の明示的positive criteriaを満たした場合のみ付与。現在のatlas入力ではこれを満たせないため、atlasはこの値を出力しない |

no-hit / unresolved のgeneは保持され、`rank_contribution` は `none` のままです（atlasはscoreを一切書かない）。これは[#184](CLI.ja.md)の `family_supported_copy_unresolved` / `subfamily_unresolved` と同じ「不明は不明のまま残す」方針と整合します。

**lineage claimのprovenance / coverage** は `interval_atlas_manifest.json` の `lineage_claim` に記録されます:

- `claim_basis` — 各値の定義（判定根拠）
- `vocabulary` — 使用可能な語彙の一覧
- `reference_species_searched` / `reference_roles` / `search_method` / `thresholds` — family-index manifestから引き継いだ検索対象・閾値（manifestが無ければ `not_supplied`）
- `provider_coverage` — orthology provider一覧、family-indexで区間内遺伝子に割当済みのreference subject数
- `no_hit_count` — 入力ごとのno-hit数と `no_detectable_homology` 遺伝子数
- `outgroup_coverage` — outgroup negative searchは現在 `not_supplied`（供給されない限り `lineage_restricted_supported` は付与されない）

## Family-first research route

`family_review_queries.tsv` は、**遺伝子IDではなくfamily / process単位**の調査クエリを出します。phenotypeを使っていないので、クエリも phenotype-free です。

```text
QTL内に CLV family 3 genes
  -> このfamilyのdevelopmental roleとmutant phenotypeを調べる

QTL内に transcription / chromatin bin 4 genes
  -> このクラスで報告されているprocessとregulatorを調べる

分類不能 2 genes
  -> conserved unknown proteinのortholog phenotype / expression / knockout文献を調べる
```

LLMを使う場合も、いきなり全候補から原因遺伝子を選ばせず、`interval atlas + phenotype -> 調べるべきfamily/processを複数提案 -> source-backed retrieval -> verified bait panel -> 既存のdeterministic mapping/ranking` と段を踏みます。

## 安全境界

- family分類そのものをcausal evidenceとして自動加点しない（全出力の `rank_contribution` は `none`）
- 「有名なfamily」であることだけでTop昇格させない
- unknown / orphan / dark candidateを自動的に下位除外しない
- no-hitを `lineage_restricted_supported` にしない（no-hitは探索上の未知であり、進化的結論ではない）
- **discovery bucketとhomologous familyを混同しない**（naming bucket単独でcopy集合を定義せず、family-to-copyのgateはhomology層を要求する。issue #184）
- broad-function binの件数をgene scoreへ直接加えない
- family assignment confidence と causal confidence を混同しない
- 複数reference種が**同一annotation chain**経由で一致しても独立証拠として水増ししない（providerは記録するがvoteさせない）
- phenotype overlay前のatlasをfingerprint付きで保存する

## 既存機能との関係

- `family-index`: リファレンスごとのfamily割当をそのまま再利用（重複計算しない）
- `pathway-interval-falsify`: atlasで見つけたfamily/process仮説を後段で反証
- `llm-pathway-*`: phenotype-only laneに加え、interval-atlas-assisted laneを別モードとして持てる
- `annotation-rescue` / `dark-candidates`: unknown bucketの `required_next_action` として接続
- `orthofinder-import` / `omamer-import` / `synteny-import`: `--orthology` の入力元

## 使い方

```bash
candidate-qtl interval-atlas \
  --candidates candidates.tsv \
  --family-index out/family_index/family_index.tsv \
  --domains domains.tsv \
  --orthology orthology.tsv \
  --annotations reference_annotations.tsv \
  --flags candidate_flags.tsv \
  --out-dir out/interval_atlas
```

optionalな表は全て省略できます。何も無ければ全遺伝子が `no_informative_annotation` / `dark_or_short_orf` bucketに入り、**1件も落ちません**。
