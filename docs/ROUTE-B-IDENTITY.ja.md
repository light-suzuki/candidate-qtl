# Route B: phenotype-free identity material

`route-b-prepare` の既存local BLAST経路で生成する調査材料です。
Wet採否の正本は [canonical membership](CANONICAL-WET-MEMBERSHIP.ja.md) を参照してください。
この文書は採否Tierを新設せず、identity artifactの意味と限界だけを定義します。

## 出力

- `route_b_identity_hits.tsv`: 閾値を通った取得済みhitをすべて保持するgene–reference表。
  accession、species、taxon、gene symbol、protein name、description、identity、coverage、
  e-value、bitscore、provider/database、raw header recordsを含みます。
- `route_b_identity_inventory.tsv`: 区間全geneの状態表。hitがないgene、protein queryのないgene、
  providerが利用不能だったgeneも消しません。
- `route_b_preparation_manifest.json`: `outputs.identity_hits` / `outputs.identity_inventory` と
  SHA-256、coverageの分母・分子、provider状態を記録します。

既存の `route_b_orthology.tsv` にもmetadata列を保持しますが、relationshipは
`best_hit_homology` のままです。厳密なorthology、trait support、単一copy確定へは昇格しません。
既存interval atlasは書き換えず、rootのcanonical atlasと生成元atlasのbyte一致を維持します。

## 名前・種の扱い

BLASTの`staxids` / `sscinames`、明示的なUniProt `OS/OX/GN`、NCBI末尾のorganism表記を読みます。
IDのsuffixや名称からspeciesを推測しません。未解決の値は`unknown`です。
同一query–referenceに複数の異なるheaderがあった場合、全recordをJSONへ保持し、
矛盾するscalar値を勝手に一つ選ばず`conflicting`とします。
同じ種のhitだけでも受理し、種数のquotaは置きません。

## coverageの意味

`identity_coverage`はBLASTまたは既存annotationのidentity/function handleがある割合です。
BLASTだけの割合は`blast_identity_coverage`で別記します。
`functional_annotation_coverage`、`source_backed_literature_coverage`、
`dark_unresolved_coverage`は別々の分子と分母を持ちます。

**literature coverageはreference source rosterの存在率です。引用の実在確認、gene–trait binding、
focal geneへの機能移転を検証した割合ではありません。** すべてrank contributionはnoneです。

provider unavailable、provider failure、malformed output、正常終了したno-threshold-hit、
protein query未供給を区別します。いずれもnegative biological evidenceには変換しません。

## 保持範囲と検証限界

既存のresidue-aligned identity/coverage gateとBLAST retrieval limitは変更していません。
保持対象は「取得済みで閾値を通った全hit」であり、DB全体の全相同配列を取得する保証ではありません。
BLAST自身の`max_target_seqs`境界ではDB順序でtieが解決されるため、DBを並べ替えても
全検索結果が一致するとは主張しません。取得後のHSP順序・重複・同点では保持結果を変えません。

ローカルのsynthetic unit/consumer testsは、metadata・状態分離・全gene inventoryの検証です。
実DBでの検索網羅性、installed CLI全体、full regression suite、候補保持性能は別の検証です。

## 実環境での配線診断（2026-09-21）

`tests.test_type1_native_wiring.NativeIdentityTests`はinstalled module CLIと実`makeblastdb`/`blastp`を使い、
1,000-gene inventory、mixed-species header、taxonomy DBなし、空reference corpus、protein欠損geneを確認します。
同一配列7hitに対して取得上限5ではDB並べ替えで取得集合が変わることを実測します。
この上限は現行検索の限界であり、全DBの同点hit保持を保証しません。
未知のpipe形式IDを共通接頭辞へ縮約する不具合と、空corpusでcoverage計算が例外になる不具合を修正しました。
合成配列による配線検査であり、生物学的な候補保持性能の測定ではありません。

仕様参照: [NCBI BLAST+ options](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
（custom outfmt、`staxids`、`sscinames`、`max_target_seqs`）。
