# Route C projection evidence contract

Issue #663 の実装では、Direct C の post-freeze projection を次の2層に分けます。

```text
route_c_direct_projection_evidence.tsv
  exact / orthology / BLAST の全評価・全hitを保存
        ↓ deterministic reducer
route_c_direct_projection.tsv
  focal geneごとの選択basisとQTL内判定
```

## 選択規則

同じ focal gene を複数の方法が支持するときは上書きせず、全行をledgerへ残します。最終basisは次の固定順です。

```text
exact focal stable ID
> explicit orthology-reference mapping
> threshold-gated BLAST
```

同じgeneへの複数方法は `corroborating_methods` です。方法ごとに異なるgeneへ着地した場合は `mapping_conflict` とし、exact行が含まれていてもDirect-C canonical supportへ自動昇格させません。候補と各証拠は削除せずreview ledgerへ残します。

## missingness

BLASTファイルが存在するだけでは、そのrecordをBLAST評価済みとみなしません。

- `not_evaluated_not_supplied`: hits artifactなし
- `not_evaluated_no_panel`: artifactはあるが当該recordのpanel行なし
- `evaluated_no_hit`: 当該recordを評価したがeligible interval hitなし
- `evaluated_with_hit`: eligible hitあり

## 再現性

atlas行、BLAST行、重複BLAST行の入力順に依存しないcanonical sortとreducerを使います。`projection`と`projection_evidence`のTSVは同じ科学入力に対してbyte-stableです。raw inputのSHA-256はmanifestへ別途保持します。

この変更はscore、BLAST threshold、Wet budget、cross-species trait-transfer規則を変更しません。
