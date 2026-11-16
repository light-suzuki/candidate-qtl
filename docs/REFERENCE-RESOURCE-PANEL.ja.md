# 比較 reference-resource panel と local-only fallback

`reference-panel` は、比較ゲノム解析で使う **resource契約**(FASTA/GFF3/proteome/DB/API)をrole別に宣言し、その可用性を監査する。これはLLMが提案する「形質関連の既知遺伝子panel」(`reference_gene_panel.tsv`)とは別物であり、両者を混同しない。

## role

| role | 用途 |
| --- | --- |
| `native` | QTL座標・候補抽出の基準genome(必須) |
| `close_anchor` | 高品質annotation・synteny用の近縁種(optional) |
| `functional_reference` | mutant/pathway/literatureが豊富な参照種(optional) |
| `outgroup` | gene-family重複・保存性確認用(optional) |

role・resource_type・taxon・識別子は汎用のまま扱い、特定作物名やgene ID体系をコードにhard-codeしない。

## 実行

~~~powershell
candidate-qtl reference-panel --panel reference_resource_panel.tsv --out-dir panel_out --json
~~~

入力TSVの必須列は `resource_id`・`role`・`resource_type`・`taxon`・`source`・`uri_or_path`。任意列は `assembly_id`・`annotation_version`・`proteome_version`・`db_version`・`external_api`・`review_status`・`notes`。テンプレートは [`templates/reference_resource_panel.tsv`](../templates/reference_resource_panel.tsv)。

## 受入条件との対応

- **native FASTA/GFF3だけで主ワークフローが完走する**: 本コマンドは `native` roleのresourceがローカルに存在することだけを要求し、他roleは無くてもよい。主ワークフロー(`case-run`)はこのpanelに依存しない。
- **追加resourceはoptional、欠落時はskip理由をmachine-readableに残す**: ローカル欠落は `availability=skipped` / `skip_reason=local_resource_missing`、リモート/外部APIは `availability=declared_external` / `skip_reason=external_resource_not_fetched_local_first` として `reference_resource_panel.tsv` とmanifestの `skip_reasons` に記録する。
- **provenanceをmanifestへ保存**: taxon、assembly/annotation/proteome/DB version、source、per-resource fingerprint(SHA256またはsize-tiered)、external API使用と送信先を `reference_resource_panel_manifest.json` に保存する。
- **作物名・gene ID体系をhard-codeしない**: enumはgenericなrole/resource_typeのみ。
- **known-gene panelと混同しない**: 別ファイル・別schema(`REFERENCE_RESOURCE_PANEL_SCHEMA`)・別column設計。

## 出力

- `reference_resource_panel.tsv`: 検証済みpanel + `availability` / `skip_reason` / `resource_sha256` / `resource_size_bytes`。
- `reference_resource_panel_manifest.json`: role/availability集計、native可用性、skip理由、external API使用、fingerprint、fingerprint mode、科学的限界。

resourceが存在することはprovenanceであり、homology・orthology・形質関連の証拠ではない。外部resourceはlocal-first契約として宣言のみで、この経路では取得しない。

## Pisum LIS dataset snapshot

`config/pisum_lis_reference_registry.yml` は、現行LIS Datastoreで確認した `Cameor.gnm1.ann1.7SZR`、`Cameor.gnm2.ann1.KT6X`、`JI2822.gnm1.ann1.D0QS`、`ZW6.gnm1.ann1.TKZX` を、gene namespaceとは別のKEY4付きdataset identityとして保持する。`freeze_reference_resource_snapshots` はrun時に各datasetの公式`MANIFEST.*.yml`と`CHECKSUM.*.md5`を個別にfull SHA-256 freezeする。direct targetの欠損は`blocked_required_resource_missing`、optional datasetの欠損はresource単位の`unavailable`/`not_evaluated`になり、別datasetの利用可否を巻き込まない。

このsnapshot層はassembly間identityを推定せず、候補を削除しない。gene mappingとmulti-reference candidate unionは`assembly-bridge` / `multi-reference-locus-union`の別契約であり、verified 1:1以外のone-to-many、PAV、alternate-only geneを保持する。
