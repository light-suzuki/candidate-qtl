# Direct C の論文報告symbolと参照配列の対応

変更分類: `CORE_INTEGRATION`。既存の `--identity-mapping` とreference bridgeを接続する任意のreview経路です。Tier、候補順位、verified supportの条件は変更しません。

論文で報告されたsymbolと登録CDS名が異なる場合、名前の類似や共出現だけでは対応を認可しません。原論文の同じgene→accession関係を読んだ担当者が、凍結原文中の引用、場所、citation、種、正確なaccessionと判断理由を明示します。機械は引用の所在・endpoint・hash・retrieval scopeを確認し、関係の意味は明記されたreviewに依拠します。これはdatabaseによる同一遺伝子証明ではありません。

## reviewの準備と実行

既存のclaim bindingを作った凍結世代に対して、次のPython APIで未承認templateを作れます。source/modelへの問い合わせは行いません。

```python
from candidate_qtl.route_c_paper_identity_mapping import prepare_paper_mapping_review

result = prepare_paper_mapping_review(
    "c-freeze/route_c_direct_rescue_frozen.tsv",
    "c-bound/route_c_direct_rescue_freeze_manifest.bound.json",
    "paper-identity-review",
)
print(result["review"])
```

`paper_identity_mapping_review.tsv` の空欄を、実際のsource reviewに基づいて埋めます。空欄のままのtemplateはconsumerに拒否されます。独立reviewを実施していない者の名や実施していないcuration methodを書いてはいけません。

| 追加列 | 内容 |
|---|---|
| evidence_type | `reviewed_paper_symbol_sequence_v1` |
| citation | この凍結recordとclaim bindingに対応する同じ論文 |
| source_scope | 実際に凍結取得済みの `abstract` または `full_text` |
| locator | 論文内の該当箇所 |
| curation_method | 実際に行った `curated_human` または `independent_agent_review` |
| mapping_relation | `reported_gene_sequence` |
| mapping_reason | 共出現ではなく、そのgeneの配列であると読める理由・留保 |
| taxonomy_record_path / taxonomy_record_sha256 | 種名表記が異なる場合だけ必要な、公式taxonomy recordのpath/SHA（両方空欄は可） |

既存mapping列も使います。sourceは `reported_gene_symbol`、targetは `reference_sequence_id`、source_namespaceは `paper_reported_gene`、mapping_basisは `curated_id_map`、rank_contributionは `none` です。`quoted_mapping` は凍結原文内の実際の引用で、symbolと正確なaccessionを含めます。引用内に2つの文字列があるだけでは関係を証明しません。複数段落を使う場合も、保存原文の連続した該当範囲として引用します。

`source_record_path` とSHAはtemplateが指す **凍結source_review_records** を保持します。別途保存したHTMLや任意の文字列ファイルへ差し替えてはいけません。必要な配列登録節が元の凍結取得範囲にない場合は、原資料の取得・抽出・reviewを新世代に記録し、その新世代を凍結してから使います。元の凍結sourceや研究回答は改変しません。

```text
route-c-direct-project --frozen <frozen.tsv> --freeze-manifest <bound.json>
  --atlas <atlas.tsv> --qtl-id <id> --out-dir <new-generation>
  --identity-mapping <reviewed-paper-edge.tsv>
  --reference-bridge <reference_sequence_bridge_manifest.json>
```

新しいmappingのmanifest contractは `route-c-paper-source-identity-mapping-1.0.0`。既存のstable ID/locus→構造化UniProt xrefのcontractとschemaはそのまま再生します。

## 保持する境界

- 原recordのreference_sequence_idが既知なら完全一致が必須です。空欄の場合だけ、同じ凍結論文のreview済みedgeが認可したtargetを投影用viewに補います。原recordの空欄とhashは保持します。
- 同種の公式translation、taxid、CDS、sequence SHA、bridge request、citation、全祖先hashを確認します。raw eligible-hitsだけではこのpaper edgeを認可しません。取得不能は未評価です。
- 原reference/focal種名と公式配列recordのorganismが一致しない場合は、bridgeの同じtaxidを持つ公式taxonomy recordが必要です。既存offline loaderの `taxonomy-<taxid>.json`（taxonId、scientificName、synonyms）をpath/SHAで封印し、3者の名前がそのscientificName/synonyms内にあることを確認します。人が作ったalias表・reviewer名・共通属だけでは認可しません。原名は書き換えず、別taxidや無関係な名前、taxonomy改ざんは拒否します。公式raw responseと正規化の出典は調査receiptに保持します。
- symbolしかなくidentity未解決だったrecordは、sourceがverified、traitがdirect target/related、directnessがD0-D2の場合に限って投影用viewのidentityを解決できます。source未確認、context-only、異種host、variant未確認などの別の不足は解消しません。原pending状態は変更しません。
- projection ledgerは正確なsequence endpointと `reviewed_paper_sequence_mapping` を記録します。最終consumerもpaper・mapping・bridge・原recordのhashを再検証します。
- このsymbol経路のWet採用には、reviewed paper edgeと同じsequence endpointの実BLAST投影が必要です。名前やaccession文字列と区間IDの偶然の一致だけでは採用しません。BLASTはreviewの投影根拠であり、verified supportへ昇格しません。D2 associationに新しいcausal perturbationの肯定を要求しません。複数copyは保持し、競合するmappingや不正なendpointは拒否します。
- 実験hostと配列の由来が異なる関係はこの経路では同一視しません。名前aliasの追加、symbolをlocusへ書き換える回避、真値に合わせた例外は行いません。

検証: `tests/test_route_c_paper_identity_mapping.py` が中立fixtureの実BLAST・CLI・members・Wet、原空欄保持、複数copy、改ざんと不正edgeを確認します。これは科学benchmark性能の主張ではありません。
