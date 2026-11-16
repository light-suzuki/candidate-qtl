# Route C claim-binding producer

Direct C の citation 解決・cue一致・verbatim quote確認だけでは、次の要素が同じ主張に属することは証明できません。

```text
gene / locus / variant
source species / accession
experiment / allele / perturbation
target phenotype
```

`route_c_claim_binding_producer` は、人が一次資料を読んで確認した結合情報を、既存の `route-c-source-claim-binding-1.0.0` sidecarとして検証・封印します。モデルのrelation labelだけでpositive supportを作る機能ではありません。

## 1. review packetを作る

```bash
python -m candidate_qtl.route_c_claim_binding_producer prepare \
  --frozen route_c_direct_rescue_frozen.tsv \
  --freeze-manifest route_c_direct_rescue_freeze_manifest.json \
  --out-dir binding_review
```

生成される `route_c_claim_binding_review.tsv` には、record roster、source-verifyで解決したcitation、accepted extraction、frozen identityが事前入力されます。各verification axis、`curation_method`（人が確認したら `curated_human`、独立agentなら `independent_agent_review`）、`reviewed_by`、判断理由を記入します。

reviewで利用できるscopeの上限は、freezeでhashを固定した `source_review_records.json` の `available_text_scope` と実際の保存sectionから決めます。`fullText` を保存済みなら、分類用の `retrieval_scope` / `classification_scope` が `title_abstract` のままでも `full_text` reviewへ配送できます。分類、引用、identity、curationの判定は変更しません。captureを持たない旧freezeは元のretrieval scopeを使い、captureを宣言した世代のhash欠損・不正なavailability・保存sectionとの不一致は拒否します。取得されなかったcitationは本文scopeを許可しません。

## 2. reviewを封印する

```bash
python -m candidate_qtl.route_c_claim_binding_producer seal \
  --frozen route_c_direct_rescue_frozen.tsv \
  --freeze-manifest route_c_direct_rescue_freeze_manifest.json \
  --review binding_review/route_c_claim_binding_review.tsv \
  --out-dir bound_freeze
```

出力:

- `route_c_source_claim_binding.tsv`
- `route_c_direct_rescue_freeze_manifest.bound.json`

元のfreeze manifestは変更しません。bound manifestは、frozen row、source-validation TSV、verified source JSON、accepted quote、review input、binding sidecarをSHA-256で結びます。下流の `route-c-direct-project` と `route-c-direct-members` にはbound manifestを渡します。

## verified_completeの条件

少なくとも次を全て満たす必要があります。

- citationがsource-verifyで解決済み
- quoteがsource-verifyのaccepted verbatim extractionと一致
- review identityがfrozen recordの同じID/symbol/locusである
- source organism、experimental target、phenotype、perturbation-result bindingが明示的にverified
- named reviewerと`curation_method=curated_human`または`independent_agent_review`（どちらも上の機械検査を同じく通す。方法はbound manifestの`curation_method_counts`に記録）
- variant-only、related-trait-only、title-onlyは現行v1契約では昇格しない

不足は`review_required`または`not_evaluated`、矛盾は`contradicted`になります。provider failureや未評価をbiological negativeへ変換しません。

このproducerは#661のRoute C固有の安全経路を実装し、Issueはclose済みです。A/B/C共通の自動semantic relation proofを実装したという意味ではありません。`curated_human`は実際に人が原典を確認した場合だけ記録し、agentの判断は`independent_agent_review`として区別します。A/C優先（2026-11-08）では、機械検査を通ったagent reviewのbindingもsupportを認めます。projectionの残受入は#663で管理します。
