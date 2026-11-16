# Route B の複数参照保持

ローカルBLASTが取得した参照のうち、identity・coverage・evalue条件を満たすものを
全件、raw subject IDとともに保持する。最高bit scoreの1件への縮約はしない。
参照集合と正規化済み出力順は、取得済み行の順序・同点・完全重複に依存しない。

**取得範囲と保持範囲は別である。** `max_target_seqs=5`は今回変更していない。
DB検索が返さなかった第6件以降やDB内同点の取得順まで解決したという主張ではない。
`retained_reference_pair_count`と`genes_with_a_usable_hit`を別々に記録する。
互換用relationship名`best_hit_homology`は局所homologyを表し、唯一のbest/orthologを意味しない。

## canonical consumer

`interval_atlas`は同一geneへの複数referenceを集合として保持する。
従来のB0単一reference欄は互換のため残すが、semantic bundleはcanonical atlasの隣の
`route_b_preparation_manifest.json`が明示したreference tableも読む。

- atlasとcorpusのSHA-256を照合し、別runや改変された資料を拒否する。
- atlasの全reference endpointを、corpusの`reference_gene_id`へexact joinする。
- `reference_material_cards`としてgene・reference・species・機能・sourceを一組ずつ保持する。
- 出典IDは既存のclosed-world source rosterへ追加する。別参照の出典を先頭参照の根拠へ付け替えない。
- 原文の出典/identity/transfer検証、supported union、Wetの規則は変更しない。
- bare atlasなどpreparation sidecarがない入力では追加資料を自動探索しない。

family-indexを使った場合も、`best_subject`と`family_subjects`の明示endpointを既存corpusへexact joinする。UniProtの`sp|accession|entry` / `tr|accession|entry`はaccessionを照合し、raw値とendpointの取得経路をcardに残す。gene symbolから参照を推測しない。assignment状態・identity・query coverage・evalue・bit score・次点との競合情報もcopy別にatlasからsemantic promptへ渡す。これらはsole orthologの確定や検証済み機能支持ではない。

semantic importは凍結済みunitのカードとfamily endpointに明記されたexact referenceだけを提案用の閉世界台帳へ通す。資料の到達だけでsupported unionへ昇格せず、source/transfer検証が未完の候補は`unverified_route_b_review`として確認対象に保持する。件数を満たすための追加BLAST・文献探索は行わない。

これで「機能不明の最高hitだけが残り、機能解析済みの次点がLLMの前で消える」損失を防ぐ。
調査材料の到達と、そのclaimがsupportを得ることは別である。複数候補をsole orthologにしない。
古い単一参照lookupを直接消費する診断経路の網羅性までは、この変更で保証しない。

## 検証

`python -m unittest tests.test_route_b_reference_materials -v`

1,000遺伝子の合成atlasから、実semantic builderのpromptまで次点の機能・出典が届き、
候補遺伝子の個数が増減しないことを検証する。`tests.test_issue823_semantic_endpoint_wiring`では架空のfamily-only資料をcanonical prepare、bundle、import、source/transfer検証、review候補freezeへ通し、copy保持とsupport未昇格も確認する。実モデル・実Wetの性能試験ではない。
