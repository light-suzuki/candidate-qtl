# Evidence weight registryの設計根拠

`candidate-qtl` のweightは、親系統の配列決定や追加実験を行う候補順を決めるための順序尺度です。因果確率、posterior probability、p値、effect sizeではありません。既知遺伝子を後から見てweightを調整する用途にも使いません。

既定値は `templates/evidence_weight_registry.tsv` に集約されています。実行時にはregistry versionとSHA-256を保存し、`score_contributions.tsv` の各行へdirectness、trait specificity、dependency rule、rationale、referenceを付加します。

## 階層

1. 現在の材料でのmutant、complementation等の直接的遺伝証拠
2. 既知QTG、独立allele、形質検証済みortholog
3. recombinant core、明示されたQTL peak/support interval
4. trait-specificな比較ゲノム、発現、domain evidence
5. generic homology、generic domain、GO presence、interval membershipはcontextのみ

ECOはevidence typeとassertion provenanceを区別して記録する枠組みを提供します。GO evidence codeもexperimental、phylogenetic、computational、author/curator statement、electronic annotationを分けています。このためtoolのhit confidenceを科学的evidence strengthとして扱いません。BLASTのe-valueが強くても、generic similarityだけなら正のweightはありません。

KnetMinerはgeneとtraitの間のevidence pathを追跡可能にする考え方を採ります。本workflowでもscoreの合計だけでなく、各contributionからsourceと根拠へ戻れることを必須にします。ただしKnetMinerのscoreを転載したものではありません。

## 独立性と重複防止

同一BLAST/domain/orthology chain、同一expression experiment、同一mapping resultから派生したcodeは独立voteにしません。rankingでは同じevidence family・directionの最大寄与だけを採用します。Top判定の独立性には、別のcitable sourceまたは明示された`independence_group`が必要です。

## 未解決点

- `homology_supported`と`expression_detected`のlegacy weightは`experimental`です。新規データは原則contextまたはtrait-specific codeへ変換します。
- 作物、形質、mapping precisionによるweight calibrationは未完です。#69の摂動・leave-one-family-out結果で感度を表示し、単一の最適値を主張しません。
- against evidenceはsupportと相殺して消さず、conflict flagと元sourceを保持します。

## 参照

- Evidence & Conclusion Ontology: https://evidenceontology.org/
- Gene Ontology evidence codes: https://geneontology.org/docs/guide-go-evidence-codes/
- Hassani-Pak et al. 2021, KnetMiner: https://doi.org/10.1111/pbi.13583
- Sequence Ontology: https://www.sequenceontology.org/
