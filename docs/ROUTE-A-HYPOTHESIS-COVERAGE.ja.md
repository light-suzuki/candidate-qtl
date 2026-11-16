# Route A causal-layer coverage contract

Issue: [#658]()

この機能はRoute Aの仮説をscoreするものではない。phenotype-onlyの初期仮説が1つの因果層へ固定化していないかを、truth・QTL区間・Route Bを見ずに監査する契約である。

## 固定ラベル

`causal_layer` は次の6値だけを許す。

- `perception`
- `signaling`
- `regulation`
- `execution`
- `alternative`
- `unclear`

ラベル数を均等にはしない。表現型から特定層へ局在するのが自然な場合、その偏りは保持する。ラベル、LLM confidence、coverage verdictはいずれも `rank_contribution=none` で、candidate scoreやverified supportを作らない。

## 実行順序

```text
phenotype-investigate
  -> phenotype-localization-packet
  -> fresh LLM/agent A
  -> phenotype-localization-freeze
  -> fresh LLM/agent B: labeled hypotheses (target 12; count is advisory)
  -> mechanism-import
  -> mechanism-coverage-audit
       pass: downstream verification
       gap: fresh LLM/agent C, at most 3 complementary hypotheses
            -> mechanism-import --proposal-round 2
            -> mechanism-coverage-audit (initial + gap table)
            -> downstream verification
```

localization payloadは次の厳密JSONで、`supporting_phrase` はphenotype文字列からのverbatim phraseでなければならない。

```json
{
  "causal_localization": "broad|upstream|signaling|regulation|execution|unknown",
  "localization_confidence": "high|medium|low",
  "supporting_phrase": "verbatim phenotype phrase"
}
```

freeze後のlocalization hashがmechanism importのprovenanceと一致しなければcoverage auditは停止する。初期portfolioは目標12行だが、行数だけで受理/拒否しない。well-formedでfamilyが一意なら実際の行数を保持し、補充は最大3行までcoverage gapに応じて要求する。同一familyは機械的に先勝ちでdeduplicateし、削除したhypothesis IDをauditへ記録する。機構の広がり（5件以上）は advisory として記録し、coverage不足は`review_required`へ送る。

補充promptには初回と共通のJSON template、使用可能なラベル、参照遺伝子・出典が不明な場合の扱いを含める。独立agentに初回の履歴やschemaの推測を要求しない。補充行のIDはimporterがround別に付け、初回IDとの衝突を避ける。

## 偏りの扱い

- `broad` / `unknown` またはconfidence `low`: `perception` / `signaling` / `regulation` / `execution` の欠落を補う。単一ラベルが60%を超える場合は、未優勢層を最大3件の範囲で補う。
- confidence `high` の局在型: 局在層への偏りを許し、隣接層1件以上と `alternative` 2件以上を要求する。
- confidence `medium` の局在型: 局在層への偏りを許すが、`unclear` 以外の別層を2件以上保持する。

これらは件数の強制均等化ではなく、coverage holeの検出である。最大3件で穴が閉じない場合は `review_required` とし、15件を超えて生成しない。

## blind境界

localization promptはphenotype packetだけを見る。補充promptへ渡してよいのは同じphenotype packet、要求ラベル、既使用family名だけである。次は渡さない。

- QTL interval / candidate IDs / interval annotation
- Route B output
- truth gene / truth family / answer paper / score
- benchmark historyまたは前回生成prose

`--candidates` はleak検査専用で、promptへ挿入されない。候補IDがpromptまたはimport済みmechanism表へ現れた場合はfail closedする。coverage auditは各mechanism表に対応する `mechanism_import_manifest.json` を必須とし、`interval_informed=false`、blind discovery claim、proposal round、出力SHA-256、先行localization freezeのSHA-256を再検証する。初回importはfamily-node portfolio成立も必須とする。

このhash照合が証明するのはimportとfreezeの対応であり、外部LLMが実際に見た文脈や生成時刻ではない。`external_generation_order_verified=false`を記録する。実LLM評価では、先にlocalizationをfreezeし、そのhash・使用prompt・生成開始時刻・raw responseをexecutor側のprovenanceにも保存する。既存payloadの後付けラベル変更をfresh generationとして扱わない。

## CLI例

```bash
candidate-qtl phenotype-localization-packet \
  --investigation run/01/phenotype_investigation.tsv --out-dir run/02

candidate-qtl phenotype-localization-freeze \
  --payload localization.raw.json \
  --investigation run/01/phenotype_investigation.tsv --out-dir run/03

candidate-qtl mechanism-import \
  --payload initial.raw.json --investigation run/01/phenotype_investigation.tsv \
  --localization-freeze-manifest run/03/phenotype_localization_freeze_manifest.json \
  --out-dir run/04 --proposal-round 1

candidate-qtl mechanism-coverage-audit \
  --mechanisms run/04/mechanism_hypotheses.tsv \
  --mechanism-import-manifest run/04/mechanism_import_manifest.json \
  --localization-freeze-manifest run/03/phenotype_localization_freeze_manifest.json \
  --investigation run/01/phenotype_investigation.tsv --out-dir run/05
```

auditが `gap_fill_requested` を返したときだけ、`run/05/mechanism_gap_fill_prompt.txt` をfresh agent Cへ渡す。補充payloadを別の`mechanism-import --proposal-round 2`へ通し、初回と同じ`--localization-freeze-manifest`を指定する。最後に`mechanism-coverage-audit`へ`--mechanisms`と`--mechanism-import-manifest`を初期・補充の順で2回渡し、さらに `--gap-request-audit run/05/mechanism_coverage_audit.json` を渡す。要求ラベル、初期表hash、localization freeze hash、再構成したgap prompt hashのどれかが一致しなければ補充表は受理されない。`--proposal-round 2`だけではinterval-informedにならない。区間を見た場合だけ`--interval-informed`を明示し、そのpayloadはblind discoveryを名乗れない。

`guided-profile audit-predict` と `guided-profile triple-route` にもこの順序がcanonical Route Aとして配線されている。LLM生成は外部実行であり、このprofile自体は実行しない。`--phenotype-localization-payload` が無ければmechanism importは開始せず、coverageが `gap_fill_requested` の間はreference resolution以降へpathを渡さない。補充を再開するときは `--gap-fill-mechanisms`、`--gap-fill-mechanism-import-manifest`、`--gap-request-audit` を3点セットで指定する。同一localizationの再開ではfreeze hashを保ち、元gap requestのauditとpromptはhash付きファイルとして保持する。
