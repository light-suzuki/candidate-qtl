# map-based cloningとの関係

## 設計原則

> **Candidate triage does not replace map-based cloning.** It reduces the amount of mapping that must be completed before biological hypotheses can be tested, while preserving mapping as an independent source of evidence and a fallback when candidate evidence is weak.

日本語で言い直すと、このworkflowは **mappingを続けながら、ある程度信頼できる区間が得られた時点で区間内候補を網羅的に調べ、強い候補があれば早期にWet validationへ進める**ための反復型map-based cloning支援系です。mappingの代替ではありません。

## なぜこの位置づけなのか

従来のmap-based cloningでは、候補が数個になるまで区間を狭めるために追加個体・組換え体探索・marker設計・世代更新が必要で、そこに時間と労力の大半が消えます。

一方、reference genome / annotation / comparative genomics が揃った現在は、数Mb〜十数Mbの区間でも区間内全遺伝子を計算機的に確認し、既知機能・orthology・expression・genetic evidenceを統合して優先順位を付けられます。

つまり削減できるのは **「候補検証を始める前に完了させておくべきmappingの量」** であって、mappingそのものではありません。

## 反復ループ

```text
        ┌───────────────────────────────────────────────┐
        │                                               │
        v                                               │
 1. mapping を十分なところまで進める                     │
        │                                               │
        v                                               │
 2. 区間内の全候補を保持して評価する                      │  候補が弱い / Wetで否定された
        │  (interval atlas / family / orthology)        │
        v                                               │
 3. よく研究された植物の知見を探索priorとして使う          │
        │  (Arabidopsis / Medicago / Lotus ...)         │
        v                                               │
 4. 強い候補があれば Wet へ進む                          │
        │  (parental sequencing / PCR / Sanger)         │
        v                                               │
 5. 判定 ──── 強い候補なし or 否定 ─────────────────────┘
        │                          (additional mapping /
        │                           regulatory / SV / PAV /
     確定                          annotation rescue へ)
```

ステップ2〜4は、ステップ1を打ち切る理由にはなりません。**candidate triageとadditional mappingは並列かつ反復的に使います。**

## 具体的な取り決め

### `Top candidate` ≠ `mapping complete`

Top候補が出たことは、mappingが完了したことを意味しません。mapping evidence（区間、組換え体、共分離）は候補証拠とは独立に保持し続け、候補が強いという理由でmappingを止めません。

### `mapping gate` 通過 ≠ fine mapping を尽くす義務

逆方向も同様です。`mapping-gate` を通ったなら、極限までfine mappingしてからでなければ候補検証へ進めない、ということはありません。区間が信頼できる段階で候補検証を始めてよく、それがこのworkflowの狙いです。

### 候補が弱い/否定された場合の標準分岐

これは失敗ではなく、設計上の標準ルートです。

- additional mapping / recombinant screening へ戻る
- regulatory / SV / PAV lane（`promoter-review`, `sv-plan`）へ広げる
- annotation rescue（`annotation-rescue`）で欠落gene modelを回収する
- `multigenome-private` の分母を見て、候補が少ないのが生物学的事実か注釈・callabilityの穴かを判定する

candidate setが空、あるいは全件 `weakened` になったときに「候補なし」で止まらず、上のどれへ戻るかを明示することが要件です。

### biological rank と action priority を混同しない

**候補の生物学的順位**と、**次に行うべき実験の優先順位**は別物です。

- 生物学的にTop1の候補でも、それを検証するassayが高コスト・低情報量なら、次の実験としては最適ではない
- 逆に、複数候補を一度に識別できるmarkerや組換え体スクリーニングは、単一候補のSangerより先に来ることがある

`next-experiment` はこの区別を前提に設計します。候補順位そのものとは独立に、

- additional mapping / recombinant screening
- parental candidate sequencing
- marker development
- regulatory / SV follow-up

を比較し、**現在残っている候補・仮説を最も多く識別できる実験**を優先します。

## 最終的な位置づけ

このworkflowは「原因遺伝子ランキング器」ではなく、**mapping・Dry解析・Wet validationの総コストと時間を減らすdecision-support workflow**です。評価するときも、Top-k精度だけでなく「原因遺伝子に到達するまでの総工数」で見るのが正しい見方になります。

関連: Weight・mapping閾値の運用方針 / [次実験選択](NEXT-EXPERIMENT.ja.md) / 共同研究者向けガイド
