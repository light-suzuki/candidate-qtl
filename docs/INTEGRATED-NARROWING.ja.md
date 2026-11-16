# 統合候補絞り込み（integrated narrowing）

`integrated-narrowing` は **candidate membership を作る層ではなく、すでに確定した候補集合の Wet review 順序を付ける層**です。原因遺伝子の確率推定でも、新しい score でもありません。

Wet集合には、verified supportと区間対応が確認されたA/B提案のreview候補が共存できます。
後者の証拠確認が未完という状態を保持し、順序付けの再実行で候補を消したりsupportへ昇格させたりしません。

production では `wet_candidate_handoff` が先に `wet_candidates.tsv` を凍結し、`integrated-narrowing` はその集合だけを受け取ります。

```text
supported_member_union（verified subset）＋ 区間に対応した適格A/B review候補
（optional Cは既存のsupport/review member gateを経由）
        ↓
Wet budget gate
        ↓
wet_candidates.tsv       # production membership の正本
        ↓
ranked / family / copy context を同じ集合へ射影
        ↓
integrated-narrowing     # ordering only
        ↓
stage_a_narrowed         # membership は不変
```

## 不変条件

- `wet_candidates.tsv` が唯一の production membership artifact
- narrowing は candidate を追加・削除・padding・Top-N truncate しない
- `dry_score` / `dry_rank` / `dry_group` / `screening_rank` を再計算しない
- `rank_contribution=none`
- optional evidence の欠落は `not_evaluated` であり negative ではない
- ordering 後の gene set が frozen Wet membership と一致しなければ fail-closed

## 入力

必須は ranked candidate table です。family/copy/convergence 等は optional です。

| 入力 | 用途 | 欠落時 |
|---|---|---|
| ranked candidates | base ordering | 実行不可 |
| family recovery | gene → recovered family | `not_evaluated` |
| convergence | family recovery state | `not_evaluated` |
| family copy discrimination | copy state | `not_evaluated` |
| family copy shortlist | first/co-first ordering | 空 |
| remaining candidate set | survival/priority context | `not_evaluated` |

渡された gene ID が ranked universe に存在しない場合は fail-closed です。production handoff はその前に、gene-bearing context table を frozen Wet membership へ射影します。

## recovery state

recovery state は ordinal ですが、causal score ではありません。

```text
strict_convergent
> a_and_b_supported
> a_only
> b_only
> none
> not_evaluated
```

`not_evaluated` は「情報なし」、`none` は「評価したが knowledge-route recovery なし」です。いずれも candidate deletion の根拠にはなりません。

## copy state と #499

`family-to-copy` の状態は3種類に分けます。

### 1. validated copy-set state

以下は homology layer 上で family/copy set が validated された状態です。

- `preferred_copy`
- `co_preferred_tie`
- `family_supported_copy_unresolved`
- `tandem_ambiguous`
- `paralog_ambiguous`
- `copy_weakened`

### 2. evaluated but unresolved context

以下は「copy 判別を実行したが、まだ決められない」状態です。

- `insufficient_discriminating_evidence`
- `family_resolution_required`

これらは単独では family を validate/promote しません。ただし Route A/B で family がすでに独立に回収されている場合、copy 未評価より不利にしてはいけないため N2 を維持します。

```text
recovered family + not_evaluated
                 → N2

recovered family + insufficient_discriminating_evidence
                 → N2

recovered family + family_resolution_required
                 → N2
```

### 3. copy-set identity exclusion

`not_in_validated_copy_set` は上の「未解決」と同義ではありません。

これは issue #184 の contract 上、**coarse family bucket には入っていても、その validated homologous copy set の member ではない**ことを表します。そのため、拒否された coarse family identity だけを根拠に N2 を継承しません。

```text
recovered coarse family
+
not_in_validated_copy_set
        ↓
その family 由来の N2 promotion はしない
        ↓
他の independent/base evidence へ fallback
```

重要なのは、これは candidate deletion ではないことです。gene は frozen Wet membership に残り、別 family hypothesis、mapping、base rank、その他の evidence でレビューできます。

## tier 規則

| Tier | 条件 | 主な reason |
|---|---|---|
| `N1` | validated copy set の `preferred_copy` | `validated_family_preferred_copy` |
| `N2` | validated family/copy の non-single state、または recovered family + copy `not_evaluated` / 上記2つの unresolved context | `validated_family_copy_unresolved` / `family_recovered_copy_*` |
| `N3` | N1/N2 でなく base `dry_group` が `top` / `reserve` / `rescue` | `base_top_reserve_rescue` |
| `N4` | その他 | `evidence_insufficient` / `base_rank_only` |

`preferred_copy` は「最初に Wet/sequencing する copy」という operational priority であり、causal declaration ではありません。

## first Wet priority

`first_wet_priority` は次の順に既存 artifact から受け取ります。

1. `family_copy_shortlist.tsv` の `first` / `co_first`
2. remaining-candidate artifact の `high` / `medium` / `low`
3. 無ければ空

これは numeric score ではありません。

## canonical Wet handoff との境界

productionの`wet_candidate_handoff`は#817の`upper_evidence_scoped_v3`を使います。累積25を初めて超えるTierだけを比較し、上位は保持、後続は全gene台帳に残します。採用/保留は元の調査担当者の理由付き判断を`original_assessment_v1`として機械照合し、別の生物学reviewerを要求しません。20未満だけで追加探索せず、25で強制切断しません。旧`evidence_scoped_v2`・版なしscopeは独立reviewを含め旧policyで再生します。orderingはこのmembershipを変更しません。詳細は[canonical Wet policy](CANONICAL-WET-MEMBERSHIP.ja.md)に従います。

以下のfillはlegacy strict/`supported_only`/明示refinementの互換規則です。

- 7 supported genesなら7件すべてを`tier1_verified_support`として保持し、空きreview枠があれば区間内の`tier2_reference_phenotype`→`tier3_family_or_domain`の順でreview coverageを追加する
- fill rowは`membership_basis=interval_review_coverage`で、verified support / score / route flag / rank contributionを新しく持たない
- `tier4_dark`はfillに使わず、reviewable候補が不足した枠は`unfilled_budget_slots`として残す
- `wet_budget`超を任意件数へ切らず、通常のWet handoffを出す。件数だけで停止・警告・理由説明を要求しない
- 固定容量を明示した`--budget-policy strict`だけは、bounded verified refinement / reviewへ進む
- `--wet-fill-mode supported_only`は診断・比較用で、canonical defaultのtiered fillとは分ける

**`integrated-narrowing`自身はpaddingをしません。** 先に凍結された`wet_candidates.tsv`がtier1だけでもtier2/3 review coverageを含んでいても、そのexact membershipを並べるだけです。eliminated/unsupported gene が ordering slot を消費しないよう、ranked table と family/copy context は先に同じ membership へ射影されます。

## 設計原則

> 評価して「まだ決められない」ことは negative evidence ではない。一方、validated copy set から外れたことを単なる missingness として family support に戻してもいけない。

この区別により、`not_evaluated != negative` と homology/copy-set validation の両方を維持します。
