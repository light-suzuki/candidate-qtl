# Canonical Route B: B0 → 初期候補採用／支持確認 → Wet

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

初期Bの種名表示差を解決する場合は、既存record-cacheへ公式UniProt taxonomyの原JSONを `taxonomy-<taxonId>.json` として別保存できます。同じ凍結endpoint自身のtaxon ID、元claimと公式scientificNameの完全一致（case/空白のみ正規化）、その原recordのscientificName・commonName・synonymsからの表示全文再現がそろった場合だけ表示差として監査します。原claim/endpointは書き換えず、taxonomy原fileのSHAと解決記録をtransfer manifestへ保存し、sessionでは通常の入力hash封印を適用します。candidate freezeとWetのprojection再構築も、このmanifestで封印したtaxonomy原fileを再照合して使います。自動network取得・推測alias・一般的な括弧除去はせず、cache不足・真の種差・namespace差や別endpointの軸借用を許可しません。これはinitial候補のidentity監査のみで、strict support/legacy方針を変更しません。

Route Bのproductionには、**新規通常の初期Wet候補採用**と、**verified supportへの厳格な昇格**という異なる入口があります。下図の通常初期候補に独立原文再読は不要です。

新規bundleは配送metadata `initial_material_projection=similarity_annotation_context_v1`を保存し、`sequence_similarity_hits`内の具体的な蛋白名・header descriptionも初期調査材料に含めます。元のhitフィールド、member・参照・種・DBの対応とsource pointerを保持したUNVERIFIED similarity annotation contextであり、`reference_function`、orthology、verified supportへは昇格しません。ID・種名・出典のみ、空欄、hypothetical/uncharacterized headerは材料ありに数えません。markerのない旧freezeは従来の初期材料計数を再生し、凍結bundle・回答は変更しません。

```text
interval genes（全gene/copy台帳）
   ↓
route-b-prepare（BLAST・参照機能・B0 phenotype-blind freeze）
   ↓
route-b-semantic-bundle（family/orthogroup等の単位に圧縮）
   ↓
LLM：参照機能の調査 → 表現型との中立的照合
  ├─ 新文献が必要：discovery保存 → --extend-bundle → 再freeze
  └─ 凍結資料内：既存unit/sourceを使う
   ↓
route-b-semantic-import（claim・論文別9項目から初期Tier・unit対応）
   ↓
新規通常：論文別初期Tierの高い方から有用5報をprimaryへ（元の全根拠は保持）
   ↓
source-verify（新規通常は同じprimary集合の書誌ID/titleと既存gene ID等を確認）
   ├─ 新規 initial_bibliography_v1
   │    → 参照/memberの対応照合・gene Tierの適用
   │    → 適格なB review候補 → 第1段Wet（独立原文再読は不要）
   │
   └─ verified support／旧 source_function_v1
        → 独立原文・機能review＋transfer verification
        → route_b_verified_candidates.tsv → supported_member_union
   ↓
A由来の候補・根拠とgene/copy単位で合流
   ↓
wet_candidate_handoff → wet_candidates.tsv
   ↓
必要な場合に境界比較・採用／保留
```

通常の文献調査は元の参照遺伝子または共有機構ごとに約20候補を調べ、各論文の初期Tier・理由を記録して高Tier順の有用5報を次段へ渡します。取り込み時の`paper_evidence_table.tsv`には各論文の観察事実・確認箇所・Tier根拠と`primary` / `audit_only`の選択状態を記録します。copy間で調査を共有し、元の全論文記録・反証・選外理由を保持します。5報への選択は既存の最強適用根拠によるgene Tierを変更せず、不足だけで適格候補を除外しません。`paper_followup_requests.json`は空欄や具体的な上位Tier手がかりを原担当に限定差し戻しするための依頼であって、実モデルの再実行ではありません。第1段候補を`supported_member_union`へ自動昇格させず、旧bundleは旧policyで再生します。

## なぜこの経路がcanonicalか

Route Bは「区間内にある遺伝子と相同遺伝子の既知機能を調べ、その事実から今回の表現型への具体的なつながりを判断する」経路です。事実確認と条件付き推論を分け、候補にするための物語を作りません。LLM/agentは資料の意味と適用を判断し、機械はID・種・出典・引用・hash・凍結member対応を照合します。照合成功は因果の証明ではありません。

- B0 inventoryはphenotypeを見る前にfreezeする。
- semantic bundleはB0の全gene identityを保持する。
- semantic responseはまずproposalです。新規通常policyでは書誌・参照/memberの条件を満たせば第1段review候補になれますが、verified supportへは別の厳格な条件が必要です。
- `route_b_verified_candidates.tsv`はsupported unionへ、適格な`route_b_review_candidates.tsv`は別入力からWetへ渡す。reviewはverified supportではない。
- group-level supportでcopyが解けない場合は複数copyを保持する。
- provider failure / responseなしは`not_evaluated`でありnegativeではない。
- family / orthogroup / domain / reference gene / literatureのいずれも無いunitへは現行の凍結prompt生成規則でrequestを出さない(#496 Phase 2)。unitと全memberはbundleに残り、状態は`not_requested_no_reviewable_material` = **not_evaluated**であって、調査済み空claimや生物学的absenceではない。未依頼unitは`route_b_semantic_units_not_requested.tsv`に記録し、manifestの`compression_observability`が`semantic_unit_count` / `requested_unit_count` / `skipped_unit_count`を別々に持つ。追加資料が得られた場合は新しいbundle/promptをfreezeし、旧世代の未依頼を調査完了へ読み替えない。

## 新規B-only reviewの条件と互換性

旧bundleは `b_review_membership_policy=source_function_v1` を固定します。出典で確認した参照機能の内容、参照gene/種と区間unitの対応、今回の表現型への具体的な条件付き推論を独立reviewで確認してから、B単独のWet reviewへ渡します。対象種の因果論文や同一表現型の実証、厳密な1:1 copy確定は不要です。Type3の`mechanistic_review_seed`を残せますが、出典IDだけ、汎用kinase/TF/GO、内部知識だけでは足りません。

原著抄録・本文やprovenance付きの具体的curated function recordを使えます。出典なしseedや未探索・取得失敗・対応不明は、元rawと`route_b_review_eligibility.tsv`、全gene台帳へ理由付きで残します。Bの資料不足はA/Cの独立根拠の拒否理由にはしません。旧版のないbundleは`legacy_grounded_review`で再生し、旧freezeを書き換えません。採否規則の正本は [Wet membership](CANONICAL-WET-MEMBERSHIP.ja.md) です。

通常の親Agent入口と新資料のfreeze・必要時の支持確認は [通常運用](AGENTIC-QTL-CASE.ja.md)、固定sessionの再開は [claim review session](CLAIM-REVIEW-SESSION.ja.md) を使います。下記guided-profileのresponse自動検出だけで、Web調査・実レビューまで完了したとは扱いません。

## 提案する主張の粒度

1 claimは、出典で個別に確認できる1つの観察事実・機能・過程に限定します。
別の遺伝子・種・表現型についての実験を一つの強い主張に合成せず、
各`source_ids`がそのclaimを支持するように分けます。ゲノム資源の論文や相同性だけを機能実験とは扱いません。
`reference_gene`はunitに供給されたreference endpointをそのまま使い、論文中のsymbolへ勝手に置き換えません。

`biological_function` / `mechanism`には報告された機能・過程を記載し、
`phenotype_relation`では観察事実と対象表現型への条件付き推論を区別します。
未検証の拡張・種間移転・copyの不確実性は`uncertainty_note` / `alternative_hypothesis`へ残します。
関連する一部の機構を示すために全表現型の説明を要求しませんが、部分的な証拠で複合主張全体を確認済みにはしません。
これは提案指示の改善であり、source・ID・種・transferの検証を省略するものではありません。

## legacy `interval-reverse-lookup` B1

`guided-profile --profile triple-route`の既存B1 overlayは互換性・診断のため当面残しますが、**production agent/LLM pathではありません**。

- role: `diagnostic_only`
- canonical extensionはB1 artifactを読みません (`legacy_b1_consumed=false`)
- `--shared-route`もlegacy B1 diagnostic overlayの入力に限定されます
- production Wet membershipへの情報はcanonical semantic laneのverified候補と適格review候補を使います。legacy B1を別の採否入力にはしません。

つまり、B1 prompt packageとsemantic laneの2本を競合する本流として扱いません。

## guided-profileからの到達

既存CLI surfaceを増やさず、同じ`guided-profile --profile triple-route --execute`を2段階で使います。

### 1. 最初のexecute

semantic responseがまだ無い場合:

```text
state = semantic_review_required
provider_state = not_evaluated_no_semantic_response
```

`b04_canonical_semantic/01_bundle/route_b_semantic_prompts.jsonl`を出して停止します。これはmissではありません。

manifestの`expected_semantic_run_dir`は:

```text
<profile-out>/b04_canonical_semantic/02_llm_run
```

を指します。

### 2. semantic responseを作る

既存`llm-run`を使う場合は、上の規定directoryへ出力します。

```text
candidate-qtl llm-run \
  --prompts <profile-out>/b04_canonical_semantic/01_bundle/route_b_semantic_prompts.jsonl \
  --out-dir <profile-out>/b04_canonical_semantic/02_llm_run \
  --provider ... \
  --model ...
```

agentが直接回答を作る場合も、同じcontractの:

```text
<profile-out>/b04_canonical_semantic/02_llm_run/llm_candidate_responses.json
```

を置きます。canonical semantic promptsを`llm-run`へ渡すと、同じdirectoryへ
`llm_run_manifest.json`が生成され、sibling bundleを含む
`bundle_sha256`、`prompt_sha256`、`response_sha256`が自動記録されます。
siblingの`route_b_semantic_bundle.json`が無いcanonical入力は実行前にfail closedします。
再実行時は3 artifactの現在byte hashと照合し、欠落・不一致・provider/modelだけの
legacy manifestはfail closedします。responseが無いfirst passは従来どおりbundle/promptを
生成して`semantic_review_required`で停止します。通常の非canonical promptsでは既存の
`llm-run` manifest互換性を維持します。

### 3. 同じguided-profileを再execute

同じ`guided-profile --profile triple-route --execute ...`を再実行すると、canonical extensionが規定位置のresponseを自動検出し:

```text
semantic import
 -> source verify
 -> transfer verify
 -> candidate freeze
 -> supported A/B union + eligible explicit review inputs
 -> Wet handoff
```

まで進みます。

この方式はprovider-neutral artifact boundaryを保ちます。Wet budgetの既定20は目安であり、候補の切断上限や調査完了条件ではありません。独立review/追加資料は既存入力として明示し、response自動検出を意味レビュー実行と混同しません。

## upstream/provider gap

canonical extensionはpartial profileをbiology missへ変換しません。

- phenotypeなし → `not_runnable_missing_phenotype`
- B0 atlasなし → `not_runnable_upstream_b0_missing`
- semantic responseなし → `semantic_review_required`
- source provider/cache不足 → source verification側の`not_evaluated`
- Route B verified candidateが0でもRoute A supported memberは保持し、handoffはdegraded coverageを明示する

## Membership boundary

最終membershipの正本は`wet_candidates.tsv`です。Route B semantic candidate freezeはA/B unionへの入力であり、独自Top-Nを作りません。

詳細は [CANONICAL-WET-MEMBERSHIP.ja.md](CANONICAL-WET-MEMBERSHIP.ja.md) を参照してください。
