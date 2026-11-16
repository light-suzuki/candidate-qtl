# A/Bの調査材料を保持する

通常の新規A/Bは `initial_bibliography_v1`を使います。元のA/Bの具体的な機能判断、文献ID・タイトルの安価な照合、凍結された参照・区間対応に基づき、第一段の**本候補**として採用します。初期Tierは元の提案時にT1〜T6/Uと理由を付けます。第一段では独立した本文再確認や厳密なTier確定を要求せず、確認状態の警告・保留理由にも使いません。本文確認は第二段で必要な場合だけ考えます。verified supportと候補採用は別属性で、Direct C、identity・配列・member/hashの厳密な確認は維持します。旧 `source_function_v1` と版のない凍結世代は元の方針で再生します。正本は[Wet採否](CANONICAL-WET-MEMBERSHIP.ja.md)です。

## Bの参照は不足時に補う

既存の使える区間geneと参照の対応・family/domain・機能注釈を先に利用する。参照DBの新規取得・構築、全geneの厳密orthology、完全なfamily-indexを毎回の必須工程にしない。空の表や区間geneと対応しない表は「参照あり」とみなさない。材料不足はgeneごとに記録し、健全な既存材料を失わない。

このworkstationのBの配列照合では、Arabidopsis thaliana（At）を基準資産とし、材料が不足するgeneに限って近縁のcrop modelを最大1種追加する。peaではMedicago truncatula（Mt）を使うAt+Mtの範囲で運用し、種や計算を次々に増やさない。これはローカル参照資産と計算量の運用方針であり、species別の候補選抜規則ではない。A/B/Cで調べる一次文献の対象種をAt/Mtに制限しない。

蛋白配列から補う必要がある場合は、`route-b-prepare --proteins` と既存のlocal protein DBを使う。DBは明示的な`--protein-db`、または`CANDIDATE_QTL_PROTEIN_DB`で指定する。自動探索先は設定済みDBと既定のローカルDBに限られ、PC内の全ディレクトリを検索しない。`--no-auto-local-db`は自動探索を無効にする。参照が使える案件に追加検索を強制しない。

```sh
# 参照材料不足時の既存入口。各pathは案件の実在する資産へ置換する。
candidate-qtl route-b-prepare --candidates INPUT/candidates.tsv \
  --proteins INPUT/candidate_proteins.fa --protein-db RESOURCES/protein_db \
  --out-dir OUTPUT/b0 --case-id CASE
```

機能テキストのcorpusと検索先の蛋白DBは別資産。corpusに配列列がなくてもBLASTは可能であり、corpusの全accessionに配列を付け直してから始める必要はない。利用可能な`--reference-function`を追加すると、対応する参照の機能・出典を調査材料として渡せる。DB・query・検索設定と出典hashを記録し、取得済み複数hitをbest-hit一つへ縮約しない。検索の取得上限は記録され、DB全体の相同配列網羅を保証しない。

`route-b-prepare`のlocal検索は、既存のmaterial判定で未解決の区間geneだけをexact IDで元FASTAから選ぶ。headerと配列のbytesを保持した`route_b_reference_query.fa`を保存し、元FASTAを変更しない。manifestの`providers.homology_search.query_subset`に元/queryのSHA-256、元・query・スキップ・配列欠損・区間外IDと件数、検索を完了したIDを記録する。区間外IDは検索へ混ぜず台帳へ残す。重複/空IDはエラーとし、全区間gene inventoryは保持する。queryが空ならBLASTとDB探索を省略して`not_evaluated_no_protein_query`とする。検索しなかったgeneや配列欠損geneを`evaluated_no_threshold_hit`にしない。

区間protein同士の自己BLASTはparalog群の補助材料であり、他種の機能同定やtransferの代替ではない。全gene inventory、identityの有無、機能材料の有無、semantic解釈の処理状態、最終Wetを別々に確認する。DB不足や検索失敗は未評価として残し、成功した既存材料によるBの継続を止めない。

Atと必要時の1種で調べてもidentity/functionが付かないgeneは、限られた知識の下でのunknownとして残る。それ自体をworkflowの不具合とはしない。実際に検索したno-hitは追加調査の低優先度台帳に残せるが、非原因の証明でも無条件削除の根拠でもない。identity handleと機能・文献のsource materialは別軸であり、配列一致だけでsemantic review完了や形質支持へ昇格させない。

## 原文と未検証材料の保持

`route-b-prepare --reference-function` は従来のTSVに加え、取得済みUniProt JSONのディレクトリを受け取る。タンパク質名、SIMILARITY、InterPro/Pfam/PROSITE、GO、keywords、参照種と出典ハッシュを、逆引き表からsemantic unitまで渡す。参照domainは区間遺伝子自身で測定したdomainとは区別する。同じaccessionの複数原本は曖昧な先勝ちにせずエラーとする。

全レコードの文献は検索用リンク候補である。phenotype出典にはDISRUPTION PHENOTYPEに直接付属するPubMedだけを記録し、`verified_source_ids` は生成しない。注釈だけでは実験的支持や因果支持にならない。unitの `investigation_material` は機能記述・注釈・文献の有無を示す調査準備情報で、支持件数でも実行完了件数でもない。材料不足を記録し、候補判断を変え得る具体的な手掛かりがある場合に追加取得を検討する。

Mt等のnative蛋白IDを持つDBでは、配列hitと機能テキストのIDが一致するとは限らない。まずversion付きnative accessionを公開providerの記録へ照合し、種・配列・原記録hashを確認する。別DBの注釈を使う場合は明示cross-reference、または同種の全長完全一致配列で対応を記録し、native IDを`reference_gene_id`とする参照機能TSVへ接続する。元providerの`source_id`・`source_url`、元accession、対応根拠・原本hashと未対応台帳を保持する。symbolや旧IDのstemから対応を推測しない。完全一致配列だけでは一意なlocusや別記録の変異体実験との同一性は確定せず、形質支持を生成しない。`putative`等の推定注釈はそのまま保持し、hypothetical/uncharacterizedは機能既知と数えない。

ここでの機械的な対応条件は、別DBの行を同じ実体の記録として結合する場合に適用する。原論文の報告geneから別種の配列proxyへ知見を橋渡しする生物学的判断は元の調査LLMが行い、通常候補では報告gene自身の配列取得やproxy側の同じ表現型報告を必須にしない。原報のgene・種・観察は保存し、実際に投影するproxyのidentity・配列と対象geneへの対応は検査する。詳細は[報告遺伝子・代理参照・対象座位](CANONICAL-WET-MEMBERSHIP.ja.md#報告遺伝子代理参照対象座位を分ける)に従う。

Atと追加モデルの取得済み対応・identity-hit表をそれぞれ合流し、既存freezeを上書きせず新しいB0を作る。`--reference-function`を繰り返してAtとnative IDの機能TSVを渡すと、同じreferenceに属する複数の出典カードもsemantic bundleへ保持される。取得済みhitを利用する場合は再BLASTを要求しない。下流の`route-b-semantic-bundle --identity-hits`には、両DBの全hitを保持した合流表を明示する。scalarなreverse lookupの最初の参照だけを全資料と解釈せず、各referenceのカードと対応台帳も確認する。素材追加でprompt hashが変わった場合、旧semantic回答を新promptの実行結果として再利用しない。

Route Aのcanonical bridgeは、同一の凍結配列と検索原本ハッシュを検証した上で、有意な部分一致を `partial_alignment_requires_domain_or_copy_review` として保持する。identityとE-valueは検索時の凍結閾値を使い、両配列の整列部分が凍結済みminimum query length以上であることを求める。全長coverageを満たさない候補は、完全長homologyや形質支持ではない。支持unionへ入れず、競合コピーを残してWetレビューへ渡す。生の検索結果が無い旧compact表から部分一致を推測しない。

1,000遺伝子は検証サイズであり上限ではない。全memberを中間表現・分割台帳に保持し、分割ごとのTop-Kを設けない。取得原本全体をコピーごとに埋め込まず、必要な注釈と原本ハッシュを渡す。既存sharderの文字数予算を超える単一unitは明示的なoversized unitであり、黙って切り捨てない。今回の回帰検証はモデルを呼ばず、1,000遺伝子の注釈受け渡し・全件分割保持とA→Wetのreview-only境界を確認する。特定ケース・遺伝子・familyの例外を用いない。

非盲目の保存データ再生は接続診断として扱う。新しい予測精度のベンチ結果とはしない。厳密な盲目性を要求しない場合にも、支持・未検証・未評価の区別と成果物ハッシュは維持する。

通常Bの機能調査では、注釈・相同機能・内部知識を探索手掛かりとして初めからWeb/論文を利用できる。既知機能を先に整理し、具体的な表現型へのつながりがなければ理由付きで非提案にできる。未探索や取得失敗とは分ける。互換性用B-only Wet reviewの`source_function_v1`では、外部資料の実機能記述と対応・条件付き推論を確認する。新規の通常第一段は上記の初期候補policyを使う。原著1〜2報のsoft budgetは未採用で、必要な矛盾・copy差の調査を禁じない。

opt-in探索索引はannotation description、GO、family/domain、hit_descriptionを手掛かりとして保持し、原unitのowner/hashへ戻れる。索引の既定値はOFFのまま。索引・共有配置で保存bytesが減っても、実際の原文読解・LLM判断の削減を実証したことにはならない。同一参照機能質問の共有とclaim/copy固有適用は [claim review session](CLAIM-REVIEW-SESSION.ja.md) を参照する。

### 初期判断用の共有資料と集計

既存の `scripts/agentic_route_b.py prepare --research-bundle ... --research-bundle-sha256 ...` は、探索索引に加え、`research/manifest.json` の `initial_material.pages` に完全な機能・表現型・pathway・GO説明・tissue/stage等の共有資料を出す。同じfunction名でも補助情報が違えば別資料とし、同一contextだけを共有する。各bindingに元unit/member/source pointerと参照・出典・evidence情報を保持する。長い単一資料は切らずに `over_bound_bytes` を記録する。索引のpreviewだけを再要約して、補助情報を捨てない。

元担当者は判断を次の形式で保存できる。これは正式なclaim回答とは別の作業記録で、候補やTierを変更しない。複数参照の既知・不明は別行に残す。`assessed_fields` は読んだ欄を記録できる場合に付け、省略を全欄読了と解釈しない。旧判断にpointerがない場合も理由と判断を保持し、`judgments_without_material_pointer` として別集計する。記録欠落だけで原判断を消したり再研究を求めたりしない。提出されたpointerが別geneを指す等の矛盾は拒否する。

```json
{"bundle_sha256":"FROZEN_BUNDLE_SHA256","assessments":[
  {"unit_id":"u1","gene_id":"g1","state":"assessed_no_specific_bridge",
   "reason":"原担当者の判断理由","used_material_pointers":["/units/0/genes/0"],
   "assessed_fields":["reference_function"]}
]}
```

```sh
python scripts/agentic_route_b.py summarize-initial \
  --bundle FROZEN_BUNDLE --bundle-sha256 FROZEN_BUNDLE_SHA256 \
  --assessments initial_assessments.json --out NEW_SUMMARY_DIR
```

機械はbundle hash、unit/gene、pointerと欄の存在を照合し、全geneについて材料有無・原判断の内訳・読了記録のない欄を別々に出す。`assessed_claim` / `assessed_concrete_lead` / `assessed_no_specific_bridge` と `assessed_insufficient_material` が混在すれば `partially_interpreted` とし、既知の判断を消さない。判断のないgeneは `not_assessed`、不足判断のみは `insufficient_only`。これらは研究者の記録の集計であり、機械による意味判定や読了認証ではない。未読欄の全件確認や全geneの機能解明は完了条件にせず、情報不足だけで閾値緩和・参照追加・全件検索を要求しない。
