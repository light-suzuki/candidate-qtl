# 構造化証拠による二段階候補選定

## 全体像

候補選定は、全QTL遺伝子への軽量解析と、上位候補だけへのorthology精査を分離する。

~~~text
QTL区間
  -> transcript-aware抽出
  -> 全候補のBLAST/DIAMOND、domain、eggNOG、発現、位置情報
  -> canonical evidence.tsv
  -> 一次順位 Top / Reserve / Rescue / Evidence-insufficient
  -> Top/Rescueの少数候補だけreciprocal BLAST
  -> refined evidence.tsv
  -> 最終順位
  -> 親間direct sequencing
~~~

LLMはこの工程に不要である。すべての採点は構造化された証拠行と保存済み重み表から再計算できる。

## transcript-aware抽出

`extract`はtranscriptごとにCDSを保持し、別isoformのCDSを連結しない。代表transcriptは次の順で選ぶ。

1. canonical/representative annotation tag
2. 完全ORF
3. CDS長
4. transcript ID

主な出力は`candidates.tsv`、`transcripts.tsv`、全isoformの`candidate_transcript_cds.fa`、代表CDSの`candidate_cds.fa`、代表proteinの`candidate_proteins.fa`である。

## canonical evidence

`evidence.tsv`の必須列は以下である。

| 列 | 意味 |
|---|---|
| gene_id | QTL内gene ID |
| evidence_code | direct_model_genetic、domain_supported等の証拠種別 |
| direction | support、against、rescue、context |
| confidence | high、medium、low |
| value | hit ID、発現値、距離など |
| source_type | BLAST、HMMER、expression、manual curation等 |
| source_id | ファイルパス、PMID、DOI、stable accession |
| description | 人が確認できる短い説明 |

`source_id` は `reference`（citable id）/ `local_file` / `label` / `none` に分類され、Top eligibility の `independent_source_count`（独立した citable source 数）に使う。分類は **allowlist 駆動**で、`reference` として扱うのは次のみ。

- **ontology/文献トークン**: `PMID:12345`・`GO:0006355`・DOI(`10.xxxx/...`) 等。
- **namespace 付き DB id で local part がその namespace の accession 文法に一致するもの**: `UniProt:Q9LZS3`・`TAIR:AT1G69120`・`RefSeq:NP_001078807`・`EnsemblPlants:...`・`NCBIProtein:XP_...`。namespace は既知でも local part が文法不一致なら citable にしない(例: `PMID:fixture`・`UniProt:candidate` は `label`)。
- **bare accession**: token 全体が特定の既知 accession 文法(UniProt・AGI・RefSeq・GenBank/EMBL/DDBJ・PMC)に大文字 token として一致するもの(`Q9LZS3`・`P93003`・`AT1G69120`・`NP_001078807`・`AAR03725`・`PMC5694491` 等)。

**未知 namespace や未定義の compact token(`gene1234`・`trait2026`・`abc_12345`・`GENE0001` 等)は citable にしない**（`label` に落とす）。旧実装の汎用文法(`英字2–4文字+数字4桁以上`)はこれらを accession と誤認していたが、除去済み。したがって自由文ラベルは `independent_source_count` を水増しできず、Top 昇格を起こさない。**同一解析チェーン(homology/domain/orthology/GO context)由来の複数 accession は 1 independent source に束ねる**(明示的 `independence_group` があればそれを優先)。したがって、1回の BLAST から得た複数 UniProt hit だけでは `independent_experimental_support` にならない。**comparative_function family の証拠のうち、citable reference を持たず解析生成物(ファイルパス = `local_file`)だけを source とする行も同じチェーンへ束ねる。** `annotation-keywords` が出す `trait_function_match` の source_id は「その best hit の注釈ファイル」であり、homology 行を作ったのと同じ blastp 1回に由来する。ここを別 source として数えると、**BLAST 1回 + 位置証拠だけで Top へ昇格**してしまう。PMID/DOI や解決済み accession を持つ trait_function_match / domain_supported は従来どおり独立 source として数える。

BLAST/DIAMONDはidentity、query coverage、E-valueを使ってconfidenceを決める。上位2 hitのbit score差が小さい場合は`paralog_ambiguous`を追加する。HMMER、eggNOG、発現、QTL peak/境界、手動で確認した既知変異体・既知QTG証拠も同じ形式へ変換する。

**形質非依存の注釈は正の形質証拠にしない**。BLAST best hit・任意のHMMER domain・GO annotationの存在は保存性/注釈品質の情報であって、その遺伝子が特定QTL形質に関わる証拠ではない。これらは別コードの `homology_context`／`domain_context`／`go_annotation_context`（direction=context、加点なし）として出力し、trait scoreへ直接加点しない。形質と結びつく支持は、承認済みtrait termに一致した `trait_function_match`、query選定理由とsourceを持つ `curated_homolog_match`、ortholog group配置の `orthogroup_supported`、遺伝学的証拠に限る。保存性はあるが機能未知の候補は、汎用注釈だけを根拠に `conserved_uncharacterized` のRescueとして残し、Topには上げない。

QTL peakのみ（`--peak-window`なし）ではexact overlap（`peak_overlap`）だけを加点する。区間長5%のheuristicは`peak_heuristic_proximal`（context、加点なし）として距離のみ記録し、正の位置証拠にしない。mapping精度で位置支持を決めるため、区間が広いほど位置加点候補が増える逆説を避ける。スコア付きのproximityが要る場合は`--peak-window`（不確実性半幅）・LOD support interval・組換えcoreのいずれかを与える(`peak-context`参照)。

## QTL peak不確実性の取込（peak-context）

位置情報は`peak-context`で段階付き証拠へ正規化する。broad intervalへの所属はcontextであり、加点しない。

- LODプロファイルTSV（chrom・position・lod）から、`--lod-drop`（既定1.5）でsupport intervalを計算する。同一染色体上の別ピークは`multiple_peaks_on_chromosome`、別染色体上の閾値超えは`linked_or_ghost_qtl_on_other_chromosome`として警告する。
- `templates/qtl_peak_metadata.tsv`で、mapping method、support interval種別、組換えcore、marker密度、source IDを明示できる。
- 遺伝子ごとに最も強い位置証拠だけを1件出力する。強さの順は `recombinant_core`（high）> `peak_overlap`（high/medium）> `peak_proximal`（medium、`--peak-window`指定時）> `support_interval`（low）。broad intervalだけの遺伝子には正の位置証拠を付けない。
- 出力`position_evidence.tsv`はcanonical evidence形式なので、そのまま`evidence-rank --evidence`や`evidence-build --manual-evidence`へ渡せる。位置の寄与は`score_contributions.tsv`で監査でき、`peak_context.tsv`に`peak_method`・`peak_confidence`・`support_interval_type`・`position_evidence_reason`を保存する。

~~~powershell
.\scripts\candidate-qtl.ps1 peak-context `
  --candidates C:\private\qtl\extract\candidates.tsv `
  --lod-profile C:\private\qtl\lod_profile.tsv --lod-drop 1.5 `
  --out-dir C:\private\qtl\peak_context --json
~~~

## 発現とtrait-ontologyの取込

発現は機能証拠であり因果証明ではない。`evidence-build --expression`はlong-form TSV（`gene_id`・`tissue`・`value`必須、`stage`・`condition`・`unit`・`replicate`・`source_id`任意）を受け取り、次を守る。

- target/backgroundの倍率は**同じ`unit`内**でのみ計算する。TPM・raw count・microarrayを互いに割らない。
- **backgroundを測定していない場合、`expression_target_specific`を付けない**。target閾値を超えれば`expression_detected`のみ出し、比較不能を`expression_comparison_unavailable`（context、加点なし）として記録する。「発現している」ことと「target特異的」は別で、サンプリング不足を正の証拠へ変換しない。`specificity_class`は`comparison_unavailable`となり、`comparison_available=no`を残す。
- `expression_summary.tsv`にunitごとの行を出し、`specificity_class`・`expressed`・`non_expression`・`comparison_available`・`conditions_observed`・`source_id`を別列で保持する。
- `--target-condition`でcontrol/treatmentなどの条件を絞れる。正規化unitはevidence descriptionにも残す。

形質語はローカルのみで正規化する。`trait-ontology-map`は`templates/trait_terms.tsv`とローカルの`templates/trait_ontology_map.tsv`（term→ontology_id）を突き合わせ、`original_term`を保持したままPTO/PO/CO IDを付ける。`match_type`（exact/synonym/manual/llm_suggested）と`approval_status`を区別し、`llm_suggested`や未承認synonymは`pending_review`のままにする。ontology未対応の語は`unmapped`として残し、ワークフローを止めない。外部APIは必須ではない。

キーワード照合に渡す語彙は`trait-terms build`で統制する。raw trait termsを4クラス（`direct_phenotype`/`organ_stage`/`biological_process`/`mechanism_hypothesis`）に分類し、`direct`/`context`のweight classを付ける。任意でontology ID、人手synonym（承認扱い）、LLM proposal（常に`pending_review`）をmergeする。`--benchmark --blocked-term <遺伝子名>`では既知原因遺伝子・locus名に一致する語をanswer-leakageとして除外し、`_audit.tsv`に理由を残す。承認済みの語だけを出力し、`annotation-keywords --keyword`へ渡す。自由入力の`--keyword`は手動モードとして維持する。

## 採点規則

既定かつauthoritativeなsource of truthはversion・根拠・hashを持つ`templates/evidence_weight_registry.tsv`である。`templates/evidence_weights.tsv`は後方互換のlegacy/custom入力例であり、既定値ではない。実行時には実際に使用したregistry/profileと重みが`weights_used.tsv`およびdecision manifestへ保存される。

- 同じ`evidence_code`が複数DBから得られても、geneごとに最大1回だけ採点する。
- confidenceで重みをhigh=1、medium=0.75、low=0.5に調整する。
- againstは減点、rescueとcontextは採点しない。
- Topには、最小scoreに加えて、強い遺伝学的証拠1種類、または独立した証拠family 2種類以上を要求する。
- 発現だけ、BLAST hitだけ、論文数だけではTopにしない。
- 全証拠が0の候補を遺伝子ID順でTopにしない。

`score_contributions.tsv`に遺伝子ごとの全加減点を保存する。重みは仮説生成とSanger順序のための設定値であり、因果性の確率ではない。

## 自動Rescueと早期分岐

以下はscoreとは別にRescueへ残す。

- CDSがない、またはORFが不完全
- 保存相同性はあるがdomain・既知機能がない
- 対象組織で特異的に発現するが他の証拠が弱い
- QTL peak近傍だが他の証拠が弱い

`candidate_path_decision.json`は、十分なcoding候補がある場合は`test_coding_candidates_first`、Top候補が不足するか機能・遺伝学的支持がない場合は`broaden_to_regulatory_sv_te_cnv`を返す。これは制御変異やSVが原因だという断定ではなく、promoter、UTR、TE、CNV、PAV、未注釈領域を早期に調べる分岐である。

## 一貫実行

~~~powershell
.\scripts\candidate-qtl.ps1 candidate-pipeline `
  --candidates C:\private\qtl\extract\candidates.tsv `
  --blast C:\private\qtl\annotation\blastp_hits_1.tsv `
  --hmmer C:\private\qtl\annotation\hmmer_hits.tbl `
  --eggnog C:\private\qtl\annotation\candidate_eggnog.emapper.annotations `
  --expression C:\private\qtl\expression_long.tsv `
  --target-tissue meristem `
  --manual-evidence C:\private\qtl\manual_evidence.tsv `
  --peak Chr1:250000 `
  --interval Chr1:100000-500000 `
  --protein-fasta C:\private\qtl\extract\candidate_proteins.fa `
  --target-proteome C:\private\reference\target_proteome.fa `
  --relative-proteome relative1=C:\private\reference\relative1.fa `
  --relative-proteome relative2=C:\private\reference\relative2.fa `
  --out-dir C:\private\qtl\candidate_pipeline `
  --top-n 7 --reserve-n 3 --num-threads 4 --json
~~~

一次実行後、`03_orthology_refinement/run_orthology_refinement.sh`を実行し、`orthology-import`でcanonical evidenceへ変換する。得られたTSVを`--refined-evidence`に渡して同じ`candidate-pipeline`を再実行すると、`05_final_rank`が作られる。
