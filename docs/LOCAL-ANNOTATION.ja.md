# QTL区間限定アノテーションとAI非依存リサーチ

## 設計原則

計算対象は全染色体や全proteomeではなく、candidate-qtl extractで確定したQTL区間の候補遺伝子だけである。全ゲノムで一度作ったBLAST/DIAMOND/HMMERデータベースは再利用してよいが、照会するqueryは候補遺伝子集合に限定する。

この設計で節約するのは、query数、結果の人手確認数、外部API検索数、primer設計数である。参照DBの構築自体は初回だけ重いことがある。

## GOの扱い

狭いQTL区間で必要なのは、まず各候補のGO、domain、ortholog、文献根拠を記録する **GO annotation** である。数遺伝子から数十遺伝子の区間でGO enrichmentを行っても背景集合・検出力・多重検定の影響が大きいため、enrichment結果を原因遺伝子の根拠には使わない。

候補ごとのGO根拠は go_evidence.tsv に以下を残す。

- GO IDとterm
- 直接注釈かorthologからの推定か
- source database、accession、URL
- 実験根拠か計算根拠か

## ローカル相同性・domain計算

protein FASTAがある場合、candidate-qtl annotate-planで候補proteinだけを抜き出し、実行scriptを作る。

~~~powershell
candidate-qtl annotate-plan --candidates C:\private\qtl\extract\candidates.tsv --protein-fasta C:\private\annotation\proteins.fa --out-dir C:\private\qtl\annotation --blastp-db ~/blast_databases/reference_proteins --hmmer-db ~/db/Pfam-A.hmm --json
~~~

生成される run_local_annotation.sh は、指定されたBLASTP、DIAMOND、HMMERを候補proteinだけへ実行する。DIAMOND DBとBLASTP DBは互換ではないため、別の引数にしている。

正式なprotein FASTAがある場合はそれを優先する。`extract`はtranscriptごとのCDSを分離し、canonical tag、完全ORF、CDS長の順で代表配列を作るが、`partial_or_unverified_orf`は正式注釈の代替とみなさない。

`extract`はGFF3方言差も吸収する。gene feature欠落・親gene名が未定義・transcript無しでCDSがgeneへ直付き・CDSだけ、といった入力ではplaceholder gene/transcriptを合成し、`is_placeholder_gene`列と`placeholder_gene_count`で明示する。`mRNA`以外のtranscript型（`transcript`/`*RNA`/`pseudogenic_transcript`等）と`pseudogene`/`ncRNA_gene`も認識し、`gene_biotype`を記録してCDSの無いpseudogene/非コードは`noncoding_or_pseudogene`とし、黙ってprotein化しない。`id_map.tsv`にgene/transcript/protein IDの対応を出し、manifestの`gff_dialect_qc`にduplicate ID・orphan CDS・親未解決率・placeholder数を残す。`--gff-qc-report`でGenomeTools/AGAT等のvalidator結果をmanifestへ添付できる。

eggNOG-mapperを導入してデータベースを配置済みなら、eggnog-data-dirを渡すと、候補proteinだけをqueryにしたGO/orthology annotation commandを同じscriptへ追加する。eggNOG DB準備は重いが、各QTLで再計算するのは候補proteinだけでよい。

## APIによる文献リサーチ（AIなしで標準実行）

Europe PMCは認証不要のread-only APIとして使う。candidate-qtl literature は検索結果をJSONへ固定し、URLと取得日時を残す。

~~~powershell
candidate-qtl literature --query "AT1G01010 Arabidopsis drought tolerance" --out C:\private\qtl\literature\AT1G01010.json --limit 10 --json
~~~

researchは候補ごとにEurope PMCを検索し、research_summary.tsvとliterature_evidence.tsvを作る。AIは不要である。候補順位へ使う主経路は、確認した文献内容を`templates/evidence.tsv`形式へ転記し、`candidate-pipeline`で他の証拠と統合する方法である。

~~~powershell
candidate-qtl research --candidates C:\private\qtl\extract\candidates.tsv --out-dir C:\private\qtl\research --species "Pisum sativum" --trait "target trait" --top-n 8 --json
candidate-qtl candidate-pipeline --candidates C:\private\qtl\extract\candidates.tsv --blast C:\private\qtl\annotation\blastp_hits_1.tsv --manual-evidence C:\private\qtl\manual_evidence.tsv --out-dir C:\private\qtl\candidate_pipeline --json
~~~

論文の件数・被引用数は採点しない。モデル植物の遺伝学的証拠、ortholog経由の遺伝学的証拠、既知QTGを別の`evidence_code`として記録し、PMID、DOI、またはstable accessionを必須にする。

research-packetは任意のAIに渡す独立promptも必要な場合だけ使う。APIのJSON・TSVはAIを使わずにそのままレビューできる。

## AIへのプロンプト運用

生成された ai_research_prompts.md を任意のAIへ候補単位で渡す。回答には必ずDOI、PMID、またはstable database URLを要求する。受け入れるのは、直接証拠、ortholog証拠、反証、未解決事項が分離された回答だけである。

AI出力をそのままrankの入力にしない。候補ledgerへsource IDを転記し、親間配列、marker、F2/F2:3、発現、機能検証で支持または除外する。
