# Candidate-first統合基盤

関連Issue: #120。`auto-triage`は既存の`02_...`から`08_...`を維持しながら、ID、provider状態、coverage、由来を共通artifactへ正規化する。

## Canonical identity

主キーは`assembly_id:native_gene_id`である。同じgene ID文字列でもassemblyが違えば別candidateとして扱う。

`entity_aliases.tsv`は次を同じcanonical candidateへ接続する。

- native gene
- transcript
- biological protein ID
- candidate protein FASTAのsequence record ID
- reference database record／ortholog／paralog

`namespace`、`assembly_id`、`annotation_release`、`resource_version`を各行に持つ。one-to-many、別providerに同名IDがある場合、top1/top2 near-tieを縮約しない。`extract`の`id_map.tsv`では`protein_id`と`sequence_record_id`を分離する。

## Provider contractとcoverage

| artifact | 内容 |
|---|---|
| `provider_runs.tsv` | provider/resource version、input/output fingerprint、run状態、件数、失敗 |
| `candidate_provider_coverage.tsv` | candidate×provider単位の実行・hit・失敗状態 |
| `workflow_state.json` | step signature、依存、output fingerprint、resume状態 |
| `resource_manifest.json` | native/reference/input resourceのversionとfingerprint |

coverageは空欄ではなく、`not_requested`、`not_run`、`complete`、`partial`、`no_hit`、`unmapped`、`ambiguous`、`provider_missing`、`input_invalid`、`execution_failed`、`parse_failed`等で表す。BLASTが成功してannotation tableだけが壊れた場合、mappingを捨てずannotation providerだけをfailedにする。

`--resume`は出力ファイルの存在だけでは再利用しない。candidate protein、reference proteome、runner設定、provider versionのsignatureと、保存済みoutput fingerprintが一致した場合だけ再利用する。resource内容が変われば同名・同サイズでも再計算する。

## Observationとlineage

`normalized_observations.tsv`では以下を分離する。

- `mapping_robustness`: similarity、RBH、one-to-many、paralog ambiguity
- `biological_annotation`: description、GO、mutant phenotype、pathway、expression
- `locus_specific_evidence`: parental variant等
- `derived_result`: trait match、screening evidence、rank

mappingと、それから取得したannotation・trait matchを独立証拠として重複加点しない。`evidence_lineage.tsv`の`derived_from` edgeで、raw mapping → reference annotation → trait match → screening evidence → rankを追跡できる。異なるreference namespaceで同じ外部ID文字列が使われてもlineageを混同しない。

## 実行モード

- `local`: local file／local executableだけで主要出力まで完走する。
- `llm-assisted`: localの決定論的結果を維持し、形質検索仮説だけを追加する。
- online resourceを使うreference-first経路は明示的なexternal API許可を必要とする。

LLMはcanonical ID、orthology、score、固定20件を決定しない。LLM responseのresumeもprompt、provider、model、endpointを含むrequest fingerprintが一致した場合だけ許可する。

## 科学的境界

provider一致はmapping robustnessであり、独立した生物学的証拠ではない。自動annotationと自動rankは`unreviewed`のscreening resultで、原因確率ではない。provider failureやno-hitのcandidateも表から削除しない。
