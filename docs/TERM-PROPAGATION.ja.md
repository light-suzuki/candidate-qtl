# リファレンスゲノム由来のtrait検索語伝播（決定論・区間限定・LLM非依存）

`propagate-terms` は、QTL区間候補ごとに、その候補の**オーソログが参照ゲノムで持つ機能注釈（GO/記述語）**を引き、検索語として候補へ伝播する。ベンチマークで示された「出来はtrait検索語の質に支配される」に対する、**LLMを使わない機械的な語の底上げ**。

素のLLMやhand-curated per-trait表に頼らず、原則リファレンスゲノム（正本 + あればモデル植物・近縁種）だけから語を導出する。区間候補とそのオーソログに限定するため軽量で、16GB級で回る。

## 入力

- `--candidates`: 区間候補 `candidates.tsv`（`gene_id`）。
- `--orthology`（任意）: `candidate_gene_id` → `reference_gene_id` の対応表。任意列 `relationship` / `source_genome` / `reference_role` / `reference_species` / `confidence`。role は #71 の panel role（native/close_anchor/functional_reference/outgroup）に対応。テンプレート [`templates/term_propagation_orthology.tsv`](../templates/term_propagation_orthology.tsv)。
- `--reference-annotation`（任意）: `reference_gene_id` → `term`。任意列 `ontology_id` / `ontology_name` / `source`。テンプレート [`templates/term_propagation_reference_annotation.tsv`](../templates/term_propagation_reference_annotation.tsv)。
- `--trait-terms`（任意）: 解析者のtrait語（1行1語 or `term` 列）。与えると各伝播語に `trait_relevance` を付ける。

## 実行

~~~powershell
candidate-qtl propagate-terms --candidates candidates.tsv `
  --orthology term_propagation_orthology.tsv `
  --reference-annotation term_propagation_reference_annotation.tsv `
  --trait-terms traits.txt --out-dir terms --json
~~~

## 出力

- `propagated_terms.tsv`: (候補, 語) 単位の全伝播レコード。provenance（`reference_gene_id` / `source_genome` / `reference_role` / `ontology_id` / `relationship` / `term_source`）、`term_confidence`、`trait_relevance`、`review_status=pending`、`match_type=ortholog_propagated`。
- `search_term_expansion.tsv`: 候補ごとに重複排除した検索語彙（source genome/role/reference gene集約、`max_confidence`）。
- `term_propagation_manifest.json`: term coverage、オーソログはあるが機能注釈が無い候補、オーソログ無し候補、source genome/role、skip理由、fingerprint、区間限定・再現性の記載。

## 設計上の保証

- **native ゲノムのみでも完走**: orthology / reference-annotation が無ければ伝播は空になり、`skip_reasons` と `native_only=true` を machine-readable に残して exit 0。
- **区間限定**: 候補集合に無いオーソログ行は無視。全ゲノム再計算・全遺伝子treeは行わない。
- **provenance と review**: 全語に source を保持し、`review_status` は常に `pending`（人承認まで昇格しない）。このレーンは **LLMを使わない**。
- **未評価と低関連を分離**: オーソログ無し（未評価）、オーソログ有り機能注釈無し、trait非一致（`unrelated`）を別々に記録。
- **再現性**: 同じ入力から `propagated_terms.tsv` / `search_term_expansion.tsv` が byte-stable。

## このレーンで解けないこと（明示）

伝播語は provenance 付きの検索語彙であって因果証拠ではない。語を良くしても、真の遺伝子の注釈欠損（#73 annotation rescue）や、支持獲得後の同点を割る識別証拠（position/expression）は別の層で扱う。
