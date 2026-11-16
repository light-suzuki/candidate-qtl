# 本番のfamily/copy判別へ既存証拠を渡す

2026-09-18 / CORE_INTEGRATION / #729 の受渡し部分。
未知原因の実QTLでも、既に持っている証拠をcanonical copy段階へ届けるための接続。
新しい評価式・原因遺伝子の断定・証拠の自動取得は追加しない。

## 利用できる入口

既存の `guided-profile --profile triple-route --execute` の通常入力に、必要なものだけ追加する。

| CLI | native family-to-copy引数 |
|---|---|
| `--copy-orthology FILE` | `orthology_arg` |
| `--copy-synteny FILE` | `synteny_arg` |
| `--copy-expression FILE` | `expression_arg` |
| `--copy-position FILE` | `position_arg` |
| `--copy-subfamily-map FILE` | `subfamily_map_arg` |
| `--copy-target FILE` | `target_arg` |

ファイル形式は既存 `family-to-copy` の形式を再利用する。B0用 `--orthology` を
copy判別用の表として勝手に再解釈しない。全6入力を揃える必要はない。
Pythonのcanonical extension / family-semantic wrapperには
`copy_evidence_args={"expression": "...", "synteny": "..."}` を渡せる。
この変更は `case-run` YAMLへの新しい入力規約や、直接の `guided_profile()` APIを改変しない。
CLI handlerが既存canonical extensionを呼ぶ接続を対象とする。

## 証拠と候補の境界

入力はA/Bの仮説生成やphenotype-free B0へ渡さず、supported-member unionからの
canonical family recovery後にだけ既存 `family_to_copy()` へ渡す。
Cあり・B回答なしの既存fallbackにも同じ接続を行うが、Cは入力の必須条件ではない。
familyが無い場合、証拠表から新しいfamilyや支持候補を作らない。

遺伝子数・種・形質名・既知の正解・親配列の有無による分岐は追加しない。
合成ベンチのN=1000を本番の上限にしない。既存の複数copy・不明・矛盾の意味を維持する。
省略した証拠は負の証拠ではない。明示指定した不正なパスは、黙って無視せず入力エラーにする。
表内のprovider failure / not_evaluated / one-to-manyはそのままnative consumerへ渡す。
それをどのように判定するかは既存の証拠規約であり、今回のhelperが再判定しない。

## 出典・版と再開

`b04_canonical_semantic/canonical_copy_inputs.json` に、供給した各軸の絶対パス・全内容SHA-256・
byte数を記録する。元表の種・系統・assembly・provider/version列は変更しない。
外側guided manifest、canonical manifest、family-semantic wrapper manifestにも同じ入力参照を記録する。
生物学的に評価可能な軸・実際の判別内容は既存 `family_copy_manifest.json` と
`family_copy_discrimination.tsv` を参照する。
`passed_to_family_to_copy` は受渡し成功であり、全行の生物学的検証成功を意味しない。

同じ出力先を再開する際は同じ入力を指定する。内容・供給する軸・パスが変わった場合は新しい出力revisionを使う。
以前のWetに入力履歴なしで証拠を後付けしたり、証拠を省略して別条件にしたりする操作は拒否する。
処理中の入力改変もWet引渡し前に検出する。これは汎用resume/追加証拠統合システムではない。
未供給時には新しいpinを作らず、従来のfamily-to-copy呼出しを維持する。

## 確認範囲

局所テストは実argparse登録、CLI handler、family-semantic wrapper、canonical extension、
実ファイルI/Oとhash処理を通す。1,247候補の合成入力で、6軸が消えず到達することも検査する。
ただし、biologyのproducer、source/transfer、family-to-copyの判定本体、Wet producerはmock。
この接続試験を、実データでのcopy識別精度やType3性能の実証とは呼ばない。

```text
python -m unittest discover -s tests -p test_canonical_copy_inputs.py -v
```

#729の完了には、既存の実family-to-copyによる曖昧性/矛盾/欠測fixtureと、clean installのCLI到達確認が残る。
#730の親間変異・SV/PAV・プロモーター接続、#478の参照資源拡充、#642のType3評価は別作業。
