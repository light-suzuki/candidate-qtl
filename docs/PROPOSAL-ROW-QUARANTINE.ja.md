# 不正な提案行の隔離 — CORE_FIX #698

2026-09-05。候補選定・生物学的閾値を変えず、一部の不正行による正常行の巻き添え停止を修正する。

- Route A: canonical入力の認識できない出典ID形式を含む仮説は、元行・元位置・理由付きで`mechanism_rejected_rows.json`へ隔離する。出典IDを創作・削除して通過させない。
- Route B: semantic claimのscope/anchorなど行スキーマに違反するものを`semantic_claim_rejected_rows.json`へ隔離する。family IDの補完やscopeの付け替えはしない。
- 正常行は従来の出典・identity・transfer gateへ通す。提案の取り込みはverified supportではない。
- manifestとCLI結果に隔離件数・import状態を出す。全部不正なら`not_evaluated_all_rows_rejected`（A）/`not_evaluated_all_claim_rows_rejected`（B）。Bはunit別にも記録する。
- packetのschema不整合、roster不整合、freeze hash不一致、Aのinterval-blind違反は従来通り全体エラー。隔離行もAの漏洩監査から外さない。

合成入力で正常行単独と不正行混在の出力一致を確認する。Bのresponse fingerprintは元回答全体が異なるため異なってよい。全不正と空の正常回答は区別する。
旧Museベンチの回答・凍結・採点は変更せず、この修正による再計算を新しいLLM性能とは呼ばない。
