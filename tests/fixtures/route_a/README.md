# route_a fixture

`tests/fixtures/route_a/` は blast-panel freeze 層
(`candidate_qtl.blast_panel`) と、その後の blast-panel-search E2E
(別担当: `candidate_qtl.blast_evidence`) が共有する固定値フィクスチャ。
配列はすべて決定的に生成され、**このファイルの方法で再現可能**。

## 配列一覧 (120 aa、すべて決定論的)

| ID | 配列 | 由来 |
|----|------|------|
| `g1` | `S` | seed `20260812` の 120 aa |
| `g2` | `P1` | `S` の 24 箇所変異 (identity 80%) |
| `g3` | `T` | `S[:48]` + seed `20260813` の無関係 72 aa |
| `g4` | `P2` | `S` の 36 箇所変異 (identity 70%) |
| `g5` | `D` | seed `20260813` の 120 aa (ネガコン、`S` と 4.2% identical) |
| `g6` | `PA_1` | seed `20260814` の 120 aa (`PA` と同一) |
| `g7` | `PA_2` | 同上 (`PA` と同一、ID 別) |

## 生成方法 (再現手順)

```python
import random
AA = "ACDEFGHIKLMNPQRSTVWY"

def runs(seq):  # 最長ホモポリマーランを返す
    out = []
    for ch in seq:
        if out and out[-1][0] == ch:
            out[-1][1] += 1
        else:
            out.append([ch, 1])
    return out

def ok(seq):  # 制約: ラン<=5 かつ 単一AA率<=0.2
    return (max(len(r) for r in runs(seq)) <= 5
            and max(seq.count(c) for c in AA) / len(seq) <= 0.2)

def gen(rng, length=120):
    while True:
        seq = "".join(rng.choice(AA) for _ in range(length))
        if ok(seq):
            return seq

def mutate(base, n, rng):
    seq = list(base)
    for pos in rng.sample(range(len(base)), n):
        seq[pos] = rng.choice([c for c in AA if c != seq[pos]])
    return "".join(seq)

rng1 = random.Random(20260812)
S  = gen(rng1)
P1 = mutate(S, 24, rng1)
P2 = mutate(S, 36, rng1)
rng2 = random.Random(20260813)
while True:
    T = S[:48] + "".join(rng2.choice(AA) for _ in range(72))
    if ok(T):
        break
D  = gen(rng2)          # seed 20260813 の先頭消費が入るため、T 生成より後に呼ぶ
PA = gen(random.Random(20260814))
```

- `S` の 24 変異位置 (0始点): `[5, 15, 28, 35, 37, 44, 47, 49, 50, 58, 61, 69, 84, 85, 87, 91, 93, 96, 100, 102, 103, 106, 109, 118]`
- `S` の 36 変異位置 (0始点): `[2, 3, 8, 11, 13, 15, 16, 17, 18, 29, 30, 35, 36, 38, 39, 41, 46, 49, 59, 61, 63, 69, 73, 76, 78, 81, 86, 88, 92, 94, 96, 97, 99, 101, 115, 118]`

## 配列文字列 (固定値)

- `S`   = `CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVNKQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET`
- `P1`  = `CDLYQDLICSGRWNGHYHCLKFRVGRYVYSYIMRQQFKMFLEVAYCGRRVYNRMNHMNRQDITTEQPNHGNKQKNLGSRNSEYNAMSFYKGDLVRCPLTIMEISFSQSDTPQCWSCIFNT`
- `P2`  = `CDTWQWLISSGKWMGNGWYLKFRVGRYVDKFIMRQTCNHELVVAWCRNRFDNRMNHMNNCDTTCEQPNHCNKQCNLISVNSPYNGGYQEKGNDRDCFSTSYMEPFSWSDAPQCWSMIFPT`
- `T`   = `CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNGKVRCFHEWQVGLEVHENPICCLIRFQLRNPFNSFRCWNLDIKNVNAFNLFTALLEGSWDHRAENPQDAWKI`
- `D`   = `ICHVYKKINVSTGPLVSTSDANAQGDYCCSHRKHDASPDIKCWRGDSGAIARTGQTVFERWNMACYLWRWDDHHIKWQESDWKYQTEDGLTAQLNRLPFYQNGPAPARFTSSMCNKTEIC`
- `PA`  = `MYTCNTWCWQTEEFKCHFAIGTHRPYTYCTNKWFDNQDNCCEEVEEAAPRMLMKADWYYRIDDDAPRQVSDLGPMAKINYNRGSIDRINHLCVWDLHQHYYWWCPFYKAFREYFAKEKQR`

## 配列プロパティ

| 配列 | length | 最長ラン | 最多AA率 | 検証用対比 |
|------|-------:|--------:|--------:|------------|
| `S` | 120 | 2 | 0.108 (N) | 全制約クリア (verified になりうる) |
| `P1` | 120 | 2 | 0.083 (R) | `S` と 80.0% identical |
| `P2` | 120 | 2 | 0.108 (N) | `S` と 70.0% identical |
| `T` | 120 | 2 | 0.092 (L) | 先頭 48 aa が `S[:48]` と 100% identical |
| `D` | 120 | 2 | 0.075 (S) | `S` と 4.2% identical (ネガコン) |
| `PA` | 120 | 2 | 0.083 (Y) | `PA_1` / `PA_2` と 100% identical |

## パネル (panel_ok.tsv) との対応

| panel_id | query_sequence | protein_accession | 期待される E2E の blast 挙動 |
|----------|----------------|-------------------|------------------------------|
| `bait1` | `S` | `Q9TES1` | `g1` 100% / `g3` 先頭48aa 100% / `g2` 80% / `g4` 70% に hit。`g5`(`D`) は no-hit ネガコン |
| `bait2` | `PA` | `P9TES2` | `g6` / `g7` に 100% hit |
| `bait3` | `S[:48]` | `Q9TES3` | 短い検証済み bait。`g1` / `g3` に 48 aa 完全一致 |
| `bait4` | 20 aa | `Q9TES4` | `short_query` で review_required → BLAST クエリにならない (miss にならない) |

> **注**: planner 案の accession `Q0TEST1` は UniProt 正準形
> (`[OPQ][0-9][A-Z0-9]{3}[0-9][A-Z0-9]?`、6〜7 文字で 5 文字目が数字) を満たさない
> (5 文字目が `T`)。正準形を緩めず fixture 側を `Q9TES1` 等に修正した。
> `Q0TEST1` / `P0FAKE` は形式エラー (`bad_accession_format`) のテスト例としてのみ使用する。

各rowは `verification_status=verified`、同一の `verification_accession`、
query sequence と一致する `verification_sequence_sha256`、検証recordのURL・取得時刻も
持つ。これらの検証証跡が欠ける／矛盾するrowは `review_required` となり、queryに昇格しない。

## freeze の期待値 (panel_ok.tsv)

- 全列 `rank_contribution=none`、列名に forbidden token なし。
- freeze 結果: `verified_count=3`、`review_required_count=1` (bait4)、
  `panel_status=verified` (min_panel_size=2)。
- 検証用 SHA-256 (upper-case のまま):
  - `S` = `68ed37701d75a3c3bd0e886f908a6c5ee54bdb273bb4100fa6d6f1e32da2c995`
  - `PA` = `37d634640abf6ee7978f6bdd1d14a432bb8ade88f7a79c8ce35b63f78f95b51e`
  - `S[:48]` = `185363c5f1fd6a514583d27a71876dd61b4f81a54e96f8e72452e5fc9e9511a7`
  - `bait4` (20 aa) = `94a84b2a599e4aa00c3d32c3d1b0d3e4d693346a2043354ee0b001f1c45a5253`

`interval_proteins.fa` は blast-panel-search が freeze 後に初めて開く
"interval protein" であり、freeze 自体 (`freeze_blast_panel`) はこれを
**開かない**ことが仕様 (`interval_accessed=false`)。
