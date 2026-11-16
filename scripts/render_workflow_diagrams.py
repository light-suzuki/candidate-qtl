#!/usr/bin/env python3
"""Render the shared workflow illustrations into light/dark SVG assets.

The diagram source of truth is this script: edit the nodes/arrows here and re-render rather
than hand-editing the SVGs.

    python scripts/render_workflow_diagrams.py
    python scripts/render_workflow_diagrams.py --check

The diagrams intentionally separate two different products of the current
implementation:

* ``case-run`` creates the deterministic rank and the original frozen
  ``stage_a/`` shortlist independently of A/B/C discovery artifacts.
* A/B/C family/copy artifacts are non-scoring inputs to optional integrated
  narrowing.  When the same ``triple-route`` run is also given ``--case`` and
  ``--resources`` and the bridge identity checks pass, ``case-run`` may create
  the separate ``stage_a_narrowed/`` Wet-priority shortlist.

This distinction is important: family/copy discovery never rewrites the native
``dry_score`` / ``dry_rank`` / ``screening_rank``.
"""

from __future__ import annotations

import argparse
import sys
from html import escape
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ASSETS = REPO / "docs" / "assets"

FONT = (
    "-apple-system, BlinkMacSystemFont, 'Hiragino Sans', 'Noto Sans JP', "
    "'Yu Gothic UI', 'Yu Gothic', Meiryo, 'Segoe UI', sans-serif"
)

LIGHT = {
    "name": "light",
    "bg": "#ffffff",
    "panel": "#f6f8fa",
    "panel_border": "#d8dee4",
    "box": "#ffffff",
    "border": "#d0d7de",
    "text": "#1f2328",
    "muted": "#59636e",
    "edge": "#8c959f",
    "io": ("#dafbe1", "#4ac26b", "#0a5c2b"),
    "step": ("#ffffff", "#d0d7de", "#1f2328"),
    "atlas": ("#ddf4ff", "#54aeff", "#0a3069"),
    "route_a": ("#faf0ff", "#c297ff", "#512a97"),
    "route_b": ("#fff4ec", "#f0a875", "#7d2f00"),
    "route_c": ("#ddf8fa", "#76d4dc", "#0f5257"),
    "narrow": ("#fff8c5", "#e3b341", "#6b4900"),
    "warn": ("#fff1f0", "#ff8182", "#82071e"),
}

DARK = {
    "name": "dark",
    "bg": "#0d1117",
    "panel": "#151b23",
    "panel_border": "#2f3742",
    "box": "#0d1117",
    "border": "#3d444d",
    "text": "#e6edf3",
    "muted": "#9198a1",
    "edge": "#6e7681",
    "io": ("#0c2d19", "#2ea043", "#3fb950"),
    "step": ("#161b22", "#3d444d", "#e6edf3"),
    "atlas": ("#0c2d6b", "#316dca", "#79c0ff"),
    "route_a": ("#2b1a4d", "#8957e5", "#d2a8ff"),
    "route_b": ("#3d1e00", "#bd561d", "#ffa657"),
    "route_c": ("#0a3b3e", "#207d84", "#56d4dd"),
    "narrow": ("#3b2300", "#9e6a03", "#e3b341"),
    "warn": ("#3c1618", "#da3633", "#ff7b72"),
}


def _char_em(ch: str) -> float:
    o = ord(ch)
    if o >= 0x2E80 or 0x2190 <= o <= 0x27BF:
        return 1.0
    if ch in "iljt.,:;'|!":
        return 0.30
    if ch.isupper():
        return 0.62
    return 0.52


def text_width(text: str, size: float) -> float:
    return sum(_char_em(ch) for ch in text) * size


class Canvas:
    def __init__(self, width: int, height: int, palette: dict) -> None:
        self.w = width
        self.h = height
        self.p = palette
        self.parts: list[str] = []
        self.warnings: list[str] = []

    def box(self, x, y, w, h, role="step", *, dashed=False, rx=10, stroke_width=1.5):
        fill, stroke, _ = self.p[role]
        dash = ' stroke-dasharray="7 5"' if dashed else ""
        self.parts.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"{dash}/>'
        )

    def panel(self, x, y, w, h, *, dashed=False, rx=14):
        dash = ' stroke-dasharray="7 5"' if dashed else ""
        self.parts.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
            f'fill="{self.p["panel"]}" stroke="{self.p["panel_border"]}" '
            f'stroke-width="1.5"{dash}/>'
        )

    def text(self, x, y, text, *, size=15, color=None, anchor="middle", weight="400", max_width=None, where=""):
        color = color or self.p["text"]
        if max_width is not None and text_width(text, size) > max_width:
            self.warnings.append(f"{where or text}: overflow: {text!r}")
        self.parts.append(
            f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" fill="{color}" text-anchor="{anchor}">{escape(text)}</text>'
        )

    def step(self, x, y, w, h, title, sub=None, sub2=None, *, role="step", dashed=False):
        self.box(x, y, w, h, role=role, dashed=dashed)
        _, _, fg = self.p[role]
        cx = x + w / 2
        avail = w - 24
        if sub2 is not None:
            self.text(cx, y + h / 2 - 13, title, size=16, weight="700", color=fg, max_width=avail, where=title)
            self.text(cx, y + h / 2 + 8, sub or "", size=12.5, color=self.p["muted"], max_width=avail, where=title)
            self.text(cx, y + h / 2 + 27, sub2, size=12.5, color=self.p["muted"], max_width=avail, where=title)
        elif sub is not None:
            self.text(cx, y + h / 2 - 4, title, size=16, weight="700", color=fg, max_width=avail, where=title)
            self.text(cx, y + h / 2 + 18, sub, size=12.5, color=self.p["muted"], max_width=avail, where=title)
        else:
            self.text(cx, y + h / 2 + 6, title, size=16, weight="700", color=fg, max_width=avail, where=title)

    def card(self, x, y, w, h, title, lines, role):
        fill, stroke, fg = self.p[role]
        self.parts.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" '
            f'fill="{self.p["box"]}" stroke="{stroke}" stroke-width="1.5"/>'
        )
        self.parts.append(
            f'<path d="M{x},{y + 12} a12,12 0 0 1 12,-12 h{w - 24} a12,12 0 0 1 12,12 '
            f'v30 h{-w} z" fill="{fill}" stroke="none"/>'
        )
        self.parts.append(
            f'<line x1="{x}" y1="{y + 42}" x2="{x + w}" y2="{y + 42}" stroke="{stroke}" stroke-width="1.5"/>'
        )
        self.text(x + w/2, y + 27, title, size=15, weight="700", color=fg, max_width=w-24, where=title)
        for i, line in enumerate(lines):
            self.text(x + 16, y + 68 + i*24, line, size=12.8, anchor="start", max_width=w-32, where=title)

    def tag(self, x, y, text, role="atlas", size=12):
        fill, stroke, fg = self.p[role]
        w = text_width(text, size) + 18
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w:.0f}" height="22" rx="11" fill="{fill}" stroke="{stroke}"/>')
        self.text(x + w/2, y + 15, text, size=size, color=fg, weight="700")

    def arrow(self, points, *, dashed=False, label=None, label_at=None, rotate_label=False):
        d = " ".join(f"{'M' if i == 0 else 'L'}{x},{y}" for i, (x, y) in enumerate(points))
        dash = ' stroke-dasharray="7 5"' if dashed else ""
        self.parts.append(
            f'<path d="{d}" fill="none" stroke="{self.p["edge"]}" stroke-width="1.7" '
            f'stroke-linejoin="round" stroke-linecap="round"{dash} marker-end="url(#ah)"/>'
        )
        if label and label_at:
            lx, ly = label_at
            transform = f' transform="rotate(-90 {lx} {ly})"' if rotate_label else ""
            self.parts.append(
                f'<text x="{lx}" y="{ly}" font-family="{FONT}" font-size="11.5" '
                f'fill="{self.p["muted"]}" text-anchor="middle"{transform}>{escape(label)}</text>'
            )

    def render(self) -> str:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" width="{self.w}" height="{self.h}" role="img">\n'
            f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,1 L9,5 L0,9 z" fill="{self.p["edge"]}"/></marker></defs>\n'
            f'<rect width="{self.w}" height="{self.h}" fill="{self.p["bg"]}"/>\n'
            + "\n".join(self.parts) + "\n</svg>\n"
        )


def phase1(palette: dict) -> Canvas:
    """Discovery and deterministic-rank lanes, kept visually separate."""
    c = Canvas(1120, 1160, palette)
    mid = 560

    c.step(350, 26, 420, 62, "FASTA + GFF3 + QTL区間 + phenotype", role="io")
    c.step(350, 122, 420, 70, "mapping adequacy gate", "marker liftover / QTLift・split / reverse は review")
    c.step(350, 226, 420, 70, "interval extraction", "transcript-aware candidate universe")
    c.step(820, 226, 250, 70, "annotation rescue", "miniprot・rankは変更しない", dashed=True)
    c.step(50, 334, 1020, 78, "interval atlas", "全遺伝子 → family / domain / orthogroup", "B0 inventory は phenotype-free", role="atlas")

    c.card(50, 456, 320, 154, "Route A ｜ phenotype-first", [
        "phenotype → 競合する機構", "source / accession / sequence 検証", "verified bait → interval projection / BLASTP"
    ], "route_a")
    c.card(400, 456, 320, 154, "Route B ｜ interval-first", [
        "atlas → family / orthology", "source-backed function は入力依存", "overlay / reverse-literature"
    ], "route_b")
    c.card(750, 456, 320, 154, "Route C ｜ retrieval（実験的）", [
        "phenotype-only retrieval → freeze", "freeze後だけ interval match", "native rank は変更しない"
    ], "route_c")

    c.step(235, 654, 650, 84, "triple-route-compare", "route disagreement / shared source / provider不足を診断", "convergence は score ではない")

    # Two separate downstream lanes.  There is deliberately no arrow from the
    # family/copy lane into deterministic ranking / stage_a.
    c.panel(45, 780, 500, 330)
    c.tag(62, 796, "A/B/C family・copy（non-scoring）", role="route_a")
    c.step(80, 836, 430, 68, "family-candidate-bridge", "(gene_id, family_id) + mapping / dependency state")
    c.step(80, 936, 430, 68, "family-to-copy", "preferred / ambiguous / unresolved を保持")
    c.step(80, 1036, 430, 54, "family/copy artifacts", "rank_contribution = none", role="narrow")

    c.panel(575, 780, 500, 330)
    c.tag(592, 796, "case-run deterministic lane", role="atlas")
    c.step(610, 836, 430, 68, "canonical evidence → deterministic rank", "A/B/C prose・convergenceで順位を上書きしない")
    c.step(610, 936, 430, 68, "stage_a/ 固定20", "screening_rank基準・元のreview予算", role="io")
    c.step(610, 1036, 430, 54, "Phase 2 review / Wetへ", "stage_a は family/copy と独立")

    c.arrow([(mid, 88), (mid, 120)])
    c.arrow([(mid, 192), (mid, 224)])
    c.arrow([(mid, 296), (mid, 332)])
    c.arrow([(770, 261), (818, 261)], dashed=True)
    c.arrow([(945, 296), (945, 332)], dashed=True)
    c.arrow([(mid, 412), (mid, 454)])
    c.arrow([(370, 58), (22, 58), (22, 532), (48, 532)], label="phenotype", label_at=(36, 298), rotate_label=True)
    c.arrow([(770, 58), (1098, 58), (1098, 532), (1072, 532)], label="phenotype", label_at=(1084, 298), rotate_label=True)
    c.arrow([(210, 610), (210, 630), (360, 630), (360, 652)])
    c.arrow([(560, 610), (560, 652)])
    c.arrow([(910, 610), (910, 630), (760, 630), (760, 652)])
    c.arrow([(430, 738), (430, 834)])
    c.arrow([(295, 904), (295, 934)])
    c.arrow([(295, 1004), (295, 1034)])

    # Deterministic lane source: same case candidate universe + canonical
    # evidence, not the route-comparison artifact. Dashed connector is a context
    # relationship, deliberately not a causal edge from comparison.
    c.arrow([(820, 412), (1054, 412), (1054, 870), (1042, 870)], dashed=True,
            label="case universe / canonical evidence", label_at=(1015, 775))
    c.arrow([(825, 904), (825, 934)])
    c.arrow([(825, 1004), (825, 1034)])

    c.text(mid, 1138, "実線 = 現行の通常データフロー ｜ 破線 = optional / context。family/copy → native rank の矢印は存在しない。",
           size=12.5, color=c.p["muted"], max_width=1030)
    return c


def phase2(palette: dict) -> Canvas:
    """Wet/review flow and the conditional family/copy narrowing handoff."""
    c = Canvas(1120, 950, palette)

    # Inputs carried forward from Phase 1.
    c.step(50, 30, 470, 72, "stage_a/ 固定20", "deterministic screening shortlist・不変", role="io")
    c.step(600, 30, 470, 72, "family/copy artifacts", "non-scoring・triple-route後の別artifact", role="narrow")

    c.panel(35, 138, 1050, 170)
    c.tag(52, 154, "① review / Wet")
    c.step(65, 196, 300, 82, "human / LLM review", "元Stage-Aのscore/集合は上書き不可")
    c.step(410, 196, 300, 82, "Primer / PCR / Sanger", "assay execution")
    c.step(755, 196, 300, 82, "normalized Wet state", "失敗・partial・negativeを区別")
    c.arrow([(285, 102), (285, 194)])
    c.arrow([(365, 237), (408, 237)])
    c.arrow([(710, 237), (753, 237)])

    c.panel(35, 344, 1050, 220)
    c.tag(52, 360, "② interpretation / narrowing", role="narrow")
    c.step(65, 408, 300, 86, "locus-hypothesis ledger", "coding / regulatory / SV / PAV")
    c.step(410, 408, 300, 86, "candidate-elimination", "強い位置・遺伝証拠でのみ除外")
    c.step(755, 408, 300, 86, "integrated-narrowing", "family/copy + remaining → N1..N4", role="narrow")
    c.arrow([(905, 278), (905, 406)])
    c.arrow([(365, 451), (408, 451)])
    c.arrow([(710, 451), (753, 451)])

    # Conditional family/copy handoff: same-run case/resources or registered
    # resources, and identity checks must pass.
    c.arrow([(835, 102), (835, 132), (1035, 132), (1035, 406)], dashed=True,
            label="条件付き: case-runがartifactを検証して消費", label_at=(880, 122))

    c.step(755, 520, 300, 58, "stage_a_narrowed/ Wet優先20", "narrowing評価時だけ生成", role="io", dashed=True)
    c.arrow([(905, 494), (905, 518)], dashed=True)

    c.panel(35, 610, 1050, 210)
    c.tag(52, 626, "③ next experiment")
    c.step(90, 682, 380, 82, "next-experiment", "narrowingが無くても実行可能")
    c.step(650, 656, 380, 58, "候補が強い", "recombinant / F2:3 / independent allele", role="io")
    c.step(650, 736, 380, 58, "弱い・未評価・棄却", "追加mapping / SV / rescue / expression")
    # Base path always reaches next-experiment from remaining candidate state.
    c.arrow([(560, 494), (560, 723), (472, 723)])
    # Narrowing is an optional additional input.
    c.arrow([(905, 578), (905, 600), (510, 600), (510, 723), (472, 723)], dashed=True,
            label="optional narrowing input", label_at=(705, 592))
    c.arrow([(470, 704), (648, 685)])
    c.arrow([(470, 742), (648, 765)])

    c.box(50, 850, 1020, 62, role="warn", dashed=True)
    c.text(560, 875, "現在の統合境界", size=14, weight="700", color=c.p["warn"][2])
    c.text(560, 897,
           "triple-route → case-run は --case + --resources を同一runに渡した時だけ。standalone homolog-linking の自動handoffは未実装（#426）。",
           size=12.2, color=c.p["muted"], max_width=980)
    return c


DIAGRAMS = {"workflow-phase1": phase1, "workflow-phase2": phase2}


def build() -> tuple[dict[Path, str], list[str]]:
    files: dict[Path, str] = {}
    warnings: list[str] = []
    for stem, fn in DIAGRAMS.items():
        for palette in (LIGHT, DARK):
            canvas = fn(palette)
            files[ASSETS / f"{stem}-{palette['name']}.svg"] = canvas.render()
            warnings.extend(f"{stem}-{palette['name']}: {warning}" for warning in canvas.warnings)
    return files, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail when committed assets differ from a fresh render")
    args = parser.parse_args()

    files, warnings = build()
    for warning in warnings:
        print(f"OVERFLOW {warning}", file=sys.stderr)

    if args.check:
        stale = [path for path, text in files.items() if not path.exists() or path.read_bytes() != text.encode("utf-8")]
        for path in stale:
            print(f"STALE {path.relative_to(REPO)}", file=sys.stderr)
        if stale or warnings:
            return 1
        print(f"ok: {len(files)} diagram assets up to date")
        return 0

    ASSETS.mkdir(parents=True, exist_ok=True)
    for path, text in files.items():
        path.write_bytes(text.encode("utf-8"))
        print(f"wrote {path.relative_to(REPO)}")
    return 1 if warnings else 0


if __name__ == "__main__":
    raise SystemExit(main())
