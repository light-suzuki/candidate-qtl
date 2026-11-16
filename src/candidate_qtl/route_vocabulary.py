"""Neutral shared phenotype-route vocabulary (production API)."""
from __future__ import annotations
import re
ROUTE_COLUMNS = ["phenotype_token", "pathway_family", "annotation_keywords"]
_STOPWORDS = frozenset("the a an and or of to in on with without at is are be increased decreased reduced increase decrease reduction higher lower more less altered abnormal plant plants phenotype trait leaf leaves stem seed seeds size number shape color colour".split())
_TOKEN = re.compile(r"[a-z][a-z0-9]{2,}")
def phenotype_tokens(text: str) -> list[str]:
    return list(dict.fromkeys(t for t in _TOKEN.findall(str(text).casefold()) if t not in _STOPWORDS))
def derive_route(tokens: list[str], route: object) -> tuple[set[str], set[str]]:
    pathways: set[str] = set(); keywords: set[str] = set()
    for row in route.to_dict("records"):
        if (token := str(row["phenotype_token"]).strip().casefold()) and token in tokens:
            pathways.add(str(row["pathway_family"]).strip())
            keywords.update(k.strip() for k in re.split(r"[|;,]+", str(row["annotation_keywords"])) if k.strip())
    return pathways, keywords
