"""Construct synthetic unversioned inputs before test prediction/hash freezes.

These helpers never modify saved runs. Historical adapter tests deliberately
exercise legacy_grounded_review; new bibliography policy requires resolved
metadata and must not inherit their source-free admission expectations.
"""
import json
from pathlib import Path


def unversioned_bundle(built):
    path = Path(built["bundle"])
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.pop("b_review_membership_policy", None)
    path.write_text(json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return built
