"""Lossless display transport for an already frozen Route B reading view.

Empty local cells and integrity strings remain in a hash-bound restore ledger.
No biological row, nonempty value, unit, gene/copy or source is selected here.
"""

from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path

CONTRACT = "agentic-route-b-compact-display-1.0.0"


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read(root, name, expected):
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("artifact escapes root")
    data = path.read_bytes()
    if digest(data) != expected:
        raise ValueError("artifact hash mismatch")
    return data


def restore(root, shard):
    root = Path(root)
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest["contract"] != CONTRACT:
        raise ValueError("unsupported compact contract")
    entry = next((x for x in manifest["shards"] if x["file"] == shard), None)
    if entry is None:
        raise ValueError("unknown compact shard")
    display = json.loads(read(root, shard, entry["sha256"]))
    ledger = json.loads(read(root, entry["ledger"], entry["ledger_sha256"]))
    original = Path(manifest["original_root"])
    read(original, "prepare_manifest.json", manifest["original_manifest_sha256"])
    materials = {x["short_id"]: x for x in manifest["materials"]}
    for short in display["materials"]:
        record = materials[short]
        b = read(root, record["file"], record["sha256"])
        if b != read(original, record["original_file"], record["sha256"]):
            raise ValueError("material changed")
    prompts = copy.deepcopy(display["prompts"])
    if len(prompts) != len(ledger["prompt_deltas"]):
        raise ValueError("prompt roster changed")
    for prompt, delta in zip(prompts, ledger["prompt_deltas"]):
        if prompt["prompt"] != delta["prompt_identity"]:
            raise ValueError("prompt identity changed")
        prompt["prompt_sha256"] = delta["prompt_sha256"]
        for key in ("system_material_id", "user_prefix_material_id"):
            if key in prompt:
                prompt[key] = materials[prompt[key]]["material_id"]
        for field, rows in prompt.get("unit_tables", {}).items():
            changes = delta["tables"][field]
            if len(rows) != len(changes):
                raise ValueError("table roster changed")
            for row, change in zip(rows, changes):
                row["material_id"] = materials[row["material_id"]]["material_id"]
                row["binding"].update(change["empty_cells"])
                row["reference_handle"] = change["reference_handle"]
    restored = {**ledger["assignment_envelope"], "prompts": prompts}
    if encoded(restored) != read(original, entry["original_file"], entry["original_sha256"]):
        raise ValueError("compact view does not restore full original assignment")
    return restored


def prepare(reading, out):
    original = Path(reading).resolve()
    out = Path(out)
    if out.exists():
        raise ValueError("output already exists")
    source = (original / "prepare_manifest.json").read_bytes()
    m = json.loads(source)
    if not m["reading_view"]["exact_roundtrip_verified"]:
        raise ValueError("original roundtrip not verified")
    assignments = []
    material_index = {}
    for s in m["reading_view"]["shards"]:
        a = json.loads(read(original, s["file"], s["sha256"]))
        assignments.append((s, a))
        for row in a["materials"]:
            data = read(original, row["file"], row["sha256"])
            if len(data) != row["bytes"]:
                raise ValueError("material size mismatch")
            material_index.setdefault(row["material_id"], (row, data))
    out.mkdir(parents=True)
    materials = []
    short_ids = {}
    for n, (mid, (row, data)) in enumerate(material_index.items()):
        short = f"m{n:04d}"
        short_ids[mid] = short
        name = f"materials/{short}.json"
        (out / name).parent.mkdir(exist_ok=True)
        (out / name).write_bytes(data)
        materials.append({"short_id": short, "material_id": mid, "file": name, "original_file": row["file"], "sha256": row["sha256"], "bytes": row["bytes"]})
    shards = []
    for n, (s, a) in enumerate(assignments):
        views = copy.deepcopy(a["prompts"])
        deltas = []
        for v in views:
            delta = {"prompt_identity": copy.deepcopy(v["prompt"]), "prompt_sha256": v.pop("prompt_sha256"), "tables": {}}
            for key in ("system_material_id", "user_prefix_material_id"):
                if key in v:
                    v[key] = short_ids[v[key]]
            for field, rows in v.get("unit_tables", {}).items():
                changes = []
                for row in rows:
                    empty = {k: x for k, x in row["binding"].items() if x == "" or x is None}
                    for key in empty:
                        del row["binding"][key]
                    changes.append({"empty_cells": empty, "reference_handle": row.pop("reference_handle")})
                    row["material_id"] = short_ids[row["material_id"]]
                delta["tables"][field] = changes
            deltas.append(delta)
        display = {
            "contract": CONTRACT,
            "case_id": a["case_id"],
            "instructions": "Read ALL assigned prompts and their linked materials. Short IDs are transport aliases only; gene, copy, prompt and unit IDs are unchanged. Each material file is complete and may be read once then reused. Omitted local cells are empty string/null, preserved exactly in the restore ledger. No nonempty scientific value or row is omitted. Follow linked system and user-prefix instructions. Answer every assigned prompt_id separately.",
            "materials": [short_ids[x["material_id"]] for x in a["materials"]],
            "prompts": views,
        }
        ledger = {"assignment_envelope": {k: v for k, v in a.items() if k != "prompts"}, "prompt_deltas": deltas}
        name = f"reading/b{n:04d}.json"
        ln = f"ledgers/b{n:04d}.json"
        for file, value in ((name, display), (ln, ledger)):
            (out / file).parent.mkdir(exist_ok=True)
            (out / file).write_bytes(encoded(value))
        shards.append(
            {
                "file": name,
                "sha256": digest(encoded(display)),
                "ledger": ln,
                "ledger_sha256": digest(encoded(ledger)),
                "original_file": s["file"],
                "original_sha256": s["sha256"],
                "prompt_ids": s["prompt_ids"],
                "unit_ids": s["unit_ids"],
                "display_bytes": len(encoded(display)),
                "material_bytes": s["material_bytes"],
            }
        )
    manifest = {
        "contract": CONTRACT,
        "original_root": str(original),
        "original_manifest_sha256": digest(source),
        "case_id": m["case_id"],
        "prompt_count": m["prompt_count"],
        "materials": materials,
        "shards": shards,
        "exact_roundtrip_verified": False,
    }
    (out / "manifest.json").write_bytes(encoded(manifest))
    for s in shards:
        restore(out, s["file"])
    manifest["exact_roundtrip_verified"] = True
    manifest["original_total_reading_bytes"] = sum(x["total_reading_bytes"] for x in m["reading_view"]["shards"])
    manifest["compact_total_reading_bytes"] = sum(x["display_bytes"] + x["material_bytes"] for x in shards)
    manifest["compact_one_read_material_bytes"] = sum(x["display_bytes"] for x in shards) + sum(x["bytes"] for x in materials)
    (out / "manifest.json").write_bytes(encoded(manifest))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reading", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    m = prepare(args.reading, args.out)
    print(json.dumps({k: m[k] for k in ("case_id", "prompt_count", "exact_roundtrip_verified", "original_total_reading_bytes", "compact_total_reading_bytes", "compact_one_read_material_bytes")}))
