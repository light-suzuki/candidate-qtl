"""Read explicit portions of frozen local material without clipping or inference."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

if __package__:
    from . import agentic_route_b_light as light
else:
    import agentic_route_b_light as light


def _verified(path_arg: str, expected: str) -> tuple[Path, bytes]:
    path = Path(path_arg).resolve()
    raw = path.read_bytes()
    if not expected or hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")
    return path, raw


def _pointer(data: object, pointer: str) -> object:
    if not isinstance(pointer, str) or (pointer and not pointer.startswith("/")):
        raise ValueError("field must be an exact JSON pointer")
    current = data
    for token in pointer.split("/")[1:]:
        if re.search(r"~(?![01])", token):
            raise ValueError("invalid JSON pointer escape")
        token = token.replace("~1", "/").replace("~0", "~")
        try:
            if isinstance(current, list):
                if not token.isdigit() or str(int(token)) != token:
                    raise ValueError("invalid array index")
                current = current[int(token)]
            elif isinstance(current, dict):
                current = current[token]
            else:
                raise ValueError("selector does not address a container")
        except (KeyError, IndexError) as exc:
            raise ValueError(f"unknown selector: {pointer}") from exc
    return current


def _finish(packet: dict, selected: list[dict], max_chars: int) -> dict:
    if max_chars < 1:
        raise ValueError("max_chars must be positive")
    size = len(json.dumps(selected, ensure_ascii=False))
    packet.update(selected_pointers=[item["pointer"] for item in selected],
                  omitted="all unselected ranges/pointers remain not_read",
                  research_state="not_evaluated_by_reader", selected_characters=size)
    if size > max_chars:
        packet.update(state="narrow_selection_required", read_state="not_read",
                      next_action="request narrower explicit pointers/line ranges or explicitly raise max_chars")
    else:
        packet.update(state="selected", read_state="selected_ranges_returned", selections=selected)
    return packet


def source_read(source_arg: str, source_sha256: str, *, pointers: list[str] | None = None,
                line_start: int | None = None, line_end: int | None = None,
                max_chars: int = 12000) -> dict:
    """Read already-saved source JSON pointers or an inclusive UTF-8 line range."""
    path, raw = _verified(source_arg, source_sha256)
    packet = {"source": {"path": str(path), "sha256": source_sha256}, "network": "not_used"}
    if pointers is not None and (line_start is not None or line_end is not None):
        raise ValueError("select JSON pointers or lines, not both")
    if pointers:
        data = json.loads(raw)
        selected = [{"pointer": pointer, "value": _pointer(data, pointer)} for pointer in pointers]
    elif line_start is not None or line_end is not None:
        lines = raw.decode("utf-8").splitlines(keepends=True)
        if (not isinstance(line_start, int) or not isinstance(line_end, int)
                or isinstance(line_start, bool) or isinstance(line_end, bool)
                or line_start < 1 or line_end < line_start or line_end > len(lines)):
            raise ValueError("invalid inclusive line range")
        selected = [{"pointer": f"lines:{line_start}-{line_end}",
                     "value": "".join(lines[line_start - 1:line_end])}]
    else:
        return {**packet, "state": "index_only", "read_state": "not_read",
                "bytes": len(raw), "next_action": "specify exact JSON pointers or inclusive line_start/line_end"}
    return _finish(packet, selected, max_chars)


def material_read(bundle_arg: str, bundle_sha256: str, light_index_arg: str,
                  light_index_sha256: str, *, batch_id: str,
                  unit_ids: list[str] | None = None, gene_ids: list[str] | None = None,
                  fields: list[str] | None = None, max_chars: int = 12000) -> dict:
    """Read only the assigned batch; fields are unit-relative JSON pointers.

    With gene selectors, bare field names address the original gene/card rows.
    Every returned pointer addresses the unchanged full-material sidecar.
    """
    bundle_path, bundle = light._load_bundle(bundle_arg, bundle_sha256)
    index_path, raw = _verified(light_index_arg, light_index_sha256)
    index = json.loads(raw)
    if (index.get("contract") not in {light.CONTRACT, light.V2_CONTRACT}
            or index.get("case_id") != bundle["case_id"]
            or index.get("bundle", {}).get("sha256") != bundle_sha256
            or Path(index.get("bundle", {}).get("path", "")).resolve() != bundle_path):
        raise ValueError("light index does not bind this bundle/case")
    units = {unit["unit_id"]: unit for unit in bundle["units"]}
    batches = index.get("batches")
    if not isinstance(batches, list) or not batches:
        raise ValueError("light index requires batches")
    seen_units, seen_batches = set(), set()
    chosen = None
    for batch in batches:
        bid, ids = batch.get("bid"), batch.get("unit_ids")
        if (not isinstance(bid, str) or not re.fullmatch(r"b[0-9]{2,}", bid)
                or bid in seen_batches or not isinstance(ids, list) or not ids
                or len(set(ids)) != len(ids) or seen_units.intersection(ids)
                or any(uid not in units for uid in ids)):
            raise ValueError("invalid light batch/unit roster")
        seen_batches.add(bid)
        seen_units.update(ids)
        expected = [gene for uid in ids for gene in units[uid]["member_interval_genes"]]
        if batch.get("gene_ids") != expected:
            raise ValueError("light batch gene roster mismatch")
        if bid == batch_id:
            chosen = batch
    if not seen_units <= {unit["unit_id"] for unit in light._requested_units(bundle)}:
        raise ValueError("index names units outside the requested roster")
    if chosen is None:
        raise ValueError("unknown assigned batch")
    directory = index_path.parent / batch_id
    path, material_raw = _verified(str(directory / "material.full.json"), chosen.get("material_full_sha256", ""))
    _, roster_raw = _verified(str(directory / "roster.json"), chosen.get("roster_sha256", ""))
    material = json.loads(material_raw)
    expected_units = [units[uid] for uid in chosen["unit_ids"]]
    if material != {"units": expected_units}:
        raise ValueError("full material differs from frozen bundle")
    expected_roster = {gid: row for unit in expected_units for gid, row in light._roster(unit).items()}
    if json.loads(roster_raw) != expected_roster:
        raise ValueError("roster differs from frozen bundle")
    if unit_ids and any(uid not in chosen["unit_ids"] for uid in unit_ids):
        raise ValueError("unknown/out-of-batch unit selector")
    if gene_ids and any(gid not in chosen["gene_ids"] for gid in gene_ids):
        raise ValueError("unknown/out-of-batch gene selector")
    targets = [(i, unit) for i, unit in enumerate(expected_units)
               if (not unit_ids or unit["unit_id"] in unit_ids)
               and (not gene_ids or set(gene_ids).intersection(unit["member_interval_genes"]))]
    if gene_ids and set(gene_ids) - {g for _, unit in targets for g in unit["member_interval_genes"]}:
        raise ValueError("gene/unit selectors conflict")
    packet = {"case_id": bundle["case_id"], "batch_id": batch_id,
              "source": {"path": str(path), "sha256": chosen["material_full_sha256"]},
              "bundle": {"path": str(bundle_path), "sha256": bundle_sha256},
              "light_index": {"path": str(index_path), "sha256": light_index_sha256}}
    if not (unit_ids or gene_ids or fields):
        return {**packet, "state": "index_only", "read_state": "not_read",
                "units": [{"unit_id": u["unit_id"], "gene_ids": u["member_interval_genes"],
                           "fields": list(u)} for u in expected_units]}
    selected = []
    for i, unit in targets:
        prefix = f"/units/{i}"
        if gene_ids:
            candidates = [(f"/{collection}/{j}", row) for collection in ("genes", "reference_material_cards")
                          for j, row in enumerate(unit.get(collection, []))
                          if row.get("gene_id", row.get("interval_gene_id")) in gene_ids]
            for field in fields or [""]:
                matches = 0
                for pointer, row in candidates:
                    if field.startswith("/"):
                        if field != pointer and not field.startswith(pointer + "/"):
                            continue
                        value, suffix = _pointer(unit, field), field
                    elif field:
                        if field not in row:
                            continue
                        value, suffix = row[field], pointer + "/" + field.replace("~", "~0").replace("/", "~1")
                    else:
                        value, suffix = row, pointer
                    selected.append({"pointer": prefix + suffix, "value": value})
                    matches += 1
                if not matches:
                    raise ValueError(f"unknown/out-of-gene field selector: {field}")
        else:
            for field in fields or [""]:
                pointer = field if not field or field.startswith("/") else "/" + field.replace("~", "~0").replace("/", "~1")
                selected.append({"pointer": prefix + pointer, "value": _pointer(unit, pointer)})
    return _finish(packet, selected, max_chars)
