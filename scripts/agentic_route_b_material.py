"""Lossless, readable Route B transport; never groups or selects biology."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import textwrap
from typing import Any


CONTRACT = "agentic-route-b-shared-material-1.0.0"
DEFAULT_MAX_MATERIAL_BYTES = 262144
DISPLAY_CONTRACT = "agentic-route-b-lossless-display-1.0.0"
DISPLAY_MAX_LINE_CHARS = 1000
_DISPLAY_CHUNKS = "__cqtl_display_string_chunks__"
_REFERENCE_FIELDS = {
    "corpus_version", "developmental_stage", "evidence_codes", "evidence_type",
    "function_transfer_scope", "go_terms", "keywords", "literature_source_ids",
    "lookup_status", "mapping_relationship", "mapping_state", "phenotype_source_ids",
    "protein_family", "protein_name", "provenance_detail", "provenance_tier",
    "retrieved_at", "source_database", "source_release", "source_url",
    "tissue_specificity", "verified_source_ids",
}
_NOTE = (
    "Read each prompt's unit and its linked material files together with the shared "
    "system and user_prefix. A material file holds the complete original fields, "
    "not a summary or verified support. Bindings preserve every gene/reference/source "
    "association and duplicate row. Reuse material already read; answer every assigned "
    "prompt_id/route_b_unit_id separately. No family-level answer or Top-K selection. "
    "An oversized_full_material assignment retains one complete unit and all full "
    "material links; read those files in chunks instead of dropping the unit."
)


def _encoded(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def prompt_digest(prompt: dict[str, Any]) -> str:
    return _digest(_encoded(prompt))


def restore_display(data: bytes, canonical: bytes) -> Any:
    """Decode only explicitly listed display paths, never marker-like source fields."""
    display = json.loads(data)
    if display.get("contract") != DISPLAY_CONTRACT or display.get("source_sha256") != _digest(canonical):
        raise ValueError("display source/contract mismatch")
    value = display["value"]

    def at(path):
        node = value
        for key in path:
            node = node[key]
        return node

    for path in display["chunked_string_paths"]:
        node = at(path)
        if not isinstance(node, dict) or set(node) != {_DISPLAY_CHUNKS}:
            raise ValueError("invalid display string chunks")
        chunks = node[_DISPLAY_CHUNKS]
        if not isinstance(chunks, list) or not all(isinstance(chunk, str) for chunk in chunks):
            raise ValueError("invalid display string chunks")
        restored = "".join(chunks)
        if path:
            at(path[:-1])[path[-1]] = restored
        else:
            value = restored
    # Child paths still use the display keys until their parents are restored.
    for renamed in reversed(display["renamed_keys"]):
        parent = at(renamed["parent_path"])
        key = "".join(renamed["key_chunks"])
        if key in parent:
            raise ValueError("display key collision")
        parent[key] = parent.pop(renamed["display_key"])
    if _encoded(value) != canonical:
        raise ValueError("display differs from canonical material")
    return value


def _display_sidecar(canonical: bytes, source_file: str) -> tuple[bytes, dict[str, Any]]:
    paths, renamed_keys = [], []

    def render(value, **options):
        text = json.dumps(value, ensure_ascii=False, **options)
        # Line readers can treat these Unicode characters as physical breaks.
        for character in ("\u0085", "\u2028", "\u2029"):
            text = text.replace(character, f"\\u{ord(character):04x}")
        return text

    def chunks(text):
        result = []
        for index in range(0, len(text), 200):
            chunk = text[index:index + 200]
            if len(render(chunk)) > 700:
                result.extend([chunk[:100], chunk[100:]])
            else:
                result.append(chunk)
        return result

    def transform(value, path):
        if isinstance(value, str) and len(render(value)) > 512:
            paths.append(path)
            return {_DISPLAY_CHUNKS: chunks(value)}
        if isinstance(value, list):
            return [transform(item, path + [index]) for index, item in enumerate(value)]
        if isinstance(value, dict):
            result = {}
            for key, item in value.items():
                displayed = key
                if len(render(key)) > 256:
                    index = 0
                    displayed = f"__display_key_{index}__"
                    while displayed in value or displayed in result:
                        index += 1
                        displayed = f"__display_key_{index}__"
                    renamed_keys.append({"parent_path": path, "display_key": displayed,
                                         "key_chunks": chunks(key)})
                result[displayed] = transform(item, path + [displayed])
            return result
        return value

    value = transform(json.loads(canonical), [])
    receipt = read_receipt_token(_digest(canonical))
    # The receipt is the last key, so it sits at the end of the displayed file.
    display = {"contract": DISPLAY_CONTRACT, "source_sha256": _digest(canonical),
               "chunked_string_paths": paths, "renamed_keys": renamed_keys, "value": value,
               "read_receipt": receipt}
    text = render(display, indent=2) + "\n"
    lines = text.splitlines()
    maximum = max(map(len, lines), default=0)
    if maximum > DISPLAY_MAX_LINE_CHARS:
        raise ValueError("display exceeds the physical line bound")
    data = text.encode("utf-8")
    restore_display(data, canonical)
    return data, {"file": "display/" + source_file, "sha256": _digest(data), "bytes": len(data),
                  "source_file": source_file, "source_sha256": _digest(canonical),
                  "line_count": len(lines), "max_line_characters": maximum,
                  "canonical_roundtrip_verified": True, "read_receipt": receipt}


def read_receipt_token(canonical_sha256: str) -> str:
    """Transport receipt shown at the end of a display file.

    Echoing it shows the worker opened that file to its end. It is a delivery
    receipt, not proof that the material was understood or correctly judged.
    """
    return "rr-" + hashlib.sha256(f"cqtl-read-receipt-v1:{canonical_sha256}".encode()).hexdigest()[:16]


def _worker_task(reading_file: str, case_id: str, parent_bundle_sha256: str = "",
                 research_ledgers: tuple[dict[str, str], ...] = ()) -> bytes:
    """Launch instructions only; the frozen assignment owns science and IDs."""
    text = (
        f"Research discovery parent_bundle_sha256: {parent_bundle_sha256 or 'request the pinned parent hash from coordinator before returning discoveries'}.\n"
        "Independent Route B assignment. The working directory and permitted read "
        "root must be the prepare output directory.\n"
        f"The canonical assignment is {reading_file}. Read display/{reading_file} "
        "completely, in contiguous line ranges when needed; its lines are at most 1000 characters. "
        "The display envelope's value retains the original object fields. Only paths listed in "
        "chunked_string_paths represent long strings as __cqtl_display_string_chunks__ arrays: "
        "join all chunks in order with NO added characters, including system/user_prefix. "
        "renamed_keys records any long original key as key_chunks. A marker-like field at an "
        "unlisted path is original data. These are lossless display annotations, never the "
        "scientific response schema. Canonical files and prompt hashes are unchanged. "
        "Its materials[].file paths are relative to the prepare output root, "
        "not the assignment's own directory. For each path read display/<materials[].file>, "
        "the bounded-line display of that linked shared-material file, "
        "including system and user_prefix, before interpreting each unit. "
        "These explicitly assigned files are permitted reads even though they "
        "are outside reading/b*/. Do not read sibling assignments, A outputs, "
        "truth or paths outside the prepare output root. No shell, edits or children.\n"
        "Follow the original system/user_prefix claim schema. Answer EVERY "
        "assigned prompt_id/route_b_unit_id separately; preserve all members, "
        "bindings and source IDs. Reuse already read shared material. Missing "
        "or unread material means not_evaluated, never no claims or a biological "
        "negative. Read oversized_full_material units in chunks without dropping them.\n"
        "Return one JSON object only, without fences or prose: "
        f'{{"case_id":{json.dumps(case_id, ensure_ascii=False)},"responses":[...]}}. '
        "Each completed response contains prompt_id, route_b_unit_id and raw_text "
        "(a JSON string encoding the original claims object, with proper escaping). "
        "Every display file ends with a read_receipt value. Each completed response also "
        "contains read_receipts: the list of read_receipt values of the assignment display "
        "and of every material display that this unit links to, copied after you have read "
        "that file. A completed response whose unit receipts are incomplete is treated as "
        "unread: with no claims it becomes not_evaluated. "
        "Use an empty claims list only for an actually reviewed unit with no claims. "
        "The original research ledger retains its inspected scope and reason. A completed "
        "initial functional assessment does not require a literature search for "
        "every unit and does not establish biological exclusion. Mere indexing, "
        "unread or uninterpretable material, and unresolved concrete leads remain "
        "not_evaluated. Do not replace the original researcher's assessment with "
        "a keyword rule or add another biological reviewer. "
        "Web research is allowed by the original prompt. If you retrieve a relevant "
        "paper outside stable_source_ids, DO NOT discard it or cite it against the "
        "old bundle. Pause claim submission and return a separate discovery JSON "
        "with case_id, coordinator-pinned parent_bundle_sha256, and discovered_sources. "
        "Each discovery has unit_id, reference_gene/reference_species (existing frozen "
        "endpoint), reported_gene/reported_gene_id/reported_species (paper gene, possibly "
        "a different plant), source_id, source_title, source_identifier_aliases, source_url, "
        "retrieved_at, source_text and a short source-based reference_function. "
        "source_text is a concise relevant retrieved passage, not a full paper. "
        "The coordinator runs route-b-semantic-bundle --extend-bundle PARENT "
        "--parent-bundle-sha256 HASH --research-discoveries DISCOVERIES --out-dir NEW, "
        "prepares the new prompts, and resumes changed units. Unchanged prompts may "
        "reuse exact responses through reuse-plan. Discovery is not a claim or support. "
        "Each unavailable unit contains prompt_id, route_b_unit_id, "
        'status="not_evaluated" and a non-empty reason. Do not invent a units/reviews '
        "envelope. The parent records executor/model provenance from native receipts; "
        "do not guess it.\n"
    )
    if research_ledgers:
        text += "".join(
            f"Pinned research ledger, permitted read: {row['file']} (relative to the prepare output root; "
            f"SHA-256 {row['sha256']}).\n"
            for row in research_ledgers
        )
        text += (
            "Research ledger transcription mode (opt-in): transcribe the original recorded judgments from those "
            "ledgers only - claims, the nine checklist answers, source IDs and member/copy correspondence - and add "
            "no new research, search or judgment of your own. A unit with no recorded judgment in the ledgers stays "
            "not_evaluated with a reason that names 台帳に判断なし; that is neither a negative result nor a reviewed "
            "no-claim answer. A ledger entry is never a substitute for the assigned display material or its "
            "read_receipt: every linked material file is still read completely as instructed above. Output "
            "discipline: a long result (about 4,000 characters or more) is written to an artifact and only its "
            "path, SHA-256, record count and line count are shown; an omitted (truncated) part stays unread until "
            "that exact line range is read again. These rules do not depend on the execution environment, local or "
            "cloud, or on the provider.\n"
        )
    text = "\n".join(textwrap.fill(line, width=DISPLAY_MAX_LINE_CHARS,
                                  break_long_words=False, break_on_hyphens=False)
                     for line in text.splitlines()) + "\n"
    if max(map(len, text.splitlines()), default=0) > DISPLAY_MAX_LINE_CHARS:
        raise ValueError("worker task exceeds the physical line bound")
    return text.encode("utf-8")


def _reference_handle(record: dict[str, Any]) -> dict[str, Any]:
    # Labels help a reader open the relevant full card, but are not identity keys.
    return {key: record[key] for key in (
        "reference_gene", "reference_gene_id", "reference_species", "source_id",
        "source_ids", "reference_source_id", "source_url", "reference_source_url",
    ) if key in record}


def _prefer_shared_material_order(
    views: list[tuple[dict[str, Any], set[str]]],
    metadata: dict[str, dict[str, Any]],
) -> list[tuple[dict[str, Any], set[str]]]:
    """Co-locate exact shared reference material; never select or merge units.

    Ignore generic system/user-prefix reuse. Prefer a repeated full reference
    card, then a repeated gene-reference record. Ordering within each group
    remains stable, and every view appears exactly once.
    """
    def relevant(mid: str) -> bool:
        return mid.startswith(("reference-material-cards-", "gene-reference-"))

    frequency = Counter(
        mid for _, refs in views for mid in refs if relevant(mid)
    )
    groups: dict[str, list[tuple[dict[str, Any], set[str]]]] = {}
    for index, item in enumerate(views):
        reused = [mid for mid in item[1] if relevant(mid) and frequency[mid] > 1]
        cards = [mid for mid in reused if mid.startswith("reference-material-cards-")]
        eligible = cards or reused
        if eligible:
            # Expected reuse is only a transport heuristic, not biology.
            group_key = max(
                eligible,
                key=lambda mid: (metadata[mid]["bytes"] * (frequency[mid] - 1), mid),
            )
        else:
            group_key = f"unshared-{index}"
        groups.setdefault(group_key, []).append(item)
    return [item for group in groups.values() for item in group]


def plan_reading_views(
    prompts: list[dict[str, Any]], case_id: str, source_sha256: str,
    max_prompts: int, max_material_bytes: int, *, prefer_shared_material: bool = False, parent_bundle_sha256: str = "",
    research_ledgers: tuple[dict[str, str], ...] = (),
) -> tuple[dict[str, bytes], list[dict[str, Any]]]:
    """Preflight every bounded assignment before prepare creates any files."""
    if max_material_bytes < 1:
        raise ValueError("max_material_bytes must be positive")
    files: dict[str, bytes] = {}
    metadata: dict[str, dict[str, Any]] = {}
    views = []

    def share(kind: str, value: Any, refs: set[str]) -> str:
        material = {"kind": kind, "value": value}
        mid = f"{kind}-{_digest(_encoded(material))}"
        name = f"shared-material/{mid}.json"
        data = _encoded(material)
        files.setdefault(name, data)
        metadata[mid] = {"material_id": mid, "file": name, "sha256": _digest(data), "bytes": len(data)}
        refs.add(mid)
        return mid

    def binding(kind: str, record: dict[str, Any], local_keys: set[str], refs: set[str]) -> dict[str, Any]:
        local = {k: v for k, v in record.items() if k in local_keys}
        material = {k: v for k, v in record.items() if k not in local_keys}
        return {"binding": local, "material_id": share(kind, material, refs), "reference_handle": _reference_handle(material)}

    for prompt in prompts:
        refs: set[str] = set()
        view: dict[str, Any] = {
            "prompt": {k: v for k, v in prompt.items() if k not in {"system", "user"}},
            "prompt_sha256": prompt_digest(prompt),
        }
        if "system" in prompt:
            view["system_material_id"] = share("system", prompt["system"], refs)
        user = prompt.get("user")
        unit = None
        if isinstance(user, str):
            # Only the canonical exporter layout is transformed. Unknown layouts
            # remain an exact inline string instead of being guessed or normalized.
            boundary = user.find("\n\n{")
            if boundary >= 0:
                begin = boundary + 2
                try:
                    parsed = json.loads(user[begin:])
                    if isinstance(parsed, dict) and json.dumps(parsed, ensure_ascii=False, sort_keys=True) == user[begin:]:
                        unit = parsed
                except (ValueError, TypeError):
                    pass
        if unit is None:
            if "user" in prompt:
                view["user_inline"] = user
        else:
            view["user_prefix_material_id"] = share("user-prefix", user[:begin], refs)
            packed_unit = dict(unit)
            tables = {}
            for field in ("genes", "reference_material_cards", "source_provenance", "literature_link_candidates"):
                records = unit.get(field)
                if not isinstance(records, list) or not all(isinstance(r, dict) for r in records):
                    continue
                packed_unit.pop(field)
                if field == "genes":
                    tables[field] = [binding("gene-reference", r, {k for k in r if not (k.startswith("reference_") or k in _REFERENCE_FIELDS)}, refs) for r in records]
                else:
                    local_keys = {"interval_gene_id"} if field == "reference_material_cards" else {"gene_id", "interval_gene_id", "unit_id"}
                    tables[field] = [binding(field.replace("_", "-"), r, local_keys, refs) for r in records]
            view["unit"] = packed_unit
            view["unit_tables"] = tables
        views.append((view, refs))

    def render(items: list[dict[str, Any]], refs: set[str], index: int) -> tuple[str, bytes, int]:
        name = f"reading/b{index:04d}/assignment.json"
        materials = [metadata[mid] for mid in sorted(refs)]
        assignment = {"contract": CONTRACT, "case_id": case_id, "source_sha256": source_sha256,
                      "instructions": _NOTE, "max_material_bytes": max_material_bytes,
                      "reading_state": "bounded", "prompts": items, "materials": materials}
        material_bytes = sum(m["bytes"] for m in materials)
        data = _encoded(assignment)
        if len(data) + material_bytes > max_material_bytes:
            assignment["reading_state"] = "oversized_full_material"
            data = _encoded(assignment)
        return name, data, len(data) + material_bytes

    def pack(ordered_views: list[tuple[dict[str, Any], set[str]]]) -> tuple[dict[str, bytes], list[dict[str, Any]]]:
        reading_files: dict[str, bytes] = {}
        shards: list[dict[str, Any]] = []
        current: list[dict[str, Any]] = []
        current_refs: set[str] = set()

        def flush() -> None:
            name, data, size = render(current, current_refs, len(shards))
            reading_files[name] = data
            task_name = str(Path(name).with_name("worker_task.txt")).replace("\\", "/")
            task_data = _worker_task(name, case_id, parent_bundle_sha256, research_ledgers)
            reading_files[task_name] = task_data
            shards.append({"file": name, "sha256": _digest(data), "assignment_bytes": len(data),
                           "worker_task": {"file": task_name, "sha256": _digest(task_data), "bytes": len(task_data)},
                           "material_bytes": size - len(data), "total_reading_bytes": size,
                           "reading_state": "bounded" if size <= max_material_bytes else "oversized_full_material",
                           "over_bound_bytes": max(0, size - max_material_bytes),
                           "prompt_ids": [v["prompt"]["prompt_id"] for v in current],
                           "unit_ids": [v["prompt"]["route_b_unit_id"] for v in current]})

        for view, refs in ordered_views:
            proposed = current + [view]
            merged_refs = current_refs | refs
            _, _, size = render(proposed, merged_refs, len(shards))
            if current and (len(proposed) > max_prompts or size > max_material_bytes):
                flush()
                current, current_refs = [], set()
            current.append(view)
            current_refs.update(refs)
        if current:
            flush()
        return reading_files, shards

    reading_files, shards = pack(views)
    if prefer_shared_material:
        preferred = _prefer_shared_material_order(views, metadata)
        if [v["prompt"]["prompt_id"] for v, _ in preferred] != [v["prompt"]["prompt_id"] for v, _ in views]:
            preferred_files, preferred_shards = pack(preferred)

            def cost(items: list[dict[str, Any]]) -> tuple[int, int, int, int]:
                return (sum(s["total_reading_bytes"] for s in items), len(items),
                        max((s["total_reading_bytes"] for s in items), default=0),
                        sum(s["reading_state"] == "oversized_full_material" for s in items))

            baseline_cost, preferred_cost = cost(shards), cost(preferred_shards)
            # Shared cards are only a heuristic: co-location can destroy other
            # sharing or create extra workers. Publish only a measured byte gain
            # with no increase in workers, largest assignment or oversized units.
            if preferred_cost[0] < baseline_cost[0] and all(a <= b for a, b in zip(preferred_cost, baseline_cost)):
                reading_files, shards = preferred_files, preferred_shards
    files.update(reading_files)
    displays = {}
    for name, data in list(files.items()):
        if name.endswith(".json"):
            display_data, descriptor = _display_sidecar(data, name)
            files[descriptor["file"]] = display_data
            displays[name] = descriptor
    refs_by_prompt = {view["prompt"]["prompt_id"]: refs for view, refs in views}
    for shard in shards:
        assignment = json.loads(files[shard["file"]])
        shard["display_view"] = {
            "contract": DISPLAY_CONTRACT, "assignment": displays[shard["file"]],
            "materials": [displays[material["file"]] for material in assignment["materials"]],
        }
        # Each answer must echo the receipts of the assignment and of every
        # material its own unit links to, so unread units become detectable.
        own = displays[shard["file"]]["read_receipt"]
        shard["read_receipts_by_prompt"] = {
            pid: sorted({own, *(displays[metadata[mid]["file"]]["read_receipt"]
                                for mid in refs_by_prompt[pid])})
            for pid in shard["prompt_ids"]
        }
    return files, shards


def restore_reading_shard(out_dir: str | Path, reading_file: str) -> list[dict[str, Any]]:
    """Verify hashes/roster, restore complete prompts, and compare canonical shards."""
    root = Path(out_dir).resolve()
    manifest = json.loads((root / "prepare_manifest.json").read_text(encoding="utf-8"))
    reading = manifest["reading_view"]
    if reading["contract"] != CONTRACT:
        raise ValueError("unsupported reading view contract")
    entry = next((s for s in reading["shards"] if s["file"] == reading_file), None)
    if entry is None:
        raise ValueError("unknown reading shard")

    def read(name: str, sha256: str) -> bytes:
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("reading material path escapes prepared directory")
        data = path.read_bytes()
        if _digest(data) != sha256:
            raise ValueError(f"reading artifact hash mismatch: {name}")
        return data

    assignment_data = read(reading_file, entry["sha256"])
    # Old prepared generations have no launch sidecar and remain replayable.
    if "worker_task" in entry:
        task = entry["worker_task"]
        task_data = read(task["file"], task["sha256"])
        # Verify recorded bytes, not today's template: historical launch text
        # must remain replayable when a later generation changes instructions.
        if len(task_data) != task["bytes"]:
            raise ValueError("worker task size mismatch")
        if "display_view" in entry and max(map(len, task_data.decode("utf-8").splitlines()), default=0) > DISPLAY_MAX_LINE_CHARS:
            raise ValueError("worker task exceeds the physical line bound")
    assignment = json.loads(assignment_data)
    if (assignment.get("contract") != CONTRACT or assignment.get("case_id") != manifest["case_id"]
            or assignment.get("source_sha256") != manifest["source_sha256"]):
        raise ValueError("reading assignment source/case mismatch")
    materials = {}
    size = len(assignment_data)
    for m in assignment["materials"]:
        data = read(m["file"], m["sha256"])
        material = json.loads(data)
        if m["material_id"] != f"{material['kind']}-{_digest(_encoded(material))}" or len(data) != m["bytes"]:
            raise ValueError("material identity/size mismatch")
        if m["material_id"] in materials:
            raise ValueError("duplicate material ID")
        materials[m["material_id"]] = material["value"]
        size += len(data)
    expected_state = "bounded" if size <= reading["max_material_bytes"] else "oversized_full_material"
    if (size != entry["total_reading_bytes"] or entry["reading_state"] != expected_state
            or assignment["reading_state"] != expected_state
            or entry["over_bound_bytes"] != max(0, size - reading["max_material_bytes"])
            or assignment["max_material_bytes"] != reading["max_material_bytes"]
            or (expected_state == "oversized_full_material" and len(assignment["prompts"]) != 1)):
        raise ValueError("reading shard byte bound mismatch")
    restored = []
    for view in assignment["prompts"]:
        prompt = dict(view["prompt"])
        if "system_material_id" in view:
            prompt["system"] = materials[view["system_material_id"]]
        if "unit" in view:
            unit = dict(view["unit"])
            for field, rows in view["unit_tables"].items():
                unit[field] = []
                for row in rows:
                    record = materials[row["material_id"]]
                    if not isinstance(record, dict) or set(record) & set(row["binding"]):
                        raise ValueError("invalid material binding")
                    if row["reference_handle"] != _reference_handle(record):
                        raise ValueError("reference/source handle mismatch")
                    unit[field].append({**record, **row["binding"]})
            prompt["user"] = materials[view["user_prefix_material_id"]] + json.dumps(unit, ensure_ascii=False, sort_keys=True)
        elif "user_inline" in view:
            prompt["user"] = view["user_inline"]
        if prompt_digest(prompt) != view["prompt_sha256"]:
            raise ValueError("restored prompt hash mismatch")
        restored.append(prompt)
    if ([p["prompt_id"] for p in restored] != entry["prompt_ids"]
            or [p["route_b_unit_id"] for p in restored] != entry["unit_ids"]
            or len(restored) > manifest["max_prompts"]):
        raise ValueError("reading shard roster/count mismatch")
    wanted = set(entry["prompt_ids"])
    canonical = {}
    for shard in manifest["shards"]:
        if wanted.intersection(shard["prompt_ids"]):
            for line in read(shard["file"], shard["sha256"]).decode("utf-8").splitlines():
                p = json.loads(line)
                canonical[p["prompt_id"]] = p
    if any(canonical.get(p["prompt_id"]) != p for p in restored):
        raise ValueError("reading view differs from frozen canonical prompt")
    display_view = entry.get("display_view")
    if display_view is not None:
        if display_view.get("contract") != DISPLAY_CONTRACT:
            raise ValueError("unsupported display view contract")
        displays = [display_view["assignment"], *display_view["materials"]]
        expected_files = {reading_file, *(material["file"] for material in assignment["materials"])}
        if len(displays) != len(expected_files) or {item["source_file"] for item in displays} != expected_files:
            raise ValueError("display material roster mismatch")
        for descriptor in displays:
            original = read(descriptor["source_file"], descriptor["source_sha256"])
            display_data = read(descriptor["file"], descriptor["sha256"])
            lines = display_data.decode("utf-8").splitlines()
            maximum = max(map(len, lines), default=0)
            if (len(display_data) != descriptor["bytes"] or len(lines) != descriptor["line_count"]
                    or maximum != descriptor["max_line_characters"] or maximum > DISPLAY_MAX_LINE_CHARS):
                raise ValueError("display size/line bound mismatch")
            restore_display(display_data, original)
    return restored
