"""Bounded display pages for existing boundary views, without scientific selection."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

try:
    from scripts.agentic_route_b_material import _display_sidecar, _encoded, restore_display
except ModuleNotFoundError:
    from agentic_route_b_material import _display_sidecar, _encoded, restore_display


CONTRACT = "agentic-boundary-reading-pages-1.0.0"
MAX_PAGE_BYTES = 12 * 1024
MAX_PAGE_LINES = 160
MAX_LINE_CHARS = 1000


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sources(manifest: dict) -> list[dict]:
    return [manifest["material_catalog"], manifest["comparison_index_artifact"], *manifest["reading_views"]]


def _page_files(source_file: str, count: int) -> list[str]:
    folder = Path("reading_pages") / Path(source_file).with_suffix("")
    return [(folder / f"page-{index:05d}.txt").as_posix() for index in range(count)]


def _task(shard: dict, sources: list[dict], manifest: dict) -> bytes:
    shared = [manifest[key]["file"] for key in ("material_catalog", "comparison_index_artifact")]
    view = next(item["file"] for item in manifest["reading_views"] if item["shard_id"] == shard["shard_id"])
    selected = [next(item for item in sources if item["source_file"] == name) for name in [view, *shared]]
    lines = [
        f"Boundary comparison: {manifest['case_id']} / {shard['shard_id']}",
        f"source_wet_sha256: {manifest['source_wet_sha256']}",
        f"comparison_context_sha256: {manifest.get('comparison_context', {}).get('sha256', 'not supplied')}",
        "Reading root is .. relative to this worker_task.md directory (the directory containing manifest.json).",
        "Resolve ALL paths below from that reading root, even when your working directory is nested elsewhere.",
        "The permitted read scope must include the whole reading root. Do not read outside it through provenance paths.",
        "Read every numbered page below in increasing order, once per source; reuse already read shared sources.",
        "Each page is at most 12 KiB / 160 lines, with at most 1000 characters per line.",
        "These are conservative transport bounds, not a guarantee for every tool. Check the actual returned last line/EOF.",
        "If output truncates, continue at the next actually returned line within that SAME page before advancing.",
        "Pages concatenate with no extra characters into a display-only JSON envelope, not the scientific response schema.",
        "Read envelope.value. Only paths listed in chunked_string_paths are string chunks to join without separators.",
        "renamed_keys records any escaped long keys. Literal marker-looking fields elsewhere remain original data.",
        "Read the complete reading view: candidate_input AND primary_paper_evidence when present, request_policy,",
        "phenotype, higher_tier_context and candidate_material. comparison_context alone is not all comparison evidence.",
        "Resolve candidate_material inline/ref bindings using the shared catalog and read the shared comparison index.",
        "Follow reader_instruction in the view; retain all candidate/paper/copy bindings and known counterevidence.",
        "Source hashes/pages prove delivery integrity, not model receipt, comprehension or completed scientific comparison.",
        "Use the existing boundary response schema and executor/reviewer role. Missing material is not biological negative.",
    ]
    if manifest.get("material_reading_projection"):
        lines += ["candidate_input rows, higher_tier_context, primary_paper_evidence, request_policy and comparison_scope",
                  "also use the SAME catalog inline/ref bindings. Resolve them without rereading a shared node.",
                  "Valid structured *_json values are decoded for reading; exact original strings remain in canonical shards."]
    for item in selected:
        lines += ["", f"Source: {item['source_file']} (SHA-256 {item['source_sha256']})",
                  f"Read all {len(item['pages'])} pages: {item['pages'][0]['file']} through {item['pages'][-1]['file']}."]
    data = ("\n".join(lines) + "\n").encode("utf-8")
    _check_bounds(data)
    return data


def _check_bounds(data: bytes) -> None:
    lines = data.decode("utf-8").splitlines()
    if (len(data) > MAX_PAGE_BYTES or len(lines) > MAX_PAGE_LINES
            or max(map(len, lines), default=0) > MAX_LINE_CHARS):
        raise ValueError("boundary reading page exceeds transport bounds")


def write_reading_pages(root: Path, manifest: dict) -> dict:
    sources = []
    for source in _sources(manifest):
        original = (root / source["file"]).read_bytes()
        canonical = _encoded(json.loads(original))
        display, _ = _display_sidecar(canonical, source["file"])
        pages, current = [], []
        size = 0
        for line in display.splitlines(keepends=True):
            if current and (size + len(line) > MAX_PAGE_BYTES or len(current) == MAX_PAGE_LINES):
                pages.append(b"".join(current))
                current, size = [], 0
            current.append(line)
            size += len(line)
        if current:
            pages.append(b"".join(current))
        descriptors = []
        for name, data in zip(_page_files(source["file"], len(pages)), pages):
            _check_bounds(data)
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            descriptors.append({"file": name, "sha256": _sha(data), "bytes": len(data),
                                "line_count": len(data.splitlines())})
        restore_display(b"".join(pages), canonical)
        sources.append({"source_file": source["file"], "source_sha256": _sha(original),
                        "display_sha256": _sha(display), "pages": descriptors})
    tasks = []
    for shard in manifest["shards"]:
        name = (Path(shard["file"]).parent / "worker_task.md").as_posix()
        data = _task(shard, sources, manifest)
        (root / name).write_bytes(data)
        tasks.append({"shard_id": shard["shard_id"], "file": name, "sha256": _sha(data)})
    return {"contract": CONTRACT, "sources": sources, "tasks": tasks}


def validate_reading_page_generation(root: Path, manifest: dict) -> None:
    """Do not reinterpret a damaged page-bearing generation as a legacy one."""
    marker = manifest.get("reading_delivery_contract")
    if marker is not None and marker != CONTRACT:
        raise ValueError("unknown boundary reading delivery contract")
    artifacts = (root / "reading_pages").exists() or any(
        (root / Path(shard["file"]).parent / "worker_task.md").exists() for shard in manifest.get("shards", []))
    delivery = manifest.get("reading_pages")
    if delivery is None:
        if marker is not None or artifacts:
            raise ValueError("boundary reading page metadata missing for page-bearing generation")
        return
    if not isinstance(delivery, dict) or any(not manifest.get(key) for key in (
            "material_catalog", "comparison_index_artifact", "reading_views")):
        raise ValueError("boundary reading pages require canonical material descriptors")


def validate_reading_pages(root: Path, manifest: dict) -> None:
    validate_reading_page_generation(root, manifest)
    delivery = manifest.get("reading_pages")
    if delivery is None:
        return  # Frozen generations before bounded pages keep their existing contract.
    if delivery.get("contract") != CONTRACT:
        raise ValueError("unknown boundary reading page contract")
    expected = _sources(manifest)
    sources = delivery.get("sources", [])
    if [item.get("source_file") for item in sources] != [item["file"] for item in expected]:
        raise ValueError("boundary reading page source roster mismatch")

    def checked(descriptor):
        path = (root / descriptor["file"]).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file():
            raise ValueError("boundary reading page/task missing or outside reading root")
        data = path.read_bytes()
        if _sha(data) != descriptor["sha256"]:
            raise ValueError("boundary reading page/task hash mismatch")
        _check_bounds(data)
        return data

    for source, original_descriptor in zip(sources, expected):
        original = (root / original_descriptor["file"]).read_bytes()
        if source["source_sha256"] != original_descriptor["sha256"] or _sha(original) != source["source_sha256"]:
            raise ValueError("boundary reading page source hash mismatch")
        pages = source.get("pages", [])
        if not pages or [p["file"] for p in pages] != _page_files(source["source_file"], len(pages)):
            raise ValueError("boundary reading page order/roster mismatch")
        restored = []
        for descriptor in pages:
            data = checked(descriptor)
            if len(data) != descriptor["bytes"] or len(data.splitlines()) != descriptor["line_count"]:
                raise ValueError("boundary reading page size mismatch")
            restored.append(data)
        display = b"".join(restored)
        if _sha(display) != source["display_sha256"]:
            raise ValueError("boundary reading display hash mismatch")
        restore_display(display, _encoded(json.loads(original)))
    tasks = delivery.get("tasks", [])
    if [task.get("shard_id") for task in tasks] != [shard["shard_id"] for shard in manifest["shards"]]:
        raise ValueError("boundary reading task roster mismatch")
    for task, shard in zip(tasks, manifest["shards"]):
        expected_file = (Path(shard["file"]).parent / "worker_task.md").as_posix()
        if task["file"] != expected_file or checked(task) != _task(shard, sources, manifest):
            raise ValueError("boundary reading task differs from source/page bindings")
