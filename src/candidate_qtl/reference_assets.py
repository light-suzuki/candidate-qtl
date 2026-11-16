"""Independent reference copies and content-bound interval extraction (#843).

Original inputs are only read. Managed files are independent, read-only copies,
not filesystem snapshots or protection against concurrent external writers.
Every registration checks full source and managed-file hashes; cache hits never
use mtime/size as a substitute for content identity.
"""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import stat
import sys
import tempfile

from .core import parse_interval, parse_seqid_map, runtime_path, sha256
from .optional_deps import require_pysam
from .reference_extract import extract_candidates


ASSET_CONTRACT = "managed-reference-1.0.0"
INTERVAL_CONTRACT = "reference-interval-cache-1.0.0"
EXTRACTION_POLICY_VERSION = "complete-gene-model-1.0.0"
PROTEIN_DB_CONTRACT = "shared-protein-db-1.0.0"
ASSET_MANIFEST = "reference_asset.json"
CACHE_MANIFEST = "interval_cache.json"
EXTRACT_MANIFEST = "extract_manifest.json"
INTERVAL_OUTPUTS = {
    "candidates_tsv": "candidates.tsv",
    "transcripts_tsv": "transcripts.tsv",
    "id_map_tsv": "id_map.tsv",
    "candidate_genomic_fasta": "candidate_genomic.fa",
    "candidate_cds_fasta": "candidate_cds.fa",
    "candidate_transcript_cds_fasta": "candidate_transcript_cds.fa",
    "candidate_proteins_fasta": "candidate_proteins.fa",
}


def _key(value: object) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                     separators=(",", ":")).encode("utf-8")).hexdigest()


def _plain_path(value: str | Path) -> Path:
    path = runtime_path(value).absolute()
    # Check before resolve(): resolving first would hide a symlinked ancestor.
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError(f"reference/cache paths must not contain symlinks: {path}")
    return path.resolve()


def _file_stat(path: Path) -> list[int]:
    _plain_path(path)
    value = path.stat()
    if not stat.S_ISREG(value.st_mode) or value.st_nlink != 1:
        raise ValueError(f"reference/cache input must be a regular, non-hardlinked file: {path}")
    return [value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns]


def _unchanged(paths: dict[str, Path], expected: dict[str, list[int]]) -> None:
    for name, path in paths.items():
        if _file_stat(path) != expected[name]:
            raise ValueError(f"reference resource changed during use: {path}")


def _capture(paths: dict[str, Path]) -> tuple[dict, dict]:
    before = {name: _file_stat(path) for name, path in paths.items()}
    records = {
        name: {"path": str(path), "size_bytes": before[name][2], "sha256": sha256(path),
               "method": "full_sha256", "strength": "strong", "mtime_ns": before[name][3]}
        for name, path in paths.items()
    }
    _unchanged(paths, before)
    return records, before


def _contents(records: dict) -> dict:
    return {name: {"sha256": value["sha256"], "size_bytes": value["size_bytes"]}
            for name, value in records.items()}


def _version(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("a non-empty reference version is required")
    return value.strip()


def _separate(output: Path, inputs: list[Path]) -> None:
    if any(output.is_relative_to(path) or path.is_relative_to(output) for path in inputs):
        raise ValueError("reference/cache output must be separate from source inputs and cache entries")


def _directory(path: Path) -> None:
    _plain_path(path)
    path.mkdir(parents=True, exist_ok=True)
    if not path.is_dir():
        raise ValueError(f"reference/cache directory required: {path}")
    _plain_path(path)


def _roster(directory: Path, names: set[str]) -> dict[str, Path]:
    _plain_path(directory)
    if not directory.is_dir() or {p.name for p in directory.iterdir()} != names:
        raise ValueError(f"incomplete or unexpected files in reference/cache entry: {directory}")
    paths = {name: directory / name for name in sorted(names)}
    for path in paths.values():
        _file_stat(path)
    return paths


def _read_json(path: Path) -> tuple[dict, str, list[int]]:
    before = _file_stat(path)
    raw = path.read_bytes()
    if _file_stat(path) != before:
        raise ValueError(f"reference manifest changed during read: {path}")
    try:
        value = json.loads(raw)
    except (ValueError, UnicodeDecodeError) as error:
        raise ValueError(f"invalid reference/cache manifest: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"reference/cache manifest must be a JSON object: {path}")
    return value, hashlib.sha256(raw).hexdigest(), before


def _write_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def _copy(source: Path, target: Path) -> None:
    _file_stat(source)
    _plain_path(target)
    # Exclusive creation cannot truncate an existing symlink or hardlink.
    with source.open("rb") as reader, target.open("xb") as writer:
        shutil.copyfileobj(reader, writer, 1024 * 1024)
        writer.flush()
        os.fsync(writer.fileno())


@contextmanager
def _staging(parent: Path, key: str):
    """Publish only complete entries; simultaneous creators fail with retryable work."""
    lock = parent / (".lock-" + key)
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise ValueError(f"reference/cache creation already in progress; retry: {key}") from error
    temporary = None
    try:
        temporary = Path(tempfile.mkdtemp(prefix=".tmp-" + key + "-", dir=parent))
        yield temporary
    finally:
        if temporary is not None and temporary.exists():
            shutil.rmtree(temporary)
        lock.rmdir()


def _publish(temporary: Path, destination: Path) -> None:
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"reference/cache entry appeared during creation; retry: {destination}")
    for path in temporary.iterdir():
        _file_stat(path)
        path.chmod(0o444)
    temporary.rename(destination)


def _load_asset(manifest_path: Path, *, expected_sha: str | None = None,
                expected_identity: dict | None = None) -> tuple[dict, dict, str, dict]:
    value, manifest_sha, manifest_stat = _read_json(manifest_path)
    identity = value.get("identity")
    if (value.get("contract") != ASSET_CONTRACT or not isinstance(identity, dict)
            or identity.get("contract") != ASSET_CONTRACT
            or value.get("key") != _key(identity) or manifest_path.parent.name != value.get("key")
            or (expected_sha is not None and manifest_sha != expected_sha)
            or (expected_identity is not None and identity != expected_identity)):
        raise ValueError("reference asset manifest identity/hash mismatch")
    gff_name = identity.get("gff_name")
    if gff_name not in {"annotation.gff3", "annotation.gff3.gz"}:
        raise ValueError("invalid managed annotation filename")
    names = {"reference.fa", "reference.fa.fai", gff_name}
    _roster(manifest_path.parent, names | {ASSET_MANIFEST})
    paths = {name: manifest_path.parent / name for name in sorted(names)}
    records, observed = _capture(paths)
    if (value.get("files") != _contents(records)
            or identity.get("sources") != {"fasta": _contents(records)["reference.fa"],
                                           "gff3": _contents(records)[gff_name]}):
        raise ValueError("managed reference file SHA mismatch")
    _unchanged({ASSET_MANIFEST: manifest_path}, {ASSET_MANIFEST: manifest_stat})
    observed[ASSET_MANIFEST] = manifest_stat
    return value, records, manifest_sha, observed


def prepare_reference_asset(fasta_arg: str, gff_arg: str, cache_dir_arg: str, *,
                            reference_version: str) -> dict:
    """Copy/index a reference once; recheck original and managed full SHA on reuse."""
    version = _version(reference_version)
    sources = {"fasta": _plain_path(fasta_arg), "gff3": _plain_path(gff_arg)}
    cache = _plain_path(cache_dir_arg)
    _separate(cache, list(sources.values()))
    source_records, source_stats = _capture(sources)
    with sources["fasta"].open("rb") as stream:
        if stream.read(2) == b"\x1f\x8b":
            raise ValueError("managed reference requires an uncompressed FASTA")
    gff_name = "annotation.gff3.gz" if sources["gff3"].suffix == ".gz" else "annotation.gff3"
    identity = {"contract": ASSET_CONTRACT, "reference_version": version,
                "sources": _contents(source_records), "gff_name": gff_name}
    key = _key(identity)
    parent = cache / "references"
    _directory(parent)
    entry = parent / key
    hit = entry.exists() or entry.is_symlink()
    if not hit:
        with _staging(parent, key) as temporary:
            _copy(sources["fasta"], temporary / "reference.fa")
            _copy(sources["gff3"], temporary / gff_name)
            copied, _ = _capture({"fasta": temporary / "reference.fa", "gff3": temporary / gff_name})
            if _contents(copied) != _contents(source_records):
                raise ValueError("reference changed while creating independent copy")
            _unchanged(sources, source_stats)
            require_pysam("extract").faidx(str(temporary / "reference.fa"))
            names = {"reference.fa", "reference.fa.fai", gff_name}
            paths = _roster(temporary, names)
            files, _ = _capture(paths)
            _write_json(temporary / ASSET_MANIFEST, {
                "contract": ASSET_CONTRACT, "key": key, "identity": identity,
                "reference_version": version, "source_fingerprints": source_records,
                "files": _contents(files),
                "index_builder": {"pysam": importlib.metadata.version("pysam")},
                "protection": "independent read-only copies; not an OS snapshot or writer lock",
            })
            _unchanged(sources, source_stats)
            _publish(temporary, entry)
    manifest_path = entry / ASSET_MANIFEST
    _, _, manifest_sha, _ = _load_asset(manifest_path, expected_identity=identity)
    _unchanged(sources, source_stats)
    return {"fasta": str(entry / "reference.fa"), "gff3": str(entry / gff_name),
            "manifest": str(manifest_path), "manifest_sha256": manifest_sha,
            "cache_hit": hit, "source_fingerprints": source_records}


def _extraction_implementation() -> dict:
    package = Path(__file__).parent
    return {
        "policy": EXTRACTION_POLICY_VERSION,
        "code_sha256": {name: sha256(package / name) for name in
                        ("reference_extract.py", "core.py", "reference_assets.py", "schemas.py")},
        "dependencies": {name: importlib.metadata.version(name) for name in ("biopython", "pandas", "pysam")},
        "python": sys.version.split()[0],
    }


def _seqid_mapping(value: str | None) -> dict[str, str]:
    if value is None:
        return {}
    raw = runtime_path(value)
    if raw.exists() or raw.is_symlink():
        path = _plain_path(raw)
        before = _file_stat(path)
        mapping = parse_seqid_map(str(path))
        _unchanged({"seqid_map": path}, {"seqid_map": before})
        return mapping
    return parse_seqid_map(value)


def _load_interval(entry: Path, identity: dict) -> tuple[dict, dict, str, dict]:
    names = set(INTERVAL_OUTPUTS.values()) | {EXTRACT_MANIFEST}
    _roster(entry, names | {CACHE_MANIFEST})
    value, manifest_sha, manifest_stat = _read_json(entry / CACHE_MANIFEST)
    if (value.get("contract") != INTERVAL_CONTRACT or value.get("identity") != identity
            or value.get("key") != _key(identity) or entry.name != value.get("key")):
        raise ValueError("interval cache identity mismatch")
    records, observed = _capture({name: entry / name for name in sorted(names)})
    if value.get("files") != _contents(records):
        raise ValueError("interval cache artifact SHA mismatch")
    extracted, _, _ = _read_json(entry / EXTRACT_MANIFEST)
    if (extracted.get("outputs") != INTERVAL_OUTPUTS
            or extracted.get("interval") != identity["interval_input"]
            or extracted.get("fasta_fingerprint", {}).get("sha256") != identity["fasta_sha256"]
            or extracted.get("gff3_fingerprint", {}).get("sha256") != identity["gff3_sha256"]
            or extracted.get("fasta_index_fingerprint", {}).get("sha256") != identity["fai_sha256"]):
        raise ValueError("interval cache extraction lineage mismatch")
    _unchanged({CACHE_MANIFEST: entry / CACHE_MANIFEST}, {CACHE_MANIFEST: manifest_stat})
    observed[CACHE_MANIFEST] = manifest_stat
    return extracted, records, manifest_sha, observed


def extract_reference_interval(asset: dict, interval_arg: str, out_arg: str, *,
                               seqid_map_arg: str | None = None,
                               interval_cache_dir_arg: str | None = None) -> dict:
    """Reuse only phenotype-free extraction, with fresh independent case outputs."""
    manifest = _plain_path(asset["manifest"])
    anchor = asset.get("manifest_sha256")
    if not isinstance(anchor, str) or len(anchor) != 64 or any(c not in "0123456789abcdef" for c in anchor):
        raise ValueError("reference asset requires a full SHA-256 manifest anchor")
    registered, files, asset_sha, asset_stats = _load_asset(manifest, expected_sha=asset["manifest_sha256"])
    fasta = manifest.parent / "reference.fa"
    gff = manifest.parent / registered["identity"]["gff_name"]
    if _plain_path(asset["fasta"]) != fasta or _plain_path(asset["gff3"]) != gff:
        raise ValueError("reference asset paths do not match its manifest")
    chrom, start, end = parse_interval(interval_arg)
    if start < 1:
        raise ValueError("QTL interval requires 1-based positive coordinates")
    interval = f"{chrom}:{start}-{end}"
    mapping = _seqid_mapping(seqid_map_arg)
    if any(any(character in value for character in "\t\r\n") for pair in mapping.items() for value in pair):
        raise ValueError("seqid mapping cannot contain tabs or newlines")
    identity = {
        "contract": INTERVAL_CONTRACT, "fasta_sha256": files["reference.fa"]["sha256"],
        "gff3_sha256": files[gff.name]["sha256"], "fai_sha256": files["reference.fa.fai"]["sha256"],
        # Raw formatting is written into candidates.tsv by the existing
        # extractor. Bind both forms to preserve byte-identical fresh output.
        "interval": interval, "interval_input": interval_arg,
        "seqid_map": dict(sorted(mapping.items())),
        "implementation": _extraction_implementation(),
    }
    key = _key(identity)
    parent = (_plain_path(interval_cache_dir_arg) if interval_cache_dir_arg is not None
              else manifest.parents[2] / "intervals")
    out = _plain_path(out_arg)
    _separate(parent, [manifest.parent])
    _separate(out, [manifest.parent, parent])
    output_names = set(INTERVAL_OUTPUTS.values()) | {EXTRACT_MANIFEST}
    if any((out / name).exists() or (out / name).is_symlink() for name in output_names):
        raise ValueError("interval extraction requires fresh output artifacts")
    _directory(parent)
    entry = parent / key
    hit = entry.exists() or entry.is_symlink()
    if not hit:
        with _staging(parent, key) as temporary:
            mapping_path = None
            if mapping:
                mapping_path = temporary / "seqid_map.tsv"
                mapping_path.write_text("".join(f"{a}\t{b}\n" for a, b in sorted(mapping.items())), encoding="utf-8")
            extracted = extract_candidates(str(fasta), str(gff), interval_arg, str(temporary),
                                           seqid_map_arg=str(mapping_path) if mapping_path else None,
                                           read_only_reference=True, full_hash=True)
            if mapping_path is not None:
                mapping_path.unlink()
            extracted.pop("manifest", None)
            extracted["outputs"] = dict(INTERVAL_OUTPUTS)
            (temporary / EXTRACT_MANIFEST).unlink()
            _write_json(temporary / EXTRACT_MANIFEST, extracted)
            records, _ = _capture(_roster(temporary, output_names))
            _write_json(temporary / CACHE_MANIFEST, {
                "contract": INTERVAL_CONTRACT, "key": key, "identity": identity, "files": _contents(records),
            })
            _unchanged({name: manifest.parent / name for name in asset_stats}, asset_stats)
            _publish(temporary, entry)
    extracted, records, cache_sha, cached_stats = _load_interval(entry, identity)
    _directory(out)
    for name in INTERVAL_OUTPUTS.values():
        _copy(entry / name, out / name)
    local_records, _ = _capture({name: out / name for name in INTERVAL_OUTPUTS.values()})
    if _contents(local_records) != {name: _contents(records)[name] for name in INTERVAL_OUTPUTS.values()}:
        raise ValueError("interval cache changed while copying case artifacts")
    _unchanged({name: entry / name for name in cached_stats}, cached_stats)
    _unchanged({name: manifest.parent / name for name in asset_stats}, asset_stats)
    # Input fingerprints are lineage; case outputs always name the case copy.
    extracted.update(interval=interval_arg, fasta=str(fasta), gff3=str(gff),
                     fasta_fingerprint=files["reference.fa"], gff3_fingerprint=files[gff.name],
                     fasta_index_fingerprint=files["reference.fa.fai"])
    extracted["outputs"] = {field: str(out / name) for field, name in INTERVAL_OUTPUTS.items()}
    extracted["reference_asset"] = {"manifest": str(manifest), "manifest_sha256": asset_sha}
    extracted["interval_cache"] = {"key": key, "hit": hit, "manifest": str(entry / CACHE_MANIFEST),
                                   "manifest_sha256": cache_sha,
                                   "extraction_manifest_sha256": records[EXTRACT_MANIFEST]["sha256"]}
    _write_json(out / EXTRACT_MANIFEST, extracted)
    return extracted | {"manifest": str(out / EXTRACT_MANIFEST)}


def _protein_files(prefix: Path) -> dict[str, Path]:
    _plain_path(prefix)
    files = {path.name: path for path in sorted(prefix.parent.iterdir())
             if path.name.startswith(prefix.name + ".")}
    if not files:
        raise ValueError("protein_db must name a formatted BLAST database prefix")
    if any(name.lower().endswith(".pal") for name in files):
        raise ValueError("shared protein DB aliases (.pal) are unsupported; specify a direct database prefix")
    volumes = [{name[:-4] for name in files if name.endswith(suffix)} for suffix in (".pin", ".phr", ".psq")]
    if not volumes[0] or not volumes[0] == volumes[1] == volumes[2]:
        raise ValueError("shared protein DB requires complete .pin/.phr/.psq files for every volume")
    for path in files.values():
        _file_stat(path)
    return files


def capture_shared_protein_db(prefix_arg: str, *, version: str) -> dict:
    """Bind a direct or multi-volume BLAST database without copying/building it."""
    version = _version(version)
    prefix = _plain_path(prefix_arg)
    paths = _protein_files(prefix)
    files, source_stats = _capture(paths)
    if set(_protein_files(prefix)) != set(paths):
        raise ValueError("shared protein DB file roster changed during capture")
    return {"contract": PROTEIN_DB_CONTRACT, "prefix": str(prefix), "version": version,
            "file_roster": sorted(files), "files": files, "full_sha256": _key(_contents(files)),
            "_source_stats": source_stats}


def validate_shared_protein_db(snapshot: dict) -> dict:
    """Reject any DB filename, metadata, or full-content change across B0 use."""
    if snapshot.get("contract") != PROTEIN_DB_CONTRACT:
        raise ValueError("invalid shared protein DB snapshot contract")
    current = capture_shared_protein_db(snapshot["prefix"], version=snapshot["version"])
    if current != snapshot:
        raise ValueError("shared protein DB changed during use (roster/stat/SHA mismatch)")
    return snapshot
