#!/usr/bin/env python3
"""Submit and inspect QTLift jobs through its local REST API.

This client intentionally stores only job metadata and inferred intervals. Genome
FASTA/GFF3 contents remain local to the QTLift server.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def request_json(base_url: str, method: str, path: str, payload: dict | None = None) -> dict:
    data = None if payload is None else json.dumps(payload).encode("utf-8")
    request = Request(
        f"{base_url.rstrip('/')}{path}",
        data=data,
        method=method,
        headers={"Content-Type": "application/json", "User-Agent": "candidate-qtl-qtlift-client/0.1"},
    )
    try:
        with urlopen(request, timeout=60) as response:
            return json.load(response)
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"QTLift HTTP {exc.code}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"Cannot reach QTLift at {base_url}: {exc.reason}") from exc


def write_json(path: str | None, payload: dict) -> None:
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text + "\n", encoding="utf-8")
    print(text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Local API client for QTLift coordinate harmonization")
    parser.add_argument("--base-url", default="http://127.0.0.1:8765")
    sub = parser.add_subparsers(dest="command", required=True)

    health = sub.add_parser("health")
    health.add_argument("--out")

    scan = sub.add_parser("scan")
    scan.add_argument("--genome-root")
    scan.add_argument("--local-library", action="store_true", help="Use QTLift's configured @this-pc genome library.")
    scan.add_argument("--out")

    submit = sub.add_parser("submit")
    submit.add_argument("--genome-root")
    submit.add_argument("--local-library", action="store_true", help="Resolve reference names from local_genomes.json via @this-pc.")
    submit.add_argument("--target-ref", required=True)
    submit.add_argument("--source-ref", required=True)
    submit.add_argument("--contig", required=True)
    submit.add_argument("--start", type=int, required=True)
    submit.add_argument("--end", type=int, required=True)
    submit.add_argument("--peak", type=int)
    submit.add_argument("--name", default="")
    submit.add_argument("--preset", default="Standard")
    submit.add_argument("--target-contig", default="")
    submit.add_argument("--all-genes", action="store_true")
    submit.add_argument("--mapping-backend", default="auto")
    submit.add_argument("--markers-json", help="JSON object mapping marker names to nucleotide sequences")
    submit.add_argument("--params-json", help="Additional QTLift params JSON object")
    submit.add_argument("--tool-paths-json", help="Explicit tool paths JSON object")
    submit.add_argument("--provider-options-json", help="QTLift provider options JSON object")
    submit.add_argument("--wait", action="store_true")
    submit.add_argument("--poll-seconds", type=float, default=2.0)
    submit.add_argument("--timeout-seconds", type=float, default=3600.0, help="Maximum wait time; <=0 disables the limit.")
    submit.add_argument("--out")

    status = sub.add_parser("status")
    status.add_argument("job_id")
    status.add_argument("--out")
    return parser


def parse_json_object(value: str | None, label: str) -> dict:
    if not value:
        return {}
    parsed = json.loads(value)
    if not isinstance(parsed, dict):
        raise ValueError(f"{label} must be a JSON object")
    return parsed


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "health":
            result = request_json(args.base_url, "GET", "/api/health")
            write_json(args.out, result)
            return 0
        if args.command == "scan":
            genome_root = "@this-pc" if args.local_library else args.genome_root
            if not genome_root:
                raise ValueError("provide --genome-root or --local-library")
            result = request_json(args.base_url, "POST", "/api/genomes/scan", {"genome_root": genome_root})
            write_json(args.out, result)
            return 0
        if args.command == "status":
            result = request_json(args.base_url, "GET", f"/api/jobs/{args.job_id}")
            write_json(args.out, result)
            return 0

        if args.start > args.end:
            raise ValueError("--start must be <= --end")
        genome_root = "@this-pc" if args.local_library else args.genome_root
        if not genome_root:
            raise ValueError("provide --genome-root or --local-library")
        library = request_json(args.base_url, "POST", "/api/genomes/scan", {"genome_root": genome_root})
        names = sorted(str(item.get("name")) for item in library.get("genomes", []))
        missing = [name for name in (args.source_ref, args.target_ref) if name not in names]
        if missing:
            raise ValueError(f"reference(s) not found in genome library/root: {', '.join(missing)}; available: {', '.join(names)}")
        payload = {
            "genome_root": genome_root,
            "target_ref": args.target_ref,
            "source_ref": args.source_ref,
            "contig": args.contig,
            "start": args.start,
            "end": args.end,
            "peak": args.peak,
            "name": args.name,
            "preset": args.preset,
            "params": parse_json_object(args.params_json, "--params-json"),
            "markers": parse_json_object(args.markers_json, "--markers-json"),
            "tool_paths": parse_json_object(args.tool_paths_json, "--tool-paths-json"),
            "mapping_backend": args.mapping_backend,
            "provider_options": parse_json_object(args.provider_options_json, "--provider-options-json"),
            "all_genes": args.all_genes,
            "target_contig": args.target_contig,
        }
        result = request_json(args.base_url, "POST", "/api/jobs", payload)
        job_id = result.get("job_id") or result.get("id")
        if args.wait and job_id:
            started = time.monotonic()
            while True:
                result = request_json(args.base_url, "GET", f"/api/jobs/{job_id}")
                status = str(result.get("status", "")).lower()
                if status in {"completed", "failed", "cancelled", "canceled"}:
                    break
                if args.timeout_seconds > 0 and time.monotonic() - started >= args.timeout_seconds:
                    raise RuntimeError(f"QTLift job {job_id} did not finish within {args.timeout_seconds:g} seconds")
                time.sleep(max(args.poll_seconds, 0.5))
        write_json(args.out, result)
        return 0
    except (ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
