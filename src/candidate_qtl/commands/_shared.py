"""Argparse helpers shared by cli.py and the feature command modules.

Kept here (importing nothing from cli) so both the monolithic legacy parser and
the migrated feature modules can reuse them without an import cycle.
"""

from __future__ import annotations

import argparse


def add_json_flag(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--json", action="store_true", help="Emit a stable JSON result to stdout.")


def add_primerblast_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--primerblast-root", required=True, help="primerblast-oss checkout root.")
    parser.add_argument("--genome", required=True, help="Reference FASTA for primerblast-oss.")
    parser.add_argument("--db", action="append", required=True, help="BLAST nucleotide DB; repeat for multi-reference screening.")
    parser.add_argument("--num-threads", type=int, default=4, help="Laptop-safe worker count for BLAST-backed primer checks.")
