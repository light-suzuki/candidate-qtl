"""Runtime access to dependencies that are not installable everywhere.

`pysam` has no native Windows build, and it used to be imported at the top of
`core`, which every other module imports. That made the whole package
unimportable on native Windows -- including the commands that only read and
write TSVs and never touch a genome. Importing it at the point of use instead
keeps the table-only half of the CLI available wherever Python runs, and turns a
bare `ModuleNotFoundError` into a message that says where to run the command.
"""

from __future__ import annotations

from types import ModuleType


# Commands whose work is sequence I/O (FASTA indexing, region extraction, VCF
# reading) and therefore cannot run without pysam. Reported by `doctor` so the
# split is discoverable rather than something you find by hitting an error.
PYSAM_DEPENDENT_COMMANDS = (
    "annotation-rescue",
    "auto-triage",
    "benchmark-reference",
    "candidate-pipeline",
    "case-run",
    "discover",
    "exon-sequencing-plan",
    "extract",
    "gff-preflight",
    "guided-profile",
)

_PYSAM_UNAVAILABLE = (
    "{operation} needs pysam, which has no native Windows build.\n"
    "  - In WSL or Linux: python -m pip install -c constraints.txt -e '.[alignment]'\n"
    "  - From Windows PowerShell: use scripts\\candidate-qtl.ps1, which forwards to WSL.\n"
    "Table-only commands (rank, evidence-rank, wet-state, hypothesis-ledger, "
    "next-experiment, ...) do not need pysam and run natively."
)


def pysam_available() -> bool:
    """Report whether the sequence-I/O half of the CLI can run here."""
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return True


def require_pysam(operation: str) -> ModuleType:
    """Return the pysam module, or explain where to run ``operation`` instead.

    Raises RuntimeError (not ImportError) so the CLI treats a missing optional
    runtime as an expected, actionable failure rather than a defect.
    """
    try:
        import pysam
    except ImportError as error:
        raise RuntimeError(_PYSAM_UNAVAILABLE.format(operation=operation)) from error
    return pysam
