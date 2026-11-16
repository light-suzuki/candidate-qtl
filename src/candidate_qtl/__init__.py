"""Candidate-first QTL validation utilities.

The package initializer is intentionally side-effect free. CLI routes import
their hardened implementations directly so importing ``candidate_qtl`` never
changes functions or shared dictionaries in other modules.
"""

# Single source of truth for the version. `pyproject.toml` declares
# `dynamic = ["version"]` and resolves it from this literal via
# `[tool.setuptools.dynamic]` (attr = "candidate_qtl.__version__"), so a
# built wheel and a manifest stamped by `from . import __version__` can
# never disagree. Keep this a plain string literal -- setuptools reads it by
# parsing the AST, not by importing the module or its dependencies.
__version__ = "0.24.0"
