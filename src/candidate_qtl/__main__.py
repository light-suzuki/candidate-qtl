from .cli import run


# `run`, not `main`: expected input errors must exit non-zero with one actionable
# line, not a traceback. Keep this in step with the `candidate-qtl` console
# script in pyproject.toml so both entry points behave identically.
raise SystemExit(run())
