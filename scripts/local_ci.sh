#!/usr/bin/env bash
# Run the GitHub Actions CI workflow's Linux jobs locally, with no billable
# minutes. Every step below mirrors a step in .github/workflows/ci.yml; when
# that file changes, change this one with it.
#
#   scripts/local_ci.sh                # every Linux job
#   scripts/local_ci.sh test lint      # only the named jobs
#   scripts/local_ci.sh --list         # what can be run
#
# Jobs: test, schema-extra, lint, build, integration-blast.
# The windows-wrapper job is PowerShell-only; run scripts/local_ci.ps1 for it.
#
# WHY THIS SCRIPT CREATES ITS OWN VIRTUALENVS
# -------------------------------------------
# The reachability gate runs `python -m candidate_qtl` as a subprocess with
# PYTHONPATH deliberately stripped (so a caller cannot smuggle in a different
# source tree). The subprocess therefore imports whatever `candidate_qtl` the
# interpreter resolves -- i.e. the *installed* package. Point it at a venv that
# holds an editable install of a DIFFERENT checkout (easy to do when several
# git worktrees share one venv) and the gate silently tests that other
# checkout's code against this one's fixtures; the canonical baseline then
# fails with "baseline artifact not produced". So each repo root gets its own
# venv, keyed by its path, under LOCAL_CI_CACHE.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE_ROOT="${LOCAL_CI_CACHE:-$HOME/.cache/candidate-qtl-local-ci}"
# Keyed by repo path so sibling worktrees never share an environment. Kept
# outside the repo: a venv on a /mnt/c drive is slow and would need gitignoring.
ENV_ID="$(printf '%s' "$REPO_ROOT" | sha256sum | cut -c1-16)"
ENV_DIR="$CACHE_ROOT/$ENV_ID"
PINNED_VENV="$ENV_DIR/pinned"     # test / build helper: pinned deps, editable
SCHEMA_VENV="$ENV_DIR/schema"     # optional [schema] extra
CLEAN_VENV="$ENV_DIR/clean"       # wheel install, no source tree on the path

ALL_JOBS=(test schema-extra lint build integration-blast)
FAILED=()
PASSED=()
SKIPPED=()

log()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
step() {
    local name="$1"; shift
    printf '  %-48s' "$name"
    if "$@" >"$ENV_DIR/step.log" 2>&1; then
        printf 'PASS\n'
        return 0
    fi
    printf 'FAIL\n'
    sed 's/^/      /' "$ENV_DIR/step.log" | tail -20
    return 1
}

ensure_venv() {  # ensure_venv <dir> <pip-args...>
    local dir="$1"; shift
    if [ ! -x "$dir/bin/python" ]; then
        python3 -m venv "$dir" || return 1
        "$dir/bin/pip" install --quiet --upgrade pip || return 1
    fi
    "$dir/bin/pip" install --quiet "$@"
}

job_test() {
    local ok=0
    ensure_venv "$PINNED_VENV" -c constraints.txt -e ".[alignment]" coverage pytest || return 1
    local py="$PINNED_VENV/bin/python"
    step "doctor --json"                     "$py" -m candidate_qtl doctor --json || ok=1
    step "coverage run (core partition)"     "$py" -m coverage run scripts/run_test_partition.py --partition core || ok=1
    step "coverage report (fail_under=80)"   "$py" -m coverage report \
        --omit="src/candidate_qtl/*benchmark*.py,src/candidate_qtl/audit_freeze.py,src/candidate_qtl/commands/benchmark.py,src/candidate_qtl/challenge_aggregation.py,src/candidate_qtl/command_log_validation.py,src/candidate_qtl/corpus_masking.py,src/candidate_qtl/knowledge_coverage.py,src/candidate_qtl/lane_coverage.py,src/candidate_qtl/query_eligibility.py" || ok=1
    step "compileall src scripts"            "$py" -m compileall -q src scripts || ok=1
    step "qtlift_client.py --help"           "$py" scripts/qtlift_client.py --help || ok=1
    step "compare_qtl_sources.py --help"     "$py" scripts/compare_qtl_sources.py --help || ok=1
    step "candidate_qtl_extensions.py --help" "$py" scripts/candidate_qtl_extensions.py --help || ok=1
    step "claim_review_session.py --help" "$py" scripts/claim_review_session.py --help || ok=1
    step "generate_cli_reference.py --check" "$py" scripts/generate_cli_reference.py --check || ok=1
    step "reachability gate (source)"        "$py" scripts/reachability_gate.py --mode source --repo-root "$REPO_ROOT" || ok=1
    step "render_workflow_diagrams.py --check" "$py" scripts/render_workflow_diagrams.py --check || ok=1
    return $ok
}

job_schema_extra() {
    local ok=0
    ensure_venv "$SCHEMA_VENV" -c constraints.txt -e ".[schema]" || return 1
    local py="$SCHEMA_VENV/bin/python"
    step "pandera backend is really available" "$py" -c \
        'from candidate_qtl import schemas; assert schemas.pandera_available(), "schema extra installed but pandera_available() is False"' || ok=1
    step "unittest tests.test_schemas"        "$py" -m unittest tests.test_schemas || ok=1
    return $ok
}

job_lint() {
    ensure_venv "$PINNED_VENV" -c constraints.txt -e . ruff || return 1
    step "ruff check src scripts tests" "$PINNED_VENV/bin/ruff" check src scripts tests
}

job_build() {
    local ok=0
    ensure_venv "$PINNED_VENV" -c constraints.txt -e . build || return 1
    rm -rf "$REPO_ROOT/dist" "$CLEAN_VENV"
    step "python -m build (sdist + wheel)" "$PINNED_VENV/bin/python" -m build || return 1
    step "clean venv install from the wheel" bash -c \
        "python3 -m venv '$CLEAN_VENV' && '$CLEAN_VENV/bin/pip' install --quiet --upgrade pip && '$CLEAN_VENV/bin/pip' install --quiet -c constraints.txt '$REPO_ROOT'/dist/*.whl" || return 1
    step "candidate-qtl --version"        "$CLEAN_VENV/bin/candidate-qtl" --version || ok=1
    step "import from site-packages"      "$CLEAN_VENV/bin/python" -c \
        'import candidate_qtl; print("import ok", candidate_qtl.__version__)' || ok=1
    step "default wheel excludes pysam"   "$CLEAN_VENV/bin/python" -c \
        'import importlib.util; assert importlib.util.find_spec("pysam") is None' || ok=1
    # The default weight registry must be packaged: ranking loads it with no
    # source tree present. This guards the vendored-data packaging.
    step "packaged weight registry loads"  "$CLEAN_VENV/bin/python" -c \
        'from candidate_qtl.ranking_policy import default_registry_path, load_weight_registry; p=default_registry_path(); assert p.exists(), p; r,_=load_weight_registry(); print("registry ok", len(r), "rows from", p)' || ok=1
    step "reachability gate (installed)"  "$CLEAN_VENV/bin/python" scripts/reachability_gate.py \
        --mode installed --python "$CLEAN_VENV/bin/python" --repo-root "$REPO_ROOT" || ok=1
    return $ok
}

job_integration_blast() {
    if ! command -v blastp >/dev/null 2>&1 || ! command -v makeblastdb >/dev/null 2>&1; then
        printf '  %-48s%s\n' "blastp/makeblastdb on PATH" "SKIP (install ncbi-blast+ to run this job)"
        return 2
    fi
    ensure_venv "$PINNED_VENV" -c constraints.txt -e ".[alignment]" || return 1
    step "real local-BLAST discover smoke" "$PINNED_VENV/bin/python" -m unittest \
        tests.test_integration_blast tests.test_blast_panel_integration tests.test_claim_review_session || return 1
    step "practical packet real-BLAST handoff" "$PINNED_VENV/bin/python" \
        -m unittest discover -s workflow -p test_practical_packet_run.py
}

run_job() {
    local job="$1" rc=0
    log "$job"
    case "$job" in
        test)              job_test || rc=$? ;;
        schema-extra)      job_schema_extra || rc=$? ;;
        lint)              job_lint || rc=$? ;;
        build)             job_build || rc=$? ;;
        integration-blast) job_integration_blast || rc=$? ;;
        *) echo "unknown job: $job" >&2; return 1 ;;
    esac
    case $rc in
        0) PASSED+=("$job") ;;
        2) SKIPPED+=("$job") ;;
        *) FAILED+=("$job") ;;
    esac
}

case "${1:-}" in
    --list|-l) printf '%s\n' "${ALL_JOBS[@]}"; exit 0 ;;
    --help|-h) sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
esac

cd "$REPO_ROOT" || exit 1
mkdir -p "$ENV_DIR"
JOBS=("$@")
[ ${#JOBS[@]} -eq 0 ] && JOBS=("${ALL_JOBS[@]}")

echo "repo:  $REPO_ROOT"
echo "envs:  $ENV_DIR"
for job in "${JOBS[@]}"; do run_job "$job"; done

log "summary"
[ ${#PASSED[@]}  -gt 0 ] && echo "  passed:  ${PASSED[*]}"
[ ${#SKIPPED[@]} -gt 0 ] && echo "  skipped: ${SKIPPED[*]}"
[ ${#FAILED[@]}  -gt 0 ] && echo "  FAILED:  ${FAILED[*]}"
echo "  (the windows-wrapper job is PowerShell-only: scripts/local_ci.ps1)"
[ ${#FAILED[@]} -eq 0 ]
