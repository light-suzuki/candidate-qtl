# The windows-wrapper job from .github/workflows/ci.yml, run locally.
#
# On Windows the Python package is not installable (pysam is Linux/macOS only)
# and the runtime is WSL, so CI validates the Windows-native entry point -- the
# PowerShell wrapper -- rather than a native package install. This script is
# that job; scripts/local_ci.sh covers every Linux job.
#
#   pwsh -File scripts/local_ci.ps1

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$wrapper = Join-Path $repoRoot 'scripts/candidate-qtl.ps1'

if (-not (Test-Path $wrapper)) {
    Write-Error "wrapper not found: $wrapper"
    exit 1
}

$tokens = $null
$errors = $null
[System.Management.Automation.Language.Parser]::ParseFile((Resolve-Path $wrapper), [ref]$tokens, [ref]$errors) | Out-Null

if ($errors.Count -gt 0) {
    $errors | ForEach-Object { Write-Error $_ }
    exit 1
}

Write-Host 'windows-wrapper: candidate-qtl.ps1 parses cleanly'
