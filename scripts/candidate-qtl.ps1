[CmdletBinding()]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$CliArgs
)

$ErrorActionPreference = "Stop"
$RepositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
if ($RepositoryRoot -notmatch '^(?<drive>[A-Za-z]):\\(?<tail>.*)$') {
    throw "The repository must be on a drive mounted by WSL: $RepositoryRoot"
}
$LinuxRoot = "/mnt/$($Matches.drive.ToLower())/" + $Matches.tail.Replace('\', '/')
$PrimerblastRoot = Join-Path (Split-Path $RepositoryRoot -Parent) "primerblast-oss"

function Quote-Bash([string]$Value) {
    return "'" + $Value.Replace("'", "'""'""'") + "'"
}

$QuotedArguments = ($CliArgs | ForEach-Object { Quote-Bash $_ }) -join " "
$EnvironmentPrefix = "PYTHONPATH=$(Quote-Bash "$LinuxRoot/src")"
if (Test-Path -LiteralPath $PrimerblastRoot) {
    $PrimerblastLinuxRoot = "/mnt/c/" + $PrimerblastRoot.Substring(3).Replace('\', '/')
    $EnvironmentPrefix = "PRIMERBLAST_OSS_ROOT=$(Quote-Bash $PrimerblastLinuxRoot) $EnvironmentPrefix"
}
$Command = "cd $(Quote-Bash $LinuxRoot) && $EnvironmentPrefix python3 -m candidate_qtl $QuotedArguments"
wsl.exe -e bash -lc $Command
exit $LASTEXITCODE
