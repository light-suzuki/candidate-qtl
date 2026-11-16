[CmdletBinding()]
param(
    [switch]$AddToUserPath
)

$ErrorActionPreference = "Stop"
$RepositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$BinDirectory = Join-Path $env:USERPROFILE ".local\bin"
$WrapperPath = Join-Path $BinDirectory "candidate-qtl.ps1"
New-Item -ItemType Directory -Path $BinDirectory -Force | Out-Null

$Wrapper = @'
& '{0}\scripts\candidate-qtl.ps1' @args
exit $LASTEXITCODE
'@ -f $RepositoryRoot
Set-Content -LiteralPath $WrapperPath -Value $Wrapper -Encoding utf8NoBOM

if ($AddToUserPath) {
    $UserPath = [Environment]::GetEnvironmentVariable("Path", "User")
    $Parts = @($UserPath -split ";" | Where-Object { $_ })
    if ($Parts -notcontains $BinDirectory) {
        [Environment]::SetEnvironmentVariable("Path", (($Parts + $BinDirectory) -join ";"), "User")
    }
    if (($env:Path -split ";") -notcontains $BinDirectory) {
        $env:Path = "$BinDirectory;$env:Path"
    }
}

Write-Output "Installed: $WrapperPath"
Write-Output "Run: candidate-qtl doctor --json"
