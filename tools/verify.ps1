# Maple可能报错后仍返回0；同时核对输出和PASS标志。
# 示例：.\tools\verify.ps1 -MaplePath 'D:\Maple\bin.X86_64_WINDOWS\cmaple.exe'
[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$MaplePath)

$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $MaplePath -PathType Leaf)) {
    throw 'MaplePath must name an existing cmaple executable.'
}
$repoPath = Split-Path -Parent $PSScriptRoot
$testFiles = @('verify_maple.mpl','verify_deep_core.mpl','verify_foundation.mpl',
               'verify_research.mpl','verify_cartesian.mpl','verify_integration.mpl')
$totalChecks = 0
Push-Location -LiteralPath $repoPath
try {
    foreach ($testFile in $testFiles) {
        $output = & $MaplePath -q (Join-Path 'tests' $testFile) 2>&1
        $exitCode = $LASTEXITCODE
        $logText = $output -join "`n"
        Write-Output $logText
        $pass = [regex]::Match($logText,'PASS:\s*(\d+)')
        if ($exitCode -ne 0 -or $logText -match '(?mi)^\s*(Error|Warning),' -or -not $pass.Success) {
            throw "Maple verification failed: $testFile (exit $exitCode)."
        }
        $totalChecks += [int]$pass.Groups[1].Value
    }
    Write-Output "All $($testFiles.Count) Maple suites passed: $totalChecks checks."
} finally {
    Pop-Location
}
