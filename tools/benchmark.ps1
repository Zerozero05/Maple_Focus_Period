# Each timed Compute runs once in a fresh Maple process. No reference is computed
# until after the cold timer stops. This script writes no persistent files.
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$MaplePath,
    [ValidateRange(1, 100)]
    [int]$Repeats = 3
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $MaplePath -PathType Leaf)) {
    throw 'MaplePath must identify an existing Maple command-line executable.'
}
$mapleExecutable = (Resolve-Path -LiteralPath $MaplePath).Path
$repoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$sourcePath = Join-Path $repoRoot 'src/focus_period.mpl'
$sourceHash = (Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash
$mapleRepoRoot = $repoRoot.Replace('\', '/')
$culture = [System.Globalization.CultureInfo]::InvariantCulture

$probe = @'
restart:
currentdir("__REPO_ROOT__"):
read "src/focus_period.mpl":
variant := "__VARIANT__":
X := y+x^2+x*y+2*y^2+x^3:
Y := -x+x*y+y^2+x*y^2:
started_cpu := time(): started_wall := time[real]():
if variant="quotient_dense" then
    sample := FocusPeriod:-Compute(X,Y,x,y,7,6,method="quotient",
                                  optimize=false,factor_orbits=false):
elif variant="direct_dense" then
    sample := FocusPeriod:-Compute(X,Y,x,y,7,6,
                                  optimize=false,factor_orbits=false):
elif variant="direct_no_factor" then
    sample := FocusPeriod:-Compute(X,Y,x,y,7,6,factor_orbits=false):
else
    sample := FocusPeriod:-Compute(X,Y,x,y,7,6):
end if:
elapsed_cpu := time()-started_cpu:
elapsed_wall := time[real]()-started_wall:
printf("COLD,%s,%.6f,%.6f\n",variant,elapsed_cpu,elapsed_wall):
printf("STRUCTURE,%s,factor=%a,u=%a,p=%a\n",variant,sample["orbit_factor"],
       sample["orbit_structure"]["u_possible_orders"],
       sample["structure"]["period_possible_orders"]):
reference := FocusPeriod:-Compute(X,Y,x,y,7,6,method="quotient",
                                  optimize=false,factor_orbits=false):
for n from 2 to 7 do
    if not evalb(simplify(expand(sample["u"][n]-reference["u"][n]),trig)=0)
       or not evalb(simplify(expand(sample["du"][n]-reference["du"][n]),trig)=0)
       or not evalb(simplify(expand(sample["g"][n]-reference["g"][n]),trig)=0) then
        error "Cross-check mismatch in u/du/g order %1",n;
    end if;
end do:
for n from 0 to 6 do
    if not evalb(simplify(expand(sample["B"][n]-reference["B"][n]),trig)=0)
       or not evalb(simplify(expand(sample["p"][n]-reference["p"][n]),trig)=0) then
        error "Cross-check mismatch in B/p order %1",n;
    end if;
end do:
printf("PASS,%s,cold timer preceded cross-checks\n",variant):
'@

function Get-Median([double[]]$Values) {
    $ordered = @($Values | Sort-Object)
    $middle = [int][Math]::Floor($ordered.Count / 2)
    if (($ordered.Count % 2) -eq 1) { return $ordered[$middle] }
    return ($ordered[$middle - 1] + $ordered[$middle]) / 2
}

$configurations = @('quotient_dense', 'direct_dense', 'direct_default', 'direct_no_factor')
$measurements = @()
Write-Output 'model=mixed_quadratic_cubic,Ng=7,Np=6,cold_process=true'
Write-Output "source_sha256=$sourceHash"
Write-Output 'repeat,position,configuration,cpu_s,wall_s,verification'
for ($repeat = 1; $repeat -le $Repeats; $repeat++) {
    $sequence = @($configurations)
    if (($repeat % 2) -eq 0) { [Array]::Reverse($sequence) }
    $offset = [int][Math]::Floor(($repeat - 1) / 2) % $sequence.Count
    for ($position = 0; $position -lt $sequence.Count; $position++) {
        $variant = $sequence[($position + $offset) % $sequence.Count]
        $code = $probe.Replace('__REPO_ROOT__', $mapleRepoRoot).Replace('__VARIANT__', $variant)
        # A new native process is launched for every iteration. Exit 0 alone is
        # insufficient: cmaple may print an Error and still return success.
        $rawOutput = @($code | & $mapleExecutable -q 2>&1)
        $exitCode = $LASTEXITCODE
        $output = ($rawOutput | ForEach-Object { [string]$_ }) -join "`n"
        $timing = [regex]::Match($output, "(?m)^COLD,$variant,([0-9.]+),([0-9.]+)\r?$")
        if ($exitCode -ne 0 -or $output -match '\bError\b' -or -not $timing.Success -or
            $output -notmatch "(?m)^PASS,$variant,cold timer preceded cross-checks\r?$") {
            throw "Maple benchmark failed for $variant (exit $exitCode):`n$output"
        }
        $cpuSeconds = [double]::Parse($timing.Groups[1].Value, $culture)
        $wallSeconds = [double]::Parse($timing.Groups[2].Value, $culture)
        $measurements += [pscustomobject]@{
            Configuration = $variant; CpuSeconds = $cpuSeconds; WallSeconds = $wallSeconds
        }
        Write-Output ([string]::Format($culture, '{0},{1},{2},{3:F6},{4:F6},PASS',
            $repeat, ($position + 1), $variant, $cpuSeconds, $wallSeconds))
        Write-Output ([regex]::Match($output, '(?m)^STRUCTURE,.*$').Value.Trim())
    }
}
if ((Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash -ne $sourceHash) {
    throw 'Source changed during measurement; results do not describe one fixed candidate.'
}
Write-Output 'configuration,repeats,median_cpu_s,median_wall_s'
foreach ($variant in $configurations) {
    $samples = @($measurements | Where-Object { $_.Configuration -eq $variant })
    $medianCpu = Get-Median -Values @($samples | ForEach-Object { $_.CpuSeconds })
    $medianWall = Get-Median -Values @($samples | ForEach-Object { $_.WallSeconds })
    Write-Output ([string]::Format($culture, '{0},{1},{2:F6},{3:F6}',
        $variant, $Repeats, $medianCpu, $medianWall))
}
Write-Output 'BENCHMARK PASS: each cold measurement used a fresh process; all u/du/g/B/p checks passed.'
