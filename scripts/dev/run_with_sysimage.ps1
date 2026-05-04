param(
    [switch]$BuildIfMissing,
    [switch]$Strict,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$JuliaArgs
)

$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
$SysimagePath = Join-Path $RepoRoot "build\\JuliaRelaxTime.dll"
$MetaPath = Join-Path $RepoRoot "build\\JuliaRelaxTime.sysimage.json"

function Get-JuliaVersion {
    $versionLine = & julia --version
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to run julia --version"
    }
    if ($versionLine -match "([0-9]+\.[0-9]+\.[0-9]+)") {
        return $Matches[1]
    }
    throw "Could not parse Julia version from: $versionLine"
}

function Ensure-Sysimage {
    if ((Test-Path $SysimagePath) -and (Test-Path $MetaPath)) {
        return
    }
    if (-not $BuildIfMissing) {
        if ($Strict) {
            throw "Sysimage or metadata missing. Re-run with -BuildIfMissing."
        }
        return
    }
    & julia --project="$RepoRoot" (Join-Path $RepoRoot "scripts\\dev\\build_sysimage.jl")
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to build sysimage"
    }
}

Ensure-Sysimage

$UseSysimage = $false
if ((Test-Path $SysimagePath) -and (Test-Path $MetaPath)) {
    $meta = Get-Content $MetaPath | ConvertFrom-Json
    $currentVersion = Get-JuliaVersion
    if ($meta.julia_version -eq $currentVersion) {
        $UseSysimage = $true
    } elseif ($Strict) {
        throw "Sysimage Julia version $($meta.julia_version) does not match current Julia $currentVersion"
    }
}

$cmd = @()
if ($UseSysimage) {
    $cmd += "--sysimage=$SysimagePath"
}
$cmd += "--project=$RepoRoot"
$cmd += $JuliaArgs

& julia @cmd
exit $LASTEXITCODE
