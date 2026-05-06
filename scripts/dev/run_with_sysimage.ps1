param(
    [switch]$BuildIfMissing,
    [ValidateSet("fallback", "strict", "rebuild")]
    [string]$MismatchPolicy = "rebuild",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$JuliaArgs
)

$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
$SysimagePath = Join-Path $RepoRoot "build\\JuliaRelaxTime.dll"
$MetaPath = Join-Path $RepoRoot "build\\JuliaRelaxTime.sysimage.json"

if ($BuildIfMissing) {
    if ($MismatchPolicy -eq "strict") {
        throw "-BuildIfMissing cannot be combined with -MismatchPolicy strict"
    }
    $MismatchPolicy = "rebuild"
}

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

function Get-PlatformFamily {
    if ($IsWindows -or $env:OS -eq "Windows_NT") {
        return "windows"
    }
    if ($IsMacOS) {
        return "macos"
    }
    if ($IsLinux) {
        return "linux"
    }
    return "unknown"
}

function Get-PlatformArch {
    $arch = [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture.ToString().ToLowerInvariant()
    switch ($arch) {
        "x64" { return "x86_64" }
        "arm64" { return "aarch64" }
        default { return $arch }
    }
}

function Build-LocalSysimage {
    & julia --project="$RepoRoot" (Join-Path $RepoRoot "scripts\\dev\\build_sysimage.jl")
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to build sysimage"
    }
}

function Get-HeadCommit {
    try {
        $head = (& git -C $RepoRoot rev-parse HEAD 2>$null)
        if ($LASTEXITCODE -eq 0 -and -not [string]::IsNullOrWhiteSpace($head)) {
            return $head.Trim()
        }
    } catch {
    }
    return $null
}

function Get-SysimageCompatibility {
    param(
        $Meta,
        [string]$CurrentVersion,
        [string]$CurrentFamily,
        [string]$CurrentArch,
        [string]$HeadCommit
    )

    if (-not $Meta.julia_version) {
        return @{ Compatible = $false; Reason = "metadata missing julia_version" }
    }
    if ($Meta.julia_version -ne $CurrentVersion) {
        return @{ Compatible = $false; Reason = "Julia version $($Meta.julia_version) does not match current Julia $CurrentVersion" }
    }
    if ($Meta.platform_family -and $Meta.platform_family -ne $CurrentFamily) {
        return @{ Compatible = $false; Reason = "platform family $($Meta.platform_family) does not match current platform $CurrentFamily" }
    }
    if ($Meta.platform_arch -and $Meta.platform_arch -ne $CurrentArch) {
        return @{ Compatible = $false; Reason = "platform arch $($Meta.platform_arch) does not match current arch $CurrentArch" }
    }
    if (-not [string]::IsNullOrWhiteSpace($HeadCommit)) {
        if (-not $Meta.git_commit) {
            return @{ Compatible = $false; Reason = "metadata missing git_commit" }
        }
        if ($Meta.git_commit -ne $HeadCommit) {
            return @{ Compatible = $false; Reason = "sysimage git commit $($Meta.git_commit) does not match current HEAD $HeadCommit" }
        }
    }
    return @{ Compatible = $true; Reason = "ok" }
}

$currentVersion = Get-JuliaVersion
$currentFamily = Get-PlatformFamily
$currentArch = Get-PlatformArch
$headCommit = Get-HeadCommit

if (-not ((Test-Path $SysimagePath) -and (Test-Path $MetaPath))) {
    switch ($MismatchPolicy) {
        "strict" {
            throw "Sysimage or metadata missing."
        }
        "rebuild" {
            Build-LocalSysimage
        }
        default {
            Write-Warning "Sysimage or metadata missing; falling back to plain julia --project=."
        }
    }
}

$UseSysimage = $false
if ((Test-Path $SysimagePath) -and (Test-Path $MetaPath)) {
    $meta = Get-Content $MetaPath | ConvertFrom-Json
    $compat = Get-SysimageCompatibility -Meta $meta -CurrentVersion $currentVersion -CurrentFamily $currentFamily -CurrentArch $currentArch -HeadCommit $headCommit
    if ($compat.Compatible) {
        $UseSysimage = $true
    } else {
        switch ($MismatchPolicy) {
            "strict" {
                throw "Incompatible sysimage: $($compat.Reason)"
            }
            "rebuild" {
                Write-Warning "Incompatible sysimage detected; rebuilding local sysimage. Reason: $($compat.Reason)"
                Build-LocalSysimage
                $meta = Get-Content $MetaPath | ConvertFrom-Json
                $compat = Get-SysimageCompatibility -Meta $meta -CurrentVersion $currentVersion -CurrentFamily $currentFamily -CurrentArch $currentArch -HeadCommit $headCommit
                if (-not $compat.Compatible) {
                    throw "Rebuilt sysimage is still incompatible: $($compat.Reason)"
                }
                $UseSysimage = $true
            }
            default {
                Write-Warning "Incompatible sysimage detected; falling back to plain julia --project=. Reason: $($compat.Reason)"
            }
        }
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
