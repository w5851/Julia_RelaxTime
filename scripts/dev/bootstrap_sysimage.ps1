param(
    [string]$Repo = "w5851/Julia_RelaxTime",
    [string]$ReleaseTag = "latest",
    [string]$OutputDir = "",
    [switch]$Overwrite,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $RepoRoot "build"
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

function Get-SysimageExtension {
    param([string]$Family)
    switch ($Family) {
        "windows" { return "dll" }
        "macos" { return "dylib" }
        "linux" { return "so" }
        default { throw "Unsupported platform family for sysimage bootstrap: $Family" }
    }
}

function Get-ArchiveFormat {
    param([string]$Family)
    if ($Family -eq "windows") {
        return "zip"
    }
    return "tar.gz"
}

function Get-AssetName {
    param(
        [string]$Family,
        [string]$Arch,
        [string]$JuliaVersion,
        [string]$ArchiveFormat
    )
    return "JuliaRelaxTime-sysimage-$Family-$Arch-julia$JuliaVersion.$ArchiveFormat"
}

$family = Get-PlatformFamily
$arch = Get-PlatformArch
$juliaVersion = Get-JuliaVersion
$archiveFormat = Get-ArchiveFormat -Family $family
$sysimageExt = Get-SysimageExtension -Family $family
$assetName = Get-AssetName -Family $family -Arch $arch -JuliaVersion $juliaVersion -ArchiveFormat $archiveFormat

if ($ReleaseTag -eq "latest") {
    $assetUrl = "https://github.com/$Repo/releases/latest/download/$assetName"
} else {
    $assetUrl = "https://github.com/$Repo/releases/download/$ReleaseTag/$assetName"
}

$sysimageTarget = Join-Path $OutputDir "JuliaRelaxTime.$sysimageExt"
$metaTarget = Join-Path $OutputDir "JuliaRelaxTime.sysimage.json"

if (((Test-Path $sysimageTarget) -or (Test-Path $metaTarget)) -and -not $Overwrite) {
    throw "Target sysimage or metadata already exists in $OutputDir. Re-run with -Overwrite to replace it."
}

if ($DryRun) {
    Write-Host "Resolved asset: $assetName"
    Write-Host "Download URL: $assetUrl"
    Write-Host "Output dir: $OutputDir"
    exit 0
}

$null = New-Item -ItemType Directory -Force -Path $OutputDir
$tempRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("JuliaRelaxTime-sysimage-" + [System.Guid]::NewGuid().ToString("N"))
$null = New-Item -ItemType Directory -Force -Path $tempRoot
$archivePath = Join-Path $tempRoot $assetName
$extractDir = Join-Path $tempRoot "extract"
$null = New-Item -ItemType Directory -Force -Path $extractDir

try {
    Invoke-WebRequest -Uri $assetUrl -OutFile $archivePath

    if ($archiveFormat -eq "zip") {
        Expand-Archive -LiteralPath $archivePath -DestinationPath $extractDir -Force
    } else {
        & tar -xzf $archivePath -C $extractDir
        if ($LASTEXITCODE -ne 0) {
            throw "Failed to extract archive: $archivePath"
        }
    }

    $sysimageFile = Get-ChildItem -Path $extractDir -Recurse -File -Filter ("JuliaRelaxTime." + $sysimageExt) | Select-Object -First 1
    $metaFile = Get-ChildItem -Path $extractDir -Recurse -File -Filter "JuliaRelaxTime.sysimage.json" | Select-Object -First 1
    if (-not $sysimageFile -or -not $metaFile) {
        throw "Release asset does not contain JuliaRelaxTime.$sysimageExt and JuliaRelaxTime.sysimage.json"
    }

    Copy-Item -LiteralPath $sysimageFile.FullName -Destination $sysimageTarget -Force
    Copy-Item -LiteralPath $metaFile.FullName -Destination $metaTarget -Force

    Write-Host "Bootstrapped sysimage: $sysimageTarget"
    Write-Host "Bootstrapped metadata: $metaTarget"
} finally {
    if (Test-Path $tempRoot) {
        Remove-Item -LiteralPath $tempRoot -Recurse -Force
    }
}
