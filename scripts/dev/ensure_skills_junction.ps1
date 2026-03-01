param(
    [ValidateSet('check', 'fix')]
    [string]$Mode = 'check',

    [string]$RepoRoot = (Get-Location).Path,

    [string]$GlobalSkillsPath = "$env:USERPROFILE\.agents\skills"
)

$ErrorActionPreference = 'Stop'

$localSkillsPath = Join-Path $RepoRoot '.agents\skills'

function Get-SkillsState {
    param(
        [string]$LocalPath,
        [string]$ExpectedTarget
    )

    if (-not (Test-Path $LocalPath)) {
        return [pscustomobject]@{
            Exists = $false
            IsJunction = $false
            Target = $null
            IsHealthy = $false
        }
    }

    $item = Get-Item -LiteralPath $LocalPath
    $isJunction = $item.LinkType -eq 'Junction'
    $target = $null

    if ($isJunction) {
        $target = ($item.Target | Select-Object -First 1)
    }

    $expectedNormalized = [System.IO.Path]::GetFullPath($ExpectedTarget)
    $targetNormalized = if ($target) { [System.IO.Path]::GetFullPath([string]$target) } else { $null }

    return [pscustomobject]@{
        Exists = $true
        IsJunction = $isJunction
        Target = $target
        IsHealthy = ($isJunction -and ($targetNormalized -eq $expectedNormalized))
    }
}

if (-not (Test-Path -LiteralPath $GlobalSkillsPath)) {
    if ($Mode -eq 'check') {
        Write-Error "Global skills path not found: $GlobalSkillsPath"
        exit 2
    }

    New-Item -ItemType Directory -Path $GlobalSkillsPath -Force | Out-Null
}

$state = Get-SkillsState -LocalPath $localSkillsPath -ExpectedTarget $GlobalSkillsPath

if ($Mode -eq 'check') {
    if ($state.IsHealthy) {
        Write-Host "OK: $localSkillsPath is a Junction -> $GlobalSkillsPath"
        exit 0
    }

    if (-not $state.Exists) {
        Write-Host "NOT_OK: $localSkillsPath does not exist"
    } elseif (-not $state.IsJunction) {
        Write-Host "NOT_OK: $localSkillsPath exists but is not a Junction"
    } else {
        Write-Host "NOT_OK: $localSkillsPath points to '$($state.Target)' instead of '$GlobalSkillsPath'"
    }

    exit 1
}

if ($state.IsHealthy) {
    Write-Host "OK: already healthy Junction -> $GlobalSkillsPath"
    exit 0
}

if (Test-Path -LiteralPath $localSkillsPath) {
    $backupPath = Join-Path (Join-Path $RepoRoot '.agents') ("skills.bak-" + (Get-Date -Format 'yyyyMMdd-HHmmss'))
    Rename-Item -LiteralPath $localSkillsPath -NewName (Split-Path $backupPath -Leaf)
    Write-Host "Backup created: $backupPath"
}

New-Item -ItemType Junction -Path $localSkillsPath -Target $GlobalSkillsPath | Out-Null

$final = Get-SkillsState -LocalPath $localSkillsPath -ExpectedTarget $GlobalSkillsPath
if (-not $final.IsHealthy) {
    Write-Error "Failed to create healthy Junction at $localSkillsPath"
    exit 3
}

Write-Host "FIXED: $localSkillsPath is now Junction -> $GlobalSkillsPath"
exit 0
