$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$transportWorkspace = 'D:\CodexCompute\Julia_RelaxTime\transport-preflight-20261023-v1'
$projectDirectory = Join-Path $transportWorkspace 'project'
$evidenceDirectory = Join-Path $transportWorkspace 'evidence'
$runtimeDirectory = Join-Path $transportWorkspace 'runtime'
$downloadDirectory = Join-Path $transportWorkspace 'downloads'
$statusPath = Join-Path $transportWorkspace 'remote_status.json'
$sourceArchive = Join-Path $transportWorkspace 'transport_preflight_source_v1.zip'
$sourceManifest = Get-Content -LiteralPath (Join-Path $transportWorkspace 'package_manifest.json') -Raw | ConvertFrom-Json
$runtimeArchive = Join-Path $downloadDirectory 'julia-1.12.5-win64.zip'
$runtimeExpectedHash = '591b2053086215260be8b50071b3cb3d73fd0ac0a46d259c5d8eee77d4ed99a7'
$runtimeExpectedBytes = 273077393
$runtimeUrl = 'https://julialang-s3.julialang.org/bin/winnt/x64/1.12/julia-1.12.5-win64.zip'
$portableJulia = Join-Path $runtimeDirectory 'julia-1.12.5\bin\julia.exe'
$remoteRunStart = [DateTime]::UtcNow
$remoteStage = 'created'

function Write-RemoteStatus([string]$stage, [bool]$finished = $false, [object]$exitCode = $null, [string]$errorMessage = '') {
    $script:remoteStage = $stage
    [ordered]@{
        host = $env:COMPUTERNAME
        identity = [Security.Principal.WindowsIdentity]::GetCurrent().Name
        workspace = $transportWorkspace
        stage = $stage
        finished = $finished
        exit_code = $exitCode
        started_utc = $remoteRunStart.ToString('o')
        updated_utc = [DateTime]::UtcNow.ToString('o')
        runtime_version = '1.12.5'
        source_archive_sha256 = $sourceManifest.zip_sha256
        provenance = 'immutable dirty-worktree source archive; no calculation commit'
        error = $errorMessage
    } | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding UTF8
}

function Test-SourceFiles {
    foreach ($entry in $sourceManifest.source_files.PSObject.Properties) {
        $filePath = Join-Path $projectDirectory $entry.Name
        if (-not (Test-Path -LiteralPath $filePath -PathType Leaf)) { throw "Source file missing: $($entry.Name)" }
        $actualHash = (Get-FileHash -LiteralPath $filePath -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($actualHash -ne $entry.Value) { throw "Source identity mismatch: $($entry.Name)" }
    }
}

try {
    foreach ($directory in @($evidenceDirectory, $runtimeDirectory, $downloadDirectory)) {
        New-Item -ItemType Directory -Path $directory -Force | Out-Null
    }
    Write-RemoteStatus 'verify_source_archive'
    $sourceHash = (Get-FileHash -LiteralPath $sourceArchive -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($sourceHash -ne $sourceManifest.zip_sha256) { throw 'Source archive SHA256 mismatch' }
    if (Test-Path -LiteralPath $projectDirectory) { throw 'Isolated project directory already exists; refuse to overwrite' }
    Expand-Archive -LiteralPath $sourceArchive -DestinationPath $projectDirectory
    Test-SourceFiles
    [ordered]@{ files_checked = $sourceManifest.files; passed = $true; source_archive_sha256 = $sourceHash } |
        ConvertTo-Json | Set-Content -LiteralPath (Join-Path $evidenceDirectory 'source_verification_before.json') -Encoding UTF8

    Write-RemoteStatus 'download_portable_runtime'
    $versions = Invoke-RestMethod -Uri 'https://julialang-s3.julialang.org/bin/versions.json' -TimeoutSec 60
    $official = @($versions.'1.12.5'.files | Where-Object { $_.url -eq $runtimeUrl })
    if ($official.Count -ne 1 -or $official[0].sha256 -ne $runtimeExpectedHash -or $official[0].size -ne $runtimeExpectedBytes) {
        throw 'Portable runtime metadata differs from official pinned release'
    }
    Invoke-WebRequest -UseBasicParsing -Uri $runtimeUrl -OutFile $runtimeArchive -TimeoutSec 1200
    $runtimeActualHash = (Get-FileHash -LiteralPath $runtimeArchive -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($runtimeActualHash -ne $runtimeExpectedHash -or (Get-Item -LiteralPath $runtimeArchive).Length -ne $runtimeExpectedBytes) {
        throw 'Portable Julia archive verification failed'
    }
    [ordered]@{ version = '1.12.5'; url = $runtimeUrl; sha256 = $runtimeActualHash; bytes = $runtimeExpectedBytes; official_versions_verified = $true } |
        ConvertTo-Json | Set-Content -LiteralPath (Join-Path $evidenceDirectory 'runtime_metadata.json') -Encoding UTF8
    Write-RemoteStatus 'extract_portable_runtime'
    Expand-Archive -LiteralPath $runtimeArchive -DestinationPath $runtimeDirectory
    if (-not (Test-Path -LiteralPath $portableJulia -PathType Leaf)) { throw 'Portable Julia executable missing after extraction' }
    $versionOutput = & $portableJulia --startup-file=no --version
    if ($LASTEXITCODE -ne 0 -or $versionOutput -notmatch '1\.12\.5') { throw "Unexpected Julia runtime: $versionOutput" }

    $env:JULIA_NUM_THREADS = '1'
    $env:OPENBLAS_NUM_THREADS = '1'
    $env:JULIA_PKG_PRECOMPILE_AUTO = '0'
    $taskTempDirectory = Join-Path $evidenceDirectory 'tmp'
    New-Item -ItemType Directory -Path $taskTempDirectory -Force | Out-Null
    $env:TEMP = $taskTempDirectory
    $env:TMP = $taskTempDirectory
    Write-RemoteStatus 'instantiate_project'
    $dependencyProcess = Start-Process -FilePath $portableJulia -ArgumentList @(
        '--startup-file=no', '--threads=1', '--project=.', '--heap-size-hint=1200M',
        '-e', '"using Pkg; Pkg.instantiate(;allow_autoprecomp=false)"'
    ) -WorkingDirectory $projectDirectory -WindowStyle Hidden -PassThru -Wait `
      -RedirectStandardOutput (Join-Path $evidenceDirectory 'instantiate.stdout.log') `
      -RedirectStandardError (Join-Path $evidenceDirectory 'instantiate.stderr.log')
    if ($dependencyProcess.ExitCode -ne 0) { throw "Pkg.instantiate exit $($dependencyProcess.ExitCode)" }
    Test-SourceFiles
    Write-RemoteStatus 'manual_workflow_integration'
    $testStart = [DateTime]::UtcNow
    $testProcess = Start-Process -FilePath $portableJulia -ArgumentList @(
        '--startup-file=no', '--threads=1', '--project=.', '--heap-size-hint=1200M',
        (Join-Path $transportWorkspace 'run_manual_preflight.jl'), $projectDirectory, $evidenceDirectory
    ) -WorkingDirectory $projectDirectory -WindowStyle Hidden -PassThru -Wait `
      -RedirectStandardOutput (Join-Path $evidenceDirectory 'manual_integration.stdout.log') `
      -RedirectStandardError (Join-Path $evidenceDirectory 'manual_integration.stderr.log')
    [ordered]@{
        exit_code = $testProcess.ExitCode
        elapsed_seconds = ([DateTime]::UtcNow - $testStart).TotalSeconds
        host = $env:COMPUTERNAME
        test_scope = 'small custom-profile workflow integration; not GitHub runner production timing'
    } | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $evidenceDirectory 'manual_process_result.json') -Encoding UTF8
    Test-SourceFiles
    [ordered]@{ files_checked = $sourceManifest.files; passed = $true; source_archive_sha256 = $sourceHash } |
        ConvertTo-Json | Set-Content -LiteralPath (Join-Path $evidenceDirectory 'source_verification_after.json') -Encoding UTF8
    if ($testProcess.ExitCode -ne 0) { throw "Manual workflow integration exit $($testProcess.ExitCode)" }
    Write-RemoteStatus 'complete' $true 0
} catch {
    Write-RemoteStatus $remoteStage $true 1 $_.Exception.Message
    $_ | Out-String | Set-Content -LiteralPath (Join-Path $transportWorkspace 'setup_error.txt') -Encoding UTF8
    exit 1
}
