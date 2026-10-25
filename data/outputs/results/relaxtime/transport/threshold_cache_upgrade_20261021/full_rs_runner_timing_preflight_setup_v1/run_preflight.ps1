param([Parameter(Mandatory=$true)][ValidateSet('setup','unit','dryrun')][string]$Stage)
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$work='D:\CodexCompute\Julia_RelaxTime\full-rs-runner-timing-preflight-20261026-v1'
$project=Join-Path $work 'project'
$julia='D:\CodexCompute\Julia_RelaxTime\transport-preflight-20261023-v1\runtime\julia-1.12.5\bin\julia.exe'
$package=Get-Content -LiteralPath (Join-Path $work 'package_manifest.json') -Raw | ConvertFrom-Json
$statusPath=Join-Path $work ($Stage+'.status.json')
$started=[DateTime]::UtcNow
$arguments=@()
$verified=$false
if(Test-Path -LiteralPath $statusPath){throw "Stage already exists: $Stage"}

function Write-Status([bool]$finished,[int]$exitCode,[string]$errorText='') {
    [ordered]@{stage=$Stage;finished=$finished;exit_code=$exitCode;error=$errorText;
        started_utc=$started.ToString('o');updated_utc=[DateTime]::UtcNow.ToString('o');
        arguments=$script:arguments;source_verified_before_and_after=$script:verified} |
        ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusPath -Encoding UTF8
}

function Test-Source {
    if((Get-FileHash -LiteralPath $package.base_source_zip -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.base_source_archive_sha256){throw 'Base source archive changed'}
    if((Get-FileHash -LiteralPath (Join-Path $work 'overlay.zip') -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.overlay_archive_sha256){throw 'Timing overlay archive changed'}
    foreach($file in $package.source_files.PSObject.Properties){
        if((Get-FileHash -LiteralPath (Join-Path $project $file.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $file.Value){throw "Source changed: $($file.Name)"}
    }
}

try {
    Write-Status $false 0
    if($env:COMPUTERNAME -ne 'W'){throw 'Unexpected host'}
    if($Stage -eq 'setup'){
        if(Test-Path -LiteralPath $project){throw 'Use a new preflight project'}
        if((Get-FileHash -LiteralPath $package.base_source_zip -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.base_source_archive_sha256){throw 'Base source archive mismatch'}
        if((Get-FileHash -LiteralPath (Join-Path $work 'overlay.zip') -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.overlay_archive_sha256){throw 'Overlay mismatch'}
        Expand-Archive -LiteralPath $package.base_source_zip -DestinationPath $project
        Expand-Archive -LiteralPath (Join-Path $work 'overlay.zip') -DestinationPath $project -Force
        Test-Source
        $verified=$true
    } else {
        if(@(Get-Process -Name julia -ErrorAction SilentlyContinue).Count -gt 0){throw 'Another Julia process is running; wait for the qualification queue'}
        Test-Source
        $env:JULIA_NUM_THREADS='1'
        $env:OPENBLAS_NUM_THREADS='1'
        $env:JULIA_PKG_PRECOMPILE_AUTO='0'
        Remove-Item Env:THRESHOLD_SOURCE_MANIFEST -ErrorAction SilentlyContinue
        $arguments=@('--startup-file=no','--threads=1',"--project=$project",'--heap-size-hint=1200M')
        if($Stage -eq 'unit'){
            $env:UNIT_FILES='relaxtime/test_transport_runner_timing.jl'
            $arguments += (Join-Path $project 'tests\unit\runtests.jl')
        } else {
            $env:TRANSPORT_RUNNER_TIMING_PATH=Join-Path $work 'runtime_timing.toml'
            $arguments += @((Join-Path $project 'benchmark\relaxtime\transport_runner_timing.jl'),
                '--mode','fixed-T-sparse-muB','--T-list','120','--muB-list','0',
                '--xi-list','-0.5,0.5','--outdir',(Join-Path $work 'dry_run'),
                '--case-name','runner_timing_preflight_20261026','--dry-run','--channel-diagnostics')
        }
        $process=Start-Process -FilePath $julia -ArgumentList $arguments -WorkingDirectory $project -WindowStyle Hidden -Wait -PassThru `
            -RedirectStandardOutput (Join-Path $work ($Stage+'.stdout.log')) -RedirectStandardError (Join-Path $work ($Stage+'.stderr.log'))
        Test-Source
        $verified=$true
        if($process.ExitCode -ne 0){throw "$Stage exited $($process.ExitCode)"}
    }
    Write-Status $true 0
    exit 0
} catch {
    Write-Status $true 1 $_.Exception.Message
    throw
}
