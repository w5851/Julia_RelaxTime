param([Parameter(Mandatory=$true)][ValidateSet('setup','selection','unit','cheap','kernel','outer','joint','production')][string]$Stage)
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$work='D:\CodexCompute\Julia_RelaxTime\full-rs-mode-b-qualification-20261025-v1'
$project=Join-Path $work 'project'
$julia='D:\CodexCompute\Julia_RelaxTime\transport-preflight-20261023-v1\runtime\julia-1.12.5\bin\julia.exe'
$python='C:\Users\82306\AppData\Local\Programs\Python\Python313\python.exe'
$manifestPath=Join-Path $work 'package_manifest.json'
$sourceArchive=Join-Path $work 'source.zip'
$statusPath=Join-Path $work ($Stage+'.status.json')
$started=[DateTime]::UtcNow
$runs=@()
if (Test-Path -LiteralPath $statusPath) { throw "Stage already started: $Stage" }
function Write-Status([bool]$finished,[int]$exitCode,[string]$message='') {
    [ordered]@{host=$env:COMPUTERNAME;stage=$Stage;finished=$finished;exit_code=$exitCode;error=$message;
        started_utc=$started.ToString('o');updated_utc=[DateTime]::UtcNow.ToString('o');runs=$script:runs} |
        ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $statusPath -Encoding UTF8
}
function Test-Source {
    if ((Get-FileHash -LiteralPath $sourceArchive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.zip_sha256) { throw 'Source archive mismatch' }
    foreach ($file in $package.source_files.PSObject.Properties) {
        if ((Get-FileHash -LiteralPath (Join-Path $project $file.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $file.Value) { throw "Source mismatch: $($file.Name)" }
    }
    $inputManifest=Join-Path $work 'input_manifest.json'
    if (Test-Path -LiteralPath $inputManifest) {
        $inputHashes=Get-Content -LiteralPath $inputManifest -Raw | ConvertFrom-Json
        foreach ($file in $inputHashes.files.PSObject.Properties) {
            if ((Get-FileHash -LiteralPath (Join-Path $work ('inputs\'+$file.Name)) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $file.Value) { throw "Input mismatch: $($file.Name)" }
        }
    }
}
function Invoke-Checked([string]$name,[string]$program,[string[]]$arguments) {
    Test-Source
    $start=[DateTime]::UtcNow
    $process=Start-Process -FilePath $program -ArgumentList $arguments -WorkingDirectory $project -WindowStyle Hidden -Wait -PassThru `
        -RedirectStandardOutput (Join-Path $work ($name+'.stdout.log')) -RedirectStandardError (Join-Path $work ($name+'.stderr.log'))
    Test-Source
    $script:runs += [ordered]@{name=$name;process_exit_code=$process.ExitCode;arguments=$arguments;
        elapsed_seconds=([DateTime]::UtcNow-$start).TotalSeconds;source_verified_before_and_after=$true}
    Write-Status $false 0
    if ($process.ExitCode -ne 0) { throw "$name exited $($process.ExitCode)" }
}
try {
    Write-Status $false 0
    if ($env:COMPUTERNAME -ne 'W') { throw 'Unexpected host' }
    $package=Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    if ($package.zip_sha256 -ne '4b30436a7d84ce7089b5f7379b69c51e0c51b2eedf8654f40f76a7bc34afa123') { throw 'Unexpected source package' }
    if (-not (Test-Path -LiteralPath $project)) {
        if ($Stage -ne 'setup') { throw 'Run setup first' }
        if ((Get-FileHash -LiteralPath $sourceArchive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.zip_sha256) { throw 'Source archive mismatch' }
        Expand-Archive -LiteralPath $sourceArchive -DestinationPath $project
    }
    Test-Source
    $env:JULIA_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:JULIA_PKG_PRECOMPILE_AUTO='0'
    $env:THRESHOLD_SOURCE_MANIFEST=$manifestPath
    if ($Stage -eq 'setup') {
        Invoke-Checked 'python_unit' $python @('-m','unittest','discover','-s',(Join-Path $project 'tests\unit\python'),'-p','test_transport_mode_b_qualification.py','-v')
    } elseif ($Stage -eq 'selection') {
        $background='D:\CodexCompute\Julia_RelaxTime\full-rs-mode-b-background-20261025-v1'
        $review=Get-Content -LiteralPath (Join-Path $background 'reviewed_result.json') -Raw | ConvertFrom-Json
        if (-not $review.mode_b_background_qualified -or $review.points -ne 909) { throw 'Mode-B background review is incomplete' }
        Invoke-Checked 'selection' $python @((Join-Path $project 'scripts\analysis\relaxtime\prepare_transport_mode_b_qualification.py'),
            (Join-Path $background 'background'),(Join-Path $work 'inputs'))
        $hashes=[ordered]@{}
        foreach ($name in @('review_points.csv','channel_plan.csv','selection.json')) {
            $hashes[$name]=(Get-FileHash -LiteralPath (Join-Path $work ('inputs\'+$name)) -Algorithm SHA256).Hash.ToLowerInvariant()
        }
        [ordered]@{files=$hashes;source_archive_sha256=$package.zip_sha256;background_review_sha256=(Get-FileHash -LiteralPath (Join-Path $background 'reviewed_result.json') -Algorithm SHA256).Hash.ToLowerInvariant()} |
            ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $work 'input_manifest.json') -Encoding UTF8
        $inputPaths=@('review_points.csv','channel_plan.csv','selection.json') | ForEach-Object { Join-Path $work ('inputs\'+$_) }
        Compress-Archive -LiteralPath $inputPaths -DestinationPath (Join-Path $work 'inputs.zip')
    } else {
        if (@(Get-Process -Name julia -ErrorAction SilentlyContinue).Count -gt 0) { throw 'Another Julia process is running' }
        $arguments=@('--startup-file=no','--threads=1',"--project=$project",'--heap-size-hint=1200M')
        if ($Stage -eq 'unit') {
            $arguments += (Join-Path $project 'tests\unit\relaxtime\test_transport_validation_selection.jl')
        } else {
            $selection=Get-Content -LiteralPath (Join-Path $work 'inputs\selection.json') -Raw | ConvertFrom-Json
            if ($Stage -in @('kernel','outer','joint','production')) {
                $cheap=Get-Content -LiteralPath (Join-Path $work 'cheap\status.json') -Raw | ConvertFrom-Json
                if (-not $cheap.layer_target_passed) { throw 'Cheap layer qualification is incomplete' }
            }
            if ($Stage -eq 'joint') {
                foreach ($layer in @('kernel','outer')) {
                    $result=Get-Content -LiteralPath (Join-Path $work ($layer+'\status.json')) -Raw | ConvertFrom-Json
                    if (-not $result.layer_target_passed) { throw "$layer qualification is incomplete" }
                }
            }
            $arguments += @((Join-Path $project 'scripts\analysis\relaxtime\run_transport_mode_b_qualification.jl'),
                (Join-Path $work 'inputs'),(Join-Path $work $Stage),$Stage)
            if ($Stage -eq 'production') { $arguments += (Join-Path $work 'joint_summary') }
        }
        Invoke-Checked $Stage $julia $arguments
        if ($Stage -eq 'joint') {
            Invoke-Checked 'joint_summary' $python @((Join-Path $project 'scripts\analysis\relaxtime\summarize_transport_joint_validation.py'),
                (Join-Path $work 'joint'),'--outdir',(Join-Path $work 'joint_summary'),'--expected-points',[string]$selection.points.Count)
        }
    }
    Test-Source
    Write-Status $true 0
    exit 0
} catch {
    Write-Status $true 1 $_.Exception.Message
    exit 1
}
