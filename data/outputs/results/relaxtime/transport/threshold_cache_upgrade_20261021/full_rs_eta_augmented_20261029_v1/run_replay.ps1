$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
$casePath = 'D:/CodexCompute/trfull-eta-augmented-20261029-v1'
$sourcePath = 'D:/CodexCompute/trfull-20261026-v1/case/source'
$etaPath = 'D:/CodexCompute/trmott-eta-fine-20261029-v1'
$juliaExe = 'D:/CodexCompute/Julia_RelaxTime/transport-preflight-20261023-v1/runtime/julia-1.12.5/bin/julia.exe'
$env:JULIA_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:JULIA_PKG_PRECOMPILE_AUTO = '0'
if (Test-Path -LiteralPath ($casePath + '/replay_outcome.json')) { throw 'replay outcome already exists' }
Get-Content -LiteralPath ($casePath + '/selection/selection.json') -Raw | ConvertFrom-Json | Select-Object points, added, recomputed, fine_points, tighter_points
$argv = @('--startup-file=no', '--threads=1', ('--project=' + $sourcePath),
    ($casePath + '/complete_transport_eta_thermo.jl'), $sourcePath, $casePath, $etaPath)
$started = [DateTime]::UtcNow
$process = Start-Process -FilePath $juliaExe -ArgumentList $argv -WorkingDirectory $sourcePath -WindowStyle Hidden -PassThru `
    -RedirectStandardOutput ($casePath + '/replay.stdout.log') -RedirectStandardError ($casePath + '/replay.stderr.log')
$handle = $process.Handle
$process.WaitForExit()
$code = $process.ExitCode
[ordered]@{status = $(if ($code -eq 0) { 'completed' } else { 'failed' }); exit_code = $code;
    started_utc = $started.ToString('o'); completed_utc = [DateTime]::UtcNow.ToString('o'); arguments = $argv} |
    ConvertTo-Json -Depth 5 | Set-Content -LiteralPath ($casePath + '/replay_outcome.json') -Encoding UTF8
Get-Content -LiteralPath ($casePath + '/replay.stdout.log') -Tail 3
Get-Content -LiteralPath ($casePath + '/replay.stderr.log') -Tail 18
exit $code
