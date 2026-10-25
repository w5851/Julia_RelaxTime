$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$work='D:\CodexCompute\Julia_RelaxTime\full-rs-runner-review-checks-20261026-v1'
$python='C:\Users\82306\AppData\Local\Programs\Python\Python313\python.exe'
$project=Join-Path $work 'project'
if($env:COMPUTERNAME -ne 'W'){throw 'Unexpected compute host'}
if(Test-Path -LiteralPath $project){throw 'Use a new test project'}
$package=Get-Content -LiteralPath (Join-Path $work 'package_manifest.json') -Raw | ConvertFrom-Json
if((Get-FileHash -LiteralPath (Join-Path $work 'source.zip') -Algorithm SHA256).Hash.ToLowerInvariant() -ne $package.source_zip_sha256){throw 'Source archive mismatch'}
Expand-Archive -LiteralPath (Join-Path $work 'source.zip') -DestinationPath $project
function Test-Source {
    foreach($file in $package.source_files.PSObject.Properties){
        if((Get-FileHash -LiteralPath (Join-Path $project $file.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $file.Value){throw "Source changed: $($file.Name)"}
    }
}
Test-Source
$process=Start-Process -FilePath $python -ArgumentList @('-B','-m','unittest','discover','-s','tests/unit/python','-p','test_transport_runner_pilot.py','-v') -WorkingDirectory $project -WindowStyle Hidden -Wait -PassThru `
    -RedirectStandardOutput (Join-Path $work 'unit.stdout.log') -RedirectStandardError (Join-Path $work 'unit.stderr.log')
Test-Source
$log=Get-Content -LiteralPath (Join-Path $work 'unit.stderr.log') -Raw
$match=[regex]::Match($log,'Ran (\d+) tests? in')
$passed=$process.ExitCode -eq 0 -and $match.Success -and [int]$match.Groups[1].Value -eq 6 -and $log -match '(?m)^OK\s*$'
$result=[ordered]@{host=$env:COMPUTERNAME;exit_code=$process.ExitCode;tests=if($match.Success){[int]$match.Groups[1].Value}else{0};passed=$passed;source_verified_before_and_after=$true;source_zip_sha256=$package.source_zip_sha256;controller_sha256=(Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLowerInvariant();julia_started=$false;updated_utc=[DateTime]::UtcNow.ToString('o')}
$result | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $work 'reviewed_result.json') -Encoding UTF8
Get-Content -LiteralPath (Join-Path $work 'unit.stderr.log')
$result | ConvertTo-Json
if(-not $passed){exit 1}
exit 0
