# Auto-throttling launcher: watches free system RAM and only starts a
# new config once there's enough headroom -- no manual Task Manager
# checking needed. Safe to just launch and walk away.
#
# Usage (from PowerShell, in this folder):
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass   # if blocked
#   .\run_all_1kyr_auto.ps1

$configs = @(
    @("spim","K1_3.02e-09"), @("spim","K2_1.19e-08"), @("spim","K3_4.69e-08"),
    @("spim","K4_1.85e-07"), @("spim","K5_7.28e-07"), @("spim","K6_2.87e-06"),
    @("spim","K7_1.13e-05"), @("space","K1_3.02e-09"), @("space","K2_1.19e-08"),
    @("space","K3_4.69e-08"), @("space","K4_1.85e-07"), @("space","K5_7.28e-07"),
    @("space","K6_2.87e-06"), @("space","K7_1.13e-05")
)

$minFreeGB       = 6    # never start a new job if it would leave less than
                         # this much RAM free -- headroom for Windows itself
$hardMaxParallel = 8    # absolute ceiling regardless of free memory, as a
                         # safety net against a bad memory reading
$checkEverySec   = 15   # how often to re-check while waiting
$warmupSec       = 90   # after starting a job, wait this long before the
                         # next decision -- gives it time to ramp up to its
                         # real memory footprint so the free-RAM reading
                         # used for the *next* decision is accurate

function Get-FreeGB {
    $os = Get-CimInstance Win32_OperatingSystem
    return [math]::Round(($os.FreePhysicalMemory * 1KB) / 1GB, 2)
}

New-Item -ItemType Directory -Force -Path logs | Out-Null
$queue = New-Object System.Collections.Generic.Queue[Object]
foreach ($cfg in $configs) { $queue.Enqueue($cfg) }

$running = @()  # PSCustomObjects: Process, Name

Write-Host "Starting. Free RAM: $(Get-FreeGB) GB. Will keep at least $minFreeGB GB free, max $hardMaxParallel concurrent."

while ($queue.Count -gt 0 -or ($running | Where-Object { -not $_.Process.HasExited }).Count -gt 0) {

    foreach ($f in ($running | Where-Object { $_.Process.HasExited })) {
        Write-Host "[$(Get-Date -Format HH:mm:ss)] Finished: $($f.Name)"
    }
    $running = @($running | Where-Object { -not $_.Process.HasExited })

    $freeGB = Get-FreeGB

    if ($queue.Count -gt 0 -and $running.Count -lt $hardMaxParallel -and $freeGB -ge $minFreeGB) {
        $cfg = $queue.Dequeue()
        $model, $kname = $cfg
        $log = "logs\${model}_${kname}_1kyr.log"
        $venvPython = ".\erodibility_env\Scripts\python.exe"
        $arg = "/c $venvPython 03_run_one_1m_1kyr.py $model $kname > `"$log`" 2>&1"
        $p = Start-Process -FilePath "cmd.exe" -ArgumentList $arg -WindowStyle Hidden -PassThru
        $running += [PSCustomObject]@{ Process = $p; Name = "$model $kname" }
        Write-Host "[$(Get-Date -Format HH:mm:ss)] Started $model $kname (PID $($p.Id)) -- free was $freeGB GB, $($running.Count - 1) already running, $($queue.Count) left"
        Start-Sleep -Seconds $warmupSec
    } else {
        Write-Host "[$(Get-Date -Format HH:mm:ss)] Waiting -- free=$freeGB GB, running=$($running.Count), queued=$($queue.Count)"
        Start-Sleep -Seconds $checkEverySec
    }
}

Write-Host "All 14 configs complete. Results -> results_1kyr_dt10\, logs -> logs\*_1kyr.log"
