param(
    [Parameter(Mandatory = $true)]
    [string] $SharedCase,
    [Parameter(Mandatory = $true)]
    [string] $LocalCase,
    [Parameter(Mandatory = $true)]
    [string] $RemoteRoot,
    [Parameter(Mandatory = $true)]
    [string] $Julia,
    [Parameter(Mandatory = $true)]
    [string] $Python,
    [Parameter(Mandatory = $true)]
    [string] $ProjectName,
    [Parameter(Mandatory = $true)]
    [string] $OutputStem,
    [Parameter(Mandatory = $true)]
    [string] $Formulation,
    [Parameter(Mandatory = $true)]
    [double] $FrequencyStart,
    [Parameter(Mandatory = $true)]
    [double] $FrequencyEnd,
    [Parameter(Mandatory = $true)]
    [int] $FrequencyIncrements,
    [ValidateSet("constants", "transient")]
    [string] $Mode = "constants",
    [Parameter(Mandatory = $true)]
    [string] $PSCADVersion,
    [ValidateScript({ $_ -eq 0 -or ($_ -ge 1024 -and $_ -le 65535) })]
    [int] $PSCADPort = 0,
    [Parameter(Mandatory = $true)]
    [ValidateRange(0, 2)]
    [int] $Verbosity,
    [Parameter(Mandatory = $true)]
    [ValidateRange(1, 2147483647)]
    [int] $TimeoutSeconds
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
# Bounds remote initialization as well as the TCP connection, independently of the run timeout.
$PSCADProbeTimeoutSeconds = 60

. (Join-Path $PSScriptRoot "lifecycle.ps1")

function Remove-DirectoryWithRetry {
    param(
        [Parameter(Mandatory = $true)]
        [string] $Path,
        [int] $Attempts = 30,
        [int] $DelayMilliseconds = 500,
        [switch] $BestEffort
    )

    for ($attempt = 1; $attempt -le $Attempts; $attempt++) {
        if (-not (Test-Path -LiteralPath $Path)) {
            return
        }
        try {
            Remove-Item -LiteralPath $Path -Recurse -Force -ErrorAction Stop
            return
        } catch {
            if ($attempt -eq $Attempts) {
                if ($BestEffort) {
                    return
                }
                throw
            }
            Start-Sleep -Milliseconds $DelayMilliseconds
        }
    }
}

function Copy-Result {
    param(
        [Parameter(Mandatory = $true)]
        [string] $Source,
        [Parameter(Mandatory = $true)]
        [string] $Destination
    )

    New-Item -ItemType Directory -Path $Destination -Force | Out-Null
    if (Test-Path -LiteralPath $Source -PathType Container) {
        Get-ChildItem -LiteralPath $Source -Force |
            Copy-Item -Destination $Destination -Recurse -Force
    }
}

function Write-NewLines {
    param(
        [Parameter(Mandatory = $true)]
        [string] $Path,
        [Parameter(Mandatory = $true)]
        [ref] $Count,
        [Parameter(Mandatory = $true)]
        [bool] $ErrorStream
    )

    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        return
    }
    $lines = @(Get-Content -LiteralPath $Path)
    for ($index = $Count.Value; $index -lt $lines.Count; $index++) {
        if ($ErrorStream) {
            [Console]::Error.WriteLine($lines[$index])
        } else {
            [Console]::Out.WriteLine($lines[$index])
        }
    }
    $Count.Value = $lines.Count
}

if (-not (Test-Path -LiteralPath $SharedCase -PathType Container)) {
    throw "PSCAD exchange directory is unavailable: $SharedCase"
}
if (-not (Test-Path -LiteralPath $Julia -PathType Leaf)) {
    throw "Julia executable is unavailable: $Julia"
}
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python executable is unavailable: $Python"
}

$sharedFull = [IO.Path]::GetFullPath($SharedCase).TrimEnd('\') + '\'
$localFull = [IO.Path]::GetFullPath($LocalCase).TrimEnd('\') + '\'
if ($sharedFull.StartsWith($localFull, [StringComparison]::OrdinalIgnoreCase) -or
    $localFull.StartsWith($sharedFull, [StringComparison]::OrdinalIgnoreCase)) {
    throw "PSCAD exchange and scratch directories must not overlap"
}
$ownerPath = Join-Path $LocalCase "owner.txt"
$pscadOwnerPath = Join-Path $LocalCase "pscad-owner.txt"
Stop-PSCADCase -LocalCase $LocalCase
if (Test-Path -LiteralPath $LocalCase) {
    Remove-DirectoryWithRetry -Path $LocalCase
}
New-Item -ItemType Directory -Path $LocalCase -Force | Out-Null
Get-ChildItem -LiteralPath $SharedCase -Force |
    Where-Object Name -ne "outputs" |
    Copy-Item -Destination $LocalCase -Recurse -Force

$toolkit = Join-Path $LocalCase "toolkit"
$runner = Join-Path $toolkit "runner.jl"
$probe = Join-Path $toolkit "probe.py"
$project = Join-Path $toolkit "Project.toml"
$projectFile = Join-Path $LocalCase "generated.pscx"
$result = Join-Path $LocalCase "result"
$sharedOutput = Join-Path $SharedCase "outputs"
$stdoutPath = Join-Path $result "stdout.txt"
$stderrPath = Join-Path $result "stderr.txt"
$probePath = Join-Path $result "pscad-probe.txt"
$probeErrorPath = Join-Path $result "pscad-probe-stderr.txt"
$exitPath = Join-Path $result "exit-code.txt"
$invokePath = Join-Path $LocalCase "invoke-runner.cmd"

foreach ($required in @($runner, $probe, $project, $projectFile)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Staged PSCAD input is missing: $required"
    }
}

New-Item -ItemType Directory -Path $result -Force | Out-Null
New-Item -ItemType Directory -Path $sharedOutput -Force | Out-Null
$commandValues = @(
    $Julia,
    $Python,
    $toolkit,
    $runner,
    $projectFile,
    $result,
    $ProjectName,
    $OutputStem,
    $RemoteRoot,
    $Formulation,
    $Mode,
    $stdoutPath,
    $stderrPath,
    $exitPath
)
foreach ($value in $commandValues) {
    if ($value.IndexOfAny(@([char] '"', [char] '%', [char] 10, [char] 13)) -ge 0) {
        throw "PSCAD runner command value contains an unsupported character: $value"
    }
}
$pscadOwned = $false
$pscadProcessId = 0
$process = $null
$exitCode = 1
$timedOut = $false
$stdoutCount = 0
$stderrCount = 0
$cleanupFailures = [Collections.Generic.List[string]]::new()
try {
    if ($PSCADPort -ne 0) {
        $listeners = @(Get-NetTCPConnection -LocalPort $PSCADPort -State Listen -ErrorAction SilentlyContinue |
            Select-Object -ExpandProperty OwningProcess -Unique)
        if ($listeners.Count -eq 0) {
            $findExeOutput = Join-Path $result "pscad-find-exe-stdout.txt"
            $findExeError = Join-Path $result "pscad-find-exe-stderr.txt"
            $findExeCode = Invoke-PSCADNative -Executable $Python -Arguments @("-c", "from mhi.common import process; print(process.find_exe('PSCAD', version='$PSCADVersion', x64=True))") -StdoutPath $findExeOutput -StderrPath $findExeError
            $exe = ([string] (Get-Content -LiteralPath $findExeOutput -Raw)).Trim()
            if ($findExeCode -ne 0 -or [string]::IsNullOrWhiteSpace($exe) -or -not (Test-Path -LiteralPath $exe -PathType Leaf)) {
                throw "Could not resolve the PSCAD $PSCADVersion executable with mhi.common"
            }
            $exe = [IO.Path]::GetFullPath($exe.Trim())
            $runId = Split-Path (Split-Path $LocalCase -Parent) -Leaf
            $taskName = "LineCableModels-PSCAD-$runId"
            $startedAt = [DateTime]::UtcNow
            try {
                $action = New-ScheduledTaskAction -Execute $exe -Argument "/startup:au /port:$PSCADPort /silence:true /splash:false"
                $principal = New-ScheduledTaskPrincipal -UserId ([Security.Principal.WindowsIdentity]::GetCurrent().Name) -LogonType Interactive -RunLevel Limited
                $settings = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -MultipleInstances IgnoreNew
                Register-ScheduledTask -TaskName $taskName -TaskPath '\' -Action $action -Principal $principal -Settings $settings -Force | Out-Null
                Start-ScheduledTask -TaskName $taskName -TaskPath '\'
                $deadline = [DateTime]::UtcNow.AddSeconds(90)
                do {
                    $listeners = @(Get-NetTCPConnection -LocalPort $PSCADPort -State Listen -ErrorAction SilentlyContinue |
                        Select-Object -ExpandProperty OwningProcess -Unique)
                    if ($listeners.Count -gt 0) { break }
                    Start-Sleep -Milliseconds 500
                } while ([DateTime]::UtcNow -lt $deadline)
                if ($listeners.Count -eq 0) {
                    # PSCAD can start without opening its automation port.
                    # Keep ownership so the outer finally can still close it.
                    $applications = @(Get-CimInstance Win32_Process -Filter "Name = 'Pscad.exe'" |
                        Where-Object {
                            $_.ExecutablePath -ieq $exe -and
                            $_.CreationDate.ToUniversalTime() -ge $startedAt -and
                            $null -ne $_.CommandLine -and
                            $_.CommandLine -match "/port:$PSCADPort(\s|$)"
                        })
                    if ($applications.Count -eq 1) {
                        $application = $applications[0]
                        $pscadProcessId = $application.ProcessId
                        Set-Content -LiteralPath $pscadOwnerPath -Value @($pscadProcessId, $exe, $PSCADPort, $application.CreationDate.ToUniversalTime().ToString('o'))
                        $pscadOwned = $true
                    }
                    throw "PSCAD did not listen on port $PSCADPort within 90 seconds. A logged-on desktop session is required for the station account (a disconnected RDP session is fine), or start PSCAD with /port:$PSCADPort."
                }
                if ($listeners.Count -ne 1) {
                    throw "Multiple processes are listening on PSCAD port $PSCADPort"
                }
                $pscadProcessId = $listeners[0]
                $application = Get-CimInstance Win32_Process -Filter "ProcessId = $pscadProcessId"
                if ($null -eq $application -or $application.ExecutablePath -ine $exe -or
                    $application.CreationDate.ToUniversalTime() -lt $startedAt -or
                    $null -eq $application.CommandLine -or
                    $application.CommandLine -notmatch "/port:$PSCADPort(\s|$)") {
                    throw "Listener on port $PSCADPort does not belong to the PSCAD desktop start"
                }
                Set-Content -LiteralPath $pscadOwnerPath -Value @($pscadProcessId, $exe, $PSCADPort, $application.CreationDate.ToUniversalTime().ToString('o'))
                $pscadOwned = $true
            } finally {
                try {
                    Unregister-ScheduledTask -TaskName $taskName -TaskPath '\' -Confirm:$false -ErrorAction Stop
                } catch {
                    $cleanupFailures.Add("Unregister scheduled task: $($_.Exception.Message)")
                }
            }
        } else {
            if ($listeners.Count -ne 1) {
                throw "Multiple processes are listening on PSCAD port $PSCADPort"
            }
            $pscadProcessId = $listeners[0]
        }
        [Console]::Out.WriteLine("PSCAD port=$PSCADPort PID=$pscadProcessId owned=$pscadOwned")
        $probeExitCode = Invoke-PSCADNative -Executable $Python -Arguments @("-u", $probe, "127.0.0.1", $PSCADPort, $PSCADProbeTimeoutSeconds) -StdoutPath $probePath -StderrPath $probeErrorPath
        if ($probeExitCode -ne 0) {
            throw "PSCAD on port $PSCADPort is not answering remote calls within $PSCADProbeTimeoutSeconds s (busy, blocked by a dialog, or its desktop session is not active); check the instance or use a port without a listener so the run starts its own PSCAD."
        }
    }
    $runnerOwned = ($PSCADPort -eq 0 -or $pscadOwned).ToString().ToLowerInvariant()
    $runnerArguments = @(
        "`"--project=$toolkit`"",
        "--startup-file=no",
        "`"$runner`"",
        "`"$projectFile`"",
        "`"$result`"",
        "`"$ProjectName`"",
        "`"$OutputStem`"",
        "`"$Formulation`"",
        "`"$FrequencyStart`"",
        "`"$FrequencyEnd`"",
        "`"$FrequencyIncrements`"",
        "`"$PSCADVersion`"",
        "`"$Verbosity`"",
        "`"$Mode`"",
        "`"$PSCADPort`"",
        "`"$runnerOwned`"",
        "`"$RemoteRoot`"",
        "`"$TimeoutSeconds`""
    )
    $invokeLines = @(
        "@echo off",
        "setlocal",
        "set `"JULIA_PYTHONCALL_EXE=$Python`"",
        "`"$Julia`" $($runnerArguments -join ' ') 1>`"$stdoutPath`" 2>`"$stderrPath`"",
        "set `"PSCAD_EXIT=%ERRORLEVEL%`"",
        ">`"$exitPath`" echo %PSCAD_EXIT%",
        "exit /b %PSCAD_EXIT%"
    )
    [IO.File]::WriteAllLines($invokePath, $invokeLines, [Text.Encoding]::ASCII)

    $process = Start-Process `
        -FilePath $env:ComSpec `
        -ArgumentList @("/d", "/s", "/c", "`"$invokePath`"") `
        -PassThru `
        -WindowStyle Hidden
    Set-Content -LiteralPath $ownerPath -Value @($process.Id, $invokePath)
    if ($Verbosity -ge 2) {
        [Console]::Out.WriteLine("Started PSCAD runner process $($process.Id)")
    }

    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    while (-not $process.HasExited) {
        if ($Verbosity -ge 2) {
            Write-NewLines -Path $stdoutPath -Count ([ref] $stdoutCount) -ErrorStream $false
            Write-NewLines -Path $stderrPath -Count ([ref] $stderrCount) -ErrorStream $true
        }
        if ($stopwatch.Elapsed.TotalSeconds -ge $TimeoutSeconds) {
            $timedOut = $true
            break
        }
        Start-Sleep -Milliseconds 500
        $process.Refresh()
    }

    if ($timedOut) {
        [Console]::Error.WriteLine(
            "PSCAD runner exceeded the configured timeout of $TimeoutSeconds seconds"
        )
        $exitCode = 124
    } else {
        $process.WaitForExit()
        if (Test-Path -LiteralPath $exitPath -PathType Leaf) {
            $recordedExit = (Get-Content -LiteralPath $exitPath -Raw).Trim()
            if (-not [int]::TryParse($recordedExit, [ref] $exitCode)) {
                [Console]::Error.WriteLine(
                    "PSCAD runner wrote an invalid exit code: $recordedExit"
                )
                $exitCode = 1
            }
        } else {
            [Console]::Error.WriteLine("PSCAD runner did not write an exit code")
            $exitCode = 1
        }
        if ($Verbosity -ge 2) {
            [Console]::Out.WriteLine(
                "PSCAD runner process $($process.Id) exited with code $exitCode"
            )
        }
    }

    if ($Verbosity -ge 2) {
        Write-NewLines -Path $stdoutPath -Count ([ref] $stdoutCount) -ErrorStream $false
        Write-NewLines -Path $stderrPath -Count ([ref] $stderrCount) -ErrorStream $true
    }
} catch {
    # Report the primary error before any cleanup can fail.
    [Console]::Error.WriteLine($_.Exception.ToString())
    $exitCode = 1
} finally {
    try {
        if ($null -ne $process -and -not $process.HasExited) {
            Stop-OwnedRunner -OwnerPath $ownerPath -OwnedRoot $LocalCase
        }
    } catch { $cleanupFailures.Add("Stop runner: $($_.Exception.Message)") }
    try {
        Stop-EMTDC -CaseDirectory $LocalCase
    } catch { $cleanupFailures.Add("Stop EMTDC: $($_.Exception.Message)") }
    try {
        if ($pscadOwned) {
            Stop-OwnedPSCAD -OwnerPath $pscadOwnerPath -GraceSeconds 10
        }
    } catch { $cleanupFailures.Add("Stop owned PSCAD: $($_.Exception.Message)") }
    try {
        Copy-Result -Source $result -Destination $sharedOutput
    } catch { $cleanupFailures.Add("Copy outputs: $($_.Exception.Message)") }
    try {
        Remove-Item -LiteralPath $ownerPath -Force -ErrorAction SilentlyContinue
    } catch { $cleanupFailures.Add("Remove runner owner: $($_.Exception.Message)") }
    foreach ($failure in $cleanupFailures) {
        [Console]::Error.WriteLine("Cleanup failure: $failure")
    }
    if ($cleanupFailures.Count -gt 0 -and $exitCode -eq 0) { $exitCode = 1 }
}

if ($exitCode -eq 0) {
    Remove-DirectoryWithRetry -Path $LocalCase -BestEffort
}
exit $exitCode
