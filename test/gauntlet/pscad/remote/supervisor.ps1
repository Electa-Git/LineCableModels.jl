param(
    [Parameter(Mandatory = $true)]
    [string] $SharedCase,
    [Parameter(Mandatory = $true)]
    [string] $LocalCase,
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
    [string] $EarthField,
    [Parameter(Mandatory = $true)]
    [int] $EarthValue,
    [Parameter(Mandatory = $true)]
    [string] $EarthReadback,
    [Parameter(Mandatory = $true)]
    [double] $FrequencyStart,
    [Parameter(Mandatory = $true)]
    [double] $FrequencyEnd,
    [Parameter(Mandatory = $true)]
    [int] $FrequencyIncrements,
    [Parameter(Mandatory = $true)]
    [string] $PSCADVersion,
    [Parameter(Mandatory = $true)]
    [ValidateRange(0, 2)]
    [int] $Verbosity,
    [Parameter(Mandatory = $true)]
    [ValidateRange(1, 2147483647)]
    [int] $TimeoutSeconds
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"

function Stop-OwnedRunner {
    param(
        [Parameter(Mandatory = $true)]
        [string] $OwnerPath,
        [Parameter(Mandatory = $true)]
        [string] $OwnedRoot
    )

    if (-not (Test-Path -LiteralPath $OwnerPath -PathType Leaf)) {
        return
    }

    $owner = @(Get-Content -LiteralPath $OwnerPath)
    if ($owner.Count -ne 2) {
        throw "Invalid PSCAD runner owner file: $OwnerPath"
    }

    $processId = 0
    if (-not [int]::TryParse($owner[0], [ref] $processId)) {
        throw "Invalid PSCAD runner PID in $OwnerPath"
    }

    $runner = [IO.Path]::GetFullPath($owner[1])
    $root = [IO.Path]::GetFullPath($OwnedRoot).TrimEnd('\') + '\'
    if (-not $runner.StartsWith($root, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to stop a process outside the PSCAD case directory: $runner"
    }

    $process = Get-CimInstance Win32_Process -Filter "ProcessId = $processId"
    if ($null -eq $process) {
        Remove-Item -LiteralPath $OwnerPath -Force
        return
    }

    if ($null -eq $process.CommandLine -or
        $process.CommandLine.IndexOf($runner, [StringComparison]::OrdinalIgnoreCase) -lt 0) {
        throw "PID $processId no longer belongs to the recorded PSCAD runner"
    }

    & taskkill.exe /PID $processId /T /F | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "Could not stop PSCAD runner process tree $processId"
    }
    Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
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
    throw "VirtioFS case directory is unavailable: $SharedCase"
}
if (-not (Test-Path -LiteralPath $Julia -PathType Leaf)) {
    throw "Julia executable is unavailable: $Julia"
}
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python executable is unavailable: $Python"
}

$ownerPath = Join-Path $LocalCase "owner.txt"
Stop-OwnedRunner -OwnerPath $ownerPath -OwnedRoot $LocalCase
if (Test-Path -LiteralPath $LocalCase) {
    Remove-Item -LiteralPath $LocalCase -Recurse -Force
}
New-Item -ItemType Directory -Path $LocalCase -Force | Out-Null
Get-ChildItem -LiteralPath $SharedCase -Force |
    Where-Object Name -ne "outputs" |
    Copy-Item -Destination $LocalCase -Recurse -Force

$toolkit = Join-Path $LocalCase "toolkit"
$runner = Join-Path $toolkit "runner.jl"
$project = Join-Path $toolkit "Project.toml"
$projectFile = Join-Path $LocalCase "generated.pscx"
$result = Join-Path $LocalCase "result"
$sharedOutput = Join-Path $SharedCase "outputs"
$stdoutPath = Join-Path $result "stdout.txt"
$stderrPath = Join-Path $result "stderr.txt"
$exitPath = Join-Path $result "exit-code.txt"
$invokePath = Join-Path $LocalCase "invoke-runner.cmd"

foreach ($required in @($runner, $project, $projectFile)) {
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
    $Formulation,
    $EarthField,
    "$EarthValue",
    $EarthReadback,
    $stdoutPath,
    $stderrPath,
    $exitPath
)
foreach ($value in $commandValues) {
    if ($value.IndexOfAny(@([char] '"', [char] '%', [char] 10, [char] 13)) -ge 0) {
        throw "PSCAD runner command value contains an unsupported character: $value"
    }
}
$runnerArguments = @(
    "`"--project=$toolkit`"",
    "--startup-file=no",
    "`"$runner`"",
    "`"$projectFile`"",
    "`"$result`"",
    "`"$ProjectName`"",
    "`"$OutputStem`"",
    "`"$Formulation`"",
    "`"$EarthField`"",
    "`"$EarthValue`"",
    "`"$EarthReadback`"",
    "`"$FrequencyStart`"",
    "`"$FrequencyEnd`"",
    "`"$FrequencyIncrements`"",
    "`"$PSCADVersion`"",
    "`"$Verbosity`""
)
$invokeLines = @(
    "@echo off",
    "setlocal",
    "set `"JULIA_PYTHONCALL_EXE=$Python`"",
    "`"$Julia`" $($runnerArguments -join ' ') 1>`"$stdoutPath`" 2>`"$stderrPath`"",
    "set `"GAUNTLET_EXIT=%ERRORLEVEL%`"",
    ">`"$exitPath`" echo %GAUNTLET_EXIT%",
    "exit /b %GAUNTLET_EXIT%"
)
[IO.File]::WriteAllLines($invokePath, $invokeLines, [Text.Encoding]::ASCII)

$process = $null
$exitCode = 1
$timedOut = $false
$stdoutCount = 0
$stderrCount = 0
try {
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
        Stop-OwnedRunner -OwnerPath $ownerPath -OwnedRoot $LocalCase
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
    [Console]::Error.WriteLine($_.Exception.ToString())
    if ($null -ne $process -and -not $process.HasExited) {
        Stop-OwnedRunner -OwnerPath $ownerPath -OwnedRoot $LocalCase
    }
    $exitCode = 1
} finally {
    Copy-Result -Source $result -Destination $sharedOutput
    Remove-Item -LiteralPath $ownerPath -Force -ErrorAction SilentlyContinue
}

if ($exitCode -eq 0) {
    Remove-Item -LiteralPath $LocalCase -Recurse -Force -ErrorAction SilentlyContinue
}
exit $exitCode
