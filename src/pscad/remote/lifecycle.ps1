# Shared process lifecycle for the supervisor and caller-side cancellation.
$ErrorActionPreference = "Stop"

# Native stderr is non-terminating in Windows PowerShell 5.1 only with Continue.
# This function uses a local error preference so callers can decide from the exit code.
function Invoke-PSCADNative {
    param([string] $Executable, [string[]] $Arguments,
        [string] $StdoutPath, [string] $StderrPath)

    New-Item -ItemType Directory -Path (Split-Path $StdoutPath -Parent) -Force | Out-Null
    $ErrorActionPreference = "Continue"
    $global:LASTEXITCODE = 1
    & $Executable @Arguments 1> $StdoutPath 2> $StderrPath
    return $global:LASTEXITCODE
}

function Stop-PSCADProcessTree {
    param([int] $ProcessId, [string] $CaseDirectory)

    $result = Join-Path $CaseDirectory "result"
    $stdout = Join-Path $result "taskkill-$ProcessId-stdout.txt"
    $stderr = Join-Path $result "taskkill-$ProcessId-stderr.txt"
    $code = 1
    try {
        $code = Invoke-PSCADNative -Executable "taskkill.exe" -Arguments @("/PID", $ProcessId, "/T", "/F") -StdoutPath $stdout -StderrPath $stderr
    } catch {
        [Console]::Error.WriteLine($_.Exception.Message)
    }
    [Console]::Error.WriteLine("taskkill PID=$ProcessId exit code=$code")
    foreach ($path in @($stdout, $stderr)) {
        if (Test-Path -LiteralPath $path -PathType Leaf) {
            foreach ($line in @(Get-Content -LiteralPath $path)) {
                [Console]::Error.WriteLine($line)
            }
        }
    }
    return $code
}

function Stop-OwnedRunner {
    param(
        [Parameter(Mandatory = $true)]
        [string] $OwnerPath,
        [Parameter(Mandatory = $true)]
        [string] $OwnedRoot,
        [int] $GraceSeconds = 10
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
    if ($processId -le 0) {
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

    $exiting = $process.ExecutionState -eq 7 -or
        [string]::IsNullOrWhiteSpace($process.ExecutablePath) -or
        [string]::IsNullOrWhiteSpace($process.CommandLine) -or
        [string]::IsNullOrWhiteSpace([string] $process.CreationDate)
    if ($exiting) {
        Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
        return
    }
    if ($process.CommandLine.IndexOf($runner, [StringComparison]::OrdinalIgnoreCase) -lt 0) {
        [Console]::Error.WriteLine("Warning: PID $processId no longer belongs to the recorded PSCAD runner; leaving the process untouched and removing its stale owner file.")
        Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
        return
    }

    $stopCode = Stop-PSCADProcessTree -ProcessId $processId -CaseDirectory $OwnedRoot
    if ($stopCode -ne 0) {
        $deadline = [DateTime]::UtcNow.AddSeconds($GraceSeconds)
        do {
            $remaining = Get-CimInstance Win32_Process -Filter "ProcessId = $processId"
            if ($null -eq $remaining) { break }
            $exiting = $remaining.ExecutionState -eq 7 -or
                [string]::IsNullOrWhiteSpace($remaining.ExecutablePath) -or
                [string]::IsNullOrWhiteSpace($remaining.CommandLine) -or
                [string]::IsNullOrWhiteSpace([string] $remaining.CreationDate)
            if (-not $exiting -and
                ($remaining.ExecutablePath -ine $process.ExecutablePath -or
                 $remaining.CreationDate -ne $process.CreationDate -or
                 $remaining.CommandLine.IndexOf($runner, [StringComparison]::OrdinalIgnoreCase) -lt 0)) {
                [Console]::Error.WriteLine("Warning: PID $processId no longer belongs to the recorded PSCAD runner; leaving the process untouched.")
                break
            }
            if ([DateTime]::UtcNow -ge $deadline) {
                if (-not $exiting) { throw "Could not stop PSCAD runner process tree $processId" }
                break
            }
            Start-Sleep -Milliseconds 500
        } while ($true)
    }
    Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
}

function Stop-EMTDC {
    param([string] $CaseDirectory, [int] $GraceSeconds = 10)

    $failures = [Collections.Generic.List[string]]::new()
    $root = [IO.Path]::GetFullPath($CaseDirectory).TrimEnd('\') + '\'
    # Simulation executables are built beneath the staged case. Never select
    # by executable name alone: restored user projects can also run.
    foreach ($simulation in @(Get-CimInstance Win32_Process)) {
        if ($null -ne $simulation.ExecutablePath -and
            $simulation.ExecutablePath.StartsWith($root, [StringComparison]::OrdinalIgnoreCase) -and
            $simulation.ExecutablePath.EndsWith('.exe', [StringComparison]::OrdinalIgnoreCase)) {
            try {
                $stopCode = Stop-PSCADProcessTree -ProcessId $simulation.ProcessId -CaseDirectory $CaseDirectory
                if ($stopCode -ne 0) {
                    $deadline = [DateTime]::UtcNow.AddSeconds($GraceSeconds)
                    do {
                        $remaining = Get-CimInstance Win32_Process -Filter "ProcessId = $($simulation.ProcessId)"
                        if ($null -eq $remaining) { break }
                        $exiting = $remaining.ExecutionState -eq 7 -or
                            [string]::IsNullOrWhiteSpace($remaining.ExecutablePath) -or
                            [string]::IsNullOrWhiteSpace($remaining.CommandLine) -or
                            [string]::IsNullOrWhiteSpace([string] $remaining.CreationDate)
                        if (-not $exiting -and
                            ($remaining.ExecutablePath -ine $simulation.ExecutablePath -or
                             $remaining.CreationDate -ne $simulation.CreationDate)) {
                            [Console]::Error.WriteLine("Warning: PID $($simulation.ProcessId) no longer belongs to the recorded EMTDC simulation; leaving the process untouched.")
                            break
                        }
                        if ([DateTime]::UtcNow -ge $deadline) {
                            if (-not $exiting) { throw "Could not stop EMTDC simulation process $($simulation.ProcessId)" }
                            break
                        }
                        Start-Sleep -Milliseconds 500
                    } while ($true)
                }
            } catch { $failures.Add($_.Exception.Message) }
        }
    }
    if ($failures.Count -gt 0) { throw ($failures -join [Environment]::NewLine) }
}

function Stop-OwnedPSCAD {
    param([string] $OwnerPath, [int] $GraceSeconds = 0)

    if (-not (Test-Path -LiteralPath $OwnerPath -PathType Leaf)) {
        return
    }
    # Only a desktop PSCAD started by this supervisor gets an owner file.
    # Include creation time to distinguish the process from a reused PID.
    $owner = @(Get-Content -LiteralPath $OwnerPath)
    $processId = 0
    $port = 0
    if ($owner.Count -ne 4 -or
        -not [int]::TryParse($owner[0], [ref] $processId) -or $processId -le 0 -or
        -not [int]::TryParse($owner[2], [ref] $port) -or $port -lt 1024 -or $port -gt 65535 -or
        [IO.Path]::GetFileName($owner[1]) -ine 'Pscad.exe') {
        throw "Invalid PSCAD application owner file: $OwnerPath"
    }
    $deadline = [DateTime]::UtcNow.AddSeconds($GraceSeconds)
    $stopFailed = $false
    do {
        $application = Get-CimInstance Win32_Process -Filter "ProcessId = $processId"
        if ($null -eq $application) {
            Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
            return
        }
        # CIM can keep a terminating process after its identity fields disappear.
        # ExecutionState 7 is terminated. Never kill an incomplete process record.
        $exiting = $application.ExecutionState -eq 7 -or
            [string]::IsNullOrWhiteSpace($application.ExecutablePath) -or
            [string]::IsNullOrWhiteSpace($application.CommandLine) -or
            [string]::IsNullOrWhiteSpace([string] $application.CreationDate)
        if (-not $exiting) {
            if ($application.ExecutablePath -ine $owner[1] -or
                $application.CreationDate.ToUniversalTime().ToString('o') -cne $owner[3]) {
                [Console]::Error.WriteLine("Warning: PID $processId no longer belongs to the recorded owned PSCAD; leaving the process untouched and removing its stale owner file.")
                Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
                return
            }
            if ($application.CommandLine -notmatch "/port:$port(\s|$)") {
                [Console]::Error.WriteLine("Warning: PID $processId does not report the recorded PSCAD port $port; leaving the process untouched and removing its stale owner file.")
                Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
                return
            }
        }
        if ([DateTime]::UtcNow -lt $deadline) {
            Start-Sleep -Milliseconds 500
            continue
        }
        if ($exiting) { break }
        if ($stopFailed) {
            throw "Could not stop owned PSCAD process tree $processId"
        }
        $stopCode = Stop-PSCADProcessTree -ProcessId $processId -CaseDirectory (Split-Path $OwnerPath -Parent)
        $stopFailed = $stopCode -ne 0
        if (-not $stopFailed) { break }
        # Allow a fresh grace period after taskkill while PSCAD exits.
        # Check identity on each poll and make no further kill attempt.
        $deadline = [DateTime]::UtcNow.AddSeconds($GraceSeconds)
    } while ($true)
    Remove-Item -LiteralPath $OwnerPath -Force -ErrorAction SilentlyContinue
}

function Stop-PSCADCase {
    param(
        [Parameter(Mandatory = $true)]
        [string] $LocalCase,
        [int] $GraceSeconds = 10
    )

    $failures = [Collections.Generic.List[string]]::new()
    try {
        Stop-OwnedRunner -OwnerPath (Join-Path $LocalCase "owner.txt") -OwnedRoot $LocalCase -GraceSeconds $GraceSeconds
    } catch { $failures.Add($_.Exception.Message) }
    try {
        Stop-EMTDC -CaseDirectory $LocalCase -GraceSeconds $GraceSeconds
    } catch { $failures.Add($_.Exception.Message) }
    try {
        Stop-OwnedPSCAD -OwnerPath (Join-Path $LocalCase "pscad-owner.txt") -GraceSeconds $GraceSeconds
    } catch { $failures.Add($_.Exception.Message) }
    if ($failures.Count -gt 0) {
        throw ($failures -join [Environment]::NewLine)
    }
}
