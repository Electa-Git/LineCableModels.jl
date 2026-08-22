# PSCAD validation gauntlet

Each file in `cases/` owns one complete validation case. The case defines the LineCableModels problem, formulation, exact frequency samples, explicit terminal order, tolerances, and assertions. Gauntlet cases retain every terminal: phase assignments must be the contiguous integers `1:n`, and both Kron and bundle reduction must be disabled. Disposable exports and diagnostics remain in `cases/.work/<case>/`.

The gauntlet has three modes. `snapshot` reads an existing reference and never loads the PSCAD code. `live` executes PSCAD without changing the reference. `record` completes the live checks before replacing the reference. A failed operation never switches modes.

| Mode | PSCAD | Network | Reads accepted artifact | May replace artifact |
|:--|:--:|:--:|:--:|:--:|
| `snapshot` | No | No | Yes | No |
| `live` | Yes | Yes | No | No |
| `record` | Yes | Yes | Yes, when one exists | Yes, with `LINECABLEMODELS_GAUNTLET_PERSIST=true` |

## Comparisons

The live calculation compares the phase-domain `Z[i,j,:]` and `Y[i,j,:]` series returned by PSCAD with the corresponding LineCableModels series. Absolute and relative RMS errors are calculated separately for every matrix term across frequency. Gauntlet does not flatten the matrices, combine matrix terms, reorder terminals, interpolate frequencies, or post-process PSCAD results.

PSCAD does not emit modal `Z` and `Y` for the supported frequency-dependent model. A modal `LineParameters` result therefore throws a clear not-implemented error instead of being compared through a derived proxy.

Each case has three tolerance groups:

- `reference` checks LineCableModels against PSCAD.
- `regression` checks the current LineCableModels result against the accepted result stored in the artifact.
- `performance` checks time, bytes, and allocations against the accepted artifact only when the Julia, operating system, architecture, thread count, and BLAS configuration are identical.

Performance from different environments is reported as non-comparable and cannot pass or fail the case.

## Live setup

Instantiate the gauntlet test environment once from the repository root:

```bash
julia --project=test/gauntlet --startup-file=no -e \
'using Pkg; Pkg.instantiate()'
```

JLD2 and BenchmarkTools belong to this environment. They are not production dependencies or weak dependencies of LineCableModels.

Copy `local.example` to the ignored `local.jl` file and set the host, shared working directory, local Windows scratch directory, Julia executable, Python executable, and transport. The shared working directory is the Windows path for `test/gauntlet/cases/.work`. The final expression in that file must return `RemoteConfig`.

`RemoteConfig` fields have these meanings:

| Field | Meaning |
|:--|:--|
| `host` | Host name accepted by SSH or the selected wrapper. |
| `shared_root` | Windows path that points to `test/gauntlet/cases/.work` through the shared filesystem. |
| `remote_root` | Local Windows scratch directory used while PSCAD runs. |
| `julia_executable` | Julia executable on the Windows worker. |
| `python_executable` | Python executable containing `mhi.pscad`. |
| `pscad_version` | Supported PSCAD version. It must be `5.1.0`. |
| `transport` | `:ssh` or a symbol implemented by a method in `local.jl`. |
| `verbosity` | `0` for errors, `1` for milestones, or `2` for streamed runner output. |
| `timeout_seconds` | Maximum duration of one remote PSCAD calculation. |

Set `LINECABLEMODELS_GAUNTLET_CONFIG` when the configuration file is stored elsewhere. Relative paths are resolved from the repository root. Live and record modes load this file; snapshot mode does not.

The PSCAD host must have PSCAD 5.1.0, a working license, Julia 1.12, PythonCall 0.9, and `mhi.pscad` 3.1.2. Instantiate the remote runner project once on that host from the repository root:

```powershell
julia --project=test/gauntlet/pscad/remote --startup-file=no -e 'using Pkg; Pkg.instantiate()'
```

Set `JULIA_PYTHONCALL_EXE` while instantiating if PythonCall does not already use the Python installation containing `mhi.pscad`:

```powershell
$env:JULIA_PYTHONCALL_EXE = "C:\Python311\python.exe"
julia --project=test/gauntlet/pscad/remote --startup-file=no -e 'using Pkg; Pkg.instantiate(); using PythonCall'
```

Authentication belongs to the installed SSH or Tailscale tooling. The gauntlet does not acquire licenses, passwords, certificates, or tunnels.

## SSH, Tailscale, and VirtioFS notes

The commands in this section configure the Linux host and Windows PSCAD worker. The gauntlet does not run them. Snapshot mode does not use any of this setup.

### Ordinary SSH

Use `transport = :ssh` when the Windows worker is reachable through the system `ssh` command. Confirm the route from the Linux host before starting a live case:

```bash
ssh pscad-host powershell.exe -NoProfile -NonInteractive -Command '$PSVersionTable.PSVersion'
```

The matching `local.jl` ends with:

```julia
RemoteConfig(
    "pscad-host",
    raw"Z:\path\to\LineCableModels\test\gauntlet\cases\.work",
    raw"C:\LineCableModelsGauntlet",
    raw"C:\path\to\julia.exe",
    raw"C:\path\to\python.exe";
    transport = :ssh,
    verbosity = 1,
    timeout_seconds = 1800,
)
```

OpenSSH owns authentication and host-key verification. The gauntlet does not accept passwords or private-key material.

### Direct Tailscale SSH

The local `ts` wrapper can expose the VM's QEMU user-network SSH port on `127.0.0.1:10022`. Recreate the current `win11-vm` entry from a Linux terminal with:

```bash
ts ssh direct set win11-vm \
    --libvirt-user-domain win11 \
    --libvirt-uri qemu:///system \
    --host-address 127.0.0.1 \
    --port 10022 \
    --guest-port 22 \
    --start-policy on-demand \
    --boot-timeout 180
```

Inspect the entry and test the connection with:

```bash
ts ssh direct show win11-vm
ts ssh win11-vm --direct -- cmd.exe /d /s /c ver
```

The `--` before the remote command is required. Without it, an argument such as PowerShell's `-NoProfile` can be parsed as an SSH option.

With `--direct`, `ts` uses the local libvirt forward at `127.0.0.1:10022`. File and command traffic does not cross the Tailscale network.

The matching ignored `local.jl` configuration is:

```julia
function remote_command(
        ::Val{:tailscale},
        config::RemoteConfig,
        powershell::AbstractString
)
    Cmd(vcat(
        [
            "ts", "ssh", config.host, "--direct",
            "-o", "ServerAliveInterval=15",
            "-o", "ServerAliveCountMax=3",
            "--",
        ],
        _powershell_argv(powershell),
    ))
end

RemoteConfig(
    "win11-vm",
    raw"Z:\Documents\KUL\LineCableModels\test\gauntlet\cases\.work",
    raw"C:\Users\Amauri\AppData\Local\LineCableModels\Gauntlet\Remote",
    raw"C:\Users\Amauri\AppData\Local\Microsoft\WindowsApps\julia.exe",
    raw"C:\Users\Amauri\AppData\Local\LineCableModels\Gauntlet\PSCAD\py3.14\Scripts\python.exe";
    pscad_version = "5.1.0",
    transport = :tailscale,
    verbosity = 2,
    timeout_seconds = 1800,
)
```

Change the executable paths if Julia or Python moves. `verbosity = 2` prints runner progress while PSCAD is active. Use `1` for major local progress messages or `0` for errors only. `timeout_seconds` bounds one PSCAD run. Authentication remains under the installed `ts` and OpenSSH configuration.

Any other wrapper uses the same narrow seam. Set `transport` to a new symbol and define `remote_command(::Val{:symbol}, config, powershell)` in `local.jl`. The method must return the complete `Cmd` and pass `_powershell_argv(powershell)` unchanged to the Windows host. There is no transport fallback.

### VirtioFS host share

The `win11` libvirt domain exposes `/home/amartins` to Windows with the mount tag `home`. Windows mounts it as `Z:`. Inspect the active domain from Linux with:

```bash
virsh -c qemu:///system dumpxml win11
```

The domain must contain shared memory and the VirtioFS device:

```xml
<memoryBacking>
  <source type='memfd'/>
  <access mode='shared'/>
</memoryBacking>

<filesystem type='mount' accessmode='passthrough'>
  <driver type='virtiofs' queue='1024'/>
  <binary path='/usr/libexec/virtiofsd'/>
  <source dir='/home/amartins'/>
  <target dir='home'/>
</filesystem>
```

Check the mount and its services from an elevated PowerShell session in Windows:

```powershell
Get-CimInstance Win32_LogicalDisk |
    Where-Object DeviceID -eq "Z:" |
    Select-Object DeviceID, DriveType, FileSystem, ProviderName, VolumeName |
    Format-List

Get-Service VirtioFsSvc, WinFsp.Launcher
sc.exe query VirtioFsDrv

Measure-Command {
    Get-ChildItem "Z:\Documents\KUL\LineCableModels" -Force | Out-Null
}
```

`Z:` should appear as a fixed local drive with no provider name. `VirtioFsSvc`, `WinFsp.Launcher`, and `VirtioFsDrv` must be running.

### Lean Windows worker setup

The shared `Z:` tree belongs to the Linux host. Defender and Windows Search should inspect local Windows storage and leave `Z:` alone. Run the following once from an elevated PowerShell session:

```powershell
$ErrorActionPreference = "Stop"
$share = "Z:\"

if (-not (Test-Path -LiteralPath $share)) {
    throw "VirtioFS share is unavailable at $share"
}

$defender = Get-MpPreference
if ($defender.ExclusionPath -notcontains $share) {
    Add-MpPreference -ExclusionPath $share
}

$searchPolicy =
    "HKLM:\SOFTWARE\Policies\Microsoft\Windows\Windows Search\PreventIndexingCertainPaths"

New-Item -Path $searchPolicy -Force | Out-Null
New-ItemProperty `
    -Path $searchPolicy `
    -Name "1" `
    -PropertyType String `
    -Value "file:///Z:\*" `
    -Force | Out-Null

& sc.exe config VirtioFsSvc `
    start= auto `
    depend= "WinFsp.Launcher/VirtioFsDrv"

if ($LASTEXITCODE -ne 0) {
    throw "Could not configure VirtioFsSvc"
}

Restart-Service WSearch -Force

Write-Host "`nDefender exclusions:"
(Get-MpPreference).ExclusionPath

Write-Host "`nWindows Search exclusions:"
Get-ItemProperty -Path $searchPolicy

Write-Host "`nVirtioFS service:"
sc.exe qc VirtioFsSvc
```

Reboot Windows after saving active work. `VirtioFsSvc` must report `AUTO_START` without `DELAYED`, with `WinFsp.Launcher` and `VirtioFsDrv` as dependencies.

Defender continues scanning local files under `C:`. It does not scan files opened directly from `Z:`. Copy an untrusted file to local storage before opening it so Defender can inspect it.

Undo the Defender and Search exclusions with:

```powershell
Remove-MpPreference -ExclusionPath "Z:\"

Remove-ItemProperty `
    -Path "HKLM:\SOFTWARE\Policies\Microsoft\Windows\Windows Search\PreventIndexingCertainPaths" `
    -Name "1"

Restart-Service WSearch -Force
```

VirtioFS for Windows is documented by the [virtio-win project](https://github.com/virtio-win/kvm-guest-drivers-windows/wiki/Virtiofs%3A-Shared-file-system) and the [virtio-fs Windows guide](https://virtio-fs.gitlab.io/howto-windows.html). Microsoft documents [Defender path exclusions](https://learn.microsoft.com/en-us/defender-endpoint/microsoft-defender-antivirus-exclusions-configure) and [Windows Search path exclusions](https://learn.microsoft.com/en-us/windows/win32/search/-search-3x-wds-extidx-csm-scoperules).

### Remote file and process lifecycle

SSH carries commands and text logs only. The gauntlet writes the generated project and fixed runner files below `cases/.work/<case>/<variant>/`. Windows sees the same directory through `Z:` and copies it once into the local `C:` scratch directory. PSCAD runs from `C:`. The Windows supervisor copies the completed outputs back to `Z:` before it exits.

The supervisor records the Julia runner PID and its exact script path in `owner.txt`. A timeout stops that PID tree with `taskkill /T`, never by executable name. A local interruption sends the same targeted cancellation command before returning control to Julia. The next run checks any remaining owner file before replacing the scratch directory. It refuses to terminate a PID whose command line does not contain the recorded case runner path.

This covers ordinary errors, `Ctrl+C`, and a hung line-constants calculation. Killing the local Julia process without allowing cleanup can leave the remote run alive until its configured timeout. Starting that case again performs the same ownership check and removes the stale process first.

At debug verbosity, the terminal shows each runner milestone as it is written. PSCAD 5.1 exposes line-constants calculation through a blocking `compile()` call. The runner therefore prints PSCAD's accumulated project messages as soon as that call returns; it cannot stream messages that PSCAD has not yet returned through MHI. Complete logs remain in `stdout.txt`, `stderr.txt`, and `pscad-console.txt`.

## Commands

Run Example 3 against PSCAD:

```bash
LINECABLEMODELS_GAUNTLET_MODE=live \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "test/gauntlet/cases/example3.jl"); include("test/runtests.jl")'
```

This requires the live setup. It writes `test/gauntlet/cases/.work/example3/` and cannot replace the snapshot.

Run every gauntlet case against PSCAD:

```bash
LINECABLEMODELS_GAUNTLET_MODE=live \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "tag:gauntlet"); include("test/runtests.jl")'
```

This requires the live setup. Each case writes its own directory below `cases/.work/`. No snapshot is replaced.

Record or replace the Example 3 snapshot:

```bash
LINECABLEMODELS_GAUNTLET_PERSIST=true \
LINECABLEMODELS_GAUNTLET_MODE=record \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "test/gauntlet/cases/example3.jl"); include("test/runtests.jl")'
```

This requires the live setup. Both environment variables are mandatory. The command replaces the ignored staging collection under `.artifacts/example3/` and updates the local tree binding in `Artifacts.toml` only after PSCAD execution, parsing, structural checks, numerical comparisons, and case tolerances succeed.

Record every gauntlet case and bind all local artifacts with one package-test run:

```bash
LINECABLEMODELS_GAUNTLET_PERSIST=true \
LINECABLEMODELS_GAUNTLET_MODE=record \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "tag:gauntlet"); include("test/runtests.jl")'
```

The staged collection contains `snapshot.jld2`, its direct SHA-256 file, and `example3.tar.gz`. The JLD2 file records the problem, a declarative formulation record, the PSCAD reference, the accepted LineCableModels result, the reference comparison, the PSCAD wall time, and the local benchmark. The artifact tree hash is handled by `Pkg.Artifacts`; there are no campaign hashes or hash-derived directories.

Artifact names are stable, such as `pscad_gauntlet_example3`. Recording a new accepted result changes the tree hash bound in `Artifacts.toml`. Git history retains the previous binding. `SNAPSHOT_FORMAT_VERSION` changes only when the JLD2 schema is incompatible. The recording timestamp remains metadata and never becomes part of the artifact name.

After uploading `example3.tar.gz`, bind its real download URL:

```bash
test -n "$GAUNTLET_ARTIFACT_URL"
julia --project=test/gauntlet --startup-file=no -e \
'using TestItemRunner; include("test/gauntlet/support.jl"); using .GauntletSupport; publish_artifact(:example3, ARGS[1])' \
"$GAUNTLET_ARTIFACT_URL"
```

Set `GAUNTLET_ARTIFACT_URL` to the real immutable archive URL first. The command calculates the archive SHA-256 and adds the download entry to `test/gauntlet/Artifacts.toml`. It fails if the case was not recorded locally or the archive is missing. Do not commit an invented URL.

Run Example 3 from its snapshot:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "test/gauntlet/cases/example3.jl"); include("test/runtests.jl")'
```

This resolves `pscad_gauntlet_example3` through `Artifacts.toml`. It does not require PSCAD, Python, SSH, Tailscale, or `local.jl`. A missing binding, unavailable download, malformed artifact, digest mismatch, or stale case definition fails immediately.

Run every snapshot case locally:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no -e \
'push!(ARGS, "tag:gauntlet"); include("test/runtests.jl")'
```

Run the same snapshot selection used by CI:

```bash
CI=true \
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no --code-coverage=@. -e \
'push!(ARGS, "tag:gauntlet"); include("test/runtests.jl")'
```

CI rejects `live` and `record` instead of changing them to `snapshot`.

## Inspecting an accepted result

Load LineCableModels before JLD2 so stored package types are available during reconstruction:

```bash
julia --project=test/gauntlet --startup-file=no
```

```julia
using LineCableModels
using DataFrames
using JLD2

snapshot = JLD2.load(
    "test/gauntlet/.artifacts/example3/snapshot.jld2",
)
comparison = snapshot["reference_comparison"]
```

The path above is the ignored staging copy created by `record`. To inspect an
installed or downloaded artifact instead, resolve the binding explicitly:

```julia
using Pkg.Artifacts

artifacts_toml = "test/gauntlet/Artifacts.toml"
hash = artifact_hash("pscad_gauntlet_example3", artifacts_toml)
hash === nothing && error("example3 is not bound in $artifacts_toml")
artifact_exists(hash) ||
    ensure_artifact_installed("pscad_gauntlet_example3", artifacts_toml)
snapshot = JLD2.load(joinpath(artifact_path(hash), "snapshot.jld2"))
comparison = snapshot["reference_comparison"]
```

The raw comparison keeps one absolute and one relative RMS matrix for Z and Y:

```julia
comparison.Z.absolute
comparison.Z.relative
comparison.Y.absolute
comparison.Y.relative
```

Use the case's absolute reference tolerances as display floors. `DataFrame` reports the relative value as `missing` when the absolute error is already at or below that floor, and `relative_status` explains that the relative metric is not meaningful there. The raw `LineParametersComparison` remains unchanged.

```julia
errors = DataFrame(
    comparison;
    zero_atol = (Z = 1.0e-6, Y = 1.0e-9),
)
sort!(errors, [:quantity, :rms_absolute], rev = [false, true])
display(errors)
```

Inspect the accepted timing record with:

```julia
benchmark = snapshot["julia_benchmark"]
display((
    julia_minimum_ms = 1.0e3 * benchmark.minimum_seconds,
    julia_median_ms = 1.0e3 * benchmark.median_seconds,
    bytes = benchmark.bytes,
    allocations = benchmark.allocations,
    samples = benchmark.samples,
    pscad_seconds = snapshot["pscad_elapsed_seconds"],
    environment = benchmark.environment,
))
```

The generated [PSCAD gauntlet developer guide](../../docs/src/gauntlet.md) runs the same inspection with a small deterministic fixture.

## Files written by live execution

The current exporter writes below `cases/.work/<case>/current/`. The directory retains the generated PSCX project, staged runner files, PSCAD matrix outputs, `stdout.txt`, `stderr.txt`, `transport-stdout.txt`, `transport-stderr.txt`, `pscad-console.txt`, and `timing.txt`. A nonzero remote exit preserves the shared files and the local Windows scratch directory for inspection.

The stored case SHA-256 binds the snapshot to the bytes of its case file. The separate snapshot SHA-256 detects damage to the JLD2 payload. Either mismatch requires review and an explicit record run; neither can trigger regeneration automatically.
