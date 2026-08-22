# PSCAD validation gauntlet

Each file in `cases/` owns one complete validation case. The case defines the LineCableModels problem, formulation, frequency samples, retained-conductor order, reduction state, tolerances, and assertions. Stable PSCAD results are stored as human-named JLD2 files in `.artifacts/`. Disposable exports and diagnostics remain in `cases/.work/<case>/`.

The gauntlet has three modes. `snapshot` reads an existing reference and never loads the PSCAD code. `live` executes PSCAD without changing the reference. `record` completes the live checks before replacing the reference. A failed operation never switches modes.

## Live setup

Copy `local.example` to the ignored `local.jl` file and set the host, shared working directory, local Windows scratch directory, Julia executable, Python executable, and transport. The shared working directory is the Windows path for `test/gauntlet/cases/.work`. The final expression in that file must return `RemoteConfig`.

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
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(test_args=["test/gauntlet/cases/example3.jl"])'
```

This requires the live setup. It writes `test/gauntlet/cases/.work/example3/` and cannot replace the snapshot.

Run every gauntlet case against PSCAD:

```bash
LINECABLEMODELS_GAUNTLET_MODE=live \
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(test_args=["tag:gauntlet"])'
```

This requires the live setup. Each case writes its own directory below `cases/.work/`. No snapshot is replaced.

Record or replace the Example 3 snapshot:

```bash
LINECABLEMODELS_GAUNTLET_MODE=record \
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(test_args=["test/gauntlet/cases/example3.jl"])'
```

This requires the live setup. It replaces `.artifacts/example3.jld2` only after PSCAD execution, parsing, structural checks, numerical comparisons, and case tolerances succeed.

Run Example 3 from its snapshot:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(test_args=["test/gauntlet/cases/example3.jl"])'
```

This does not require PSCAD, Python, SSH, Tailscale, or `local.jl`. A missing, malformed, or stale snapshot fails immediately.

Run every snapshot case locally:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(test_args=["tag:gauntlet"])'
```

Run the same snapshot selection used by CI:

```bash
CI=true \
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=. --startup-file=no -e \
'using Pkg; Pkg.test(coverage=true, test_args=["tag:gauntlet"])'
```

CI rejects `live` and `record` instead of changing them to `snapshot`.

## Files written by live execution

Current and legacy exports use separate subdirectories under `cases/.work/<case>/`. Each subdirectory retains the generated PSCX project, staged runner files, PSCAD matrix outputs, `stdout.txt`, `stderr.txt`, `transport-stdout.txt`, `transport-stderr.txt`, `pscad-console.txt`, and `timing.txt`. A nonzero remote exit preserves the shared files and the local Windows scratch directory for inspection.

Snapshots contain one `Engine.LineParameters` object and the minimum structural and execution metadata. The stored SHA-256 binds the snapshot to the bytes of its case file. A mismatch requires an explicit `record` run.
