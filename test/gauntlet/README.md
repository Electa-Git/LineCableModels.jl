# Validation gauntlet

Each file in `cases/` owns one complete validation case. The case declares its reference backend, separate reference and LineCableModels problems and formulations, exact frequency samples, explicit terminal order, tolerances, and assertions. Both calculations use the ordinary `ProblemDefinition`, `Formulation`, and `compute!` grammar. Gauntlet cases retain every terminal: phase assignments must be the contiguous integers `1:n`, and both Kron and bundle reduction must be disabled. Disposable exports and diagnostics remain in `cases/.work/<backend>/<case>/`.

The gauntlet has three modes. `snapshot` reads an existing reference and never loads host configuration, starts a process, initializes Python, or contacts an external backend. `live` executes each case's declared backend without changing the reference. `record` completes the live checks before replacing the reference. A failed operation never switches modes.

| Mode | Live backend | Network | Reads accepted artifact | May replace artifact |
|:--|:--:|:--:|:--:|:--:|
| `snapshot` | No | No | Yes | No |
| `live` | Yes | Yes | No | No |
| `record` | Yes | Yes | No | Yes, with persistence enabled; an existing version also requires force |

## Comparisons

The live calculation compares the phase-domain `Z[i,j,:]` and `Y[i,j,:]` series returned by the reference formulation with the corresponding LineCableModels series. Absolute and relative RMS errors are calculated separately for every matrix term across frequency. Gauntlet does not flatten the matrices, combine matrix terms, reorder terminals, interpolate frequencies, or post-process external results.

The current PSCAD frequency-dependent model does not emit modal `Z` and `Y`. A modal `LineParameters` result therefore throws a clear not-implemented error instead of being compared through a derived proxy.

Each case has three metric groups:

- `reference` reports LineCableModels against PSCAD. It never passes or fails a test and never blocks recording.
- `regression` checks the current LineCableModels result against the accepted result stored in the artifact.
- `performance` checks time, bytes, and allocations against the accepted artifact only when the Julia, operating system, architecture, thread count, and BLAS configuration are identical.

Performance from different environments is reported as non-comparable and cannot pass or fail the case.

Interactive runs display the reference RMS table. Relative values whose absolute RMS is already below the case's display floor appear as `missing`; the underlying comparison remains unchanged.

## Problem and formulation selection

The case constructs the physical problem once, then creates separate problem and formulation objects for PSCAD and the owned EMT solver:

```julia
reference_problem = LineParametersProblem(
    problem.system;
    temperature = problem.temperature,
    earth_props = problem.earth_props,
    frequencies = problem.frequencies,
)

reference_formulation = Formulation(
    :pscad;
    earth_impedance = EarthImpedance.Wedepohl(),
)

formulation = Formulation(
    :EMT;
    earth_impedance = EarthImpedance.Pollaczek(),
    earth_admittance = EarthAdmittance.IdealGround(),
    insulation_admittance = InsulationAdmittance.Lossless(),
    options = (
        kron_reduction = false,
        reduce_bundle = false,
        ideal_transposition = false,
    ),
)
```

`Formulation(:pscad)` accepts these explicit earth-impedance method objects:

| Model | Method object | PSCAD field |
|:--|:--|:--|
| Overhead | `EarthImpedance.Deri()` | `EarthForm2` |
| Overhead | `EarthImpedance.DirectNumericalIntegration(:overhead)` | `EarthForm2` |
| Underground | `EarthImpedance.Wedepohl()` | `EarthForm` |
| Underground | `EarthImpedance.DirectNumericalIntegration(:underground)` | `EarthForm` |
| Underground | `EarthImpedance.Saad()` | `EarthForm` |
| Aerial/underground mutual | `EarthImpedance.Ametani()` | `EarthForm3` |
| Aerial/underground mutual | `EarthImpedance.Lucca()` | `EarthForm3` |

The PSCAD formulation uses `NativeEarthAdmittance()` and `NativeInsulationAdmittance()` for the calculations owned by PSCAD's line-data model. The EMT formulation selects its own earth and insulation admittance methods independently. No validator forces those choices to resemble each other. A deliberately mismatched pair runs normally and the RMS result exposes the consequence.

`Deri`, `Wedepohl`, `Saad`, `Ametani`, `Lucca`, and direct numerical integration are Engine-owned method concepts. They can be selected by external backends. Calling them through `Formulation(:EMT)` produces a clear not-implemented error because the owned EMT kernels do not currently implement those methods.

Live and record modes call the same public operation for both solvers:

```julia
reference = compute!(reference_problem, reference_formulation; options)
candidate = compute!(problem, formulation; options)
comparison = compare(reference, candidate)
```

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
| `timeout_seconds` | Maximum duration of one remote PSCAD calculation. |

Runner verbosity belongs to the ordinary Engine execution options, not the host configuration. A case can stream PSCAD milestones without changing local solver logging:

```julia
options = ComputeOptions(verbosity = (default = 0, PSCAD = 2))
outcome = run_case(case; options)
```

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
    timeout_seconds = 1800,
)
```

OpenSSH owns authentication and host-key verification. The gauntlet does not accept passwords or private-key material.

### Direct Tailscale SSH

The local `ts` wrapper can expose the VM's QEMU user-network SSH port on `127.0.0.1:10022`. Create an entry for the PSCAD worker from a Linux terminal with:

```bash
ts ssh direct set pscad-worker \
    --libvirt-user-domain pscad-worker \
    --libvirt-uri qemu:///system \
    --host-address 127.0.0.1 \
    --port 10022 \
    --guest-port 22 \
    --start-policy on-demand \
    --boot-timeout 180
```

Inspect the entry and test the connection with:

```bash
ts ssh direct show pscad-worker
ts ssh pscad-worker --direct -- cmd.exe /d /s /c ver
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
    "pscad-worker",
    raw"Z:\test\gauntlet\cases\.work",
    raw"C:\LineCableModelsGauntlet",
    raw"C:\path\to\julia.exe",
    raw"C:\path\to\python.exe";
    pscad_version = "5.1.0",
    transport = :tailscale,
    timeout_seconds = 1800,
)
```

Change the executable paths if Julia or Python moves. `timeout_seconds` bounds one PSCAD run. Set the `PSCAD` entry in `ComputeOptions.verbosity` to `2` to stream runner progress, `1` for milestones, or `0` for warnings only. Authentication remains under the installed `ts` and OpenSSH configuration.

Any other wrapper uses the same narrow seam. Set `transport` to a new symbol and define `remote_command(::Val{:symbol}, config, powershell)` in `local.jl`. The method must return the complete `Cmd` and pass `_powershell_argv(powershell)` unchanged to the Windows host. There is no transport fallback.

### VirtioFS host share

The `pscad-worker` libvirt domain exposes the repository to Windows with the mount tag `linecablemodels`. Windows mounts it as `Z:`. Inspect the active domain from Linux with:

```bash
virsh -c qemu:///system dumpxml pscad-worker
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
  <source dir='/path/to/LineCableModels'/>
  <target dir='linecablemodels'/>
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
    Get-ChildItem "Z:\" -Force | Out-Null
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

SSH carries commands and text logs only. The PSCAD formulation writes the generated project and fixed runner files below `cases/.work/pscad/<case>/reference/`. Windows sees the same directory through `Z:` and copies it once into the local `C:` scratch directory. PSCAD runs from `C:`. The Windows supervisor copies the completed outputs back to `Z:` before it exits.

The supervisor records the Julia runner PID and its exact script path in `owner.txt`. A timeout stops that PID tree with `taskkill /T`, never by executable name. A local interruption sends the same targeted cancellation command before returning control to Julia. The next run checks any remaining owner file before replacing the scratch directory. It refuses to terminate a PID whose command line does not contain the recorded case runner path.

This covers ordinary errors, `Ctrl+C`, and a hung line-constants calculation. Killing the local Julia process without allowing cleanup can leave the remote run alive until its configured timeout. Starting that case again performs the same ownership check and removes the stale process first.

At debug verbosity, the terminal shows each runner milestone as it is written. PSCAD 5.1 exposes line-constants calculation through a blocking `compile()` call. The runner therefore prints PSCAD's accumulated project messages as soon as that call returns; it cannot stream messages that PSCAD has not yet returned through MHI. Complete logs remain in `stdout.txt`, `stderr.txt`, and `pscad-console.txt`.

## Commands

The dedicated runner contains the native TestItemRunner selection:

```julia
@run_package_tests(filter=ti -> :gauntlet in ti.tags, verbose=true)
```

TestItemRunner discovers every case file. There is no file list to update.

Two optional boolean controls apply to the complete runner. Both accept only `true` or `false` and default to `false`.

| Variable | Effect |
|:--|:--|
| `LINECABLEMODELS_GAUNTLET_CLEANUP` | Removes only `test/gauntlet/cases/.work/` after every selected case succeeds. Failed runs retain their working files and diagnostics. Snapshots and backend archives are untouched. |
| `LINECABLEMODELS_GAUNTLET_FORCE` | Makes `record` rebuild the current Gauntlet version without loading its accepted snapshots, then replace its staged backend collections and `Artifacts.toml` bindings. It has no effect in `snapshot` or `live` mode. |

Record mode checks for a staged collection or artifact binding with the same backend name and Gauntlet version before it starts any case. The runner stops and lists every collision unless `LINECABLEMODELS_GAUNTLET_FORCE=true` is set. A forced run deletes only the current version's ignored staging directories and does not use the bound collection for regression or performance checks. It replaces the artifact bindings only after all cases finish successfully.

Run every case against its declared live backend:

```bash
LINECABLEMODELS_GAUNTLET_MODE=live \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

This requires the setup for every selected backend. Each case writes below `cases/.work/<backend>/<case>/`. No accepted collection is replaced.

Record every case and build one collection per backend:

```bash
LINECABLEMODELS_GAUNTLET_PERSIST=true \
LINECABLEMODELS_GAUNTLET_MODE=record \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

Both environment variables are mandatory. Record mode is accepted only through this runner. Each case first writes its snapshot into ignored staging. TestItemRunner must finish every selected case successfully before the runner creates or binds any backend archive. A reference comparison remains diagnostic and cannot block recording. Snapshot mode applies regression and comparable performance limits against the accepted collection.

Replace an existing collection with the same name and version, then remove the disposable working files after the runner exits:

```bash
LINECABLEMODELS_GAUNTLET_CLEANUP=true \
LINECABLEMODELS_GAUNTLET_FORCE=true \
LINECABLEMODELS_GAUNTLET_PERSIST=true \
LINECABLEMODELS_GAUNTLET_MODE=record \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

Use cleanup without force when the collection does not already exist. Use force without cleanup when the generated PSCAD projects and logs should remain available for inspection.

The source-owned version is `GAUNTLET_VERSION = v"1.0.0"`. Staging is grouped by backend:

```text
test/gauntlet/.artifacts/
└── pscad/
    └── v1.0.0/
        ├── cases/
        │   ├── benchmark_525kV_1600mm2_bipole_pscad/
        │   │   ├── snapshot.jld2
        │   │   └── snapshot.sha256
        │   └── ...
        ├── report.jld2
        ├── report.tsv
        ├── report.sha256
        └── benchmarks-pscad-v1.0.0.tar.gz
```

Every case must declare one lowercase backend symbol. One backend archive contains all snapshots recorded for that backend and version. `Pkg.Artifacts` binds the complete collection under a name such as `gauntlet_pscad_v1_0_0`. A snapshot stores the backend and Gauntlet version along with both problems, declarative formulation records, the external reference, the accepted LineCableModels result, comparisons, and timings. Finalization recomputes each stored comparison and writes one aggregate row per case to `report.jld2` and `report.tsv`.

The Gauntlet version is independent of the package version. Bump the patch slot for a corrected rerun with unchanged stored fields and KPIs. Bump the minor slot when adding retained data, cases, or a backend without invalidating existing readers. Bump the major slot when changing snapshot fields, KPI meaning, or loading behavior. Every change to an already published immutable collection requires a new version.

Use one immutable GitHub release for the Gauntlet version and one asset per backend:

```bash
REPO="Electa-Git/LineCableModels.jl"
BACKEND="pscad"
VERSION="$(julia --project=test/gauntlet --startup-file=no -e \
    'include("test/gauntlet/artifacts.jl"); using .GauntletArtifacts; print(GAUNTLET_VERSION)')"
TAG="gauntlet-v${VERSION}"
ASSET="benchmarks-${BACKEND}-v${VERSION}.tar.gz"
REPO_URL="https://github.com/$REPO"
GAUNTLET_ARTIFACT_URL="$REPO_URL/releases/download/$TAG/$ASSET"
```

Enable [GitHub release immutability](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases). Create `gauntlet-v1.0.0` as a draft, upload all backend assets, then publish it. Published immutable assets cannot be replaced or deleted individually. If an accepted publication is wrong, leave it intact and publish the correction under the next Gauntlet version.

Commit the reviewed case files and `Artifacts.toml` before creating the release tag. Create and push the annotated tag from that exact commit:

```bash
git tag -a "$TAG" -m "Gauntlet v${VERSION}"
git push origin "$TAG"
```

Create the draft and upload the recorded PSCAD collection with:

```bash
gh release create "$TAG" \
    "test/gauntlet/.artifacts/${BACKEND}/v${VERSION}/${ASSET}" \
    --repo "$REPO" \
    --verify-tag \
    --draft \
    --title "Gauntlet v${VERSION}" \
    --notes "Accepted validation references for Gauntlet v${VERSION}."
```

Upload only backend archives such as `benchmarks-pscad-v1.0.0.tar.gz`. Each archive already contains its aggregate report and report digest. Upload every backend archive to the same draft before publishing it. Publish the draft only after checking its tag target, asset list, archive SHA-256, and extracted report.

After uploading a backend archive, bind its real download URL:

```bash
test -n "$GAUNTLET_ARTIFACT_URL"
julia --project=test/gauntlet --startup-file=no -e \
'include("test/gauntlet/artifacts.jl"); using .GauntletArtifacts; publish_artifact(Symbol(ARGS[1]), ARGS[2])' \
"$BACKEND" "$GAUNTLET_ARTIFACT_URL"
```

The command calculates the archive SHA-256 and adds the download entry to `test/gauntlet/Artifacts.toml`. It fails if the backend collection was not recorded locally or its archive is missing. Do not commit an invented URL.

Run every case from the accepted backend collections:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

Each case resolves the collection named by its backend and `GAUNTLET_VERSION`, then loads its human-readable path below `cases/`. Snapshot mode does not require PSCAD, Python, SSH, Tailscale, or `local.jl`. A missing binding, missing case, unavailable download, malformed snapshot, digest mismatch, backend mismatch, version mismatch, or stale case definition fails immediately.

Run the same snapshot selection used by CI:

```bash
CI=true \
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no --code-coverage=@. \
test/gauntlet/runtests.jl
```

CI rejects `live` and `record` instead of changing them to `snapshot`.

## Aggregate report

Display one row per recorded PSCAD case from the staged collection:

```bash
julia --project=test/gauntlet --startup-file=no \
test/gauntlet/report.jl pscad
```

This command does not execute a case or contact the reference backend. It verifies every snapshot and digest, requires a successful recorded reference exit status, checks frequencies and terminal dimensions, and recomputes each stored `LineParametersBenchmark` from the stored reference and accepted results.

The report retains the raw maximum relative RMS for Z and Y. Its display-safe relative columns exclude terms whose absolute RMS is at or below `Z = 1e-6` Ω/m or `Y = 1e-9` S/m. The absolute RMS, raw relative RMS, and complete per-term comparison remain unchanged.

For the accepted PSCAD Gauntlet v1.0.0 checkpoint:

| Evidence | Value |
|:--|:--|
| Cases | 7 |
| Frequency samples per case | 101 |
| Snapshot replay | 65/65 assertions passed |
| Artifact tree SHA-1 | `ec9bd7850a744dad14a9547ad4c0e0555363459a` |
| Archive SHA-256 | `02972c090884e53fa5ea6b165516a0b715a86ebdde36520f17e523dad8046bf3` |

The accepted status means that the reference runs completed, the stored evidence passed the report integrity checks, and the current code reproduced every accepted LineCableModels result within the case-owned regression and performance limits. External-reference RMS values remain reported measurements rather than pass/fail criteria.

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
    "test/gauntlet/.artifacts/pscad/v1.0.0/cases/" *
    "benchmark_525kV_1600mm2_bipole_pscad/snapshot.jld2",
)
comparison = snapshot["reference_comparison"]
```

The path above is the ignored staging copy created by `record`. To inspect an
installed or downloaded artifact instead, resolve the binding explicitly:

```julia
using Pkg.Artifacts

artifacts_toml = "test/gauntlet/Artifacts.toml"
artifact = "gauntlet_pscad_v1_0_0"
hash = artifact_hash(artifact, artifacts_toml)
hash === nothing && error("$artifact is not bound in $artifacts_toml")
artifact_exists(hash) ||
    ensure_artifact_installed(artifact, artifacts_toml)
snapshot = JLD2.load(joinpath(
    artifact_path(hash),
    "cases",
    "benchmark_525kV_1600mm2_bipole_pscad",
    "snapshot.jld2",
))
comparison = snapshot["reference_comparison"]
```

The raw comparison keeps one absolute and one relative RMS matrix for Z and Y:

```julia
comparison.Z.absolute
comparison.Z.relative
comparison.Y.absolute
comparison.Y.relative
```

Use the case's absolute reference tolerances as display floors. `DataFrame` reports the relative value as `missing` when the absolute error is already at or below that floor, and `relative_status` explains that the relative metric is not meaningful there. The raw `LineParametersBenchmark` remains unchanged.

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
    reference_seconds = snapshot["reference_execution"].elapsed_seconds,
    environment = benchmark.environment,
))
```

The generated [PSCAD gauntlet developer guide](../../docs/src/gauntlet.md) runs the same inspection with a small deterministic fixture.

## Files written by live execution

The PSCAD formulation writes below `cases/.work/pscad/<case>/reference/`. The directory retains the generated PSCX project, staged runner files, PSCAD matrix outputs, `stdout.txt`, `stderr.txt`, `transport-stdout.txt`, `transport-stderr.txt`, `pscad-console.txt`, and `timing.txt`. A nonzero remote exit preserves the shared files and the local Windows scratch directory for inspection.

The stored case SHA-256 binds the snapshot to the bytes of its case file. The separate snapshot SHA-256 detects damage to the JLD2 payload. Either mismatch requires review and an explicit record run; neither can trigger regeneration automatically.
