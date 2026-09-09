# Approved executor images

These are **image recipes and entry checks**, not running-isolation certification.
Managed scientific launch, agent-death cleanup, remote preparation/jobs, terminal
lease ownership, the bounded authorized relay and reusable browser component are
implemented. Finite native fixtures verify the terminal transport; production
terminal availability still requires effective container acceptance. No image
is built, pulled or started by registration or public-page rendering. See the
[runtime ledger](../../RUNTIME_PLATFORM_PROGRESS.md) for current evidence.

The three final targets share one guarded Julia entry contract:

| Target / approved profile ID | Kind | Installed environment |
|---|---|---|
| `line-parameters` | scientific | LineCableModels and the existing line adapters |
| `power-flow` | scientific | Pinned PowerImpedance and the existing OHL/UGC adapters |
| `julia-terminal` | terminal | Julia stdlibs and lightweight execution core |

Build from the repository root on a provisioned build machine. These are explicit
operator actions and may download the pinned dependencies; the runtime never
performs them. Choose the actual engine on that host:

```sh
podman build --target line-parameters -f playground/worker/containers/Containerfile -t localhost/lcm-line-parameters:dev .
podman build --target power-flow -f playground/worker/containers/Containerfile -t localhost/lcm-power-flow:dev .
podman build --target julia-terminal -f playground/worker/containers/Containerfile -t localhost/lcm-julia-terminal:dev .
```

On a Docker host use the same targets and Containerfile with `docker build`.
Provision approved images through the operator's chosen registry/transfer
procedure. Registration requires the repository reference **and manifest digest**
present in the local image's `RepoDigests`; a mutable tag or config ID is not an
accepted substitute. No publishing, registry credentials or registry choice is
automatically authorized by these recipes. The profile fingerprint is that
approved manifest digest, not the native-source fingerprint.

The active project is a real directory at
`/opt/lcm/source/playground/worker/profiles/active`, at the same depth as the
original profile. This preserves relative Manifest paths and shared source
includes. A short symlinked active project does not preserve Julia's lexical
dependency resolution and is deliberately avoided. Image construction copies
only declared source/lockfiles, never the whole checkout, runtime state, user home
or deployment credentials. Package compilation is not model preparation.

## One policy, two separate kinds of evidence

`runtime/src/ContainerPolicy.jl` renders the only scientific/terminal create
policy. Host prerequisites and approved local image identity are checked before
creation, and all configuration is checked again on the stopped resource. A
durable intent survives partial acquisition. The driver must release it through
receipt-verified cleanup; this helper never starts or automatically retries it.

The common fixed Julia command loads `container-guard.jl` before a scientific
script or the interactive REPL. The guard checks actual process identity,
capabilities, seccomp/no-new-privileges, cgroup limits, network interfaces, mounts,
tmpfs capacity and process limits. Missing or incompatible evidence exits 78 with
a fixed diagnostic before evaluating user code. These checks use
[Linux cgroup v2](https://docs.kernel.org/6.12/admin-guide/cgroup-v2.html) and
[mountinfo](https://man7.org/linux/man-pages/man5/proc_pid_mountinfo.5.html), not a
ready flag inferred from a successful CLI return.

For terminal acquisitions, the shared create policy supplies a fresh
`LCM_TERMINAL_READY` receipt UUID. Only after successful isolation verification
does the guard load `terminal-ready.jl`, which emits the corresponding marker
from Julia's REPL initialization hook. The session owner strips that startup
prefix and rechecks physical policy/liveness before admitting input. Existing
images must be rebuilt and approved by their new digest to supply this hook;
no automatic image pull or upgrade is performed.

Required policy: UID/GID 1000; dropped capabilities; no network; private PID, IPC
and cgroup namespaces; read-only image; no host bind mounts, sockets or devices;
finite CPU/memory/PID limits; zero additional swap; finite noexec/nosuid/nodev
scratch; no engine log file or automatic restart. The 64 KiB shared-memory
allowance is included in the scratch budget. File descriptors are capped at
1024; core files and POSIX message queues are disabled. The mutable Julia depot
and home live in disposable scratch; the installed depot is read-only. The REPL
is not a package-installation service.

Docker and Podman inspection differences are normalized only inside the shared
adapter. In particular, Docker's private PID mode is empty (not `private`), while
Podman expands dropped capabilities and normalizes tmpfs/rlimit fields. Podman's
`notmpcopyup` prevents copying image `/tmp` content into new scratch. See the
[Docker namespace definition](https://github.com/moby/moby/blob/master/api/types/container/hostconfig.go)
and [Podman create options](https://docs.podman.io/en/stable/markdown/podman-create.1.html).

Neither inspection proves that an arbitrary image is trustworthy or makes a
shared Linux kernel an absolute untrusted-code boundary. Images and engines are
operator-approved. Lease authority, interruption, agent death and exact cleanup
still belong to the physical owner; the guard does not replace them.

## Current verification boundary

`runtime/test/container_image_layout.jl` reconstructs the recipe's COPY layout,
checks local dependency paths, compares every image command with the shared policy
and imports the core with actual Julia from each relocated environment.
`worker/core/test/container_entry.jl` proves missing/host isolation denies entry
before evaluation. Unit fixtures cover effective-limit failures and both engine
inspection schemas. These are not OCI image-build or running-isolation passes.

The optional `runtime/test/container_stopped_policy.jl` gate requires
`LCM_TEST_STOPPED_CONTAINER_IMAGE` set to an already cached digest-pinned test
image. It creates and inspects a stopped resource for each process kind, then
removes only its exact receipt-verified IDs. It never starts a process, so it can
check configuration compatibility even on a host missing delegated CPU. It does
not bypass production image admission: an arbitrary base image is used only in
this explicitly create-only test. Failed test diagnostics stay private in the
reported temporary directory.

Actual image builds and running kernel enforcement require a host with all
mandatory limits. The local IT-managed Podman-only development host lacks
delegated CPU. The separately provisioned owned Kubuntu host supplies those
controls: both rootless engines have passed the full kernel/PTY and managed
scientific/terminal matrices, including graceful stop and forced agent death.
Separate-computer presentation/workbench/terminal/S3 runs also pass on each.
The [physical acceptance procedure](../../runtime/PHYSICAL_ACCEPTANCE.md)
and [ledger](../../RUNTIME_PLATFORM_PROGRESS.md) retain the distinct evidence
and exact cleanup results. These private rehearsals are not a public deployment.
Do not remove the CPU
requirement or substitute process affinity on an unsupported host.
