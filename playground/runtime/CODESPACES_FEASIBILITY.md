# Codespaces / Docker feasibility

Historical scope note: the user subsequently restored full Docker and Podman
integration on an owned Kubuntu host. See the runtime execution ledger. The
feasibility evidence below remains unchanged; the shared Julia executable-path
finding has now been corrected in the recipe and both-engine policy.

Scope: disposable Docker proof of concept only. Full Docker integration and
managed-agent acceptance remain deferred; no production runtime policy changed.

The 2026-09-08 probe uses the existing `worker/core/src/ContainerIsolation.jl`
verbatim inside the repository's digest-pinned Julia 1.12.7 base image. No
application checkout, credentials, host mounts, scientific packages, broker,
gateway or worker service is deployed.

## Result: Docker feasibility passed, managed-agent hosting unavailable

The fresh 2-core / 8-GiB West Europe Codespace ran genuine Docker/Moby
`29.7.2-2`, Linux `6.8.0-1052-azure`, cgroup v2 with the `cgroupfs` driver,
and Julia `1.12.7`. The default Codespaces environment already supplied Docker;
no custom devcontainer or Docker installation was necessary.

Both plain and TTY-attached Julia processes passed the unchanged production
kernel-isolation verifier with these effective limits:

| Resource | Observed |
|---|---|
| CPU | 0.5 logical CPU; throttling observed during a bounded busy loop |
| RAM / extra swap | 512 MiB / zero |
| Tasks, including threads | 64 |
| Writable scratch, including shared memory | 8 MiB |
| UID / GID | 1000 / 1000 |

The verifier also checked dropped capabilities, seccomp, no-new-privileges,
read-only root and cgroups, isolated networking, allowed mounts and process
limits. Scratch writes worked; a root write was denied; no Docker socket was
visible inside either container. Both Julia processes returned `sum(1:100) = 5050`.

A negative-control container with no CPU quota exited **78** with
`cpu_limit_unverified`, without reaching the user-code marker. Docker-level
termination of a running test container worked. Label-verified removal of every
owned container restored the initially empty container inventory.

## Findings to carry into deferred integration

1. **No user-systemd manager.** Codespace PID 1 is `docker-init`; a direct user
   D-Bus manager query fails. The current managed-agent service/incarnation and
   post-stop recovery contract therefore cannot run unchanged here. This pass
   does not install another supervisor or weaken that contract.
2. **Julia entrypoint mismatch.** The pinned base image has Julia at
   `/usr/local/julia/bin/julia`; the current container recipe and policy target
   `/usr/local/bin/julia`, which was absent. The first launch failed at OCI exec.
   Only this temporary probe was corrected. Reconcile the recipe and policy
   together during Docker integration, then test the built production images.
3. **Attach input lifetime matters.** Inheriting EOF from the SSH command yielded
   no captured TTY output despite a successful process exit. The probe now keeps
   its stdin pipe open until process completion; both attached streams pass.
   This is fixture evidence, not a claim that the application's terminal driver
   has a defect or has completed Docker acceptance.

The probe explicitly invokes Docker with a small fixed policy. It reuses the
actual kernel verifier but does **not** run the production policy generator,
registration, managed resource owner, browser terminal, lease recovery or
scientific image builds. Memory/PID exhaustion, output flooding, REPL interaction,
proxy deployment and a remote NATS topology are not certified by this result.

## Evidence and reproduction

- [Measured report](test/codespaces_feasibility_2026-09-08.json).
- [Standalone operator probe](test/codespaces_probe.py); not part of automatic CI.
- Verifier SHA-256:
  `bfd87b9f0516412ed6c1581b519d6a33f3408fe0d1f891f77fe4dde565f272aa`.
- Probe SHA-256:
  `90ff896c82c943c568f2b2086022d117fb9ecaf4dbbb130bcb726e0f63270587`.
- Initial entrypoint failure, TTY diagnostic and final report are also retained
  locally in `/tmp/lcm-codespaces-feasibility.2eQpZb9R/`.

For another explicitly approved disposable Codespace, copy the probe and
`worker/core/src/ContainerIsolation.jl` into the same private temporary directory.
Explicitly pull the `IMAGE` digest declared in the probe, then invoke
`python3 /absolute/temporary/directory/codespaces_probe.py`. The probe does not
pull images, install dependencies or start application services. It refuses to
run without the Codespaces environment marker, uses an explicit local engine
socket and an empty Docker client configuration, and removes only its own
randomly named, label-verified containers. Retain `report.json` before deleting
the disposable Codespace. A future kernel/image may legitimately fail these
checks; do not weaken them to match this historical result.

The pre-existing `silver-umbrella-wvpwjqp5wxpc946j` Codespace was deleted with
explicit user authorization after GitHub reported no uncommitted or unpushed
changes. The separate test Codespace was also removed after exporting results;
the remote checkout remained clean. No local Docker installation, local Podman
change, source push, application deployment or production runtime change occurred.
