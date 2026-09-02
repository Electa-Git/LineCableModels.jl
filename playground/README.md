# LineCableModels playground

This directory contains three independently runnable programs and one shared
wire contract. Quarto owns static authoring, Bonito owns live browser controls,
NATS/JetStream carries typed jobs, and a separate worker owns LineCableModels,
PowerImpedance, and solver dependencies. The publisher starts no numerical
workers and remains usable when NATS or all workers are offline.

The enforced process and cache boundaries are recorded in
[`ARCHITECTURE.md`](ARCHITECTURE.md). The end-to-end requirement and test
matrix is in [`VERIFICATION.md`](VERIFICATION.md).

## Bootstrap

From the repository root, run:

```sh
./playground/bootstrap.sh
```

The idempotent bootstrap:

- requires Julia 1.12 or newer;
- downloads the pinned Quarto CLI into the user's XDG data directory;
- verifies the official archive SHA-256 before extracting it;
- instantiates the Julia environment;
- installs the single `lcm` command;
- exposes it through `~/.local/bin`; and
- renders the initial static site.

It uses no `sudo` and does not modify shell startup files. Linux x86-64 and
ARM64 are supported. Rerunning it reuses a valid installation.

If `~/.local/bin` is not already on `PATH`, the bootstrap prints the exact
export line to add. The working-tree source remains live, so editing the Julia
package does not require reinstalling the command.

## Author and build

The home is authored in `index.qmd`; `_quarto.yml` owns site configuration and
`assets/theme.scss` owns the visual treatment. Render it with:

```sh
lcm playground build
```

The generated `_site/` directory is disposable and ignored by Git.

## Reuse a page template

Open `/templates/` in the running site to inspect the rendered gallery. Every
entry is backed by the corresponding file in `templates/`; its **Source**
control shows the complete authoring source. Copy the closest `.qmd` file,
change its title and contents, and retain the `lc-layout-*` classes you need.

The gallery currently covers a full canvas, balanced columns, a feature with a
sidebar, a top band over two columns, a 2×2 dashboard, media with narrative,
and the intended responsive Bonito viewport boundary.

## Reuse a live widget

Open `/widgets/` in the running site for a catalog of session-local Bonito
controls. Quarto owns the surrounding explanation and layout; each Bonito app
is mounted through a responsive, same-origin iframe:

```markdown
{{< bonito route="/widgets/control-panel" title="Cable controls" height="34rem" >}}
```

The local shortcode is implemented in `_extensions/bonito/`. Reusable widget
factories and their registered routes live in `src/widgets.jl`. The examples
cover sliders, numeric spinners, checkboxes, dropdowns, text fields, actions,
progress feedback, a read-only console sink, persistent tabs, namespaced
toolbars, and a composed engineering control panel. The asynchronous job panel
adds worker discovery, Run/Cancel/Retry, durable lifecycle state, retained
previous results, and broker messages without submitting work during page
construction. Layout primitives are consumed directly from `BonitoWidgets`;
the console accepts observable Julia, worker, or broker messages and starts no
execution machinery of its own.

## Build a complete workbench

Open `/workbenches/` for the architectural boundary and
`/workbenches/template` for the complete standalone demonstration. Unlike a
small live widget, a workbench is not embedded in a Quarto iframe: it owns the
browser viewport, persistent views, navigation, inspector, output dock, and
status bar.

Reusable structural types and rendering live under `src/workbench/`. A local
application module under `src/workbenches/` defines its session state,
composes those primitives, and implements typed `handle!` methods. The
template imports no broker or scientific package and performs no work merely
because its page was opened.

## Start

```sh
lcm playground start
```

`start` renders the site, registers the generated static files beside Bonito's
live routes, starts the server, and opens the
default browser. Configuration can come from `HOST`, `PORT`, and `PROXY_URL`,
or explicit options:

```sh
lcm playground start --host 127.0.0.1 --port 8080
lcm playground start --no-open
lcm playground start --no-render
```

Use `--no-render` only when `_site/index.html` already exists. Press `Ctrl+C`
to stop the server cleanly.

## Start the broker and worker

For native development, start a JetStream-enabled NATS server, initialize the
two streams with administrator credentials, then start the worker:

```sh
lcm nats init
lcm worker start
lcm nats status
```

`lcm nats` administers the playground's JetStream streams; it does not replace
or launch the external NATS server. All commands read `NATS_CONNECT_URL`;
administration may use the separate `NATS_ADMIN_URL`. The worker also accepts
`--nats`, `--id`, and `--capacity`. It connects outbound and exposes no Julia
listening port.

For an isolated local deployment:

```sh
cp playground/deploy/.env.example playground/deploy/.env
# Replace every placeholder in playground/deploy/.env, then:
lcm container resolve
lcm container start
```

`start` automatically selects a real Docker Engine on a Docker-first machine
or native Podman on a Podman host. The choice is also explicit and repeatable:

```sh
lcm container start --runtime docker
lcm container start --runtime podman
```

Set `LCM_CONTAINER_RUNTIME=auto|docker|podman` to choose a persistent default.
A `docker` command implemented by Podman's compatibility shim is recognized as
Podman; explicitly requesting Docker in that situation fails with a useful
diagnostic instead of pretending a Docker daemon exists.

The remaining lifecycle commands use the same resolved runtime and fixed
Compose project name:

```sh
lcm container status
lcm container logs --follow --tail 100 worker publisher
lcm container stop
```

`lcm container stop --volumes` also removes persistent deployment data and is
therefore intentionally opt-in. Add `--remote` to any lifecycle command to use
the mTLS/S3 profile under `deploy/remote`.

The base profile always applies memory and process-count limits. On Docker, or
on a Podman host whose user service has a delegated CPU cgroup controller, add
the CPU-quota override:

```sh
lcm container start --cpu-limits
```

Keeping CPU quotas in an explicit override lets the same base profile run under
rootless Podman installations that cannot enforce CPU cgroups; it does not
silently pretend a quota was applied.

Publisher, worker, and administrator credentials are distinct. Never expose
the example credentials or administrator URL to browser JavaScript.

## Distributed lifecycle test

Run the isolated source suites with:

```sh
julia --startup-file=no --project=playground/protocol \
  playground/protocol/test/runtests.jl
julia --startup-file=no --project=playground \
  playground/test/runtests.jl
julia --startup-file=no --project=playground/worker \
  playground/worker/test/runtests.jl
LCM_TEST_POWERIMPEDANCE=1 julia --startup-file=no \
  --project=playground/worker playground/worker/test/runtests.jl
```

The integration harness starts an isolated official NATS container with
ephemeral credentials and port, initializes JetStream, launches disposable
workers, exercises worker and broker failure, verifies redelivery,
deduplication, reconnection, reattachment, and subject permissions, and removes
all test resources afterward:

```sh
./playground/test/integration/run.sh
./playground/test/integration/run.sh --tls
./playground/test/integration/run-artifacts.sh
./playground/test/integration/run-artifacts.sh --tls
./playground/test/integration/graceful-shutdown.sh
```

Set `CONTAINER_RUNTIME=podman` for the low-level integration harnesses when the
container command is not named `docker`. The artifact harness starts a pinned MinIO server, creates separate
write-only worker and read-only publisher identities, uploads through the
worker store, and downloads both a full object and a byte range through the
publisher gateway.

After starting the local Compose profile, verify the complete container path
(native-auth broker, supervised scientific executor, shared-filesystem
artifacts, same-origin retrieval, and result-cache reuse) with:

```sh
./playground/test/integration/run-local-stack-smoke.sh \
  playground/deploy/.env
```

After starting the remote Compose profile, verify the complete container path
(mutual-TLS broker, supervised scientific executor, S3 artifact storage,
same-origin retrieval, and result-cache reuse) with:

```sh
LCM_PUBLISHER_PORT=18080 \
  ./playground/test/integration/run-remote-stack-smoke.sh \
  playground/deploy/remote/.env
```

## Remote workers

The remote profile in [`deploy/remote/`](deploy/remote/) requires mutual TLS
for every publisher, worker, and administrator connection while retaining the
same subject-level permissions. Workers connect outbound to NATS and expose no
Julia listener. The included certificate generator is strictly for a local
rehearsal; production certificates must come from the deployment's PKI.

This follows NATS's documented separation between transport encryption and
subject authorization: the TLS block secures client connections, while each
NATS identity remains limited to its explicit publish and subscribe subjects.
See the upstream [TLS](https://docs.nats.io/learn/security/encryption) and
[authorization](https://docs.nats.io/learn/security/authorization) guidance.
