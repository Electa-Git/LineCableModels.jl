# Playground execution architecture

The playground has one strict boundary: Quarto and Bonito publish; a separate
worker executes registered scientific operations. Typed JSON messages cross a
NATS/JetStream broker. The browser never talks to NATS.

```text
Quarto + Bonito publisher
          |
          | typed JSON messages
          v
     NATS + JetStream
          |
          v
 LineCableModels worker
          |
          v
 supervised Julia executor
```

The publisher, broker, worker daemon, and executor are independently runnable.
The publisher starts and remains usable when the broker, workers, or optional
scientific packages are absent.

## Workbench boundary

The browser-hosted engineering workbench is a publisher mode, not an execution
environment. `WorkbenchUI` owns the semantic shell, persistent view mounting,
intrinsic panel interactions, shared visual language, and typed action
dispatch. A concrete workbench module owns its session-local state, views, and
`handle!` methods.

The reusable workbench module imports neither NATS nor scientific packages. A
domain application may receive a publisher-side broker capability explicitly,
but that adapter remains outside the structural shell. The standalone template
at `/workbenches/template` deliberately receives no broker client and performs
no numerical work.

Quarto documents and links to complete workbench routes. It does not wrap a
full workbench in an iframe: the application shell must own viewport geometry,
focus, splitters, and persistent rendering surfaces directly.

## Command-surface styling

`assets/toolbar.css` owns toolbar control geometry and the complete normal,
hovered, active, busy, disabled, and focus states. Both `Toolbar` and `Ribbon`
include this stylesheet themselves; neither depends on `widgets.jl` or the
gallery shell for its controls. `assets/ribbon.css` only arranges/sizes those
controls and styles ribbon-owned tabs, groups, and overflow surfaces. It must
not override a composed control's foreground/background state pair.

Colors remain semantic references to `assets/brand.css`. The host provides the
shared palette and theme preference; changing it updates the mounted controls
through CSS inheritance, without recreating controls or callbacks. The same
contract applies to ordinary groups, quick access, and overflow groups.

## Presentation boundary

The reusable presentation surface is a Quarto format, not a scientific
application and not a second workbench. Quarto and Pandoc own authoring and
document compilation. Reveal owns only deck navigation, URL state, overview,
speaker notes, fullscreen, and print hooks. The `lcm-deck` format disables
Reveal's layout, centering, scaling, and transitions; LCM layout primitives are
the sole owners of projection geometry.

Slides are laid out at their real browser pixel dimensions inside a contained
16:9 stage. CSS transforms must never scale an ancestor of an active live
iframe, canvas, or WebGL surface. Reveal's overview may transform slide
thumbnails only while live children are non-interactive and visually replaced
by inert placeholders; leaving overview restores the untransformed persistent
frame. Other projector ratios receive deterministic letterboxing. Overflow is
an authoring error rather than a reason to shrink an entire slide until it
becomes unreadable.

Live presentation content is mounted through the existing same-origin Bonito
iframe boundary. The deck page does not import Bonito, NATS, scientific
packages, or broker credentials. Live frames remain mounted while navigating;
slide-entry, slide-leave, viewport-settling, viewport-settled, and print-mode
messages let a child pause activity and resize only after geometry stabilizes.
Missing publishers, brokers, or workers do not prevent static navigation.

Speaker previews and print output never instantiate a second live application.
They display an explicit placeholder containing the live-view title and a
clickable playground URL. Version 1 deliberately does not generate Makie
snapshots or capture mutable browser state for PDF output.

The format grammar is intentionally narrow: a flat sequence of slides, named
LCM layouts, standard Quarto notes, and the validated `bonito` shortcode. Raw
per-slide styles, executable page scripts, nested Reveal slides, and arbitrary
external iframe sources are outside the contract. Presentation CSS consumes
the same `assets/brand.css` palette as the website and workbench; a second deck
palette is forbidden.

Published text behavior is owned by `assets/published-text.css`, loaded by both
Quarto formats. It hides the browsing caret without disabling mouse selection,
copying, focus, or text entry. Page/template/deck schemas must not duplicate
that policy. It is deliberately absent from Bonito widget and workbench
documents, whose components retain their own caret and selection semantics.

Reveal is retained only if the hostile presentation specimen proves exact
pointer coordinates, persistent live-frame identity, stable fullscreen and
multi-resolution sizing, focus-safe navigation, inert presenter previews, and
one-page print placeholders. If it fails, the Quarto grammar, semantic markup,
layouts, and live boundary remain unchanged and only the deck controller is
replaced by a custom HTML adapter.

## UI toolkit boundary

`src/toolkit/Toolkit.jl` is the shared presentation vocabulary below gallery
pages and workbenches. Its concrete controls, fields, forms, dialogs, notices,
toasts, disclosures, property grids, data tables, and viewport frames own their
DOM behavior and component-scoped CSS once. Gallery routes instantiate those
exact types; workbench leaves compose the same types directly.

`assets/brand.css` remains the only palette authority and
`assets/control-contract.css` remains the cross-surface native-control
authority. Toolkit styles consume those tokens without redefining them. X-ray
metadata is implemented beside each owned component. `SecretInput` is the
deliberate exception to ordinary reactive controls: its value does not exist in
Julia state or diagnostic metadata and crosses the browser boundary only when
the containing form is explicitly submitted.

## Enforced invariants

1. `playground/Project.toml` never depends on LineCableModels,
   PowerImpedance, PowerModels, or solver packages.
2. Page construction never submits work.
3. Bonito callbacks never synchronously wait for calculations.
4. The browser never connects directly to NATS or receives broker credentials.
5. Messages contain plain data and registered operation names, never Julia
   code, functions, expressions, or `Serialization` payloads.
6. Workers acknowledge durable jobs only after storing a terminal result.
7. The last successful result remains visible while newer work is pending.
8. Power-flow preparation is explicit and cached independently from impedance
   evaluation.
9. Missing capabilities produce a finite `unavailable` state, not a spinner.
10. Arbitrary REPL execution is a separate sandboxed feature and is not part of
    ordinary calculations.
11. Gallery and workbench compositions share toolkit implementations and
    semantic CSS contracts; gallery-only replicas are forbidden.
12. Reveal never owns presentation layout or scales live content.
13. Presentation and print pages never execute numerical work.
14. Every printable live viewport declares a public playground URL.
15. Presenter previews never mount a duplicate live iframe.

These rules are checked by `playground/test/architecture.jl`.

## Durable and transient traffic

JetStream stores job requests and terminal results. Core NATS carries progress,
logs, cancellation, capability announcements, and heartbeats:

```text
lcm.jobs.v1.<priority>.<operation>
lcm.results.v1.<job_id>
lcm.events.v1.<job_id>
lcm.control.v1.cancel.<job_id>
lcm.workers.v1.heartbeat.<worker_id>
lcm.workers.v1.capabilities.<worker_id>
```

Workers use pull consumers with explicit acknowledgement. They do not use
`JetChannel.take!`, because that helper acknowledges a message before the
scientific operation completes.

If a worker reaches NATS before an administrator has created the streams, it
remains alive and retries the bounded setup operation. Compose additionally
orders workers after a successful `runtime-init` service. Native and remote
startup therefore do not depend on a lucky process launch order.

`high` and `normal` jobs use distinct durable consumers. Workers poll every
high-priority operation before normal-priority work, while round-robin
scheduling prevents one operation from monopolizing its priority class.

## Cache boundary

Completed deterministic results are content-addressed and persistent. Prepared
resources—operating points, linearizations, factorizations, or network
models—are worker-local, leased, bounded, and reconstructible from immutable
input data. No Julia object crosses the broker.

Power-flow work is split into two named operations:

```text
powerflow.prepare(specification) -> prepared_resource_key
impedance.evaluate(prepared_resource_key, specification_hash, parameters)
    -> impedance result
```

The resource key is an opaque hint. The immutable specification remains the
authority and allows another worker to reconstruct lost state.

## Execution and failure boundary

Scientific operations run in supervised, persistent Julia executors. Each
worker-capacity slot owns one executor and reuses it across jobs. Cooperative
operations check a cancellation/deadline token between stages. A cancellation,
request deadline, operation timeout, broken protocol frame, or executor crash
terminates the affected executor and the daemon replaces it without losing its
broker connection.

Warning-level engine logger events are framed separately from protocol output,
streamed to the job console, bounded and deduplicated, and stored in the
terminal result. Detailed exception traces stay in worker logs; browser results
receive only a category, concise message, and diagnostic identifier.

## Artifact boundary

Inline JSON is capped at 64 KiB inside a wire-message ceiling of 256 KiB.
Larger results are stored by digest. The native profile uses a local
filesystem; remote and container profiles use MinIO or another S3-compatible
service. Workers have write-only object credentials and the publisher has
read-only credentials. Browser clients retrieve artifacts only through
`/artifacts/sha256/<digest>`, which supports `GET`, `HEAD`, and single byte
ranges without exposing storage credentials.

## Container resource boundary

All supplied profiles use read-only container roots, dropped capabilities,
`no-new-privileges`, bounded process counts, and memory limits. CPU quotas live
in `compose.cpu-limits.yml` so deployments with a delegated CPU cgroup can
enforce them without making the base profile unusable on rootless Podman hosts
that lack that controller. A future arbitrary-code session must require this
CPU override (or an equivalent orchestrator quota) in addition to its separate
identity, namespace, filesystem, network, and lease restrictions.

The complete requirement-to-evidence matrix and repeatable verification
commands are maintained in [`VERIFICATION.md`](VERIFICATION.md).
