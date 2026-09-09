# End-to-end verification

This document lists repeatable checks for the retained publisher/worker stack
and Runtime Platform v1. The current v1 release status and retained evidence are
in [RUNTIME_PLATFORM_PROGRESS.md](RUNTIME_PLATFORM_PROGRESS.md); the existence
of a harness is not a passing release gate. Source inspection, finite test
drivers and effective deployment isolation are distinct evidence.

## Publisher, legacy worker and presentation compatibility

The table below preserves the original stack's compatibility checks. Its
unassigned job path is not exposed by the authenticated v1 gateway. The new
runtime has separate ownership, profile, terminal and acceptance contracts
described below and in [ARCHITECTURE.md](ARCHITECTURE.md).

| Plan section | Implemented boundary | Authoritative evidence |
|---|---|---|
| 1. Architectural invariants | Publisher contains no engine packages or execution path; page construction is inert; browser assets contain no broker credentials; results remain visible while work is pending. | `ARCHITECTURE.md`, `test/architecture.jl`, `test/jobhandle.jl`, offline route/shutdown harness. |
| 2. Separate environments | Publisher, protocol, and worker have independent `Project.toml` and `Manifest.toml` files. Only the worker declares LineCableModels and PowerImpedance. | `Project.toml`, `protocol/Project.toml`, `worker/Project.toml`, architecture dependency assertions. |
| 3. Message contract | Versioned JSON-only requests, ordered events, results, failures, artifact references, and heartbeats have bounded fields and payload sizes. | `protocol/src/`, 43 protocol assertions, four golden JSON payloads. |
| 4. NATS subjects and durability | Jobs and results use JetStream; events, logs, cancellation, heartbeats, and capabilities use Core NATS. Pull consumers ACK only after terminal storage/publication. | `src/broker/Subjects.jl`, `worker/src/Consumer.jl`, `src/broker/RuntimeAdmin.jl`, lifecycle redelivery and duplicate tests. |
| 5. Non-blocking broker client | Connection supervision, reconnect, discovery, submit, cancel, retry, and reattach run outside page construction. Offline submission fails immediately. | `src/broker/BrokerClient.jl`, `test/broker_lifecycle.jl`, offline publisher harness. |
| 6. Generic UI | Shared toolkit primitives compose identically in gallery routes and workbench leaves; one `JobPanel` owns Run/Cancel/Retry, dirty state, progress, console, cache status, retained result, and a persistent result callback. | `src/toolkit/`, `src/widgets/JobControls.jl`, `src/broker/JobHandle.jl`, `test/toolkit.jl`, `test/visual_contracts.jl`, browser visual-contract audit, `test/jobhandle.jl`. |
| 7. Diagnostic worker lifecycle | Echo, delay, executor delay, progress, warning, and requested failure operations exercise the full transport before science. | `worker/src/operations/diagnostics.jl`, worker tests, broker lifecycle integration. |
| 8. Closed operation registry | Every operation declares validation, schema, timeout, capability, cache policy, and execution mode. Unknown/eval operations are rejected. | `worker/src/OperationRegistry.jl`, registry and architecture tests. |
| 9. Result and prepared caches | Persistent content-addressed results are separate from bounded, leased, single-flight worker-local preparations. | `worker/src/Cache.jl`, cache tests, local/remote repeated-request smokes. |
| 10. Real engine operations | Geometry, cable constants, frequency sweeps, and line scans call actual LineCableModels constructors and compute APIs. | `worker/src/operations/linecablemodels.jl`, 11 direct package-parity assertions. |
| 11. Prepare/evaluate split | `powerflow.prepare` owns network, power flow, and linearization; `impedance.evaluate` reuses or reconstructs it and varies only declared passive inputs. | `worker/src/operations/powerimpedance.jl`, nine heavy PowerImpedance assertions. |
| 12. Cancellation and crash isolation | Cooperative checks, request deadlines, per-operation timeouts, and hard executor replacement coexist with JetStream redelivery after daemon death. | `worker/src/Executor.jl`, worker cancellation/deadline tests, lifecycle worker/executor kill tests. |
| 13. Results and artifacts | Small results remain inline; larger results use local or S3 content-addressed storage behind a same-origin GET/HEAD/range gateway. | `worker/src/Artifacts.jl`, `src/artifacts.jl`, artifact tests, TLS MinIO role-isolation harness, both stack smokes. |
| 14. Local, remote, and container commands | The same `lcm worker start` command runs natively or in a locked image; `lcm container` resolves real Docker versus Podman and drives either local or remote Compose profile; publisher, NATS administration, and worker remain separate processes behind one CLI. | `lcm`, `src/container_runtime.jl`, both OCI-compatible Dockerfiles, local/remote Compose profiles, CLI resolver tests, local and mTLS/S3 stack smokes. |
| 15. Authentication and authorization | Publisher, worker, administrator, artifact writer, and artifact reader use distinct identities and least-privilege subjects/policies. Remote traffic requires verified mTLS/TLS. | `deploy/nats.conf`, `deploy/remote/nats-tls.conf`, MinIO policies, authorization, mTLS, and TLS artifact tests. |
| 16. Scientific subjects exclude REPL evaluation | No eval/repl operation or general code payload exists on scientific subjects. The v1 private terminal uses a separate leased process and authorized byte stream; admission requires current verified host/image isolation. | Architecture source scan, `runtime/CONFIGURATION.md`, terminal process/relay/browser tests, both physical engine matrices and the v1 implementation ledger. |
| 17. Test plan | Protocol, resilience, job semantics, scientific parity, native lifecycle, container, remote TLS, authorization, artifacts, and graceful shutdown all have executable harnesses. | Commands below. |
| 18. Delivery sequence | Transport and diagnostics precede engine adapters; supervised execution precedes PowerImpedance; security/deployment profiles are separate checks. | Repository layering and retained per-gate evidence; see the v1 ledger for incomplete release gates. |
| 19. Presentation ownership | Quarto compiles, Reveal orchestrates, LCM layouts own real-pixel geometry, and Bonito remains isolated in same-origin live frames. | `ARCHITECTURE.md`, `_extensions/lcm-deck/`, presentation architecture assertions. |
| 20. Presentation lifecycle | Live frames remain mounted; resize/fullscreen changes settle before child notification; focused controls retain keyboard ownership; overview substitutes inert placeholders and restores exact geometry without replacing sessions. | `presentations/specimen.qmd`, `test/integration/presentation_browser.mjs`. |
| 21. Presentation failure and print | Static navigation works without a broker; presenter/print surfaces use public playground links and never duplicate or print live applications. | Presentation browser checks and generated print-DOM assertions. |
| 22. Presentation UX regression | Actual menu overview/selection/exit and PDF entry points; stage and slot overflow at five resolutions; cached theme and live theme parity; visible loading/readiness; laser keyboard ownership and iframe coordinates. | `test/integration/presentation_browser.mjs`; real Chrome direct-print and preview-print PDFs plus CLI PDF, checked for nine pages and two placeholders. |
| 23. Ribbon/toolbar theme regression | Repeated light/dark changes from the publisher update mounted controls in the gallery, iframe, standalone and workbench; active icon/text/background pairs, hover, busy, disabled, all sizes, quick access, overflow and callbacks retain their shared contract without gallery CSS. | `test/ribbon.jl`, `test/integration/ribbon_fixture.jl`, `test/integration/ribbon_theme_browser.mjs`, `bash test/integration/run-ribbon.sh`. |
| 24. Interaction-state regression | Template and workbench navigation share semantic selection styles. Real pointer and keyboard round trips cover both themes, expanded/collapsed navigation, rail tooltips, disabled items, ribbon and persistent tabs, dock tabs, segmented choices, selectable data rows and toolbar actions. No simulated hover classes are allowed in owned component sources. | `test/visual_contracts.jl`, `test/integration/interaction_state_contract.mjs`, `bash test/integration/run-ribbon.sh`. |

The equation explanation gate is included in `run-presentation.sh`:
`math_notes_filter.mjs` checks real Pandoc AST output and invalid authoring;
`math_notes_browser.mjs` checks both palettes, exact equation geometry, term
switching, keyboard and selection ownership, resizing and all four stage edges,
overview, printing, MathJax rerendering, teardown and dependency failure/recovery.
Static receiver/print modes are also checked by the presentation browser suite.

The same gate runs `published_shell_browser.mjs` against every published sidebar
destination and the presentation authoring guides, in both themes. It checks
natural document scrolling, complete navigation/footer contours, absence of
nested page scroll areas and horizontal overflow, and compact-menu open/close
at five widths. `incremental_lists_browser.mjs` checks the starter's nested list
and specimen's two-list sequence using native Reveal controls, forward/backward
steps, stable reserved geometry, both themes and all-visible print output.

The presentation gate covers the supplied hostile specimen, not arbitrary
overfull authored content. It runs in isolated Chromium; physical multi-monitor
window movement and the OS print dialog remain manual checks. Overview uses
Reveal's thumbnail strip, with the current slide centered and arrow navigation.
PDF preview is a static reload; only overview and the browser's direct print
lifecycle retain an existing live session.

## Repeatable commands

Runtime Platform v1 has one aggregate entry point, using the existing isolated
harnesses and retaining a private per-gate log and TSV report:

```sh
bash playground/test/integration/run-runtime-platform.sh all
# Focused groups: unit, browser, transport, scientific, host. `list` is read-only.
# `physical` is a separate opt-in group requiring approved installed images.
CONTAINER_RUNTIME=docker bash playground/test/integration/run-runtime-platform.sh transport
```

The default test engine is Podman. `host` checks native, genuine Podman and
genuine Docker independently through the existing CLI; a Docker-named Podman
shim is not Docker evidence. Transport/scientific fixtures reuse the inspected
local command prefix and filtered environment for every engine action, including
cleanup; inherited remote contexts cannot redirect a later fixture command.
Exit 1 means a failed automated gate; exit 2 means
unavailable host prerequisites. A selected-group pass does not certify effective
executor isolation or a physical second-computer deployment. Follow
[physical acceptance](runtime/PHYSICAL_ACCEPTANCE.md) for both engine matrices
and the separate-computer consumer gate. Consult
[the implementation ledger](RUNTIME_PLATFORM_PROGRESS.md) for recorded evidence;
neither this command nor a private rehearsal is a public deployment.

`test/runtime_conformance.jl` derives the new runtime/scientific component inventory
from actual `Bonito.jsrender` and X-ray inspection dispatch, expands abstract
families and requires normal-constructor fixtures. It checks owned CSS, metadata
redaction and repeated session rendering with diagnostics enabled/disabled.
Missing render/inspection/fixture coverage fails. Existing pre-v1 broker widgets
retain their separate legacy tests. Browser theme/interaction coverage remains
in the real Bonito, scientific UI, ribbon and X-ray harnesses; serialization tests
are not represented as visual browser passes.

The X-ray CSS preview gate is `bash playground/test/integration/run-xray.sh`.
It starts an isolated browser profile and mock Julia host, without NATS or
numerical workers. It exercises click-only selection, conditional CSS and
instance isolation (including duplicated portable stylesheets), both themes,
live binding updates without editor replacement, validation, source-aware
export, reset, teardown and read-only host policy, plus the real workbench.
The nested-tree fixture covers component-only, recursive, hidden-descendant,
and global recovery, sibling isolation, invalid-only drafts, and original-to-
proposed comparisons that survive selection and theme changes.
`playground/test/xray_preview.jl` covers the shared Julia editor contract.

The publisher unit suite also checks that disabled X-ray does not call inspection
hooks and that the build-time UI workload leaves listener/cleanup/upload state
unchanged without importing numerical engines. The offline publisher gate uses
actual GETs, checks every published index alias and the existing widget probes,
then verifies graceful SIGINT shutdown. Cold-render measurements and the latest
whole-suite results are recorded separately in the runtime implementation ledger.

For startup failures, `LCM_RUNTIME_TRACE_COMPILE=1` enables retained compiler
timings in `runtime/test/run-broker.sh` (full/terminal/scientific) and the legacy
`test/integration/run.sh` worker fixture. The latter uses the same Julia module
and arguments directly for tracing; run without this flag for actual CLI
acceptance. Failed legacy workers retain private diagnostic logs before normal
resource cleanup. Browser terminal failures retain screenshots, owner-filtered
assignment/inventory state and exceptions; no transport deadline or automatic
input replay is enabled by diagnostics.

From the repository root:

```sh
julia --startup-file=no --project=playground/protocol \
  playground/protocol/test/runtests.jl
julia --startup-file=no --project=playground \
  playground/test/runtests.jl
julia --startup-file=no --project=playground/worker \
  playground/worker/test/runtests.jl
LCM_TEST_POWERIMPEDANCE=1 julia --startup-file=no \
  --project=playground/worker playground/worker/test/runtests.jl

./playground/test/integration/run.sh
./playground/test/integration/run.sh --tls
./playground/test/integration/run-artifacts.sh --tls
./playground/test/integration/graceful-shutdown.sh

lcm presentation check playground/presentations/specimen.qmd
lcm presentation build playground/presentations/specimen.qmd
lcm presentation export playground/presentations/specimen.qmd --pdf
./playground/test/integration/run-presentation.sh
bash playground/test/integration/run-ribbon.sh
bash playground/test/integration/run-xray.sh
```

The Compose profiles add deployable-artifact coverage:

```sh
lcm container resolve
lcm container start
lcm container status
./playground/test/integration/run-local-stack-smoke.sh playground/deploy/.env

./playground/deploy/remote/generate-dev-certs.sh
lcm container start --remote
lcm container status --remote
./playground/test/integration/run-remote-stack-smoke.sh playground/deploy/remote/.env
```

Run the same legacy Compose commands with `--runtime docker` and
`--runtime podman` to pin a host explicitly. Their `--cpu-limits` option composes
the optional legacy quota overlay. This is **not** the Runtime Platform v1
executor policy: v1 requires its declared effective limits and refuses an
unsupported host; private terminals never use an unconfined native fallback.

The remote smoke uses separate publisher and worker containers/network
identities and the same outbound-only worker connection used on another
computer. Moving that worker to a physical host changes certificates and
environment values, not application code or protocol.
