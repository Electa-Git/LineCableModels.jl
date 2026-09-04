# End-to-end verification

This matrix audits the implementation against the complete execution plan.
It distinguishes source inspection from runtime evidence so a narrow unit test
cannot be mistaken for proof of distributed behavior.

## Requirement matrix

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
| 16. Stateful REPL deferred safely | No eval/repl operation or general code payload exists. A future REPL requires a separate leased sandbox and identity. | Architecture source scan and documented container boundary. |
| 17. Test plan | Protocol, resilience, job semantics, scientific parity, native lifecycle, container, remote TLS, authorization, artifacts, and graceful shutdown all have executable harnesses. | Commands below. |
| 18. Delivery sequence | Transport and diagnostics precede engine adapters; supervised execution precedes PowerImpedance; security/deployment profiles are last. | Repository layering and the complete passing verification set. |

## Repeatable commands

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

Run the same commands with `--runtime docker` and `--runtime podman` to pin a
host explicitly. `--cpu-limits` composes the optional quota overlay; it remains
off by default because some rootless Podman services lack delegated CPU
cgroups.

The remote smoke uses separate publisher and worker containers/network
identities and the same outbound-only worker connection used on another
computer. Moving that worker to a physical host changes certificates and
environment values, not application code or protocol.
