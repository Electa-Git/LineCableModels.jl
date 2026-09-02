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
