# Remote TLS deployment

Remote publishers and workers use the same NATS protocol as the local profile,
but client connections require mutually authenticated TLS. The server verifies
the client certificate and NATS still applies the publisher, worker, and
administrator subject permissions independently.

For a local TLS rehearsal only, generate disposable 30-day certificates:

```sh
./generate-dev-certs.sh
cp .env.example .env
# Replace every NATS and MinIO secret in .env.
lcm container start --remote
```

Runtime selection is automatic. It can be pinned on either kind of host:

```sh
lcm container start --remote --runtime docker
lcm container start --remote --runtime podman
```

Memory and process-count limits are part of the base profile. Where the
container host supports CPU cgroup delegation, apply CPU quotas as well:

```sh
lcm container start --remote --cpu-limits
```

The separate override is required for portability: some rootless Podman user
services do not receive the CPU controller and must reject, rather than fake,
an unconditional CPU quota.

Production deployments must replace the generated files with certificates from
the organization's PKI, keep private keys outside Git, restrict key-file
permissions, and expose port 4222 only through an appropriate firewall or
private network. Clients use `tls://` and provide `NATS_TLS_CA_PATH`,
`NATS_TLS_CERT_PATH`, and `NATS_TLS_KEY_PATH`. The Julia publisher and worker
read those variables directly through NATS.jl; neither opens an inbound Julia
port. The connection URL remains `tls://` from configuration through the TLS
handshake. If the network address is not present in the server certificate,
set `NATS_TLS_SERVER_NAME` to a DNS identity that is present; hostname and CA
verification remain mandatory. The repository pins an auditable NATS.jl 0.1.1
transport patch in `playground/vendor/NATS`; see its `LCM_PATCHES.md`.

The rehearsal profile also starts a TLS-protected MinIO server and creates two
least-privilege storage identities: the worker can write content-addressed
objects and the publisher can read them through its same-origin artifact
gateway. The MinIO API binds to `127.0.0.1` by default, so running this example
does not expose object storage to the network.

For a worker on another computer, use an organization-managed S3-compatible
endpoint with a certificate for its real DNS name. Configure both processes
with:

```text
LCM_ARTIFACT_BACKEND=s3
LCM_S3_ENDPOINT=https://objects.example.org
LCM_S3_BUCKET=lcm-artifacts
LCM_S3_PREFIX=linecablemodels
AWS_ACCESS_KEY_ID=...
AWS_SECRET_ACCESS_KEY=...
AWS_REGION=us-east-1
```

Give the worker identity `PutObject` only for the configured prefix and the
publisher identity `GetObject` only. Do not expose the MinIO root identity to
either process. The included MinIO service is a local rehearsal; publishing it
through a firewall, load balancer, or reverse proxy is an infrastructure
administrator decision and is intentionally not enabled by this repository.
