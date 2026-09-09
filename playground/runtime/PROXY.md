# Private trusted-user gateway

`Caddyfile.example` and `proxy.example.toml` are operator templates, not an
installed deployment. Use an approved host/domain and the owned `lcm runtime`
gateway, never the legacy `lcm playground start` listener. A private Caddy 2.10.2
rehearsal on the owned Kubuntu host has passed the Docker terminal/two-user gate;
this does not install or certify a public deployment. Consult the
[acceptance ledger](../RUNTIME_PLATFORM_PROGRESS.md) for each engine and scenario.

The template uses Caddy's HTTPS site handling and
[hashed Basic authentication](https://caddyserver.com/docs/caddyfile/directives/basic_auth).
It provisions two example identities: `operator` (runtime administrator) and
`researcher`. Replace/add accounts in the server-owned configuration, keeping
administrator names aligned with the runtime file. This is trusted-user access,
not public signup, SSO or a hostile-code sandbox.

## Configuration and validation

1. Copy both templates to private operator configuration; adjust all paths,
   the runtime origin and the approved hostname. Set `LCM_PUBLIC_HOST` to that
   hostname without a scheme. HTTPS is required; do not publish port 8080 or
   change its loopback listener to `0.0.0.0`.
2. Generate an independent high-entropy proxy key and store it in the runtime's
   mode-0600 `proxy_key_file`. Supply the identical value as `LCM_PROXY_KEY` only
   to the Caddy service through its private credential/environment mechanism.
   Keep it out of command arguments, shell history, source control and logs.
3. Run `caddy hash-password` interactively for each account. Supply the resulting
   hashes as `LCM_OPERATOR_PASSWORD_HASH` / `LCM_RESEARCHER_PASSWORD_HASH` through
   the same private service configuration. Never supply plaintext passwords in
   the Caddyfile. Do not reuse broker or artifact credentials.
4. Build the published site, then validate the copied runtime configuration.
   Validate the copied Caddyfile in the intended service environment:

   ```sh
   ./playground/lcm playground build
   ./playground/lcm runtime check --config /absolute/private/runtime.toml
   caddy validate --config /absolute/private/Caddyfile --adapter caddyfile
   ```

   Caddy installation, certificate provisioning, service installation and public
   exposure require operator approval. Do not paste adapted configuration into
   ordinary logs: environment expansion can materialize the private proxy key.

Start the runtime and approved Caddy service only after validation. The runtime
command is `./playground/lcm runtime start --config /absolute/private/runtime.toml`.
Worker configuration/provisioning follows [CONFIGURATION.md](CONFIGURATION.md);
authentication does not approve a worker, prepare a model or enable a terminal.

## Boundary and deployment checks

The two Caddy `handle` branches are exclusive. Private runtime/application paths
require authentication; public documents, capability discovery and shared assets
do not. Authenticated upstream requests replace the principal/key headers with
server-owned values. The public branch removes them entirely. Both remove Basic
credentials and any asserted UI-host key before forwarding. This uses Caddy's
[request matchers](https://caddyserver.com/docs/caddyfile/matchers) and
[upstream header controls](https://caddyserver.com/docs/caddyfile/directives/reverse_proxy#headers).
No URI rewriting or separate WebSocket route bypass is introduced.

The gateway independently verifies the actual proxy peer, private key, owner,
Origin and mutation marker. It does not trust `X-Forwarded-For` as a principal.
Keep the proxy and gateway local to each other for this profile; a remote agent
still connects outbound to NATS, not to a public Julia listener.

Before exposure, verify on the approved host:

- Anonymous home/deck/asset requests succeed without a login prompt or allocation.
- Private APIs and WebSocket upgrades reject missing/wrong credentials.
- `researcher` cannot enroll workers or attach to the operator's run/terminal.
- Forged principal, proxy-key and UI-host-key headers grant no authority.
- Foreign-Origin mutations/upgrades are rejected after authentication.
- Direct loopback requests without the proxy key cannot access private routes;
  port 8080 is not reachable from another machine.
- Stopping a run, disconnecting the broker and restarting the proxy preserve the
  documented ownership/cleanup behavior, not volatile Julia session memory.

The opt-in [physical acceptance group](PHYSICAL_ACCEPTANCE.md) runs the shipped
template with private loopback HTTPS and verified test certificates. It checks
two-user HTTP/WebSocket ownership and graceful/crash cleanup. Every deployment
still needs the hostname, certificate and exposure checks above; unit tests alone
are not evidence that a particular proxy has run.
